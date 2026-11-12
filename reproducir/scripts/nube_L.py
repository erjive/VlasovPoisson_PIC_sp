"""La nube fuera de equilibrio con momento angular distribuido.

Dato inicial: F_0(Q, J) = exp(-sin^2(Q/2)/s_Q^2) J^2 exp(-J^2/s_J^2), s_Q = s_J = 0.1, por
exp(-(L - L_0)^2/s_L^2), L_0 = 2, en el isócrono, con autogravedad. Con L igual para todas las
partículas es la nube de vlasov-poisson_PIC (reproducir/scripts/bulto.py), que deja un núcleo denso
en la órbita circular cuando la masa pasa de un 2 % de la central. El modelo reducido de ese
repositorio (reproducir/scripts/reducido.py, docs/modelo_reducido) predice que una dispersión de L
del orden del 1 % deshace el núcleo. Aquí se comprueba con este código.

Las partículas son las del modelo reducido: la retícula de puntos medios en (Q, J) recortada donde
el peso baja de 1e-7 de su máximo, repetida en N_L puntos medios de |L - L_0| <= ancho s_L. Entran
con state = checkpoint. Las medidas son las del modelo reducido (reducido.medidas_plano), tomadas en
el plano de zeta = sqrt(2J) exp(iQ) con las variables de ángulo-acción del potencial de cada
instante y la L de cada partícula. Las de bulto.py en el plano (r, p_r) solo valen con L igual: con
L distribuida cada partícula oscila alrededor de su propio radio R(L).

La retícula en L. Las N_L capas vuelven a estar en fase en t_rec = 2 pi/(|dOmega/dL| dL), con dL el
paso de la retícula, y esa recurrencia se nota desde mucho antes: desde t_rec - 70/s_L para una
coherencia de 1e-3. Las corridas del modelo reducido del 2026-11-11 (N_L = 4, 6, 10 y 20 puntos en
|L - L_0| <= 2 s_L, t_rec = 10417 para s_L = 0.02 y 0.04) están afectadas a partir de unas 20 tau_1:
la suma libre de las capas, que para una gaussiana continua sigue cayendo, llega a 0.4 (s_L = 0.01),
0.9 (0.02) y 0.7 (0.04) al final de la corrida. coherencia_libre() lo calcula. Las corridas de
esta serie usan |L - L_0| <= 2.5 s_L y N_L = 8, 12, 20 y 30, con las que la suma libre no se aparta
de la continua más de 7e-4 en toda la corrida. Las que terminan en v repiten las retículas del
2026-11-11, para comparar con aquellas corridas del modelo partícula a partícula.

    python3 nube_L.py preparar                  condiciones iniciales y .par
    python3 nube_L.py tiempo [--solo ...]       corrida corta de cada caso y costo estimado
    python3 nube_L.py correr [--solo ...]       lo que falte, una corrida detrás de otra, 4 hilos
    python3 nube_L.py analizar [--solo ...]     mapa (Q, J) de cada instantánea y medidas
    python3 nube_L.py tanda <nombre>            corre y analiza una tanda (TANDAS), en orden
    python3 nube_L.py consistencia              C021 frente a la corrida B021 de vlasov-poisson_PIC
    python3 nube_L.py tabla                     medidas frente al modelo reducido
    python3 nube_L.py reticulas                 coherencia libre de la retícula en L de cada corrida

Necesita el repositorio vlasov-poisson_PIC junto a este (o en la variable VP_PIC_FIJO): de él salen
la retícula, las medidas y las corridas con las que se compara. Todo queda en exe/rep/12_nube/.
"""
import os, sys, subprocess, time, argparse, warnings, hashlib, datetime
for _k, _v in (('OMP_NUM_THREADS', '4'), ('OMP_PLACES', 'cores'), ('OMP_PROC_BIND', 'close')):
    os.environ.setdefault(_k, _v)
import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.abspath(os.path.join(AQUI, '..', '..'))
FIJO = os.environ.get('VP_PIC_FIJO', os.path.join(os.path.dirname(RAIZ), 'vlasov-poisson_PIC'))
if not os.path.isdir(os.path.join(FIJO, 'reproducir', 'scripts')):
    sys.exit(f'no encuentro vlasov-poisson_PIC en {FIJO}; defina VP_PIC_FIJO')
sys.path.insert(0, os.path.join(RAIZ, 'tools'))
sys.path.insert(0, os.path.join(FIJO, 'reproducir', 'scripts'))
warnings.filterwarnings('ignore', category=RuntimeWarning)
import bulto as B
import reducido as R
from demo_eta import envolvente
from aa_numerico_L import MapaAA, phi_iso, aa_isocrono
from raw_io import CorridaRaw

EXE = os.path.join(RAIZ, 'exe')
GRUPO = '12_nube'
BASE = os.path.join(EXE, 'rep', GRUPO)
PARDIR = os.path.join(RAIZ, 'reproducir', 'corridas', GRUPO)
FIJO_BULTO = B.BASE                             # corridas de L fija: vlasov-poisson_PIC/exe/demo_eta/bulto
B.BASE = BASE                                   # las medidas de bulto.py, sobre las corridas de esta serie
L0, SJ, SQ = 2.0, 0.1, 0.1
TFIN, DT, SALIDA = 10000, 0.05, 200
TAU1, KAPPA = B.TAU1, R.KAPPA
TS = [1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18, 20, 25, 30]

# nombre: masa del gas, retícula completa (N_J, N_Q), puntos, anchura s_L y semiancho (en s_L) de la
# retícula en L, corrida del modelo reducido y corrida de L fija con las que se compara, y cambios
# en el .par.
def _c(m, nl, sl, red=None, fijo=None, ancho=2.5, nj=200, nq=50, **cambios):
    return dict(m=m, nj=nj, nq=nq, nl=nl, sl=sl, ancho=ancho, red=red, fijo=fijo, cambios=cambios)


CORRIDAS = {
    # Consistencia entre los dos códigos: las partículas de B021, todas con L = L_0.
    'C021':        _c(0.021, 1, 0.0, 'R021', 'B021', nj=400, nq=100),
    # Barrido en s_L con la retícula de 200 x 50 en (Q, J) (150 x 14 tras el recorte).
    'M030_s000':   _c(0.03, 1,  0.0,   'rL0_0.03'),
    'M030_s005':   _c(0.03, 8,  0.005, 'rL005b_0.03'),
    'M030_s010':   _c(0.03, 12, 0.01,  'rL010b_0.03'),
    'M030_s020':   _c(0.03, 20, 0.02,  'rL020b_0.03'),
    'M030_s040':   _c(0.03, 30, 0.04,  'rL040b_0.03'),
    'M010_s000':   _c(0.01, 1,  0.0,   'rL0_0.01'),
    'M010_s020':   _c(0.01, 20, 0.02,  'rL020b_0.01'),
    # Convergencia de s_L = 0.02: más puntos en L, y mitad del paso radial con el mismo paso de tiempo.
    'M030_s020n':  _c(0.03, 30, 0.02),
    'M030_s020r':  _c(0.03, 20, 0.02, dr=0.05, courant=2.0),
    # Las retículas en L del modelo reducido del 2026-11-11 (afectadas por la recurrencia; ver arriba).
    'M030_s005v':  _c(0.03, 4,  0.005, 'rL005_0.03', ancho=2.0),
    'M030_s010v':  _c(0.03, 6,  0.01,  'rL010_0.03', ancho=2.0),
    'M030_s020v':  _c(0.03, 10, 0.02,  'rL020_0.03', ancho=2.0),
    'M030_s040v':  _c(0.03, 20, 0.04,  'rL040_0.03', ancho=2.0),
    'M010_s020v':  _c(0.01, 10, 0.02,  'rL020_0.01', ancho=2.0),
}
# Tandas, en el orden en que conviene correrlas. Costo medido con `tiempo` (corridas) y estimado (análisis).
TANDAS = {
    'primera':      ['M030_s000', 'M030_s020', 'M030_s010'],
    'resto':        ['M030_s040', 'M030_s005', 'M010_s000', 'M010_s020'],
    'convergencia': ['M030_s020n', 'M030_s020r'],
    'modelo':       ['M030_s020v', 'M030_s010v', 'M030_s040v', 'M030_s005v', 'M010_s020v'],
}


def ruta(*p):
    return os.path.join(BASE, *p)


def isocrono_rp(Q, J, L):
    """(Q, J, L) -> (r, p_r) en el isócrono (G = M = b = 1). Es la inversión de bulto.dato_isocrono
    (Binney y Tremaine, sección 3.5.2) con L por partícula: s = 2 + c (1 - e cos eta),
    r^2 = s (s - 2), Q = eta - e c/(c + 1) sin eta, c = (J + c_L)^2 - 1, c_L = (L + sqrt(L^2 + 4))/2."""
    E = -0.5/(J + 0.5*(L + np.sqrt(L*L + 4)))**2
    c = -0.5/E - 1.0
    e = np.sqrt(np.clip(1 - (L*L/c)*(1 + 1/c), 0, None))
    k = e*c/(c + 1)
    eta = Q.copy()
    for _ in range(60):
        eta = eta - (eta - k*np.sin(eta) - Q)/(1 - k*np.cos(eta))
    s_ = 2 + c*(1 - e*np.cos(eta))
    r = np.sqrt(s_*(s_ - 2))
    p = np.sqrt(np.clip(2*(E + 1/(1 + np.sqrt(1 + r*r)) - 0.5*L*L/(r*r)), 0, None))*np.where(np.sin(eta) >= 0, 1.0, -1.0)
    return r, p


def particulas(n):
    """Q, J, L, f y (r, p_r) de las partículas de la corrida n. f es el valor de la DF que lee el
    código; la masa de cada partícula es proporcional a f L."""
    c = CORRIDAS[n]
    if c['nl'] == 1:
        Q, J, f = R.reticula(c['nj'], c['nq'], sj=SJ, sq=SQ)
        L = np.full_like(Q, L0)
    else:
        Q, J, L, w = R.reticula_L(c['nj'], c['nq'], c['nl'], c['sl'], L0, c['ancho'], sj=SJ, sq=SQ)
        f = w/L
    r, p = isocrono_rp(Q, J, L)
    return Q, J, L, f, r, p


def capas_L(n):
    """Valores de L de la retícula de la corrida n y su peso (masa relativa de cada capa)."""
    c = CORRIDAS[n]
    if c['nl'] == 1:
        return np.array([L0]), np.array([1.0])
    Ln = L0 + c['ancho']*c['sl']*(2*(np.arange(c['nl']) + 0.5)/c['nl'] - 1.0)
    return Ln, np.exp(-(Ln - L0)**2/c['sl']**2)*Ln


def coherencia_libre(Ln, w, t, J=0.0):
    """|sum_l w_l exp(-i Omega(J, L_l) t)|/sum_l w_l: lo que queda de la suma de las capas en la mezcla
    libre. Para una gaussiana continua es exp(-(dOmega/dL s_L t)^2/4); con N_L capas vuelve a 1 en t_rec."""
    om = (J + R.c_de_L(Ln))**-3
    return np.abs((w[None, :]*np.exp(-1j*om[None, :]*np.asarray(t, float)[:, None])).sum(1))/w.sum()


def reticulas():
    """Para cada corrida con L distribuida: paso de la retícula, t_rec y la mayor diferencia, hasta el
    final de la corrida, entre la coherencia libre de la retícula y la de la gaussiana continua."""
    t = np.linspace(0, TFIN, 4001); h = 1e-5
    dOdL = abs((R.c_de_L(L0 + h))**-3 - (R.c_de_L(L0 - h))**-3)/(2*h)
    print(f'|dOmega/dL| = {dOdL:.4f} en J = 0; fin de las corridas en t/tau_1 = {TFIN/TAU1:.1f}')
    for n, c in CORRIDAS.items():
        if c['nl'] == 1:
            continue
        Ln, w = capas_L(n); dL = Ln[1] - Ln[0]
        fino = L0 + 6*c['sl']*(2*(np.arange(600) + 0.5)/600 - 1.0); wf = np.exp(-(fino - L0)**2/c['sl']**2)*fino
        dif = [np.abs(coherencia_libre(Ln, w, t, J) - coherencia_libre(fino, wf, t, J)) for J in (0.0, SJ)]
        tarde = t >= 20*TAU1
        print(f'{n:11} s_L = {c["sl"]:5.3f}  N_L = {c["nl"]:2d} en ±{c["ancho"]:g} s_L  dL = {dL:.5f}  t_rec = {2*np.pi/dOdL/dL/TAU1:5.1f} tau_1'
              f'  max |S - S_gauss|: {max(d.max() for d in dif):.1e} (hasta 20 tau_1: {max(d[~tarde].max() for d in dif):.1e})')


def forma(n):
    """Filas en J y columnas en Q que quedan tras el recorte, y puntos en L."""
    c = CORRIDAS[n]
    filas, columnas = B.recorte_reticula(c['nj'], c['nq'], 1e-7, SJ, SQ)
    return len(filas), len(columnas), c['nl']


def preparar():
    os.makedirs(ruta('ic'), exist_ok=True); os.makedirs(PARDIR, exist_ok=True)
    for n, c in CORRIDAS.items():
        Q, J, L, f, r, p = particulas(n)
        nrc, npc, nlc = forma(n)
        assert len(r) == nrc*npc*nlc
        # El mapa analítico del código devuelve (Q, J) del dato: comprueba la inversión con L por partícula.
        Qc, Jc, _ = aa_isocrono(r, p, L)
        eq = np.abs(np.angle(np.exp(1j*(Qc - Q)))).max(); ej = np.abs(Jc - J).max()
        m = f*L; nef = m.sum()**2/(m**2).sum()
        ic = ruta('ic', f'{n}.dat')
        np.savetxt(ic, np.column_stack([r, p, L, f]), fmt='%.17e')
        nota = ''
        if c['fijo']:                           # las mismas partículas que la corrida de L fija
            g = os.path.join(FIJO_BULTO, 'ic', f'{c["fijo"]}.dat')
            if os.path.exists(g):
                a = np.loadtxt(g)
                nota = (f'; frente a {c["fijo"]}: max |dr| = {np.abs(a[:, 0] - r).max():.1e}, max |dp| = {np.abs(a[:, 1] - p).max():.1e},'
                        f' max |df/f| = {np.abs(a[:, 2]/f - 1).max():.1e}')
        sl = c['sl'] if c['nl'] > 1 else 0.25   # con una sola L el intervalo solo fija dlc, que se cancela
        v = [('state', 'checkpoint'), ('CheckPointfile', f'rep/{GRUPO}/ic/{n}.dat'),
             ('Nrc', nrc), ('Npc', npc), ('Nlc', nlc), ('a0', c['m']),
             ('l0', L0), ('sl', c['sl'] if c['nl'] > 1 else 0.5), ('sr', SJ), ('sp', SQ),
             ('lminc', L0 - (c['ancho'] if c['nl'] > 1 else 2.0)*sl), ('lmaxc', L0 + (c['ancho'] if c['nl'] > 1 else 2.0)*sl),
             ('rminc', 3.0), ('rmaxc', 14.0), ('pminc', -0.25), ('pmaxc', 0.25),
             ('dr', 0.1), ('rmin', 0.0), ('rmax', 20.0), ('pmax', 2.0), ('courant', 1.0),
             ('Nt', int(round(TFIN/DT))), ('integrator', 'yoshida4'), ('bsplineorder', 1),
             ('BGtype', 'Isochrone'), ('autointeraction', '.true.'),
             ('spatial_output', SALIDA), ('field_output', SALIDA), ('time_output', 100*SALIDA),
             ('output_format', 'raw'), ('directory', f'rep/{GRUPO}/{n}')]
        v = [(k, c['cambios'].get(k, x)) for k, x in v]
        assert not set(c['cambios']) - {k for k, _ in v}
        cab = [f'# Nube {n}: F_0(Q, J) con s_Q = s_J = {SJ} por exp(-(L - {L0:g})^2/s_L^2), s_L = {c["sl"]:g},',
               f'# M_gas = {c["m"]:g}, retícula {c["nj"]} x {c["nq"]} recortada ({nrc} x {npc}) por N_L = {nlc}'
               + (f' en |L - {L0:g}| <= {c["ancho"]:g} s_L:' if nlc > 1 else ':'),
               f'# {len(r)} partículas, N_ef = {nef:.0f}. Generado por reproducir/scripts/nube_L.py.']
        if c['cambios']:
            cab.append('# Cambios: ' + ', '.join(f'{k} = {x}' for k, x in c['cambios'].items()) + '.')
        open(os.path.join(PARDIR, f'{n}.par'), 'w').write('\n'.join(cab + [f'{k:<16} = {x}' for k, x in v]) + '\n')
        print(f'{n:11} {len(r):6d} partículas ({nrc} x {npc} x {nlc}), N_ef = {nef:6.0f}; ida y vuelta del mapa:'
              f' |dQ| <= {eq:.1e}, |dJ| <= {ej:.1e}{nota}')


def _ejecutar(n, extra=(), log=None):
    env = dict(os.environ, OMP_NUM_THREADS='4', OMP_PLACES='cores', OMP_PROC_BIND='close')
    t0 = time.time()
    r = subprocess.run(['./VP_PIC', os.path.join(PARDIR, f'{n}.par'), *extra], cwd=EXE,
                       stdout=open(log, 'w'), stderr=subprocess.STDOUT, env=env)
    return r.returncode, time.time() - t0


def paso_de_tiempo(log):
    for l in open(log):
        if 'Time step fixed at size' in l:
            return float(l.split(':')[1])
    return float('nan')


def tiempo(solo=None, pasos=2000):
    """Corrida corta de cada caso, con las salidas de la corrida completa (una instantánea cada SALIDA
    pasos), para estimar el costo de la serie."""
    total = 0.0
    for n in CORRIDAS:
        if solo and n not in solo:
            continue
        log = ruta(f'tiempo_{n}.log')
        cod, s = _ejecutar(n, [f'Nt={pasos}', f'directory=rep/{GRUPO}/tiempo_{n}', f'time_output={pasos}'], log)
        npart = int(np.prod(forma(n))); est = s*TFIN/DT/pasos
        hecho = os.path.exists(ruta(f'{n}.ok'))
        total += 0.0 if hecho else est
        print(f'{"OK" if cod == 0 else "FALLO":5} {n:11} {npart:6d} partículas: {s:6.1f} s en {pasos} pasos, dt = {paso_de_tiempo(log):g}'
              f' -> {est/60:6.1f} min la corrida completa{" (ya corrida)" if hecho else ""}', flush=True)
        subprocess.run(['rm', '-rf', ruta(f'tiempo_{n}')])
    print(f'total de lo que falta: {total/60:.0f} min ({total/3600:.1f} h), sin el análisis')


def metadatos(n, par):
    """Qué produjo la corrida: sha256 del ejecutable y del .par, y commit del repositorio."""
    sha = lambda f: hashlib.sha256(open(f, 'rb').read()).hexdigest()
    git = lambda *a: subprocess.run(['git', *a], cwd=RAIZ, capture_output=True, text=True).stdout.strip()
    src = git('log', '-1', '--format=%h', '--', 'src/') + ('+cambios' if git('status', '--porcelain', 'src/') else '')
    with open(os.path.join(PARDIR, 'METADATOS.txt'), 'a') as fo:
        fo.write(f'{n:11} {datetime.datetime.now().isoformat(timespec="minutes")}  par {sha(par)[:16]}'
                 f'  VP_PIC {sha(os.path.join(EXE, "VP_PIC"))[:16]}  src {src}  HEAD {git("rev-parse", "--short", "HEAD")}\n')


def correr(solo=None):
    for n in CORRIDAS:
        if (solo and n not in solo) or os.path.exists(ruta(f'{n}.ok')):
            continue
        cod, s = _ejecutar(n, log=ruta(f'{n}.log'))
        dt = paso_de_tiempo(ruta(f'{n}.log'))
        print(f'{"OK" if cod == 0 else "FALLO":5} {n} ({s:.0f} s), dt = {dt:g}', flush=True)
        if cod == 0:
            assert abs(dt - DT) < 1e-12, f'el código eligió dt = {dt}, no {DT}'
            open(ruta(f'{n}.ok'), 'w').write('')
            metadatos(n, os.path.join(PARDIR, f'{n}.par'))


class MapaL(MapaAA):
    """MapaAA con el radio de la órbita circular calculado una vez por cada valor distinto de L. Las
    partículas de esta serie tienen N_L valores de L, y la bisección por partícula era la cuarta parte
    del costo del mapa. El resultado es el mismo bit a bit que el de MapaAA."""

    def r_circular(self, L, it=100):
        L = np.asarray(L, float)
        v, inv = np.unique(L, return_inverse=True)
        return MapaAA.r_circular(self, v, it)[inv].reshape(L.shape)

    def _raiz(self, E, L, a, b, it=60):
        """60 bisecciones en vez de 80: el intervalo inicial mide menos de 80 y tras 60 mitades queda por
        debajo de la separación entre dos números de máquina, así que las 20 restantes no cambian nada
        (comprobado bit a bit sobre 42000 partículas en tres instantáneas)."""
        return MapaAA._raiz(self, E, L, a, b, it)


def _mapear(tarea):
    raw, indices = tarea
    c = CorridaRaw(raw)
    rg = c.r; phib = phi_iso(rg)
    sal = []
    for i in indices:
        s = c[i]
        r, p, L = s['r_part'], s['p_part'], s['l_part']
        ps = s['potential'] - phib
        Q, J, _ = MapaL(rg, ps)(r, p, L)        # potencial de este instante, L de cada partícula
        Qd, Jd, _ = aa_isocrono(r, p, L)        # isócrono desnudo
        sal.append((s['time'], r, p, Q, J, ps, Qd, Jd, s['total_energy']))
    c.cerrar()
    return sal


def analizar(solo=None, procesos=4):
    """Escribe <n>/bulto.npz con las mismas claves que bulto.analizar de vlasov-poisson_PIC, más L:
    t, hk (mapa del potencial del instante), hk_desnudo, phi_self, r, w (masa relativa, f L), L y,
    en precisión simple, R, P, Q, J de cada partícula en cada instantánea."""
    from multiprocessing import Pool
    for n in CORRIDAS:
        if (solo and n not in solo) or os.path.exists(ruta(n, 'bulto.npz')) or not os.path.exists(ruta(f'{n}.ok')):
            continue
        t0 = time.time()
        raw = ruta(n, 'vlasov_output.raw')
        c = CorridaRaw(raw); nt = len(c); rg = c.r; s0 = c[0]; c.cerrar()
        w, L = s0['fl'], s0['l_part']
        grupos = [(raw, list(range(i, min(i + 25, nt)))) for i in range(0, nt, 25)]
        t, R_, P_, Q_, J_, ps, hk, hkd, en = ([] for _ in range(9))
        malos = 0
        with Pool(procesos) as pool:
            for b in pool.imap(_mapear, grupos):
                for ti, r, p, q, j, pot, qd, jd, e in b:
                    t.append(ti); ps.append(pot); en.append(e)
                    malos = max(malos, int((~np.isfinite(q) | ~np.isfinite(j)).sum()))
                    hk.append([np.nansum(w*np.exp(-1j*m*q)) for m in range(6)])
                    hkd.append([np.nansum(w*np.exp(-1j*m*qd)) for m in range(6)])
                    R_.append(r.astype(np.float32)); P_.append(p.astype(np.float32))
                    Q_.append(q.astype(np.float32)); J_.append(j.astype(np.float32))
        t, ps, hk, hkd = np.array(t), np.array(ps), np.array(hk), np.array(hkd)
        hk = hk/hk[:, :1].real
        hkd = hkd/hkd[:, :1].real
        np.savez(ruta(n, 'bulto.npz'), t=t, hk=hk, hk_desnudo=hkd, phi_self=ps, r=rg, w=w, L=L, energia=np.array(en),
                 R=np.array(R_), P=np.array(P_), Q=np.array(Q_), J=np.array(J_))
        print(f'  analizado {n}: {len(t)} instantáneas, {len(w)} partículas ({time.time()-t0:.0f} s);'
              f' sin mapa en alguna instantánea: {malos}; |dE/E| <= {np.abs(np.array(en)/en[0] - 1).max():.1e}', flush=True)


def plano(npz, sal):
    """Medidas en el plano de zeta (reducido.medidas_plano) de un bulto.npz; se guardan en sal."""
    if os.path.exists(sal):
        return np.load(sal)
    d = np.load(npz)
    t, w = d['t'], np.asarray(d['w'], float)
    Q, J = d['Q'], d['J']
    ok = np.isfinite(Q).all(0) & np.isfinite(J).all(0)
    X = np.empty((len(t), int(ok.sum())), np.float32); Y = np.empty_like(X)
    for k in range(len(t)):
        a = np.sqrt(2.0*J[k][ok].astype(float)); q = Q[k][ok].astype(float)
        X[k], Y[k] = a*np.cos(q), a*np.sin(q)
    M = R.medidas_plano(t, X, Y, w[ok], SJ, TAU1, KAPPA)
    np.savez(sal, fuera=1.0 - w[ok].sum()/w.sum(), **M)
    return np.load(sal)


def medidas(n):
    return plano(ruta(n, 'bulto.npz'), ruta(n, 'medidas.npz'))


def por_capas(n):
    """Primer armónico de cada capa de L, h_l(t) = sum_{i en l} m_i exp(-i Q_i)/sum_{i en l} m_i, del PIC.
    Devuelve t, |h_1| total, la media de |h_l| sobre las capas (pesada con su masa) y su cociente, la
    coherencia entre capas: 1 si todas oscilan en fase, ~N_ef^(-1/2) con fases al azar y menos si se
    abren en abanico. Se guarda en <n>/capas.npz."""
    sal = ruta(n, 'capas.npz')
    if not os.path.exists(sal):
        d = np.load(ruta(n, 'bulto.npz'))
        t, w, L = d['t'], np.asarray(d['w'], float), d['L']; Q = d['Q']
        v, inv = np.unique(L, return_inverse=True)
        W = np.bincount(inv, weights=w)
        h = np.empty((len(t), len(v)), complex)
        for k in range(len(t)):
            e = np.exp(-1j*Q[k].astype(float)); e[~np.isfinite(e)] = 0.0
            h[k] = (np.bincount(inv, weights=w*e.real) + 1j*np.bincount(inv, weights=w*e.imag))/W
        np.savez(sal, t=t, L=v, W=W, h=h)
    d = np.load(sal)
    return _coherencia(d['t'], d['h'], d['W'])


def por_capas_reducido(nombre):
    """Lo mismo para una corrida del modelo reducido."""
    d = np.load(os.path.join(R.DATOS, nombre + '.npz'))
    t, m = d['t'], d['m']
    L = d['L'] if 'L' in d.files and len(d['L']) == len(m) else np.full(len(m), L0)   # corridas de L igual
    v, inv = np.unique(L, return_inverse=True)
    W = np.bincount(inv, weights=m)
    h = np.empty((len(t), len(v)), complex)
    for k in range(len(t)):
        z = d['X'][k].astype(float) + 1j*d['Y'][k].astype(float)
        e = np.conj(z)/(np.abs(z) + 1e-300)
        h[k] = (np.bincount(inv, weights=m*e.real) + 1j*np.bincount(inv, weights=m*e.imag))/W
    return _coherencia(t, h, W)


def _coherencia(t, h, W):
    total = np.abs((h*W).sum(1))/W.sum(); capa = (np.abs(h)*W).sum(1)/W.sum()
    return t, total, capa, total/capa


def numeros(M):
    """Los números de la tabla a partir de las medidas en el plano de zeta."""
    t = M['t']; u = t >= 20*TAU1; v = (t >= 10*TAU1) & (t <= 20*TAU1)
    env = envolvente(M['ac'], t)
    z = (t >= 2*TAU1) & (t <= 30*TAU1)
    k = np.argmin(np.where(z, env, np.inf))
    return dict(h1=np.abs(M['h1'][v]).mean(), h1t=np.abs(M['h1'][u]).mean(), rho_c=M['rho_c'][u].mean(), m02=M['m02'][u].mean(),
                mc=float(M['mnuc']), dr=M['dr'][u].mean(), ac=[env[np.argmin(np.abs(t - x*TAU1))] for x in TS],
                acmin=env[k], tmin=t[k]/TAU1, actarde=env[t >= 25*TAU1].mean())


def campo(n):
    """Media en 10 <= t/tau_1 <= 20 de la desviación D del potencial del gas (bulto.campo_particulas)."""
    t, D = B.campo_particulas(n)
    return D[(t >= 10*TAU1) & (t <= 20*TAU1)].mean()


def elipse(nu):
    """De bulto.nucleo (elipse en el plano (r, p_r)): amplitud de la oscilación radial del núcleo."""
    t = nu['t']; k = int(round(45/(t[1] - t[0])))
    media = np.array([nu['rc'][max(0, i - k):i + k + 1].mean() for i in range(len(t))])
    return envolvente(np.abs(nu['rc'] - media), t)


def consistencia(n='C021'):
    """La corrida n de este código frente a la de L fija con las mismas partículas."""
    f = CORRIDAS[n]['fijo']
    a = np.load(ruta(n, 'bulto.npz')); b = np.load(os.path.join(FIJO_BULTO, f, 'bulto.npz'))
    t = a['t']; assert np.allclose(t, b['t'], atol=1e-6)
    print(f'{n} (VlasovPoisson_PIC_sp) frente a {f} (vlasov-poisson_PIC); tau_1 = {TAU1:.2f}')
    print(f'  dato inicial: max |dr| = {np.abs(a["R"][0] - b["R"][0]).max():.1e}, max |dp| = {np.abs(a["P"][0] - b["P"][0]).max():.1e}'
          f' (precisión simple); masa: max |dm/m| = {np.abs((a["w"]/a["w"].sum())/(b["w"]/b["w"].sum()) - 1).max():.1e}')
    # Mismo mapa: el de L por partícula sobre la primera y la última instantánea de la corrida de L fija.
    for k in (0, len(t) - 1):
        m = MapaAA(b['r'], b['phi_self'][k])
        Q, J, _ = m(b['R'][k].astype(float), b['P'][k].astype(float), np.full(b['R'].shape[1], L0))
        print(f'  mapa con L por partícula sobre {f}, t = {t[k]:.0f}: max |dJ| = {np.abs(J - b["J"][k]).max():.1e},'
              f' max |dzeta| = {np.abs(np.sqrt(2*J)*np.exp(1j*Q) - np.sqrt(2*b["J"][k].astype(float))*np.exp(1j*b["Q"][k].astype(float))).max():.1e}')
    ha, hb = np.abs(a['hk'][:, 1]), np.abs(b['hk'][:, 1])
    print('  |h_1|(t), curva a curva: mayor diferencia relativa al máximo de |h_1| en cada tramo de t/tau_1')
    for x, y in ((0, 1), (1, 2), (2, 5), (5, 10), (10, 20), (20, 31)):
        z = (t >= x*TAU1) & (t <= y*TAU1)
        print(f'    [{x:2d}, {y:2d}]: {np.abs(ha - hb)[z].max()/hb[z].max():.2e}   medias {ha[z].mean():.4e} {hb[z].mean():.4e}')
    da = np.sqrt(np.mean((a['R'].astype(float) - b['R'].astype(float))**2, axis=1))
    print('  distancia rms entre las posiciones de las mismas partículas en t/tau_1 = 0.1, 0.5, 1, 2, 5, 10, 30: '
          + ' '.join(f'{da[np.argmin(np.abs(t - x*TAU1))]:.1e}' for x in (0.1, 0.5, 1, 2, 5, 10, 30)))
    Ma, Mb = medidas(n), plano(os.path.join(FIJO_BULTO, f, 'bulto.npz'), ruta(f'{f}_fijo_medidas.npz'))
    Na, Nb = numeros(Ma), numeros(Mb)
    B.BASE = FIJO_BULTO; nub = B.nucleo(f); tb, Db = B.campo_particulas(f); mub = B.escalas(f)[2]
    B.BASE = BASE; nua = B.nucleo(n); ta, Da = B.campo_particulas(n); mua = B.escalas(n)[2]
    v = (t >= 10*TAU1) & (t <= 20*TAU1); u = t >= 20*TAU1
    ea, eb = elipse(nua), elipse(nub)
    filas = [('mu', mua, mub), ('|h_1|, 10-20 tau_1', Na['h1'], Nb['h1']), ('|h_1|, >= 20 tau_1', Na['h1t'], Nb['h1t']),
             ('D, 10-20 tau_1', Da[v].mean(), Db[v].mean()),
             ('m02 (plano zeta)', Na['m02'], Nb['m02']), ('rho_c (plano zeta)', Na['rho_c'], Nb['rho_c']),
             ('m_c, disco en zeta', Na['mc'], Nb['mc']), ('m_c, elipse en (r, p_r)', float(nua['mnuc']), float(nub['mnuc'])),
             ('ancho dr del núcleo, zeta', Na['dr'], Nb['dr']), ('ancho dr, elipse', nua['dr'][u].mean(), nub['dr'][u].mean()),
             ('a_c mínima, zeta', Na['acmin'], Nb['acmin']), ('a_c en >= 25 tau_1, zeta', Na['actarde'], Nb['actarde']),
             ('a_c en >= 25 tau_1, elipse', ea[t >= 25*TAU1].mean(), eb[t >= 25*TAU1].mean())]
    print(f'  {"medida":28} {n:>11} {f:>11}   cociente')
    for q, x, y in filas:
        print(f'  {q:28} {x:11.4e} {y:11.4e}   {x/y:6.3f}')
    print(f'  mínimo de a_c (zeta) en t/tau_1 = {Na["tmin"]:.1f} y {Nb["tmin"]:.1f}')
    print('  a_c(t) en t/tau_1 =', TS)
    for q, x in ((f'{n} zeta', Na['ac']), (f'{f} zeta', Nb['ac']),
                 (f'{n} elipse', [ea[np.argmin(np.abs(t - x*TAU1))] for x in TS]),
                 (f'{f} elipse', [eb[np.argmin(np.abs(t - x*TAU1))] for x in TS])):
        print(f'    {q:14} ' + ' '.join(f'{y:.3f}' for y in x))


def tabla():
    """Medidas de cada corrida frente al modelo reducido: exe/rep/12_nube/tabla.txt."""
    v = lambda t: (t >= 10*TAU1) & (t <= 20*TAU1)
    u = lambda t: t >= 20*TAU1
    lineas = ['Nube con L distribuida: PIC (VlasovPoisson_PIC_sp) y modelo reducido. sigma_L = s_L/sqrt(2).',
              '|h_1|: media en 10-20 tau_1 y en t >= 20 tau_1. capa: media sobre las capas de L del |h_1| de cada una, y coh: cociente',
              'entre el |h_1| total y esa media (1: capas en fase), las dos en 10-20 y en t >= 20 tau_1. rho_c, m02 y dr: medias en',
              't >= 20 tau_1; m_c: media en t >= 20 tau_1. a_c: mínimo de la envolvente en 2-30 tau_1 (y cuándo) y media en t >= 25 tau_1.',
              'D: media en 10-20 tau_1.',
              f'{"corrida":11} {"M_gas":>5} {"s_L":>5} {"N_L":>3} {"":8} {"|h_1|":>8} {"tarde":>8} {"capa":>8} {"tarde":>8} {"coh":>5} {"tarde":>5} '
              f'{"rho_c":>5} {"m02":>5} {"m_c":>5} {"dr":>5} {"a_c min":>7} {"t/tau1":>6} {"a_c":>5} {"D":>8} {"mu":>4}']

    def fila(nombre, N, capas):
        t, tot, capa, coh = capas
        return (f'{nombre} {N["h1"]:8.2e} {N["h1t"]:8.2e} {capa[v(t)].mean():8.2e} {capa[u(t)].mean():8.2e} {coh[v(t)].mean():5.2f} {coh[u(t)].mean():5.2f} '
                f'{N["rho_c"]:5.2f} {N["m02"]:5.3f} {N["mc"]:5.3f} {N["dr"]:5.3f} {N["acmin"]:7.3f} {N["tmin"]:6.1f} {N["actarde"]:5.3f}')

    for n, c in CORRIDAS.items():
        if not os.path.exists(ruta(n, 'bulto.npz')):
            continue
        lineas.append(fila(f'{n:11} {c["m"]:5g} {c["sl"]:5g} {c["nl"]:3d} {"PIC":8}', numeros(medidas(n)), por_capas(n))
                      + f' {campo(n):8.2e} {B.escalas(n)[2]:4.2f}')
        if c['red'] and os.path.exists(os.path.join(R.DATOS, c['red'] + '.npz')):
            lineas.append(fila(f'{"":11} {"":5} {"":5} {"":3} {"reducido":8}', numeros(R.diagnosticos(c['red'])), por_capas_reducido(c['red'])))
    open(ruta('tabla.txt'), 'w').write('\n'.join(lineas) + '\n')
    print('\n'.join(lineas))


def tanda(nombre):
    """Corre y analiza las corridas de una tanda, una detrás de otra, y rehace la tabla."""
    for n in TANDAS[nombre]:
        correr([n]); analizar([n])
    tabla()


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('paso', choices=['preparar', 'tiempo', 'correr', 'analizar', 'tanda', 'consistencia', 'tabla', 'reticulas'])
    ap.add_argument('nombre', nargs='?', choices=list(TANDAS), help='tanda (solo para el paso tanda)')
    ap.add_argument('--solo', nargs='+', help='solo estas corridas')
    a = ap.parse_args()
    if a.solo and set(a.solo) - set(CORRIDAS):
        sys.exit('no existen: ' + ', '.join(sorted(set(a.solo) - set(CORRIDAS))))
    if a.paso == 'tanda':
        tanda(a.nombre or sys.exit('falta el nombre de la tanda: ' + ', '.join(TANDAS)))
    elif a.paso in ('preparar', 'consistencia', 'tabla', 'reticulas'):
        globals()[a.paso]()
    else:
        globals()[a.paso](a.solo)
