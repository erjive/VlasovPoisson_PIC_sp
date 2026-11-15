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

    python3 nube_L.py preparar [--solo ...]     condiciones iniciales y .par
    python3 nube_L.py tiempo [--solo ...]       corrida corta de cada caso y costo estimado
    python3 nube_L.py correr [--solo ...]       lo que falte, una corrida detrás de otra, 4 hilos
    python3 nube_L.py analizar [--solo ...]     mapa (Q, J) de cada instantánea
    python3 nube_L.py medir [--solo ...]        medidas de cada corrida, en archivos pequeños
    python3 nube_L.py tanda <nombre>            corre, analiza y mide una tanda (TANDAS), en orden
    python3 nube_L.py lista [tandas o corridas] nombres de las corridas (para los guiones del cluster)
    python3 nube_L.py prueba                    la cadena completa con una corrida corta, y su costo
    python3 nube_L.py modelo [corridas]         series del modelo reducido (serie:CPUs), o sus corridas
    python3 nube_L.py consistencia              C021 frente a la corrida B021 de vlasov-poisson_PIC
    python3 nube_L.py tabla                     medidas frente al modelo reducido
    python3 nube_L.py reticulas                 coherencia libre de la retícula en L de cada corrida
    python3 nube_L.py figura                    figura del artículo (docs/articulo/figuras/bulto_L.pdf)

Necesita el repositorio vlasov-poisson_PIC junto a este (o en la variable VP_PIC_FIJO): de él salen
la retícula, las medidas y las corridas con las que se compara. Todo queda en exe/rep/12_nube/.
`tabla` solo lee los archivos pequeños que deja `medir` (medidas.npz, capas.npz, campo.npz,
resumen.npz), así que las corridas se pueden hacer en otra máquina y traer solo esos archivos
(reproducir/cluster/). Las corridas usan VP_HILOS hilos (4 por omisión) y el análisis VP_PROCESOS
procesos (4).
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
M_MIX = 0.2/(np.e*np.sqrt(np.pi))               # la mayor m_c de un gas mezclado a acción fija (artículo, Sección 6C)

# nombre: masa del gas, retícula completa (N_J, N_Q), puntos, anchura s_L y semiancho (en s_L) de la
# retícula en L, corrida del modelo reducido y corrida de L fija con las que se compara, y cambios
# en el .par.
def _c(m, nl, sl, red=None, fijo=None, ancho=2.5, nj=200, nq=50, perfil=None, **cambios):
    return dict(m=m, nj=nj, nq=nq, nl=nl, sl=sl, ancho=ancho, red=red, fijo=fijo, perfil=perfil, cambios=cambios)


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
    # Frontera en el plano de mu y sigma_L (2026-11-13): la misma serie a M_gas/M = 0.02 y 0.05. Según la
    # estimación B2 del relevo (cizalla en L frente a la rotación propia del núcleo, con la m_c del PIC), el
    # núcleo se pierde entre s_L = 0.007 y 0.014 a 0.02, y entre 0.064 y 0.13 a 0.05.
    'M020_s000':   _c(0.02, 1,  0.0,    'rL0_0.02'),
    'M020_s0025':  _c(0.02, 6,  0.0025, 'rL0025b_0.02'),
    'M020_s005':   _c(0.02, 8,  0.005,  'rL005b_0.02'),
    'M020_s010':   _c(0.02, 12, 0.01,   'rL010b_0.02'),
    'M020_s020':   _c(0.02, 20, 0.02,   'rL020b_0.02'),
    'M050_s000':   _c(0.05, 1,  0.0,    'rL0_0.05'),
    'M050_s020':   _c(0.05, 20, 0.02,   'rL020b_0.05'),
    'M050_s040':   _c(0.05, 30, 0.04,   'rL040b_0.05'),
    # Con s_L >= 0.08 las capas de L menor recurren antes (giran más rápido): hacen falta N_L = 70 y 120 para que
    # la suma libre no se aparte de la gaussiana más de 7e-4 en toda la corrida (con 60 y 90, 7e-3 y 4e-2).
    'M050_s080':   _c(0.05, 70, 0.08,   'rL080b_0.05'),
    'M050_s120':   _c(0.05, 120, 0.12),      # sin modelo: con 252 000 partículas la suma directa tardaría ~1 día
    # 2026-11-15: entre s_L = 0.04 (el núcleo sigue, m_c = 0.15) y 0.08 (ya no), para ubicar dónde desaparece con
    # M_gas/M = 0.05. N_L = 50 cumple la cota de 7e-4 de `reticulas` (con 45, 1.1e-3).
    'M050_s060':   _c(0.05, 50, 0.06,   'rL060b_0.05'),
    # 2026-11-15: otras formas del dato en J con L igual (vlasov-poisson_PIC/reproducir/scripts/formas.py). La masa da el
    # mu de la nube del artículo con M_ref = 0.02, 0.03 y 0.05 (el número del nombre es 1000 M_ref). p2 es la nube del
    # artículo en la misma retícula de 400 x 100; p4 es J^4 exp(-J^2/s_J^2) y h05 es (J - 0.05)_+^2 exp(-J^2/s_J^2), cuyo
    # soporte no llega a la órbita circular. Las de h05 llegan a 30 tau_1 = 11 800.
    'F_p2_020':    _c(0.02,    1, 0.0, 'R020',       nj=400, nq=100),
    'F_p2_030':    _c(0.03,    1, 0.0, 'R030',       nj=400, nq=100),
    'F_p2_050':    _c(0.05,    1, 0.0, 'R050',       nj=400, nq=100),
    'F_p4_020':    _c(0.01849, 1, 0.0, 'rF_p4_020',  nj=400, nq=100, perfil=('p', 4.0)),
    'F_p4_030':    _c(0.02774, 1, 0.0, 'rF_p4_030',  nj=400, nq=100, perfil=('p', 4.0)),
    'F_p4_050':    _c(0.04623, 1, 0.0, 'rF_p4_050',  nj=400, nq=100, perfil=('p', 4.0)),
    'F_h05_020':   _c(0.01301, 1, 0.0, 'rF_h05_020', nj=400, nq=100, perfil=('hueco', 0.05), Nt=236000),
    'F_h05_030':   _c(0.01952, 1, 0.0, 'rF_h05_030', nj=400, nq=100, perfil=('hueco', 0.05), Nt=236000),
    'F_h05_050':   _c(0.03253, 1, 0.0, 'rF_h05_050', nj=400, nq=100, perfil=('hueco', 0.05), Nt=236000),
    # Segunda tanda (2026-11-15): formas que separan las predicciones de mu y de Gamma = g M/(|Omega_0'| J_med^(3/2))
    # (formas.py; el número del nombre es 1000 Gamma). El anillo exp(-(J - 0.15)^2/0.03^2) tiene la banda estrecha (mu
    # grande) y está lejos de la órbita circular (Gamma chico); la meseta es constante en 0.05 <= J <= 0.25. t_fin = 30
    # tau_1: 21 200 (anillo), 13 600 (meseta) y 11 800 (p0). El anillo va en 800 filas para que N_ef se parezca al de p2.
    'G_anillo_150': _c(0.02117, 1, 0.0, 'rG_anillo_150', nj=800, nq=100, perfil=('anillo', (0.15, 0.03)), Nt=424000),
    'G_anillo_300': _c(0.04234, 1, 0.0, 'rG_anillo_300', nj=800, nq=100, perfil=('anillo', (0.15, 0.03)), Nt=424000),
    'G_anillo_450': _c(0.06351, 1, 0.0, 'rG_anillo_450', nj=800, nq=100, perfil=('anillo', (0.15, 0.03)), Nt=424000),
    'G_meseta_150': _c(0.02109, 1, 0.0, 'rG_meseta_150', nj=400, nq=100, perfil=('meseta', (0.05, 0.25)), Nt=272000),
    'G_meseta_450': _c(0.06327, 1, 0.0, 'rG_meseta_450', nj=400, nq=100, perfil=('meseta', (0.05, 0.25)), Nt=272000),
    'G_p0_150':     _c(0.00376, 1, 0.0, 'rG_p0_150',     nj=400, nq=100, perfil=('p', 0.0), Nt=236000),
}
# Tandas, en el orden en que conviene correrlas. Costo medido con `tiempo` (corridas) y estimado (análisis).
TANDAS = {
    # lo esencial: L igual, s_L = 0.02 y 0.01 con la retícula buena, y las dos mismas con las partículas del modelo
    'primera':      ['M030_s000', 'M030_s020', 'M030_s010', 'M030_s020v', 'M030_s010v'],
    'resto':        ['M030_s040', 'M030_s005', 'M010_s000', 'M010_s020', 'M030_s040v', 'M030_s005v', 'M010_s020v'],
    'convergencia': ['M030_s020n', 'M030_s020r'],
    # La frontera a otras masas, las más largas primero (la cola no adelanta trabajos). La corrida de 252 000
    # partículas va aparte porque su análisis necesita más memoria (todo.sh la envía con 32 GB).
    'frontera':        ['M050_s080', 'M050_s040', 'M020_s020', 'M050_s020', 'M020_s010', 'M020_s005', 'M020_s0025',
                        'M020_s000', 'M050_s000'],
    'frontera_grande': ['M050_s120'],
    'frontera2':       ['M050_s060'],
    'formas':          ['F_h05_050', 'F_h05_030', 'F_h05_020', 'F_p4_050', 'F_p4_030', 'F_p4_020',
                        'F_p2_050', 'F_p2_030', 'F_p2_020'],
    'formas2':         ['G_anillo_450', 'G_anillo_300', 'G_anillo_150', 'G_meseta_450', 'G_meseta_150', 'G_p0_150'],
}
# Series del modelo reducido (reducido.py) que acompañan a las corridas, con las CPUs de su trabajo de Slurm, en el
# orden en que todo.sh las envía: la más larga primero.
SERIES_MODELO = (('Lb', 16), ('control', 16), ('Lf5g', 32), ('Lf5', 16), ('Lf2', 16), ('Lf5b', 32))


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
        Q, J, f = R.reticula(c['nj'], c['nq'], sj=SJ, sq=SQ, perfil=c['perfil'])
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
    if c['perfil'] is not None:
        Q, J, _ = R.reticula(c['nj'], c['nq'], sj=SJ, sq=SQ, perfil=c['perfil'])
        return len(np.unique(J)), len(np.unique(Q)), 1
    filas, columnas = B.recorte_reticula(c['nj'], c['nq'], 1e-7, SJ, SQ)
    return len(filas), len(columnas), c['nl']


def preparar(solo=None):
    os.makedirs(ruta('ic'), exist_ok=True); os.makedirs(PARDIR, exist_ok=True)
    for n, c in CORRIDAS.items():
        if solo and n not in solo:
            continue
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
        cab = [f'# Nube {n}: F_0(Q, J) con s_Q = s_J = {SJ} por exp(-(L - {L0:g})^2/s_L^2), s_L = {c["sl"]:g},'
               + (f' factor en J {c["perfil"][0]} ' + ' '.join(f'{v:g}' for v in np.atleast_1d(c['perfil'][1]))
                  + ' (reducido.perfil_J),' if c['perfil'] else ''),
               f'# M_gas = {c["m"]:g}, retícula {c["nj"]} x {c["nq"]} recortada ({nrc} x {npc}) por N_L = {nlc}'
               + (f' en |L - {L0:g}| <= {c["ancho"]:g} s_L:' if nlc > 1 else ':'),
               f'# {len(r)} partículas, N_ef = {nef:.0f}. Generado por reproducir/scripts/nube_L.py.']
        if c['cambios']:
            cab.append('# Cambios: ' + ', '.join(f'{k} = {x}' for k, x in c['cambios'].items()) + '.')
        open(os.path.join(PARDIR, f'{n}.par'), 'w').write('\n'.join(cab + [f'{k:<16} = {x}' for k, x in v]) + '\n')
        print(f'{n:11} {len(r):6d} partículas ({nrc} x {npc} x {nlc}), N_ef = {nef:6.0f}; ida y vuelta del mapa:'
              f' |dQ| <= {eq:.1e}, |dJ| <= {ej:.1e}{nota}')


def _ejecutar(n, extra=(), log=None):
    env = dict(os.environ, OMP_NUM_THREADS=os.environ.get('VP_HILOS', '4'), OMP_PLACES='cores', OMP_PROC_BIND='close')
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


def metadatos(n, par, segundos):
    """Qué produjo la corrida: sha256 del ejecutable y del .par, commit del repositorio, máquina y
    tiempo. Se guarda en <n>.meta junto a la corrida. En el portátil se añade además una línea a
    reproducir/corridas/12_nube/METADATOS.txt, que está versionado; dentro de un trabajo de Slurm
    no, para que el repositorio del cluster no quede con cambios que estorben a `git pull`."""
    import socket
    sha = lambda f: hashlib.sha256(open(f, 'rb').read()).hexdigest()
    git = lambda *a: subprocess.run(['git', *a], cwd=RAIZ, capture_output=True, text=True).stdout.strip()
    src = git('log', '-1', '--format=%h', '--', 'src/') + ('+cambios' if git('status', '--porcelain', 'src/') else '')
    trabajo = os.environ.get('SLURM_JOB_ID')
    open(ruta(f'{n}.meta'), 'w').write(
        f'corrida    {n}\nfecha      {datetime.datetime.now().isoformat(timespec="seconds")}\n'
        f'VP_PIC     {sha(os.path.join(EXE, "VP_PIC"))}\npar        {sha(par)}  {os.path.relpath(par, RAIZ)}\n'
        f'HEAD       {git("rev-parse", "--short", "HEAD")}\nsrc/       {src}\n'
        f'maquina    {socket.gethostname()}, {os.environ.get("VP_HILOS", "4")} hilos, '
        f'trabajo {trabajo or "sin Slurm"}, {segundos:.0f} s\n')
    if not trabajo:
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
            metadatos(n, os.path.join(PARDIR, f'{n}.par'), s)


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


def analizar(solo=None, procesos=None):
    """Escribe <n>/bulto.npz con las mismas claves que bulto.analizar de vlasov-poisson_PIC, más L:
    t, hk (mapa del potencial del instante), hk_desnudo, phi_self, r, w (masa relativa, f L), L y,
    en precisión simple, R, P, Q, J de cada partícula en cada instantánea."""
    from multiprocessing import Pool
    procesos = procesos or int(os.environ.get('VP_PROCESOS', '4'))
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


def tau_de(n):
    """tau_1 de la corrida n: el de su forma del dato en J (reducido.banda_perfil) o el de la nube del artículo."""
    pf = CORRIDAS[n]['perfil']
    if pf is None or n.startswith('prueba_'):   # la prueba de la cadena acorta TAU1 para que sus ventanas caigan en ella
        return TAU1
    return R.banda_perfil(pf, SJ)[3]


def plano(npz, sal, tau=TAU1):
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
    M = R.medidas_plano(t, X, Y, w[ok], SJ, tau, KAPPA)
    np.savez(sal, fuera=1.0 - w[ok].sum()/w.sum(), **M)
    return np.load(sal)


def medidas(n):
    return plano(ruta(n, 'bulto.npz'), ruta(n, 'medidas.npz'), tau_de(n))


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
    """Lo mismo para una corrida del modelo reducido (reducido.capas); None si no está."""
    d = R.capas(nombre)
    return None if d is None else _coherencia(d['t'], d['h'], d['W'])


def _coherencia(t, h, W):
    total = np.abs((h*W).sum(1))/W.sum(); capa = (np.abs(h)*W).sum(1)/W.sum()
    return t, total, capa, total/capa


def resumen(n):
    """mu (de la oscilación del potencial en la primera vuelta radial, bulto.escalas) y el mayor error
    relativo de la energía de la corrida n. Se guarda en <n>/resumen.npz."""
    sal = ruta(n, 'resumen.npz')
    if not os.path.exists(sal):
        if CORRIDAS[n]['perfil'] is None:
            A, wb, mu = B.escalas(n)
        else:                                   # como bulto.escalas, con la banda de la forma del dato en J
            omin, omax, djs, _ = R.banda_perfil(CORRIDAS[n]['perfil'], SJ)
            t, rg, ps, media, fl, esc = B.potencial(n)
            v = t <= 2*np.pi/omin; z = (rg >= B.RVENT[0]) & (rg <= B.RVENT[1])
            A = 0.5*(ps[v][:, z].max(0) - ps[v][:, z].min(0)).max()
            wb = np.sqrt((omax - omin)/djs*A); mu = wb/(omax - omin)
        e = np.load(ruta(n, 'bulto.npz'))['energia']
        np.savez(sal, A=A, wb=wb, mu=mu, dE=np.abs(e/e[0] - 1).max())
    return np.load(sal)


def medir(solo=None):
    """Calcula y guarda las medidas de cada corrida analizada: medidas.npz (plano de zeta), capas.npz
    (por capas de L), campo.npz (desviación del potencial), resumen.npz (mu y energía) y, con una
    sola L, nucleo.npz (la elipse de bulto.py). Son los únicos archivos que lee `tabla`."""
    for n, c in CORRIDAS.items():
        if (solo and n not in solo) or not os.path.exists(ruta(n, 'bulto.npz')) or os.path.exists(ruta(n, 'resumen.npz')):
            continue
        t0 = time.time()
        medidas(n); por_capas(n); B.campo_particulas(n)
        if c['nl'] == 1:
            B.nucleo(n)
        r = resumen(n)
        print(f'  medido {n}: mu = {float(r["mu"]):.3f}, |dE/E| <= {float(r["dE"]):.1e} ({time.time()-t0:.0f} s)', flush=True)


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
        if capas is None:
            pc = f'{"-":>8} {"-":>8} {"-":>5} {"-":>5}'
        else:
            t, tot, capa, coh = capas
            pc = f'{capa[v(t)].mean():8.2e} {capa[u(t)].mean():8.2e} {coh[v(t)].mean():5.2f} {coh[u(t)].mean():5.2f}'
        return (f'{nombre} {N["h1"]:8.2e} {N["h1t"]:8.2e} {pc} '
                f'{N["rho_c"]:5.2f} {N["m02"]:5.3f} {N["mc"]:5.3f} {N["dr"]:5.3f} {N["acmin"]:7.3f} {N["tmin"]:6.1f} {N["actarde"]:5.3f}')

    for n, c in CORRIDAS.items():
        if not os.path.exists(ruta(n, 'resumen.npz')):
            continue
        lineas.append(fila(f'{n:11} {c["m"]:5g} {c["sl"]:5g} {c["nl"]:3d} {"PIC":8}', numeros(medidas(n)), por_capas(n))
                      + f' {campo(n):8.2e} {float(resumen(n)["mu"]):4.2f}')
        if c['red'] and os.path.exists(os.path.join(R.DATOS, c['red'] + '_diag.npz')):
            lineas.append(fila(f'{"":11} {"":5} {"":5} {"":3} {"reducido":8}', numeros(R.diagnosticos(c['red'])), por_capas_reducido(c['red'])))
    open(ruta('tabla.txt'), 'w').write('\n'.join(lineas) + '\n')
    print('\n'.join(lineas))


def tanda(nombre):
    """Corre, analiza y mide las corridas de una tanda, una detrás de otra, y rehace la tabla."""
    for n in TANDAS[nombre]:
        correr([n]); analizar([n]); medir([n])
    tabla()


def prueba(n='M030_s010v', pasos=2000):
    """Prueba corta de la cadena completa en esta máquina: la corrida n recortada a `pasos` pasos (una
    instantánea cada SALIDA pasos), analizada y medida como las de producción, en
    exe/rep/12_nube/prueba_<n>/, que se borra al final. Comprueba que el código elige dt = DT, que
    la energía se conserva, que todas las partículas tienen mapa de ángulo-acción y que las medidas
    son finitas, y estima lo que tardarían la corrida y su análisis completos. Devuelve el número
    de fallas."""
    global TAU1
    if not os.path.exists(ruta('ic', f'{n}.dat')):
        preparar()
    pn = f'prueba_{n}'
    basura = [ruta(pn), ruta(f'{pn}.ok'), ruta(f'{pn}.log')]
    subprocess.run(['rm', '-rf', *basura])
    CORRIDAS[pn] = dict(CORRIDAS[n], red=None, fijo=None)
    tau, fallas = TAU1, []
    try:
        cod, s = _ejecutar(n, [f'Nt={pasos}', f'directory=rep/{GRUPO}/{pn}', f'time_output={pasos}'], ruta(f'{pn}.log'))
        dt = paso_de_tiempo(ruta(f'{pn}.log'))
        npart = int(np.prod(forma(n)))
        print(f'  corrida: {npart} partículas con {CORRIDAS[n]["nl"]} valores de L, {pasos} pasos en {s:.1f} s con '
              f'{os.environ.get("VP_HILOS", "4")} hilos, dt = {dt:g}; la corrida completa tardaría {s*TFIN/DT/pasos/60:.0f} min', flush=True)
        if cod != 0:
            fallas.append(f'VP_PIC terminó con el código {cod}')
            print(''.join(open(ruta(f'{pn}.log')).readlines()[-8:]))
        elif abs(dt - DT) > 1e-12:
            fallas.append(f'el código eligió dt = {dt}, no {DT}')
        if not fallas:
            open(ruta(f'{pn}.ok'), 'w').write('')
            t0 = time.time(); analizar([pn]); ta = time.time() - t0
            nsal = pasos//SALIDA + 1
            print(f'  análisis: {ta:.1f} s para {nsal} instantáneas con {os.environ.get("VP_PROCESOS", "4")} procesos; '
                  f'el de la corrida completa tardaría {ta*(TFIN/DT/SALIDA + 1)/nsal/60:.0f} min', flush=True)
            TAU1 = pasos*DT/31                  # que las ventanas de las medidas, hasta 30 tau_1, caigan en la prueba
            medir([pn])
            d = np.load(ruta(pn, 'bulto.npz'))
            if not (np.isfinite(d['Q']).all() and np.isfinite(d['J']).all()):
                fallas.append('hay partículas sin mapa de ángulo-acción')
            r = resumen(pn); dE, mu = float(r['dE']), float(r['mu'])
            if not dE < 1e-4:
                fallas.append(f'la energía cambia {dE:.1e}')
            N = numeros(medidas(pn)); t, tot, capa, coh = por_capas(pn)
            vals = [N[k] for k in ('h1', 'h1t', 'rho_c', 'm02', 'mc', 'dr', 'acmin', 'actarde')] + [campo(pn), mu]
            if not (np.all(np.isfinite(vals)) and np.isfinite(coh).all()):
                fallas.append('hay medidas que no son finitas')
            print(f'  medidas: |h_1| al final {tot[-1]:.3f}, coherencia entre capas {coh[-1]:.4f}, m02 {N["m02"]:.3f}, '
                  f'mu {mu:.3f}, |dE/E| {dE:.1e}', flush=True)
    finally:
        TAU1 = tau
        del CORRIDAS[pn]
        subprocess.run(['rm', '-rf', *basura])
    print('PASA la prueba de la cadena' if not fallas else 'FALLA la prueba de la cadena: ' + '; '.join(fallas))
    return len(fallas)


# Espacio de fases de la Sección 8: la serie con M_gas/M = 0.03 en t/tau_1 = 3, 10 y al final de la corrida.
FASE_CORRIDAS = ('M030_s000', 'M030_s010', 'M030_s020', 'M030_s040')
FASE_TIEMPOS = (3, 10)


def fase(solo=None):
    """Extrae de <n>/bulto.npz, en t/tau_1 = FASE_TIEMPOS y en la última instantánea, el plano (r, p_r) y las
    variables de ángulo-acción (Q, J) del potencial del instante, calculadas con la L de cada partícula, junto con
    su masa relativa w y su L, a <n>/fase.npz (unos MB). Corre donde están los bulto.npz (el cluster); fase.npz se
    trae aparte con rsync. Rehace los fase.npz que no tienen Q y J."""
    for n in (solo or FASE_CORRIDAS):
        sal = ruta(n, 'fase.npz')
        if os.path.exists(sal) and 'Q' in np.load(sal).files:
            print(f'  {n}: ya está')
            continue
        if not os.path.exists(ruta(n, 'bulto.npz')):
            print(f'  {n}: no hay bulto.npz')
            continue
        d = np.load(ruta(n, 'bulto.npz'))
        t = d['t']; idx = [int(np.argmin(np.abs(t - x*TAU1))) for x in FASE_TIEMPOS] + [len(t) - 1]
        np.savez(sal, t=t[idx], R=d['R'][idx], P=d['P'][idx], Q=d['Q'][idx], J=d['J'][idx], w=d['w'], L=d['L'])
        print(f'  {n}: fase.npz con t/tau_1 = {", ".join(f"{x:.1f}" for x in t[idx]/TAU1)}, {len(d["w"])} partículas', flush=True)


def figura_fase():
    """Figura de la Sección 8: densidad de masa del gas en el plano de la amplitud compleja de la oscilación radial,
    zeta = sqrt(2J) exp(i(Q - kappa_0 t)), con (Q, J) del potencial del instante y la L de cada partícula, para
    M_gas/M = 0.03 y cuatro anchuras en L (filas) en t/tau_1 = 3, 10 y al final (columnas). En este plano la órbita
    circular de cada capa de L es el origen: la dispersión R'(L) sigma_L de los radios de las órbitas circulares, que
    en el plano (r, p_r) emborrona el núcleo (+-0.42 en r con s_L = 0.04, más que el núcleo), no aparece, y queda el
    desfase de las capas. Es el plano en que se mide m_c. La densidad es la masa en cada celda de una malla de 0.016,
    suavizada con un núcleo gaussiano de una celda y dividida por la masa total y el área de la celda, en escala
    logarítmica común. Escribe docs/articulo/figuras/bulto_L_fase.pdf en vlasov-poisson_PIC."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.ndimage import gaussian_filter
    plt.rcParams.update({'font.size': 9})
    v = (-0.8, 0.8, -0.8, 0.8, 100, 100)

    def densidad(x, y, w):
        H, xe, ye = np.histogram2d(x, y, bins=(v[4], v[5]), range=((v[0], v[1]), (v[2], v[3])), weights=w)
        return gaussian_filter(H, 1.0)/(w.sum()*(xe[1] - xe[0])*(ye[1] - ye[0]))

    fig = plt.figure(figsize=(6.2, 7.6), layout='constrained')
    gs = fig.add_gridspec(len(FASE_CORRIDAS), 3)
    tope = None
    for i, n in enumerate(FASE_CORRIDAS):
        d = np.load(ruta(n, 'fase.npz')); t, Q, J, w = d['t'], d['Q'].astype(float), d['J'].astype(float), d['w'].astype(float)
        sl = CORRIDAS[n]['sl']
        for k in range(len(t)):
            ok = np.isfinite(Q[k]) & np.isfinite(J[k])
            a_ = np.sqrt(2.0*np.maximum(J[k][ok], 0.0)); th = Q[k][ok] - KAPPA*t[k]
            D = densidad(a_*np.cos(th), a_*np.sin(th), w[ok])
            if tope is None:
                tope = np.log10(D.max())             # la mayor densidad de la primera fila en t = 3 tau_1
            a = fig.add_subplot(gs[i, k])
            im = a.imshow(np.log10(np.maximum(D, 1e-30)).T, origin='lower', extent=v[:4], aspect='equal', cmap='magma_r',
                          vmin=tope - 4, vmax=tope, interpolation='nearest', rasterized=True)
            a.set_xticks([-0.5, 0, 0.5]); a.set_yticks([-0.5, 0, 0.5])
            if i == 0:
                a.set_title(f'$t/\\tau_1={t[k]/TAU1:.1f}$', fontsize=8)
            if i == len(FASE_CORRIDAS) - 1:
                a.set_xlabel(r'$\sqrt{2J}\cos\vartheta$')
            else:
                a.set_xticklabels([])
            if k:
                a.set_yticklabels([])
            else:
                a.set_ylabel(r'$\sqrt{2J}\sin\vartheta$')
                a.text(0.96, 0.95, f'$\\sigma_L/L_0={100*sl/np.sqrt(2)/L0:.2f}\\%$', transform=a.transAxes, fontsize=7,
                       va='top', ha='right')
    cb = fig.colorbar(im, ax=fig.axes, shrink=0.5, pad=0.02)
    cb.set_label(r'$\log_{10}$ of the mass density in the $\zeta$-plane', fontsize=8); cb.ax.tick_params(labelsize=7)
    # El anillo de baja densidad que rodea al núcleo: mínimo de la densidad promediada en el ángulo (anillos de 0.02
    # entre |zeta| = 0.1 y 0.4) frente al radio (2 J_res)^{1/2} de la resonancia con la oscilación del núcleo,
    # J_res = (g m_c M_gas/(sqrt(2) |Omega_0'|))^{2/3}, con la m_c de la corrida.
    h = 1e-5
    dOdJ = abs((h + R.c_de_L(L0))**-3 - (-h + R.c_de_L(L0))**-3)/(2*h)
    bordes = np.linspace(0, 0.5, 26); cen = 0.5*(bordes[1:] + bordes[:-1]); area = np.pi*np.diff(bordes**2)
    for n in FASE_CORRIDAS:
        d = np.load(ruta(n, 'fase.npz')); mc = numeros(medidas(n))['mc']
        jres = (R.G*mc*CORRIDAS[n]['m']/(np.sqrt(2)*dOdJ))**(2/3)
        txt = []
        for k in range(1, len(d['t'])):
            J_ = d['J'][k].astype(float); ok = np.isfinite(J_)
            dens = np.histogram(np.sqrt(2*np.maximum(J_[ok], 0)), bordes, weights=d['w'][ok])[0]/area/d['w'].sum()
            z = (cen > 0.1) & (cen < 0.4)
            txt.append(f't/tau_1 = {d["t"][k]/TAU1:.1f}: mínimo en |zeta| = {cen[np.argmin(np.where(z, dens, np.inf))]:.2f}')
        print(f'  {n}: m_c = {mc:.3f}, (2 J_res)^(1/2) = {np.sqrt(2*jres):.3f}; ' + '; '.join(txt))
    pdf = os.path.join(FIJO, 'docs', 'articulo', 'figuras', 'bulto_L_fase.pdf')
    fig.savefig(pdf, dpi=300); plt.close(fig)
    print('escrita', pdf)


def series_L():
    """Por masa del gas, las corridas medidas de las retículas sin recurrencia (sin sufijo v, n ni r), ordenadas por
    s_L, con los números del PIC (P), mu, D y, si existen, los del modelo con las mismas partículas (M)."""
    out = {}
    for n, c in CORRIDAS.items():
        if not (n.startswith('M') and n[-1].isdigit() and os.path.exists(ruta(n, 'resumen.npz'))):
            continue
        f = dict(n=n, sl=c['sl'], nl=c['nl'], x=100*c['sl']/(np.sqrt(2)*L0), P=numeros(medidas(n)), mu=float(resumen(n)['mu']),
                 D=campo(n))
        if c['red'] and os.path.exists(os.path.join(R.DATOS, c['red'] + '_diag.npz')):
            f['M'] = numeros(R.diagnosticos(c['red']))
        out.setdefault(c['m'], []).append(f)
    return {m: sorted(F, key=lambda f: f['sl']) for m, F in out.items()}


def frontera():
    """Dónde se pierde el núcleo para cada masa: el s_L en que el exceso m_c - m_mix cae a la mitad de su valor con
    L igual (interpolación lineal en s_L), en el PIC y en el modelo, y la estimación B2 del relevo: el cociente entre
    la dispersión |dOmega/dL| sigma_L de las frecuencias de las capas de L y la rotación g m_c M_gas/rho que el
    núcleo induce en su borde, rho = (0.1 s_J)^{1/2}, con la m_c del PIC con L igual."""
    h = 1e-5
    dOdL = abs((R.c_de_L(L0 + h))**-3 - (R.c_de_L(L0 - h))**-3)/(2*h)
    rho = np.sqrt(0.1*SJ)
    print(f'|dOmega/dL| = {dOdL:.4f} en J = 0, g = {R.G:.5f}, rho = {rho:.3f}, m_mix = {M_MIX:.4f}')
    for m, F in sorted(series_L().items()):
        if len(F) < 2 or F[0]['sl'] != 0:
            continue
        print(f'\nM_gas/M = {m:g} (mu = {F[0]["mu"]:.2f}): s_L, sigma_L/L_0 (%), N_L | m_c, |h_1|, masa con J < 0.2 s_J, D: PIC (modelo)')
        for f in F:
            P, Mm = f['P'], f.get('M')
            par = lambda k: f'{P[k]:.3g}' + (f' ({Mm[k]:.3g})' if Mm else '')
            print(f'   {f["sl"]:6.4f} {f["x"]:5.2f} {f["nl"]:3d} | m_c {par("mc")}  |h_1| {par("h1")}  m02 {par("m02")}  D {f["D"]:.2e}')
        mitad = {}
        if F[0]['P']['mc'] < 1.5*M_MIX:
            print('   con L igual no hay núcleo (m_c compatible con m_mix): no hay frontera que medir')
            continue
        for clave, nombre in (('P', 'PIC'), ('M', 'modelo')):
            G = [f for f in F if clave in f]
            s = np.array([f['sl'] for f in G]); e = np.array([f[clave]['mc'] for f in G]) - M_MIX
            k = np.nonzero(e < 0.5*e[0])[0] if len(G) > 1 else []
            if len(k) == 0:
                print(f'   {nombre}: el exceso no baja a la mitad en las corridas medidas')
                continue
            k = k[0]; mitad[clave] = s[k - 1] + (e[k - 1] - 0.5*e[0])/(e[k - 1] - e[k])*(s[k] - s[k - 1])
            print(f'   {nombre}: la mitad del exceso en s_L = {mitad[clave]:.4f} (sigma_L/L_0 = {100*mitad[clave]/np.sqrt(2)/L0:.2f} %)')
        w = R.G*F[0]['P']['mc']*m/rho
        cociente = lambda sl: dOdL*sl/np.sqrt(2)/w
        print(f'   B2: g m_c M_gas/rho = {w:.2e}; cociente 0.5 en s_L = {0.5*w*np.sqrt(2)/dOdL:.4f} y 1 en s_L = {w*np.sqrt(2)/dOdL:.4f}'
              + (f'; en la mitad del PIC vale {cociente(mitad["P"]):.2f}' if 'P' in mitad else ''))


def tabla_L():
    """Filas de la tabla de la Sección 8 del artículo (docs/articulo/tab_L.tex en vlasov-poisson_PIC, macro \\TabL):
    PIC y modelo para cada corrida de series_L()."""
    sci = B._sci
    filas = []
    for m, F in sorted(series_L().items()):
        for i, f in enumerate(F):
            P, Mm = f['P'], f.get('M')
            mod = lambda k, fmt: (fmt(Mm[k]) if Mm else '--')
            filas.append((f'${m:g}$ & ${f["mu"]:.2f}$' if i == 0 else ' & ') + f' & ${f["x"]:.2f}$ & ${f["nl"]}$ & '
                         f'{sci(P["h1"])} & {mod("h1", sci)} & ${P["m02"]:.3f}$ & {mod("m02", lambda x: f"${x:.3f}$")} & '
                         f'${P["mc"]:.3f}$ & {mod("mc", lambda x: f"${x:.3f}$")} \\\\')
        filas.append('\\hline')
    filas = filas[:-1]
    sal = os.path.join(FIJO, 'docs', 'articulo', 'tab_L.tex')
    open(sal, 'w').write('% Escrito por VlasovPoisson_PIC_sp/reproducir/scripts/nube_L.py articulo; no editar a mano.\n'
                         '\\newcommand{\\TabL}{%\n' + '%\n'.join(filas) + '%\n}\n')
    print('escrita', sal)


def figura():
    """Figura de la Sección 8 del artículo: |h_1| (media en 10-20 tau_1), masa del núcleo m_c y masa con J < 0.2 s_J
    (medias en t >= 20 tau_1) frente a sigma_L/L_0, del PIC (líneas continuas) y del modelo con las mismas partículas
    (discontinuas), con un tono de azul por masa del gas (de claro a oscuro al crecer la masa). Escribe
    docs/articulo/figuras/bulto_L.pdf en vlasov-poisson_PIC y una copia bulto_L.png en exe/rep/12_nube/."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    sup, tinta, tinta2, reja = '#ffffff', '#0b0b0b', '#52514e', '#e4e3df'
    plt.rcParams.update({'font.size': 9, 'axes.edgecolor': tinta2, 'axes.labelcolor': tinta2, 'xtick.color': tinta2,
                         'ytick.color': tinta2, 'text.color': tinta, 'axes.linewidth': 0.6, 'figure.facecolor': sup,
                         'axes.facecolor': sup, 'savefig.facecolor': sup})
    datos = series_L()
    masas = sorted(datos)
    tonos = dict(zip(masas, plt.cm.Blues(np.linspace(0.45, 0.95, len(masas)))))
    fig, ax = plt.subplots(1, 3, figsize=(10, 3.4))
    paneles = (('h1', r'$|h_1|/h_0$, mean over $10\leq t/\tau_1\leq 20$', True, None),
               ('mc', r'mass fraction $m_c$ of the core', False, (M_MIX, r'$m_{mix}$')),
               ('m02', r'mass fraction with $J < 0.2\,\sigma_J$', False, None))
    for a, (k, titulo, log, cota) in zip(ax, paneles):
        for m in masas:
            F = datos[m]
            a.plot([f['x'] for f in F], [f['P'][k] for f in F], color=tonos[m], lw=1.3, marker='o', ms=5, mec=sup, mew=0.7)
            G = [f for f in F if 'M' in f]
            if G:
                a.plot([f['x'] for f in G], [f['M'][k] for f in G], color=tonos[m], lw=1.1, ls=(0, (4, 2.5)), marker='D', ms=4,
                       mec=sup, mew=0.7)
        if cota:
            a.axhline(cota[0], color=tinta2, lw=0.8, ls=':', zorder=0)
            a.annotate(cota[1], (0.45, cota[0]), xycoords=('axes fraction', 'data'), xytext=(0, 3), textcoords='offset points',
                       ha='left', fontsize=8, color=tinta2)
        if log:
            a.set_yscale('log')
        else:
            a.set_ylim(0, None)
        a.set_xlim(-0.1, None)
        a.grid(True, axis='y', color=reja, lw=0.6); a.set_axisbelow(True)
        for lado in ('top', 'right'):
            a.spines[lado].set_visible(False)
        a.tick_params(length=3, width=0.6); a.tick_params(which='minor', length=2, width=0.4)
        a.set_title(titulo, loc='left', fontsize=9.5, color=tinta, pad=6)
        a.set_xlabel(r'$\sigma_L/L_0$ (percent)')
    h = [Line2D([], [], color=tonos[m], lw=1.6, label=rf'$M_{{gas}}/M={m:g}$ ($\mu={datos[m][0]["mu"]:.2f}$)') for m in masas]
    h += [Line2D([], [], color=tinta2, lw=1.3, marker='o', ms=5, mec=sup, label='PIC'),
          Line2D([], [], color=tinta2, lw=1.1, ls=(0, (4, 2.5)), marker='D', ms=4, mec=sup, label='epicyclic model')]
    fig.tight_layout(rect=(0, 0, 1, 0.84), w_pad=1.8)
    fig.legend(handles=h, loc='upper center', ncol=(len(h) + 1)//2, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, 1.0),
               handlelength=2.6, columnspacing=1.6)
    pdf = os.path.join(FIJO, 'docs', 'articulo', 'figuras', 'bulto_L.pdf')
    fig.savefig(pdf); fig.savefig(ruta('bulto_L.png'), dpi=110); plt.close(fig)
    print('escrita', pdf, 'con las masas', ', '.join(f'{m:g} ({len(datos[m])} corridas)' for m in masas))


def articulo():
    """Tabla y figura de la Sección 8 del artículo."""
    tabla_L(); figura()


def modelo(nombres):
    """Series del modelo reducido que acompañan a las corridas, como serie:CPUs; con "corridas", los nombres de
    todas sus corridas (para traer.sh)."""
    if nombres == ['corridas']:
        print(' '.join(n for s, _ in SERIES_MODELO for n, *_ in R.SERIES[s]))
    else:
        print(' '.join(f'{s}:{c}' for s, c in SERIES_MODELO))


def lista(nombres):
    """Nombres de las corridas de las tandas o corridas dadas (sin argumentos: todas las tandas)."""
    sal = []
    for x in nombres or list(TANDAS):
        if x not in TANDAS and x not in CORRIDAS:
            sys.exit(f'no existe la tanda ni la corrida {x}')
        sal += TANDAS.get(x, [x])
    print(' '.join(dict.fromkeys(sal)))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('paso', choices=['preparar', 'tiempo', 'correr', 'analizar', 'medir', 'tanda', 'lista', 'prueba',
                                     'modelo', 'consistencia', 'tabla', 'reticulas', 'figura', 'frontera', 'articulo',
                                     'fase', 'figura_fase'])
    ap.add_argument('nombres', nargs='*', help='tanda (paso tanda), o tandas y corridas (paso lista)')
    ap.add_argument('--solo', nargs='+', help='solo estas corridas')
    a = ap.parse_args()
    if a.solo and set(a.solo) - set(CORRIDAS):
        sys.exit('no existen: ' + ', '.join(sorted(set(a.solo) - set(CORRIDAS))))
    if a.paso == 'tanda':
        if len(a.nombres) != 1 or a.nombres[0] not in TANDAS:
            sys.exit('uso: nube_L.py tanda <' + '|'.join(TANDAS) + '>')
        tanda(a.nombres[0])
    elif a.paso == 'lista':
        lista(a.nombres)
    elif a.paso == 'prueba':
        sys.exit(prueba())
    elif a.paso == 'modelo':
        modelo(a.nombres)
    elif a.paso in ('consistencia', 'tabla', 'reticulas', 'figura', 'frontera', 'articulo', 'figura_fase'):
        globals()[a.paso]()
    else:
        globals()[a.paso](a.solo)
