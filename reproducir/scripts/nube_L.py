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
    python3 nube_L.py medir_capas [--solo ...]  masa del núcleo capa por capa de las corridas que no la tienen
    python3 nube_L.py medir_colas [--solo ...]  armónicos y perfiles junto a la órbita circular de las que no los tienen
    python3 nube_L.py tanda <nombre>            corre, analiza y mide una tanda (TANDAS), en orden
    python3 nube_L.py lista [tandas o corridas] nombres de las corridas (para los guiones del cluster)
    python3 nube_L.py prueba [corridas]         la cadena completa con una corrida corta, y su costo
    python3 nube_L.py modelo [corridas]         series del modelo reducido (serie:CPUs), o sus corridas
    python3 nube_L.py consistencia [corridas]   C021 frente a la corrida B021 de vlasov-poisson_PIC; CJ05, CJ20 y
                                                CK04 frente a J05_0.01, J20_0.08 y K_0.04 (otra s_J y Kepler)
    python3 nube_L.py tabla                     medidas frente al modelo reducido
    python3 nube_L.py reticulas                 coherencia libre de la retícula en L de cada corrida
    python3 nube_L.py figura                    figura del artículo (docs/articulo/figuras/bulto_L.pdf)
    python3 nube_L.py frontera                  dónde se pierde el núcleo en cada familia, y los cocientes B2 y chi
    python3 nube_L.py figura_chi                exceso de m_c frente a sigma_L y frente a chi (bulto_L_chi.pdf)

Necesita el repositorio vlasov-poisson_PIC junto a este (o en la variable VP_PIC_FIJO): de él salen
la retícula, las medidas y las corridas con las que se compara. Todo queda en exe/rep/12_nube/.
`tabla` solo lee los archivos pequeños que deja `medir` (medidas.npz, capas.npz, campo.npz,
resumen.npz), así que las corridas se pueden hacer en otra máquina y traer solo esos archivos
(reproducir/cluster/). Las corridas usan VP_HILOS hilos (4 por omisión) y el análisis VP_PROCESOS
procesos (4).
"""
import os, re, sys, subprocess, time, argparse, warnings, hashlib, datetime
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
import formas as FO
from demo_eta import envolvente
from aa_numerico_L import MapaAA, phi_iso, phi_esfera, aa_isocrono, aa_kepler
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
# en el .par. sj es la anchura s_J del dato y pot el potencial externo: 'isocrono' o 'kepler' (la masa
# puntual; en el código, BGtype = sphere, que coincide con -1/r en r >= 1, como en vlasov-poisson_PIC).
def _c(m, nl, sl, red=None, fijo=None, ancho=2.5, nj=200, nq=50, perfil=None, sj=SJ, pot='isocrono', **cambios):
    return dict(m=m, nj=nj, nq=nq, nl=nl, sl=sl, ancho=ancho, red=red, fijo=fijo, perfil=perfil, sj=sj, pot=pot,
                cambios=cambios)


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
# 2026-11-16: la ley con dos parámetros, Gamma y chi (predicción registrada antes de correr en
# vlasov-poisson_PIC/docs/articulo/prediction_L_scaling.md). chi = |dOmega/dL| sigma_L sqrt(J_med)/(g M_gas) (chi_de)
# compara la dispersión de las frecuencias de las capas de L con la rotación que la nube induce a su distancia de la
# órbita circular. En las series de M_gas/M = 0.02, 0.03 y 0.05 el exceso de m_c sobre m_mix cae a la mitad en
# chi = 0.20, 0.27 y 0.28 (`frontera`). Cada familia nueva lleva L igual y cuatro anchuras con chi = 0.15, 0.25, 0.35
# y 0.5 (s_L redondeado), con el menor N_L que cumple la cota de 7e-4 de `reticulas`, redondeado hacia arriba. Todas
# en la retícula de 200 x 50 en (Q, J) de las series anteriores (2100 partículas por capa) y hasta 30 tau_1.
#   M070, M100   dos masas mayores en el isócrono con s_J = 0.1 (Gamma = 0.80 y 1.14): ¿sigue chi_mitad cerca de 0.27?
#   J05M010      s_J = 0.05, M_gas/M = 0.01 (Gamma = 0.32; tau_1 = 589, t_fin = 17 700)
#   J20M080      s_J = 0.2, M_gas/M = 0.08 (Gamma = 0.32; tau_1 = 198, t_fin = 6000)
#   KM040        Kepler, s_J = 0.1, M_gas/M = 0.04 (Gamma = 0.34; tau_1 = 161, t_fin = 5000)
# Las tres últimas tienen el Gamma de la serie de M_gas/M = 0.03 (0.34), pero otra J_med u otro potencial: si la ley
# vale, la mitad cae en el mismo chi y en sigma_L/L_0 muy distintos (0.2-0.3, 0.9-1.1 y 0.4-0.5 % frente a 0.58 %).
# Con L igual son las nubes J05_0.01, J20_0.08 y K_0.04 de vlasov-poisson_PIC (bulto.py).
# La tabla de las familias y los nombres de sus corridas están en reducido.py (LEY, nombre_ley), de donde los toma
# también el modelo con las mismas partículas (series ley_<familia>).
LEY = R.LEY
nombre_ley = R.nombre_ley
for _f, (_m, _sj, _pot, _dtm, _tf, _S) in LEY.items():
    _op = dict(sj=_sj, pot=_pot, **({} if _tf == TFIN else {'Nt': int(round(_tf/DT))}))
    for _sl, _nl in ((0.0, 1),) + _S:
        _n, _r = nombre_ley(_f, _sl)
        CORRIDAS[_n] = _c(_m, _nl, _sl, _r, **_op)
# 2026-11-16: leyes de decaimiento de la mezcla y del campo (vlasov-poisson_PIC/reproducir/scripts/colas.py), en el
# régimen en que el gas se mezcla, como función de la autogravedad y de la forma de la FD. Retícula de 3200 x 100
# recortada (unas 62 000 partículas con la nube del artículo, N_ef ~ 7000), L igual, hasta 30 tau_1. El número del
# nombre es 1000 Gamma, con Gamma = g M/(|Omega_0'| J_med^(3/2)): la nube del artículo (p2) con Gamma = 0.025, 0.05,
# 0.10 y 0.15 (mu = 0.13 a 0.32), y las formas p4, h05, anillo y meseta con Gamma = 0.05 y 0.10. Con las retículas de
# 400 filas los armónicos n >= 2 y el potencial quedan dominados por el ruido de partículas desde unas 8 tau_1.
COLAS = {  # forma: (perfil, pasos si no son los de la serie, valores de Gamma)
    'p2':     (None, None, (0.025, 0.05, 0.10, 0.15)),
    'p4':     (('p', 4.0), None, (0.05, 0.10)),
    'h05':    (('hueco', 0.05), 236000, (0.05, 0.10)),
    'anillo': (('anillo', (0.15, 0.03)), 424000, (0.05, 0.10)),
    'meseta': (('meseta', (0.05, 0.25)), 272000, (0.05, 0.10)),
}


def masa_de_gamma(perfil, gam):
    """La masa del gas con la que la forma dada tiene el Gamma dado (cinco decimales)."""
    return round(gam*abs(FO.OMEGA_P)*FO.J_mediana(perfil)**1.5/R.G, 5)


for _f, (_pf, _nt, _G) in COLAS.items():
    for _g in _G:
        CORRIDAS[f'H_{_f}_{int(round(_g*1000)):03d}'] = _c(masa_de_gamma(_pf, _g), 1, 0.0, nj=3200, nq=100, perfil=_pf,
                                                         **({} if _nt is None else {'Nt': _nt}))
# La nube del artículo con Gamma = 0.10 en una retícula con el doble de filas y de columnas, para la convergencia.
CORRIDAS['H2_p2_100'] = _c(masa_de_gamma(None, 0.10), 1, 0.0, nj=6400, nq=200)
# Consistencia de las familias nuevas con el código de L fija: las nubes J05_0.01, J20_0.08 y K_0.04 de
# vlasov-poisson_PIC en su retícula de 400 x 25 (aquí recortada donde el peso baja de 1e-7 de su máximo), con L
# igual, hasta t = 2000. Comprueban el dato inicial, el fondo del código y el mapa de ángulo-acción con otra s_J y
# con la masa puntual (paso `consistencia`).
CORRIDAS['CJ05'] = _c(0.01, 1, 0.0, fijo='J05_0.01', nj=400, nq=25, sj=0.05, Nt=40000)
CORRIDAS['CJ20'] = _c(0.08, 1, 0.0, fijo='J20_0.08', nj=800, nq=25, sj=0.2, Nt=40000)
CORRIDAS['CK04'] = _c(0.04, 1, 0.0, fijo='K_0.04', nj=400, nq=25, pot='kepler', Nt=40000)
# Las medidas de bulto.py (banda, escalas, nucleo) leen s_J y el potencial externo de sus propias tablas.
for _n, _x in CORRIDAS.items():
    if _x['sj'] != SJ:
        B.CAMBIOS[_n] = dict(B.CAMBIOS.get(_n, {}), sr=str(_x['sj']))
    if _x['pot'] != 'isocrono':
        B.FONDO[_n] = 'puntual'
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
    # La ley con dos parámetros (2026-11-16), las más largas primero. La de 189 000 partículas va aparte, con 32 GB.
    'ley':             [nombre_ley(f, s)[0] for f in ('M100', 'M070', 'J20M080', 'KM040', 'J05M010')
                        for s in [x for x, _ in LEY[f][5]][::-1] + [0.0] if (f, s) != ('M100', 0.103)],
    'ley_grande':      ['M100_s103'],
    'consistencia_ley': ['CJ05', 'CJ20', 'CK04'],
    # Leyes de decaimiento (2026-11-16): las de más pasos primero. La de 6400 x 200 va aparte, con 32 GB.
    'colas':           ['H_anillo_100', 'H_anillo_050', 'H_meseta_100', 'H_meseta_050', 'H_h05_100', 'H_h05_050',
                        'H_p2_150', 'H_p2_100', 'H_p2_050', 'H_p2_025', 'H_p4_100', 'H_p4_050'],
    'colas_fina':      ['H2_p2_100'],
}
# Series del modelo reducido (reducido.py) que acompañan a las corridas, con las CPUs de su trabajo de Slurm, en el
# orden en que todo.sh las envía: la más larga primero.
SERIES_MODELO = (('Lb', 16), ('control', 16), ('Lf5g', 32), ('Lf5', 16), ('Lf2', 16), ('Lf5b', 32),
                 ('ley_M100', 32), ('ley_M070', 32), ('ley_J20M080', 16), ('ley_KM040', 16), ('ley_J05M010', 16))


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


def kepler_rp(Q, J, L):
    """(Q, J, L) -> (r, p_r) en el potencial de la masa puntual (GM = 1): a = (J + L)^2, e^2 = 1 - L^2/a,
    Q = eta - e sin eta, r = a (1 - e cos eta), p_r = sqrt(a) e sin eta/r. Es bulto.dato_kepler con L por partícula."""
    a = (J + L)**2
    e = np.sqrt(np.clip(1 - L*L/a, 0, None))
    eta = Q.copy()
    for _ in range(60):
        eta = eta - (eta - e*np.sin(eta) - Q)/(1 - e*np.cos(eta))
    r = a*(1 - e*np.cos(eta))
    return r, np.sqrt(a)*e*np.sin(eta)/r


def particulas(n):
    """Q, J, L, f y (r, p_r) de las partículas de la corrida n. f es el valor de la DF que lee el
    código; la masa de cada partícula es proporcional a f L."""
    c = CORRIDAS[n]
    if c['nl'] == 1:
        Q, J, f = R.reticula(c['nj'], c['nq'], sj=c['sj'], sq=SQ, perfil=c['perfil'])
        L = np.full_like(Q, L0)
    else:
        Q, J, L, w = R.reticula_L(c['nj'], c['nq'], c['nl'], c['sl'], L0, c['ancho'], sj=c['sj'], sq=SQ)
        f = w/L
    r, p = (isocrono_rp if c['pot'] == 'isocrono' else kepler_rp)(Q, J, L)
    return Q, J, L, f, r, p


def tfin_de(n):
    """Tiempo final de la corrida n: el de la serie, o el de su número de pasos si el .par lo cambia."""
    return CORRIDAS[n]['cambios'].get('Nt', int(round(TFIN/DT)))*DT


def kappa_de(n):
    """Frecuencia radial de la órbita circular del potencial externo de la corrida n."""
    return KAPPA if CORRIDAS[n]['pot'] == 'isocrono' else R.fondo(CORRIDAS[n]['pot'], L0)['kappa']


def capas_L(n):
    """Valores de L de la retícula de la corrida n y su peso (masa relativa de cada capa)."""
    c = CORRIDAS[n]
    if c['nl'] == 1:
        return np.array([L0]), np.array([1.0])
    Ln = L0 + c['ancho']*c['sl']*(2*(np.arange(c['nl']) + 0.5)/c['nl'] - 1.0)
    return Ln, np.exp(-(Ln - L0)**2/c['sl']**2)*Ln


def coherencia_libre(Ln, w, t, J=0.0, pot='isocrono'):
    """|sum_l w_l exp(-i Omega(J, L_l) t)|/sum_l w_l: lo que queda de la suma de las capas en la mezcla
    libre. Para una gaussiana continua es exp(-(dOmega/dL s_L t)^2/4); con N_L capas vuelve a 1 en t_rec."""
    om = (J + R.c_de_L(Ln, pot))**-3
    return np.abs((w[None, :]*np.exp(-1j*om[None, :]*np.asarray(t, float)[:, None])).sum(1))/w.sum()


def desfase_reticula(sl, nl, ancho, tfin, tau, pot='isocrono', sj=SJ):
    """Mayor diferencia, hasta tfin y hasta 20 tau, entre la coherencia libre de la retícula de nl puntos medios
    en |L - L_0| <= ancho s_L y la de la gaussiana continua (600 puntos en +-6 s_L), en J = 0 y en J = s_J."""
    t = np.linspace(0, tfin, 4001)
    Ln = L0 + ancho*sl*(2*(np.arange(nl) + 0.5)/nl - 1.0); w = np.exp(-(Ln - L0)**2/sl**2)*Ln
    fino = L0 + 6*sl*(2*(np.arange(600) + 0.5)/600 - 1.0); wf = np.exp(-(fino - L0)**2/sl**2)*fino
    dif = [np.abs(coherencia_libre(Ln, w, t, J, pot) - coherencia_libre(fino, wf, t, J, pot)) for J in (0.0, sj)]
    tarde = t >= 20*tau
    return max(d.max() for d in dif), max(d[~tarde].max() for d in dif)


def reticulas():
    """Para cada corrida con L distribuida: paso de la retícula, t_rec y la mayor diferencia, hasta el
    final de la corrida, entre la coherencia libre de la retícula y la de la gaussiana continua."""
    h = 1e-5
    dOdL = abs((R.c_de_L(L0 + h))**-3 - (R.c_de_L(L0 - h))**-3)/(2*h)
    print(f'|dOmega/dL| = {dOdL:.4f} en J = 0; fin de las corridas en t/tau_1 = {TFIN/TAU1:.1f}')
    for n, c in CORRIDAS.items():
        if c['nl'] == 1:
            continue
        Ln, w = capas_L(n); dL = Ln[1] - Ln[0]
        pot, tau, tfin = c['pot'], tau_de(n), tfin_de(n)
        dO = dOdL if pot == 'isocrono' else abs((R.c_de_L(L0 + h, pot))**-3 - (R.c_de_L(L0 - h, pot))**-3)/(2*h)
        todo, antes = desfase_reticula(c['sl'], c['nl'], c['ancho'], tfin, tau, pot, c['sj'])
        otra = '' if (pot, c['sj'], tfin) == ('isocrono', SJ, TFIN) else f'  [{pot}, s_J = {c["sj"]:g}, |dOmega/dL| = {dO:.4f}, fin en {tfin/tau:.1f} tau_1]'
        print(f'{n:11} s_L = {c["sl"]:5.3f}  N_L = {c["nl"]:2d} en ±{c["ancho"]:g} s_L  dL = {dL:.5f}  t_rec = {2*np.pi/dO/dL/tau:5.1f} tau_1'
              f'  max |S - S_gauss|: {todo:.1e} (hasta 20 tau_1: {antes:.1e}){otra}')


def forma(n):
    """Filas en J y columnas en Q que quedan tras el recorte, y puntos en L."""
    c = CORRIDAS[n]
    if c['perfil'] is not None:
        Q, J, _ = R.reticula(c['nj'], c['nq'], sj=c['sj'], sq=SQ, perfil=c['perfil'])
        return len(np.unique(J)), len(np.unique(Q)), 1
    filas, columnas = B.recorte_reticula(c['nj'], c['nq'], 1e-7, c['sj'], SQ)
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
        Qc, Jc, _ = (aa_isocrono if c['pot'] == 'isocrono' else aa_kepler)(r, p, L)
        eq = np.abs(np.angle(np.exp(1j*(Qc - Q)))).max(); ej = np.abs(Jc - J).max()
        m = f*L; nef = m.sum()**2/(m**2).sum()
        ic = ruta('ic', f'{n}.dat')
        np.savetxt(ic, np.column_stack([r, p, L, f]), fmt='%.17e')
        nota = ''
        if c['fijo']:                           # las mismas partículas que la corrida de L fija
            g = os.path.join(FIJO_BULTO, 'ic', f'{c["fijo"]}.dat')
            if os.path.exists(g):
                a = np.loadtxt(g)
                if len(a) == len(r):
                    nota = (f'; frente a {c["fijo"]}: max |dr| = {np.abs(a[:, 0] - r).max():.1e}, max |dp| = {np.abs(a[:, 1] - p).max():.1e},'
                            f' max |df/f| = {np.abs(a[:, 2]/f - 1).max():.1e}')
                else:                           # la retícula de L fija está sin recortar: se comparan las partículas comunes
                    d = np.abs(a[:, None, 0] - r[None, :]) + np.abs(a[:, None, 1] - p[None, :])
                    k = d.argmin(0)
                    nota = (f'; frente a {c["fijo"]} ({len(a)} partículas, sin recortar): las {len(r)} de aquí están allá con'
                            f' max |dr| + |dp| = {d.min(0).max():.1e} y max |df/f| = {np.abs(a[k, 2]/f/(a[k, 2]/f).mean() - 1).max():.1e}')
        sl = c['sl'] if c['nl'] > 1 else 0.25   # con una sola L el intervalo solo fija dlc, que se cancela
        v = [('state', 'checkpoint'), ('CheckPointfile', f'rep/{GRUPO}/ic/{n}.dat'),
             ('Nrc', nrc), ('Npc', npc), ('Nlc', nlc), ('a0', c['m']),
             ('l0', L0), ('sl', c['sl'] if c['nl'] > 1 else 0.5), ('sr', c['sj']), ('sp', SQ),
             ('lminc', L0 - (c['ancho'] if c['nl'] > 1 else 2.0)*sl), ('lmaxc', L0 + (c['ancho'] if c['nl'] > 1 else 2.0)*sl),
             ('rminc', 3.0), ('rmaxc', 14.0), ('pminc', -0.25), ('pmaxc', 0.25),
             ('dr', 0.1), ('rmin', 0.0), ('rmax', 20.0), ('pmax', 2.0), ('courant', 1.0),
             ('Nt', int(round(TFIN/DT))), ('integrator', 'yoshida4'), ('bsplineorder', 1),
             ('BGtype', 'Isochrone' if c['pot'] == 'isocrono' else 'sphere'), ('autointeraction', '.true.'),
             ('spatial_output', SALIDA), ('field_output', SALIDA), ('time_output', 100*SALIDA),
             ('output_format', 'raw'), ('directory', f'rep/{GRUPO}/{n}')]
        v = [(k, c['cambios'].get(k, x)) for k, x in v]
        assert not set(c['cambios']) - {k for k, _ in v}
        ancho_j = f's_Q = s_J = {SJ}' if c['sj'] == SJ else f's_Q = {SQ}, s_J = {c["sj"]:g}'
        if c['pot'] != 'isocrono':
            ancho_j += ', en las variables de ángulo-acción de la masa puntual (BGtype = sphere, r >= 1),'
        cab = [f'# Nube {n}: F_0(Q, J) con {ancho_j} por exp(-(L - {L0:g})^2/s_L^2), s_L = {c["sl"]:g},'
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
        npart = int(np.prod(forma(n))); est = s*tfin_de(n)/DT/pasos
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
    raw, indices, pot = tarea
    c = CorridaRaw(raw)
    rg = c.r
    if pot == 'isocrono':
        phib = phi_iso(rg)
    else:                                       # BGtype = sphere: la masa puntual en r >= 1
        phib = phi_esfera(rg)
    sal = []
    for i in indices:
        s = c[i]
        r, p, L = s['r_part'], s['p_part'], s['l_part']
        ps = s['potential'] - phib
        if pot == 'isocrono':
            Q, J, _ = MapaL(rg, ps)(r, p, L)    # potencial de este instante, L de cada partícula
            Qd, Jd, _ = aa_isocrono(r, p, L)    # isócrono desnudo
        else:
            ps = np.where(rg < 1, ps[np.argmin(np.abs(rg - 1))], ps)   # las órbitas no bajan de r = 1 (como bulto.py)
            Q, J, _ = MapaL(rg, ps, fondo=pot)(r, p, L)
            Qd, Jd, _ = aa_kepler(r, p, L)      # masa puntual desnuda
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
        grupos = [(raw, list(range(i, min(i + 25, nt))), CORRIDAS[n]['pot']) for i in range(0, nt, 25)]
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
    """tau_1 de la corrida n: el de su forma del dato en J (reducido.banda_perfil), el de su anchura s_J y su
    potencial externo (reducido.tau1) o el de la nube del artículo."""
    c = CORRIDAS[n]; pf = c['perfil']
    if n.startswith('prueba_'):                 # la prueba de la cadena acorta TAU1 para que sus ventanas caigan en ella
        return TAU1
    if pf is None:
        return TAU1 if (c['sj'], c['pot']) == (SJ, 'isocrono') else R.tau1(c['sj'], c['pot'])
    return R.banda_perfil(pf, SJ)[3]


def plano(npz, sal, tau=TAU1, sj=SJ, kappa=KAPPA):
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
    M = R.medidas_plano(t, X, Y, w[ok], sj, tau, kappa)
    np.savez(sal, fuera=1.0 - w[ok].sum()/w.sum(), **M)
    return np.load(sal)


def medidas(n):
    return plano(ruta(n, 'bulto.npz'), ruta(n, 'medidas.npz'), tau_de(n), CORRIDAS[n]['sj'], kappa_de(n))


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


def nucleo_por_capas(n):
    """Masa del núcleo de cada capa de L: la fracción de la masa de la capa dentro del disco de acción 0.05 s_J
    centrado en el máximo de densidad de esa misma capa, media en t >= 20 tau_1 (reducido.medidas_plano sobre las
    partículas de la capa). La m_c de la tabla usa un solo disco para todo el gas, así que baja tanto si las capas
    pierden su núcleo como si lo conservan pero quedan desfasadas; esta medida solo baja en el primer caso. Devuelve
    L, la masa W de cada capa y su m_c; se guarda en <n>/nucleo_capas.npz. Necesita <n>/bulto.npz."""
    sal = ruta(n, 'nucleo_capas.npz')
    if not os.path.exists(sal):
        d = np.load(ruta(n, 'bulto.npz'))
        t, w, L = d['t'], np.asarray(d['w'], float), d['L']; Q, J = d['Q'], d['J']
        ok = np.isfinite(Q).all(0) & np.isfinite(J).all(0)
        v, inv = np.unique(L, return_inverse=True); inv = inv.ravel()
        mc = np.empty(len(v)); W = np.empty(len(v))
        for l in range(len(v)):
            s = np.nonzero(ok & (inv == l))[0]
            a = np.sqrt(2.0*J[:, s].astype(float)); q = Q[:, s].astype(float)
            M = R.medidas_plano(t, (a*np.cos(q)).astype(np.float32), (a*np.sin(q)).astype(np.float32), w[s],
                                CORRIDAS[n]['sj'], tau_de(n), kappa_de(n))
            mc[l], W[l] = float(M['mnuc']), w[s].sum()
        np.savez(sal, L=v, W=W, mc=mc)
    d = np.load(sal)
    return d['L'], d['W'], d['mc']


def mc_capas(n):
    """Media sobre las capas de L, pesada con su masa, de la masa del núcleo de cada capa; None si no se ha medido."""
    if not os.path.exists(ruta(n, 'nucleo_capas.npz')):
        return None
    _, W, mc = nucleo_por_capas(n)
    return float((W*mc).sum()/W.sum())


def medir_capas(solo=None):
    """nucleo_por_capas de las corridas con L distribuida que tienen bulto.npz y aún no la tienen (las anteriores
    al 2026-11-16 se midieron sin ella; sus bulto.npz están donde se corrieron)."""
    for n, c in CORRIDAS.items():
        if (solo and n not in solo) or c['nl'] == 1 or os.path.exists(ruta(n, 'nucleo_capas.npz')):
            continue
        if not os.path.exists(ruta(n, 'bulto.npz')):
            print(f'  {n}: no hay bulto.npz')
            continue
        t0 = time.time(); nucleo_por_capas(n)
        print(f'  {n}: m_c por capa {mc_capas(n):.4f} ({time.time()-t0:.0f} s)', flush=True)


def armonicos(n):
    """h_k(t)/h_0, k = 0..5, de la corrida n (mapa del potencial de cada instante), copiados de bulto.npz a un
    archivo pequeño, <n>/armonicos.npz, para traerlos del cluster."""
    sal = ruta(n, 'armonicos.npz')
    if not os.path.exists(sal):
        d = np.load(ruta(n, 'bulto.npz'))
        np.savez(sal, t=d['t'], hk=d['hk'])
    return np.load(sal)


def perfiles(n, tiempos=(3, 5, 10, 15, 20, 25), nb=40, xmax=0.4):
    """La FD junto a la órbita circular en t/tau_1 = tiempos: en nb franjas de J entre 0 y xmax s_J, la masa por
    unidad de J (f0, la FD promediada en el ángulo salvo un factor) y el módulo de sus tres primeros armónicos,
    a_k = |sum_i m_i exp(-i k Q_i)|/dJ, todo por unidad de masa total. Cada valor es la media de las cinco
    instantáneas alrededor del tiempo pedido. La cola de h_k la fija el exponente de a_k en J -> 0 (Sección 3 del
    artículo): a_k ~ J^q da t^-(q+1). Se guarda en <n>/perfiles.npz. Necesita <n>/bulto.npz."""
    sal = ruta(n, 'perfiles.npz')
    if not os.path.exists(sal):
        d = np.load(ruta(n, 'bulto.npz'))
        t, w = d['t'], np.asarray(d['w'], float); Q, J = d['Q'], d['J']
        sj, tau = CORRIDAS[n]['sj'], tau_de(n)
        x = np.linspace(0.0, xmax, nb + 1); bordes = x*sj; dJ = np.diff(bordes); W = w.sum()
        f0 = np.zeros((len(tiempos), nb)); a = np.zeros((len(tiempos), 3, nb)); tt = np.zeros(len(tiempos))
        for i, T in enumerate(tiempos):
            k0 = int(np.argmin(np.abs(t - T*tau))); ks = [k for k in range(k0 - 2, k0 + 3) if 0 <= k < len(t)]
            tt[i] = t[k0]
            for k in ks:
                j = J[k].astype(float); q = Q[k].astype(float)
                ok = np.isfinite(j) & np.isfinite(q) & (j < bordes[-1])
                b = np.clip(np.searchsorted(bordes, j[ok], side='right') - 1, 0, nb - 1)
                f0[i] += np.bincount(b, weights=w[ok], minlength=nb)/len(ks)
                for m in (1, 2, 3):
                    z = w[ok]*np.exp(-1j*m*q[ok])
                    a[i, m - 1] += np.abs(np.bincount(b, weights=z.real, minlength=nb) + 1j*np.bincount(b, weights=z.imag, minlength=nb))/len(ks)
        np.savez(sal, x=x, t=tt, tau=tau, f0=f0/dJ/W, a=a/dJ/W)
    return np.load(sal)


def medir_colas(solo=None):
    """armonicos y perfiles de las corridas con una sola L que tienen bulto.npz y aún no los tienen (las anteriores
    al 2026-11-16 se midieron sin ellos; sus bulto.npz están donde se corrieron)."""
    for n, c in CORRIDAS.items():
        if (solo and n not in solo) or c['nl'] != 1 or not os.path.exists(ruta(n, 'bulto.npz')):
            continue
        if not (os.path.exists(ruta(n, 'armonicos.npz')) and os.path.exists(ruta(n, 'perfiles.npz'))):
            armonicos(n); perfiles(n)
            print(f'  {n}: armonicos.npz y perfiles.npz', flush=True)


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
            B.nucleo(n); armonicos(n); perfiles(n)
        else:
            nucleo_por_capas(n)
        r = resumen(n)
        print(f'  medido {n}: mu = {float(r["mu"]):.3f}, |dE/E| <= {float(r["dE"]):.1e} ({time.time()-t0:.0f} s)', flush=True)


def numeros(M, tau=None):
    """Los números de la tabla a partir de las medidas en el plano de zeta. Las ventanas van en unidades de tau
    (por omisión el tau_1 de la nube del artículo, también para las formas del dato en J: así está la tabla desde
    el 2026-11-15; las familias con otra s_J u otro potencial pasan su propio tau_1)."""
    tau = TAU1 if tau is None else tau
    t = M['t']; u = t >= 20*tau; v = (t >= 10*tau) & (t <= 20*tau)
    env = envolvente(M['ac'], t)
    z = (t >= 2*tau) & (t <= 30*tau)
    k = np.argmin(np.where(z, env, np.inf))
    return dict(h1=np.abs(M['h1'][v]).mean(), h1t=np.abs(M['h1'][u]).mean(), rho_c=M['rho_c'][u].mean(), m02=M['m02'][u].mean(),
                mc=float(M['mnuc']), dr=M['dr'][u].mean(), ac=[env[np.argmin(np.abs(t - x*tau))] for x in TS],
                acmin=env[k], tmin=t[k]/tau, actarde=env[t >= 25*tau].mean())


def tau_tabla(n):
    """tau_1 con el que la tabla abre sus ventanas: el de la corrida si es de una familia con otra s_J u otro
    potencial; si no, el de la nube del artículo (ver numeros)."""
    c = CORRIDAS[n]
    return TAU1 if (c['sj'], c['pot']) == (SJ, 'isocrono') else tau_de(n)


def campo(n):
    """Media en 10 <= t/tau_1 <= 20 de la desviación D del potencial del gas (bulto.campo_particulas)."""
    t, D = B.campo_particulas(n); tau = tau_tabla(n)
    return D[(t >= 10*tau) & (t <= 20*tau)].mean()


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


def consistencia_familia(n):
    """La corrida n (CJ05, CJ20, CK04: otra s_J o la masa puntual, con L igual) frente a la de L fija de
    vlasov-poisson_PIC con la misma nube. Las partículas no son idénticas: aquí la retícula está recortada donde el
    peso baja de 1e-7 de su máximo y allá no, así que se comparan curvas y no partículas: |h_1|(t) con el mapa del
    potencial del instante y con el del potencial externo, la masa con J < 0.2 s_J y el potencial del gas."""
    c = CORRIDAS[n]; f = c['fijo']; sj, tau = c['sj'], tau_de(n)
    a = np.load(ruta(n, 'bulto.npz')); b = np.load(os.path.join(FIJO_BULTO, f, 'bulto.npz'))
    k = min(len(a['t']), len(b['t'])); t = a['t'][:k]
    assert np.allclose(t, b['t'][:k], atol=1e-6)
    wa, wb = np.asarray(a['w'], float), np.asarray(b['w'], float)
    print(f'{n} (VlasovPoisson_PIC_sp, {len(wa)} partículas) frente a {f} (vlasov-poisson_PIC, {len(wb)}): {c["pot"]}, s_J = {sj:g}, '
          f'tau_1 = {tau:.1f}, hasta t = {t[-1]:.0f} = {t[-1]/tau:.1f} tau_1')
    print(f'  masa de la retícula completa que queda fuera de las {len(wa)} partículas de más peso: '
          f'{1 - np.sort(wb/wb.sum())[::-1][:len(wa)].sum():.1e}')
    za = (a['r'] >= B.RVENT[0]) & (a['r'] <= B.RVENT[1]); zb = (b['r'] >= B.RVENT[0]) & (b['r'] <= B.RVENT[1])
    assert np.allclose(a['r'][za], b['r'][zb])
    m02 = lambda d, w: np.array([w[d['J'][i] < 0.2*sj].sum() for i in range(k)])/w.sum()
    ma, mb = m02(a, wa), m02(b, wb)
    filas = (('|h_1|, mapa del instante', np.abs(a['hk'][:k, 1]), np.abs(b['hk'][:k, 1])),
             ('|h_1|, mapa del potencial externo', np.abs(a['hk_desnudo'][:k, 1]), np.abs(b['hk_desnudo'][:k, 1])),
             ('masa con J < 0.2 s_J', ma, mb))
    print('  mayor diferencia en cada tramo de t/tau_1, relativa al máximo de la curva de L fija en el tramo (y su media):')
    tramos = [(x, y) for x, y in ((0, 1), (1, 2), (2, 5), (5, 10), (10, 20)) if x*tau < t[-1]]
    for q, x, y in filas:
        print(f'    {q:34} ' + '  '.join(f'[{u}, {v}]: {np.abs(x - y)[z].max()/y[z].max():.1e} ({y[z].mean():.3e})'
                                        for u, v in tramos for z in [(t >= u*tau) & (t <= v*tau)]))
    dp = np.abs(a['phi_self'][:k][:, za] - b['phi_self'][:k][:, zb]).max(1)/np.abs(b['phi_self'][:k][:, zb]).max(1)
    print('    potencial del gas en 3 <= r <= 15, max |dPhi|/max |Phi|:  '
          + '  '.join(f'[{u}, {v}]: {dp[(t >= u*tau) & (t <= v*tau)].max():.1e}' for u, v in tramos))
    B.BASE = FIJO_BULTO; Ab = B.escalas(f)
    B.BASE = BASE; Aa = B.escalas(n)
    print(f'  amplitud A del potencial en la primera vuelta y mu: {Aa[0]:.5e} y {Aa[2]:.4f} frente a {Ab[0]:.5e} y {Ab[2]:.4f}')
    ea, eb = a['energia'], None
    print(f'  |dE/E| de {n}: {np.abs(ea/ea[0] - 1).max():.1e}')


def tabla():
    """Medidas de cada corrida frente al modelo reducido: exe/rep/12_nube/tabla.txt."""
    lineas = ['Nube con L distribuida: PIC (VlasovPoisson_PIC_sp) y modelo reducido. sigma_L = s_L/sqrt(2).',
              '|h_1|: media en 10-20 tau_1 y en t >= 20 tau_1. capa: media sobre las capas de L del |h_1| de cada una, y coh: cociente',
              'entre el |h_1| total y esa media (1: capas en fase), las dos en 10-20 y en t >= 20 tau_1. rho_c, m02 y dr: medias en',
              't >= 20 tau_1; m_c: media en t >= 20 tau_1. a_c: mínimo de la envolvente en 2-30 tau_1 (y cuándo) y media en t >= 25 tau_1.',
              'D: media en 10-20 tau_1.',
              f'{"corrida":11} {"M_gas":>5} {"s_L":>5} {"N_L":>3} {"":8} {"|h_1|":>8} {"tarde":>8} {"capa":>8} {"tarde":>8} {"coh":>5} {"tarde":>5} '
              f'{"rho_c":>5} {"m02":>5} {"m_c":>5} {"dr":>5} {"a_c min":>7} {"t/tau1":>6} {"a_c":>5} {"D":>8} {"mu":>4}']

    def fila(nombre, N, capas, tau):
        v = lambda t: (t >= 10*tau) & (t <= 20*tau)
        u = lambda t: t >= 20*tau
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
        tau = tau_tabla(n)
        lineas.append(fila(f'{n:11} {c["m"]:5g} {c["sl"]:5g} {c["nl"]:3d} {"PIC":8}', numeros(medidas(n), tau), por_capas(n), tau)
                      + f' {campo(n):8.2e} {float(resumen(n)["mu"]):4.2f}')
        if c['red'] and os.path.exists(os.path.join(R.DATOS, c['red'] + '_diag.npz')):
            lineas.append(fila(f'{"":11} {"":5} {"":5} {"":3} {"reducido":8}', numeros(R.diagnosticos(c['red']), tau),
                               por_capas_reducido(c['red']), tau))
    open(ruta('tabla.txt'), 'w').write('\n'.join(lineas) + '\n')
    print('\n'.join(lineas))


def tanda(nombre):
    """Corre, analiza y mide las corridas de una tanda, una detrás de otra, y rehace la tabla."""
    for n in TANDAS[nombre]:
        correr([n]); analizar([n]); medir([n])
    tabla()


def prueba(n='M030_s010v', pasos=2000):
    # (con nombre de corrida: python3 nube_L.py prueba <corrida>)
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
              f'{os.environ.get("VP_HILOS", "4")} hilos, dt = {dt:g}; la corrida completa tardaría {s*tfin_de(n)/DT/pasos/60:.0f} min', flush=True)
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
                  f'el de la corrida completa tardaría {ta*(tfin_de(n)/DT/SALIDA + 1)/nsal/60:.0f} min', flush=True)
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


def familias_L():
    """Por familia (potencial, s_J, masa del gas), las corridas medidas de las retículas sin recurrencia (nombres
    <familia>_s<número>, sin sufijo v, n ni r), ordenadas por s_L, con los números del PIC (P), mu, D, la coherencia
    entre capas (media en 10-20 tau_1; 1 con L igual) y, si existen, los del modelo con las mismas partículas (M)."""
    out = {}
    for n, c in CORRIDAS.items():
        if not (re.fullmatch(r'[A-Z0-9]+_s\d+', n) and c['perfil'] is None and os.path.exists(ruta(n, 'resumen.npz'))):
            continue
        tau = tau_tabla(n)
        tc, _, _, coh = por_capas(n)
        f = dict(n=n, sl=c['sl'], nl=c['nl'], x=100*c['sl']/(np.sqrt(2)*L0), P=numeros(medidas(n), tau), mu=float(resumen(n)['mu']),
                 D=campo(n), coh=coh[(tc >= 10*tau) & (tc <= 20*tau)].mean(),
                 mcc=(numeros(medidas(n), tau)['mc'] if c['nl'] == 1 else mc_capas(n)))
        if c['red'] and os.path.exists(os.path.join(R.DATOS, c['red'] + '_diag.npz')):
            f['M'] = numeros(R.diagnosticos(c['red']), tau)
        out.setdefault((c['pot'], c['sj'], c['m']), []).append(f)
    return {k: sorted(F, key=lambda f: f['sl']) for k, F in out.items()}


MASAS_SERIE = (0.01, 0.02, 0.03, 0.05)          # las de la Tabla y la primera figura de la Sección 8 del artículo


def series_L():
    """familias_L() para las series de la Sección 8 del artículo (isócrono, s_J = 0.1, MASAS_SERIE), por masa del
    gas. Las familias de la prueba de la ley (LEY) van en `frontera`, `figura_chi` y `tabla_chi`."""
    return {m: F for (pot, sj, m), F in familias_L().items() if (pot, sj) == ('isocrono', SJ) and m in MASAS_SERIE}


def J_mediana(sj=SJ):
    """Mediana de la acción del dato (pesada con la masa) en la retícula de 400 x 100, como formas.J_mediana de
    vlasov-poisson_PIC: (2 J_med)^(1/2) es la distancia de la nube a la órbita circular en el plano de zeta."""
    Q, J, w = R.reticula(400, 100, sj=sj, sq=SQ)
    o = np.argsort(J); c = np.cumsum(w[o])/w.sum()
    return J[o][np.searchsorted(c, 0.5)]


def constantes(pot='isocrono', sj=SJ):
    """Constantes del modelo epicíclico para la familia (potencial, s_J), con L_0 = 2: acoplamiento g, |Omega_0'| =
    |dOmega/dJ| y |dOmega/dL| en la órbita circular, y la mediana J_med de la acción del dato."""
    F = R.fondo(pot, L0); h = 1e-5
    dOdL = abs(R.c_de_L(L0 + h, pot)**-3 - R.c_de_L(L0 - h, pot)**-3)/(2*h)
    return dict(g=F['g'], dOdJ=3.0*F['c']**-4, dOdL=dOdL, jmed=J_mediana(sj))


def gamma_de(m, pot='isocrono', sj=SJ):
    """Gamma = g M_gas/(|Omega_0'| J_med^(3/2)), Eq. (Gamma) del artículo."""
    k = constantes(pot, sj)
    return k['g']*m/(k['dOdJ']*k['jmed']**1.5)


def chi_de(m, sl, pot='isocrono', sj=SJ):
    """chi = |dOmega/dL| sigma_L sqrt(J_med)/(g M_gas), con sigma_L = s_L/sqrt(2): la dispersión de las frecuencias de
    las capas de L frente a la velocidad angular g M_gas/sqrt(J_med) de la rotación que la nube induce a su distancia
    sqrt(J_med) de la órbita circular. Es el tiempo de transporte sqrt(J_med)/(g M_gas) entre el tiempo de desfase de
    las capas, como 1/Gamma lo es entre el tiempo de la cizalla en J. chi = c'(L_0) sigma_L/(Gamma J_med)."""
    k = constantes(pot, sj)
    return k['dOdL']*sl/np.sqrt(2)*np.sqrt(k['jmed'])/(k['g']*m)


def sl_de_chi(chi, m, pot='isocrono', sj=SJ):
    """El s_L con el que la familia (potencial, s_J, masa) tiene el chi dado."""
    return chi/chi_de(m, 1.0, pot, sj)


def _mitad(x, y, nivel=0.5):
    """Primer cruce hacia abajo de y por nivel*y[0] (y[0] es el valor con L igual), interpolando linealmente en x."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    k = np.nonzero(y < nivel*y[0])[0]
    if len(x) < 2 or len(k) == 0:
        return None
    k = k[0]
    return x[k - 1] + (y[k - 1] - nivel*y[0])/(y[k - 1] - y[k])*(x[k] - x[k - 1])


def frontera(verboso=True):
    """Dónde se pierde el núcleo en cada familia (potencial, s_J, masa): el s_L en que el exceso m_c - m_mix cae a
    la mitad de su valor con L igual (interpolación lineal en s_L), en el PIC y en el modelo, y dos cocientes en ese
    punto. B2 (relevo del 2026-11-12): la dispersión |dOmega/dL| sigma_L de las frecuencias de las capas de L entre
    la rotación g m_c M_gas/rho que el núcleo induce en su borde, rho = (0.1 s_J)^{1/2}, con la m_c del PIC con L
    igual. chi (2026-11-15): la misma dispersión entre la rotación g M_gas/sqrt(J_med) que la nube entera induce a
    su distancia de la órbita circular (chi_de). También el s_L en que la coherencia entre capas (media en 10-20
    tau_1; es una cota superior, depende de N_L) cae a 0.5. Devuelve, por familia, un diccionario con Gamma, las
    curvas (chi, exceso relativo) del PIC y del modelo y los valores de la mitad."""
    p = print if verboso else (lambda *a, **k: None)
    res = {}
    for (pot, sj, m), F in sorted(familias_L().items()):
        if len(F) < 2 or F[0]['sl'] != 0:
            continue
        k = constantes(pot, sj); dOdL = k['dOdL']; rho = np.sqrt(0.1*sj)
        nube = k['g']*m/np.sqrt(k['jmed'])
        p(f'\n{pot}, s_J = {sj:g}: |dOmega/dL| = {dOdL:.4f} en J = 0, g = {k["g"]:.5f}, rho = {rho:.3f}, m_mix = {M_MIX:.4f}, '
          f'J_med = {k["jmed"]:.4f}')
        p(f'M_gas/M = {m:g} (mu = {F[0]["mu"]:.2f}, Gamma = {gamma_de(m, pot, sj):.3f}): s_L, sigma_L/L_0 (%), N_L, chi | m_c, |h_1|, '
          'masa con J < 0.2 s_J: PIC (modelo) | D, coherencia')
        for f in F:
            P, Mm = f['P'], f.get('M')
            par = lambda q: f'{P[q]:.3g}' + (f' ({Mm[q]:.3g})' if Mm else '')
            p(f'   {f["sl"]:6.4f} {f["x"]:5.2f} {f["nl"]:3d} {chi_de(m, f["sl"], pot, sj):5.3f} | m_c {par("mc")}  |h_1| {par("h1")}  '
              f'm02 {par("m02")} | D {f["D"]:.2e}  coh {f["coh"]:.2f}' + ('' if f['mcc'] is None else f'  m_c por capa {f["mcc"]:.3f}'))
        if F[0]['P']['mc'] < 1.5*M_MIX:
            p('   con L igual no hay núcleo (m_c compatible con m_mix): no hay frontera que medir')
            continue
        r = dict(gamma=gamma_de(m, pot, sj), mu=F[0]['mu'], mc0=F[0]['P']['mc'])
        w = k['g']*F[0]['P']['mc']*m/rho
        for clave, nombre in (('P', 'PIC'), ('M', 'modelo')):
            G = [f for f in F if clave in f]
            s = np.array([f['sl'] for f in G]); e = np.array([f[clave]['mc'] for f in G]) - M_MIX
            r[clave] = (np.array([chi_de(m, x, pot, sj) for x in s]), e/e[0], 100*s/np.sqrt(2)/L0)
            sm = _mitad(s, e)
            if sm is None:
                p(f'   {nombre}: el exceso no baja a la mitad en las corridas medidas')
                continue
            r['mitad_' + clave] = sm
            p(f'   {nombre}: la mitad del exceso en s_L = {sm:.4f} (sigma_L/L_0 = {100*sm/np.sqrt(2)/L0:.2f} %); '
              f'chi = {chi_de(m, sm, pot, sj):.3f}')
        p(f'   B2: g m_c M_gas/rho = {w:.2e}; cociente 0.5 en s_L = {0.5*w*np.sqrt(2)/dOdL:.4f} y 1 en s_L = {w*np.sqrt(2)/dOdL:.4f}'
          + (f'; en la mitad del PIC vale {dOdL*r["mitad_P"]/np.sqrt(2)/w:.2f}' if 'mitad_P' in r else ''))
        p(f'   chi: g M_gas/sqrt(J_med) = {nube:.2e}; chi = 0.25 en s_L = {sl_de_chi(0.25, m, pot, sj):.4f}')
        sc = _mitad([f['sl'] for f in F], [f['coh'] for f in F])
        if sc is not None:
            r['mitad_coh'] = sc
            p(f'   coherencia entre capas 0.5 en s_L = {sc:.4f} (sigma_L/L_0 = {100*sc/np.sqrt(2)/L0:.2f} %): chi = '
              f'{chi_de(m, sc, pot, sj):.3f}, cociente B2 = {dOdL*sc/np.sqrt(2)/w:.2f}')
        res[(pot, sj, m)] = r
    return res


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


# Familias con otra s_J u otro potencial en la figura de chi: color y marcador (paleta categórica validada de la guía
# de gráficas; Kepler lleva el violeta que ya tiene en la figura de las formas) y rótulo.
ESTILO_FAM = {('isocrono', 0.05): ('#eb6834', 's', r'$\sigma_J=0.05$'),
              ('isocrono', 0.2):  ('#1baf7a', '^', r'$\sigma_J=0.2$'),
              ('kepler', 0.1):    ('#6250d6', 'X', 'Kepler')}


def figura_chi():
    """Figura de la Sección 8D del artículo: el exceso de la masa del núcleo sobre m_mix, relativo a su valor con L
    igual, frente a sigma_L/L_0 (izquierda) y frente a chi (derecha), del PIC (líneas continuas) y del modelo con las
    mismas partículas (discontinuas), para las familias en las que hay núcleo con L igual. La familia del artículo
    (isócrono, s_J = 0.1) va en tonos de azul, de claro a oscuro al crecer la masa; las de otra s_J u otro potencial,
    en los colores de ESTILO_FAM. Escribe docs/articulo/figuras/bulto_L_chi.pdf en vlasov-poisson_PIC, una copia
    bulto_L_chi.png en exe/rep/12_nube/, e imprime los valores de chi en la mitad del exceso."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    sup, tinta, tinta2, reja = '#ffffff', '#0b0b0b', '#52514e', '#e4e3df'
    plt.rcParams.update({'font.size': 9, 'axes.edgecolor': tinta2, 'axes.labelcolor': tinta2, 'xtick.color': tinta2,
                         'ytick.color': tinta2, 'text.color': tinta, 'axes.linewidth': 0.6, 'figure.facecolor': sup,
                         'axes.facecolor': sup, 'savefig.facecolor': sup})
    res = frontera(verboso=False)
    ref = sorted(k for k in res if k[:2] == ('isocrono', SJ))
    otras = sorted(k for k in res if k[:2] != ('isocrono', SJ))
    tonos = dict(zip(ref, plt.cm.Blues(np.linspace(0.45, 0.95, max(len(ref), 2)))))
    estilo = {k: (tonos[k], 'o', rf'$M_{{gas}}/M={k[2]:g}$ ($\Gamma={res[k]["gamma"]:.2f}$)') for k in ref}
    for k in otras:
        c, mk, rot = ESTILO_FAM[k[:2]]
        estilo[k] = (c, mk, rf'{rot}, $M_{{gas}}/M={k[2]:g}$ ($\Gamma={res[k]["gamma"]:.2f}$)')
    fig, ax = plt.subplots(1, 2, figsize=(7.6, 3.5))
    for a, col, xlab in ((ax[0], 2, r'$\sigma_L/L_0$ (percent)'), (ax[1], 0, r'$\chi$')):
        for k in ref + otras:
            c, mk, _ = estilo[k]
            a.plot(res[k]['P'][col], res[k]['P'][1], color=c, lw=1.3, marker=mk, ms=5, mec=sup, mew=0.7)
            if 'M' in res[k] and len(res[k]['M'][0]) > 1:
                a.plot(res[k]['M'][col], res[k]['M'][1], color=c, lw=1.1, ls=(0, (4, 2.5)), marker='D', ms=3.5, mec=sup, mew=0.6)
        a.axhline(0.5, color=tinta2, lw=0.8, ls=':', zorder=0)
        a.axhline(0.0, color=tinta2, lw=0.6, zorder=0)
        a.set_xlim(left=-0.03*a.get_xlim()[1]); a.set_ylim(-0.25, 1.1)
        a.grid(True, axis='y', color=reja, lw=0.6); a.set_axisbelow(True)
        for lado in ('top', 'right'):
            a.spines[lado].set_visible(False)
        a.tick_params(length=3, width=0.6)
        a.set_xlabel(xlab)
    ax[0].set_title(r'excess of the core mass, against $\sigma_L/L_0$', loc='left', fontsize=9.5, color=tinta, pad=6)
    ax[1].set_title(r'against $\chi$', loc='left', fontsize=9.5, color=tinta, pad=6)
    ax[0].set_ylabel(r'$(m_c-m_{mix})/(m_c-m_{mix})_{\sigma_L=0}$')
    h = [Line2D([], [], color=estilo[k][0], lw=1.6, marker=estilo[k][1], ms=5, mec=sup, label=estilo[k][2]) for k in ref + otras]
    filas = (len(h) + 2)//3
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.065*filas - 0.02), w_pad=1.6)
    fig.legend(handles=h, loc='upper center', ncol=3, frameon=False, fontsize=8.2, bbox_to_anchor=(0.5, 1.0),
               handlelength=2.4, columnspacing=1.4)
    ax[1].legend(handles=[Line2D([], [], color=tinta2, lw=1.3, marker='o', ms=5, mec=sup, label='PIC'),
                          Line2D([], [], color=tinta2, lw=1.1, ls=(0, (4, 2.5)), marker='D', ms=3.5, mec=sup, label='epicyclic model')],
                 loc='upper right', frameon=False, fontsize=8.2, handlelength=2.6)
    pdf = os.path.join(FIJO, 'docs', 'articulo', 'figuras', 'bulto_L_chi.pdf')
    fig.savefig(pdf); fig.savefig(ruta('bulto_L_chi.png'), dpi=110); plt.close(fig)
    print('escrita', pdf)
    print('chi en la mitad del exceso de m_c (PIC; modelo), cociente B2 en ese punto y chi con la coherencia 0.5:')
    for k in ref + otras:
        r = res[k]; pot, sj, m = k
        kc = constantes(pot, sj); w = kc['g']*r['mc0']*m/np.sqrt(0.1*sj)
        txt = [f'{pot}, s_J = {sj:g}, M_gas/M = {m:g} (Gamma = {r["gamma"]:.3f})']
        if 'mitad_P' in r:
            txt.append(f'sigma_L/L_0 = {100*r["mitad_P"]/np.sqrt(2)/L0:.2f} %, chi = {chi_de(m, r["mitad_P"], pot, sj):.3f}, '
                       f'B2 = {kc["dOdL"]*r["mitad_P"]/np.sqrt(2)/w:.2f}')
        if 'mitad_M' in r:
            txt.append(f'modelo: {100*r["mitad_M"]/np.sqrt(2)/L0:.2f} %, chi = {chi_de(m, r["mitad_M"], pot, sj):.3f}')
        if 'mitad_coh' in r:
            txt.append(f'coherencia 0.5: chi = {chi_de(m, r["mitad_coh"], pot, sj):.3f}, B2 = {kc["dOdL"]*r["mitad_coh"]/np.sqrt(2)/w:.2f}')
        print('  ' + '; '.join(txt))
    con = [k for k in ref + otras if 'mitad_P' in res[k]]
    print('  para el artículo (dos decimales): chi en la mitad, PIC: ' + ', '.join(f'{chi_de(k[2], res[k]["mitad_P"], k[0], k[1]):.2f}' for k in con)
          + '; modelo: ' + ', '.join(f'{chi_de(k[2], res[k]["mitad_M"], k[0], k[1]):.2f}' for k in con if 'mitad_M' in res[k])
          + '; sigma_L de la mitad/(Gamma J_med): '
          + ', '.join(f'{res[k]["mitad_P"]/np.sqrt(2)/(res[k]["gamma"]*constantes(k[0], k[1])["jmed"]):.2f}' for k in con))
    print('  chi de las corridas de cada familia: ' + '; '.join(f'{k[2]:g}: ' + ' '.join(f'{x:.2f}' for x in res[k]['P'][0]) for k in ref + otras))
    xs = [100*res[k]['mitad_P']/np.sqrt(2)/L0 for k in con]; cs = [chi_de(k[2], res[k]['mitad_P'], k[0], k[1]) for k in con]
    cm = [chi_de(k[2], res[k]['mitad_M'], k[0], k[1]) for k in con if 'mitad_M' in res[k]]
    print(f'  entre las {len(con)} familias: sigma_L/L_0 de la mitad de {min(xs):.2f} a {max(xs):.2f} % (factor {max(xs)/min(xs):.1f}); '
          f'chi de {min(cs):.2f} a {max(cs):.2f} (factor {max(cs)/min(cs):.2f}); modelo, chi de {min(cm):.2f} a {max(cm):.2f}')
    dif = []                                    # m_c por capa frente a m_c, en las corridas con chi <= 0.26 y con 0.26 < chi <= 0.36
    for (pot, sj, m), F in familias_L().items():
        for f in F:
            if f['nl'] > 1 and f['mcc'] is not None:
                dif.append((chi_de(m, f['sl'], pot, sj), f['mcc']/f['P']['mc'] - 1))
    if dif:
        a = [d for c, d in dif if c <= 0.26]; b = [d for c, d in dif if 0.26 < c <= 0.36]; c5 = [d for c, d in dif if c > 0.36]
        print(f'  m_c por capa frente a m_c ({len(dif)} corridas): con chi <= 0.26 difiere entre {100*min(a):+.0f} y {100*max(a):+.0f} %; '
              f'con chi = 0.35, entre {100*min(b):+.0f} y {100*max(b):+.0f} %; con chi = 0.5, entre {100*min(c5):+.0f} y {100*max(c5):+.0f} %')
    print('  última corrida con núcleo (exceso > 0.1 del de L igual) y primera sin él, en chi: '
          + '; '.join(f'{k[2]:g}: {max([x for x, e in zip(res[k]["P"][0], res[k]["P"][1]) if e > 0.1]):.2f} y '
                      f'{min([x for x, e in zip(res[k]["P"][0], res[k]["P"][1]) if e <= 0.1], default=float("nan")):.2f}' for k in ref + otras))


def tabla_chi():
    """Filas de la tabla de la Sección 8D del artículo (docs/articulo/tab_chi.tex en vlasov-poisson_PIC, macro
    \\TabChi): para cada familia con núcleo, Gamma, la anchura sigma_L/L_0 en que el exceso de m_c cae a la mitad y
    su chi, del PIC y del modelo, y el cociente con la rotación del núcleo (B2) en la mitad del PIC."""
    res = frontera(verboso=False)
    orden = sorted(res, key=lambda k: (k[0] != 'isocrono', k[1] != SJ, k[1], k[2]))
    filas, antes = [], None
    for k in orden:
        pot, sj, m = k; r = res[k]; kc = constantes(pot, sj)
        w = kc['g']*r['mc0']*m/np.sqrt(0.1*sj)
        x = lambda s: f'${100*s/np.sqrt(2)/L0:.2f}$'
        c = lambda s: f'${chi_de(m, s, pot, sj):.2f}$'
        if antes is not None and k[:2] != antes:
            filas.append('\\hline')
        antes = k[:2]
        filas.append(f'{"isochrone" if pot == "isocrono" else "Kepler"} & ${sj:g}$ & ${m:g}$ & ${r["gamma"]:.2f}$ & '
                     f'{x(r["mitad_P"])} & {x(r["mitad_M"])} & {c(r["mitad_P"])} & {c(r["mitad_M"])} & '
                     f'${kc["dOdL"]*r["mitad_P"]/np.sqrt(2)/w:.2f}$ \\\\')
    sal = os.path.join(FIJO, 'docs', 'articulo', 'tab_chi.tex')
    open(sal, 'w').write('% Escrito por VlasovPoisson_PIC_sp/reproducir/scripts/nube_L.py articulo; no editar a mano.\n'
                         '\\newcommand{\\TabChi}{%\n' + '%\n'.join(filas) + '%\n}\n')
    print('escrita', sal)


def articulo():
    """Tablas y figuras de la Sección 8 del artículo."""
    tabla_L(); figura(); tabla_chi(); figura_chi()


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
    ap.add_argument('paso', choices=['preparar', 'tiempo', 'correr', 'analizar', 'medir', 'medir_capas', 'medir_colas', 'tanda', 'lista',
                                     'prueba', 'modelo', 'consistencia', 'tabla', 'reticulas', 'figura', 'frontera', 'articulo',
                                     'fase', 'figura_fase', 'figura_chi'])
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
        sys.exit(sum(prueba(n) for n in (a.nombres or ['M030_s010v'])))
    elif a.paso == 'consistencia':
        for n in (a.nombres or ['C021']):
            consistencia(n) if n == 'C021' else consistencia_familia(n)
    elif a.paso == 'modelo':
        modelo(a.nombres)
    elif a.paso in ('tabla', 'reticulas', 'figura', 'frontera', 'articulo', 'figura_fase', 'figura_chi'):
        globals()[a.paso]()
    else:
        globals()[a.paso](a.solo)
