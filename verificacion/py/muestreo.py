# I6. Estados de muestreo aa_random y aa_halton (dftype = gauss), sin avanzar
# (Nt = 0), en una caja de (r, p, L) que contiene toda la distribucion.
#  (a) aa_random: masa total a0; todas las particulas con la misma masa
#      (f L constante); las medias pesadas por masa de L y de J coinciden
#      con las exactas dentro de 5 desviaciones del muestreo, y |h_k(0)| con
#      el del estado de cuadratura aa_quad dentro de 5/sqrt(N). Un factor
#      equivocado en la medida (el peso L, o dr dp = dQ dJ) desplazaria las
#      medias mucho mas que eso.
#        <L> = Int L^2 C(L) dL / Int L C(L) dL,   C = exp(-(L-l0)^2/sl^2)
#        <J> = 2 sr / sqrt(pi)                    (perfil J^2 exp(-J^2/sr^2))
#  (b) aa_random: la misma semilla repite la muestra bit a bit, otra semilla
#      da otra, y con seed = 0 la semilla usada queda en params_usados.par.
#  (c) aa_halton: cada nodo esta dentro de su celda, la masa es a0 y
#      |h_k(0)| coincide con el de aa_quad al 1 %.
# Uso: muestreo.py <VP_PIC> <dir_trabajo>
import numpy as np, h5py, subprocess, sys, os, re
exe, wd = sys.argv[1], sys.argv[2]
os.makedirs(wd, exist_ok=True); os.chdir(wd)
P = dict(rminc=1.0, rmaxc=10.0, pminc=-0.6, pmaxc=0.6, lminc=1.6, lmaxc=2.4, a0=1e-3, l0=2.0, sr=0.1, sl=0.2)
plantilla = """dr = 0.1
Nrc = {nr}
Npc = {np}
Nlc = {nl}
rmin = 0.0
rmax = 30.0
rminc = {rminc}
rmaxc = {rmaxc}
pminc = {pminc}
pmaxc = {pmaxc}
lminc = {lminc}
lmaxc = {lmaxc}
pmax = 2.0
Nt = 0
time_output = 1
spatial_output = 1
field_output = 1
output_format = hdf5
a0 = {a0}
l0 = {l0}
sr = {sr}
sp = 0.1
sl = {sl}
state = {state}
dftype = gauss
seed = {seed}
cutoff = 0.0
bsplineorder = 1
integrator = yoshida4
BGtype = Isochrone
autointeraction = .false.
directory = {d}
"""
def corre(d, state, nr, np_, nl, seed=0):
    open(d+'.par', 'w').write(plantilla.format(d=d, state=state, nr=nr, np=np_, nl=nl, seed=seed, **P))
    subprocess.run([exe, d+'.par'], stdout=open(d+'.log', 'w'), stderr=subprocess.STDOUT, check=True)
    g = h5py.File(d+'/vlasov_output.h5', 'r')['step_0000000000']
    return dict(r=g['r_part'][:], p=g['p_part'][:], L=g['l_part'][:], fl=g['fl'][:], hk=np.loadtxt(d+'/hk1.tl')[1:])
def accion(s):
    E = 0.5*s['p']**2+0.5*s['L']**2/s['r']**2-1/(1+np.sqrt(1+s['r']**2))
    return 1/np.sqrt(-2*E)-0.5*(s['L']+np.sqrt(s['L']**2+4))
ok = True
quad = corre('quad', 'aa_quad', 200, 40, 16)

# (a) y (b) aa_random
N = 20000
a = corre('azar_1', 'aa_random', 50, 40, 10, seed=20261004)
masa = 8*np.pi**2*a['fl'].sum()*(P['rmaxc']-P['rminc'])/50*(P['pmaxc']-P['pminc'])/40*(P['lmaxc']-P['lminc'])/10
iguales = np.ptp(a['fl'])/a['fl'].mean()
x = np.linspace(P['lminc'], P['lmaxc'], 200001); C = np.exp(-(x-P['l0'])**2/P['sl']**2)
Lex = np.trapezoid(x*x*C, x)/np.trapezoid(x*C, x); Jex = 2*P['sr']/np.sqrt(np.pi)
J = accion(a)
zL = (a['L'].mean()-Lex)/(a['L'].std()/np.sqrt(N)); zJ = (J.mean()-Jex)/(J.std()/np.sqrt(N))
eh = np.abs(a['hk']/quad['hk']-1).max()
print('  aa_random N=%d: masa/a0-1 = %.1e, dispersion de f L = %.1e' % (N, masa/P['a0']-1, iguales))
print('     <L> = %.6f (exacto %.6f, %+.2f desviaciones)   <J> = %.6f (exacto %.6f, %+.2f desviaciones)' % (a['L'].mean(), Lex, zL, J.mean(), Jex, zJ))
print('     max |h_k(0)/h_k(aa_quad) - 1| = %.2e  (5/sqrt(N) = %.2e)' % (eh, 5/np.sqrt(N)))
ok &= abs(masa/P['a0']-1) < 1e-12 and iguales < 1e-12 and abs(zL) < 5 and abs(zJ) < 5 and eh < 5/np.sqrt(N)
b = corre('azar_2', 'aa_random', 50, 40, 10, seed=20261004)
c = corre('azar_3', 'aa_random', 50, 40, 10, seed=7)
d = corre('azar_4', 'aa_random', 50, 40, 10, seed=0)
usada = int(re.search(r'^seed\s*=\s*(\d+)', open('azar_4/params_usados.par').read(), re.M).group(1))
e = corre('azar_5', 'aa_random', 50, 40, 10, seed=usada)
rep = all(a[k].tobytes() == b[k].tobytes() for k in ('r', 'p', 'L', 'fl'))
otra = a['r'].tobytes() != c['r'].tobytes()
reloj = usada != 0 and all(d[k].tobytes() == e[k].tobytes() for k in ('r', 'p', 'L', 'fl'))
print('     misma semilla: %s;  otra semilla: %s;  seed = 0 uso %d y repetirla da lo mismo: %s' % (
      'identica' if rep else 'DISTINTA', 'distinta' if otra else 'IDENTICA', usada, 'si' if reloj else 'NO'))
ok &= rep and otra and reloj

# (c) aa_halton
nr, np_, nl = 200, 40, 16
h = corre('halton', 'aa_halton', nr, np_, nl); g = corre('malla', 'aa', nr, np_, nl)
drc, dpc, dlc = (P['rmaxc']-P['rminc'])/nr, (P['pmaxc']-P['pminc'])/np_, (P['lmaxc']-P['lminc'])/nl
masa = 8*np.pi**2*h['fl'].sum()*drc*dpc*dlc
# Cada nodo se desplaza dentro de su celda: los indices de celda de los nodos
# que quedan (los no ligados se quitan) caen en la caja y no se repiten.
ir = np.floor((h['r']-P['rminc'])/drc).astype(int); ip = np.floor((h['p']-P['pminc'])/dpc).astype(int)
il = np.floor((h['L']-P['lminc'])/dlc).astype(int)
dentro = ir.min() >= 0 and ir.max() < nr and ip.min() >= 0 and ip.max() < np_ and il.min() >= 0 and il.max() < nl
unicos = len(np.unique(ir*np_*nl+ip*nl+il)) == len(ir)
ehh = np.abs(h['hk']/quad['hk']-1).max(); ehg = np.abs(g['hk']/quad['hk']-1).max()
print('  aa_halton %dx%dx%d: %d nodos con masa (aa: %d), uno por celda: %s, masa/a0-1 = %.1e' % (nr, np_, nl, len(h['r']), len(g['r']), 'si' if dentro and unicos else 'NO', masa/P['a0']-1))
print('     max |h_k(0)/h_k(aa_quad) - 1| = %.2e  (aa, sin desplazar: %.2e)' % (ehh, ehg))
ok &= abs(masa/P['a0']-1) < 1e-12 and dentro and unicos and ehh < 1e-2
print(('PASA' if ok else 'FALLA')+' I6 estados de muestreo aa_random y aa_halton')
