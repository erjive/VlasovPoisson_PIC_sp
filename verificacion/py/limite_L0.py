# I4. Limite L -> 0 en el isocrono (r0 = 1.5, p0 = 0, t = 50).
#  (a) L = 0: orbita radial que cruza el centro (reflexion de main.f90)
#      frente a RK4 en x = +-r; debe converger con orden 2 (leapfrog).
#  (b) L = 1e-3 con yoshida4: el codigo debe elegir un paso que resuelva el
#      pericentro, dt = courant r_p^2/L (cota de set_timestep), y con
#      courant = 0.25 la energia debe conservarse fuera del pericentro
#      (r > 0.1) a mejor que 1e-6. Sin la cota, con dt = 0.01, el error era
#      197 (defecto N1).
# Las dos partes usan archivos de particulas distintos: el paso es comun a
# todas las particulas, y la de L = 1e-3 cambiaria el de la orbita radial.
# Uso: limite_L0.py <VP_PIC> <dir_trabajo>
import numpy as np, h5py, subprocess, sys, os, re
exe, wd = sys.argv[1], sys.argv[2]
os.makedirs(wd, exist_ok=True); os.chdir(wd)
with open('ic_a.dat', 'w') as fo:
    for L in (0.0, 1.0):
        fo.write('%.17e %.17e %.17e %.17e\n' % (1.5, 0.0, L, 1.0))
Lb = 1e-3
with open('ic_b.dat', 'w') as fo:
    fo.write('%.17e %.17e %.17e %.17e\n' % (1.5, 0.0, Lb, 1.0))
plantilla = """dr = 1.0
Nrc = {n}
Npc = 1
Nlc = 1
rmin = 0.0
rmax = 50.0
rminc = 0.0
rmaxc = 1.0
pminc = 0.0
pmaxc = 1.0
lminc = 0.0
lmaxc = 1.0
a0 = 1.0
state = checkpoint
CheckPointfile = {ic}
courant = {c}
pmax = {pm}
Nt = {nt}
time_output = 100000000
spatial_output = {so}
field_output = {so}
output_format = hdf5
bsplineorder = 1
integrator = {integ}
BGtype = Isochrone
autointeraction = .false.
directory = {d}
"""
def lanza(d, **kw):
    open(d+'.par', 'w').write(plantilla.format(d=d, **kw))
    out = subprocess.run([exe, d+'.par'], capture_output=True, text=True, check=True).stdout
    open(d+'.log', 'w').write(out)
    return float(re.search(r'Time step fixed at size:\s*(\S+)', out).group(1))
def corre(integ, dt):
    nt = round(50/dt); d = '%s_%g' % (integ, dt)
    dtc = lanza(d, n=2, ic='ic_a.dat', c=0.5, pm=0.5/dt, nt=nt, so=nt//500, integ=integ)
    assert abs(dtc-dt) < 1e-12*dt, 'el codigo no uso el paso pedido'
    return h5py.File(d+'/vlasov_output.h5', 'r')
def acc(x):
    s = np.sqrt(1+x*x); return -x/(s*(1+s)**2)
x, v, h = 1.5, 0.0, 5e-4
for _ in range(round(50/h)):
    k1x, k1v = v, acc(x); k2x, k2v = v+h/2*k1v, acc(x+h/2*k1x); k3x, k3v = v+h/2*k2v, acc(x+h/2*k2x); k4x, k4v = v+h*k3v, acc(x+h*k3x)
    x, v = x+h/6*(k1x+2*k2x+2*k3x+k4x), v+h/6*(k1v+2*k2v+2*k3v+k4v)
e = []
for dt in (0.02, 0.01):
    f = corre('leapfrog', dt); k = sorted(y for y in f if y.startswith('step_'))[-1]
    e.append(abs(f[k]['r_part'][0]-abs(x)))
orden = np.log2(e[0]/e[1])
print('  L=0, leapfrog: |r-r_ref| = %.2e, %.2e  orden %.2f' % (e[0], e[1], orden))
# (b) pericentro exacto de la orbita: L^2/(2 r^2) + Phi(r) = E, por biseccion.
Phi = lambda r: -1/(1+np.sqrt(1+r*r))
E0 = 0.5*Lb**2/1.5**2+Phi(1.5); lo, hi = 1e-9, 1.5
for _ in range(200):
    m = np.sqrt(lo*hi)
    if 0.5*Lb**2/m**2+Phi(m) > E0: lo = m
    else: hi = m
cour = 0.25; dt_esp = cour*lo**2/Lb
kw = dict(n=1, ic='ic_b.dat', c=cour, pm=cour/0.01, integ='yoshida4')
dt = lanza('yoshida4_L1e-3', nt=0, so=1, **kw)                      # el paso que elige el codigo
nt = round(50/dt); so = nt//500
lanza('yoshida4_L1e-3', nt=(nt//so)*so, so=so, **kw)
f = h5py.File('yoshida4_L1e-3/vlasov_output.h5', 'r'); st = sorted(y for y in f if y.startswith('step_'))
r = np.array([f[k]['r_part'][0] for k in st]); p = np.array([f[k]['p_part'][0] for k in st])
E = 0.5*p**2+0.5*Lb**2/r**2+Phi(r)
dE = (np.abs(E-E[0])/abs(E[0]))[r > 0.1].max()
print('  L=1e-3, yoshida4, courant=%g: dt = %.6e (courant r_p^2/L = %.6e, sin la cota 1.0e-02), max|dE/E0| en r > 0.1 = %.2e' % (cour, dt, dt_esp, dE))
print(('PASA' if abs(orden-2) < 0.15 else 'FALLA')+' I4a orbita radial L=0 por el centro')
print(('PASA' if abs(dt/dt_esp-1) < 1e-6 and dE < 1e-6 else 'FALLA')+' I4b paso y energia con L pequeno')
