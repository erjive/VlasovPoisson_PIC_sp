# I4. Limite L -> 0 con paso fijo en el isocrono (r0 = 1.5, p0 = 0, t = 50).
#  (a) L = 0: orbita radial que cruza el centro (reflexion de main.f90)
#      frente a RK4 en x = +-r; debe converger con orden 2 (leapfrog).
#  (b) L = 1e-3 con yoshida4, dt = 0.01: la energia debe conservarse (< 1e-6).
#      Hoy falla: el paso no resuelve el pericentro r_p ~ L (defecto N1).
# Uso: limite_L0.py <VP_PIC> <dir_trabajo>
import numpy as np, h5py, subprocess, sys, os
exe, wd = sys.argv[1], sys.argv[2]
os.makedirs(wd, exist_ok=True); os.chdir(wd)
with open('ic.dat', 'w') as fo:
    for L in (0.0, 1e-3, 1.0):
        fo.write('%.17e %.17e %.17e %.17e\n' % (1.5, 0.0, L, 1.0))
plantilla = """dr = 1.0
Nrc = 3
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
CheckPointfile = ic.dat
courant = 0.5
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
def corre(integ, dt):
    nt = round(50/dt); d = '%s_%g' % (integ, dt)
    open(d+'.par', 'w').write(plantilla.format(pm=0.5/dt, nt=nt, so=nt//500, integ=integ, d=d))
    subprocess.run([exe, d+'.par'], stdout=open(d+'.log', 'w'), stderr=subprocess.STDOUT, check=True)
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
f = corre('yoshida4', 0.01); st = sorted(y for y in f if y.startswith('step_'))
Phi = lambda r: -1/(1+np.sqrt(1+r*r))
E = np.array([0.5*f[k]['p_part'][1]**2+0.5*f[k]['l_part'][1]**2/f[k]['r_part'][1]**2+Phi(f[k]['r_part'][1]) for k in st])
dE = np.abs(E-E[0]).max()/abs(E[0])
print('  L=1e-3, yoshida4, dt=0.01: max|dE/E0| = %.2e' % dE)
print(('PASA' if abs(orden-2) < 0.15 else 'FALLA')+' I4a orbita radial L=0 por el centro')
print(('PASA' if dE < 1e-6 else 'FALLA')+' I4b energia con L pequeno a paso fijo')
