# I1. Orden de convergencia en dt de euler, leapfrog y yoshida4 frente al
# integrador "analytic" (exacto) en el isocrono sin autogravedad.
# Uso: integradores.py <VP_PIC> <dir_trabajo>
import numpy as np, h5py, subprocess, sys, os
exe, wd = sys.argv[1], sys.argv[2]
os.makedirs(wd, exist_ok=True); os.chdir(wd)
rng = np.random.default_rng(1); rows = []
Phi = lambda r: -1/(1+np.sqrt(1+r*r))
F = lambda r, L: -r/(np.sqrt(1+r*r)*(1+np.sqrt(1+r*r))**2)+L*L/r**3
while len(rows) < 60:
    L = rng.choice([0.3, 1.0, 2.0]); r = rng.uniform(0.3, 4); p = rng.uniform(-0.6, 0.6)
    if Phi(r)+L*L/(2*r*r)+p*p/2 < -0.02 and abs(F(r, L)) < 1.0:
        rows.append((r, p, L, 1.0))
np.savetxt('ic.dat', np.array(rows), fmt='%.17e')
plantilla = """dr = 1.0
Nrc = 60
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
dt_switch = fix
Nt = {nt}
time_output = 100000000
spatial_output = {so}
field_output = {nt}
output_format = hdf5
bsplineorder = 1
integrator = {integ}
BGtype = Isochrone
autointeraction = .false.
directory = {d}
"""
def corre(integ, dt):
    nt = round(20/dt); d = '%s_%g' % (integ, dt)
    open(d+'.par', 'w').write(plantilla.format(pm=0.5/dt, nt=nt, so=nt//20, integ=integ, d=d))
    subprocess.run([exe, d+'.par'], stdout=open(d+'.log', 'w'), stderr=subprocess.STDOUT, check=True)
    return h5py.File(d+'/vlasov_output.h5', 'r')
def ultimo(h): return h[sorted(k for k in h if k.startswith('step_'))[-1]]
dts = [0.05, 0.025, 0.0125]; esperado = {'euler': 1, 'leapfrog': 2, 'yoshida4': 4}
ok = True
for integ in esperado:
    err = []
    for dt in dts:
        a = ultimo(corre('analytic', dt)); b = ultimo(corre(integ, dt))
        err.append(np.abs(a['r_part'][:]-b['r_part'][:]).max())
    orden = np.log2(err[1]/err[2])
    bien = abs(orden-esperado[integ]) < 0.15*esperado[integ]
    ok &= bien
    print('  %-9s max|r-r_exacto| = %s   orden observado %.2f (esperado %d)' % (integ, ' '.join('%.2e' % e for e in err), orden, esperado[integ]))
print(('PASA' if ok else 'FALLA')+' I1 orden temporal de los integradores')
