# I2. Colapso frio de una esfera uniforme (R=1, M=1) con autogravedad y sin
# fondo. Antes del cruce de cascaras cada cascara sigue la cicloide
# r = r0 cos^2(th), t = (th + sin th cos th)/sqrt(2 M/R^3) del continuo.
# Se mide el error relativo mediano en r(t=0.8) para 0.1 < r0 < 0.9 (lejos
# del centro, donde pesa el termino centrifugo, y del borde discontinuo) con
# N = 200 y 400 cascaras. Un esquema consistente con el continuo converge
# mas rapido que 1/N; un sesgo O(1/N) (defecto D1: autogravedad de cada
# cascara quitada) da un error ~ 1/N.
# El paso lo elige el codigo: con L = 1e-4 lo fija la cota del pericentro de
# set_timestep (5.0e-5; sin ella seria courant dr/pmax = 2.5e-3), asi que el
# numero de pasos hasta t = 0.8 se calcula con el paso que informa una
# corrida de cero pasos.
# Uso: colapso.py <VP_PIC> <dir_trabajo>
import numpy as np, h5py, subprocess, sys, os, re
exe, wd = sys.argv[1], sys.argv[2]
os.makedirs(wd, exist_ok=True); os.chdir(wd)
plantilla = """dr = 0.01
Nrc = {N}
Npc = 1
Nlc = 1
rmin = 0.0
rmax = 1.5
rminc = 0.0
rmaxc = 1.0
pminc = -1.0
pmaxc = 1.0
lminc = 0.0
lmaxc = 2.0e-4
a0 = 1.0
state = checkpoint
CheckPointfile = ic_{N}.dat
courant = 0.5
pmax = 2.0
dt_switch = fix
Nt = {Nt}
time_output = 100000000
spatial_output = {so}
field_output = {so}
output_format = hdf5
bsplineorder = 1
integrator = yoshida4
BGtype = null
autointeraction = .true.
directory = colapso_{N}
"""
def theta(t):
    a, b = 0.0, np.pi/2
    g = lambda x: x+np.sin(x)*np.cos(x)-np.sqrt(2)*t
    for _ in range(200):
        c = 0.5*(a+b)
        if g(a)*g(c) <= 0: b = c
        else: a = c
    return 0.5*(a+b)
med = {}
for N in (200, 400):
    L = 1e-4
    with open('ic_%d.dat' % N, 'w') as fo:
        for k in range(1, N+1):
            a = (k-1)/N; b = k/N
            fo.write('%.17e %.17e %.17e %.17e\n' % (0.5*(a+b), 0.0, L, (b**3-a**3)/L))
    open('colapso_%d.par' % N, 'w').write(plantilla.format(N=N, Nt=0, so=1))
    out = subprocess.run([exe, 'colapso_%d.par' % N], capture_output=True, text=True, check=True).stdout
    dt = float(re.search(r'Time step fixed at size:\s*(\S+)', out).group(1))
    Nt = round(0.8/dt)
    open('colapso_%d.par' % N, 'w').write(plantilla.format(N=N, Nt=Nt, so=Nt))
    subprocess.run([exe, 'colapso_%d.par' % N], stdout=open('colapso_%d.log' % N, 'w'), stderr=subprocess.STDOUT, check=True)
    h = h5py.File('colapso_%d/vlasov_output.h5' % N, 'r')
    st = sorted(k for k in h if k.startswith('step_'))
    r0 = h[st[0]]['r_part'][:]; t = h[st[-1]].attrs['time']; r = h[st[-1]]['r_part'][:]
    rref = r0*np.cos(theta(t))**2
    sel = (r0 > 0.1) & (r0 < 0.9)
    e = np.abs(r[sel]-rref[sel])/rref[sel]
    med[N] = np.median(e)
    print('  N=%d  dt=%.2e  t=%.2f  error relativo en r: mediana %.2e, maximo %.2e' % (N, dt, t, med[N], e.max()))
orden = np.log2(med[200]/med[400])
ok = med[400] < 1e-4 and orden > 1.5
print('  orden en 1/N: %.2f' % orden)
print(('PASA' if ok else 'FALLA')+' I2 colapso frio: autogravedad frente al continuo')
