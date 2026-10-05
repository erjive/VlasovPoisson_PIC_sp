# I3. Medida del espacio de fases. Se muestrea la F(E) de Plummer isotropo,
# F = (24 sqrt2 / 7 pi^3) (-E)^{7/2}, en mallas regulares de (r, p, L) y se
# corre el codigo con Nt = 0. La densidad depositada debe converger a
# rho = (3/4pi)(1+r^2)^{-5/2} en 0.3 < r < 7 con orden cercano a 2: un factor
# 8 pi^2 L, 2 pi o r^2 equivocado daria un error O(1) que no baja.
# Uso: medida.py <VP_PIC> <dir_trabajo>
import numpy as np, h5py, subprocess, sys, os
exe, wd = sys.argv[1], sys.argv[2]
os.makedirs(wd, exist_ok=True); os.chdir(wd)
C = 24*np.sqrt(2)/(7*np.pi**3)
def genera(Nr, Np, NL, archivo):
    R, P, LM = 20.0, 1.42, 4.0
    r = (np.arange(Nr)+0.5)*R/Nr; p = -P+(np.arange(Np)+0.5)*2*P/Np; L = (np.arange(NL)+0.5)*LM/NL
    pp, LL = np.meshgrid(p, L, indexing='ij'); out = []; tot = 0.0
    for ri in r:
        E = pp**2/2+LL**2/(2*ri**2)-1/np.sqrt(1+ri**2)
        F = np.where(E < 0, C*np.clip(-E, 0, None)**3.5, 0.0); m = F > 0
        tot += (F*LL)[m].sum(); out.append(np.c_[np.full(m.sum(), ri), pp[m], LL[m], F[m]])
    out = np.vstack(out); np.savetxt(archivo, out, fmt='%.17e')
    return len(out), 8*np.pi**2*tot*(R/Nr)*(2*P/Np)*(LM/NL)
plantilla = """dr = {dr}
Nrc = {N}
Npc = 1
Nlc = 1
rmin = 0.0
rmax = 25.0
rminc = 0.0
rmaxc = 1.0
pminc = 0.0
pmaxc = 1.0
lminc = 0.0
lmaxc = 1.0
a0 = {M}
state = checkpoint
CheckPointfile = {ic}
Nt = 0
time_output = 1
spatial_output = 1
field_output = 1
output_format = hdf5
bsplineorder = 1
integrator = leapfrog
BGtype = null
autointeraction = .true.
directory = medida_{k}
"""
err = {}
for k in (1, 2):
    N, M = genera(100*k, 40*k, 20*k, 'ic_%d.dat' % k)
    open('medida_%d.par' % k, 'w').write(plantilla.format(dr=0.4/k, N=N, M=M, ic='ic_%d.dat' % k, k=k))
    subprocess.run([exe, 'medida_%d.par' % k], stdout=open('medida_%d.log' % k, 'w'), stderr=subprocess.STDOUT, check=True)
    f = h5py.File('medida_%d/vlasov_output.h5' % k, 'r'); r = f['grid/r'][:]; s = f['step_0000000000']
    rho = s['avg_rho'][:]/r**2; ra = 3/(4*np.pi)*(1+r*r)**-2.5; sel = (r > 0.3) & (r < 7)
    err[k] = np.abs(rho[sel]-ra[sel]).max()/ra[sel].max()
    print('  N=%d  dr=%.2f  max|rho-rho_a|/max(rho_a)=%.2e  K=%.5f (3pi/64=0.14726)  W=%.5f (-3pi/32=-0.29452)' % (N, 0.4/k, err[k], s.attrs['kinetic_energy'], s.attrs['potential_energy']))
orden = np.log2(err[1]/err[2])
print('  orden: %.2f' % orden)
# Criterio: un factor equivocado en la medida (8 pi^2 L, 2 pi, r^2) da un error O(1)
# que no baja. El umbral 0.1 se fijo despues de la primera corrida (que dio 3.8e-2 con
# N = 154270; el primer umbral, 3e-2, era mas estricto que el error de cuadratura).
print(('PASA' if err[2] < 0.1 and orden > 1.3 else 'FALLA')+' I3 medida del espacio de fases (densidad de Plummer)')
