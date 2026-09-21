# Referencia de U3: integral exacta (50 digitos) del interpolante lineal de
# los valores nodales que uso poisson_rk, y comparacion con su M y Phi.
import mpmath as mp, sys
mp.mp.dps = 50
rows = [list(map(mp.mpf, l.split())) for l in open(sys.argv[1])]
r = [x[0] for x in rows]; rho = [x[1] for x in rows]; pot = [x[2] for x in rows]; F = [x[3] for x in rows]
Nr = len(r)
M = [mp.mpf(0)]*Nr; P = [mp.mpf(0)]*Nr
M[0] = 4*mp.pi/3*rho[0]*r[0]**3; P[0] = 2*mp.pi/3*rho[0]*r[0]**2
for i in range(1, Nr):
    a, b = r[i-1], r[i]; s = (rho[i]-rho[i-1])/(b-a); al = rho[i-1]-s*a
    c0 = M[i-1]-4*mp.pi*(al*a**3/3+s*a**4/4)
    M[i] = c0+4*mp.pi*(al*b**3/3+s*b**4/4)
    P[i] = P[i-1]+c0*(1/a-1/b)+4*mp.pi*(al*(b**2-a**2)/6+s*(b**3-a**3)/12)
shift = P[-1]+M[-1]/r[-1]; P = [p-shift for p in P]
eF = max(abs(F[i]+M[i]/r[i]**2)/abs(M[i]/r[i]**2) for i in range(Nr))
eP = max(abs(pot[i]-P[i])/abs(P[i]) for i in range(Nr))
tol = float(sys.argv[2])
print('  Nr=%d  max rel |F-F_exacta|=%s  max rel |Phi-Phi_exacta|=%s  (tolerancia %g)' % (Nr, mp.nstr(eF, 3), mp.nstr(eP, 3), tol))
print(('PASA' if max(eF, eP) < tol else 'FALLA') + ' U3 integrador de Poisson exacto para su interpolante (Nr=%d)' % Nr)
