"""h_k exacto para el estado aa_quad sin autogravedad en el isócrono, y comparación con la
salida del código.

Sin autogravedad J y L se conservan y Q = Q0 + omega(J,L) t, con omega = 1/(J+c)^3 y
c = (L + sqrt(L^2+4))/2. Con la distribución F = A_f(Q) B_f(J) C_f(L) del estado aa_quad
y la función de prueba n del código (analysish.f90),

    h_k(t) = a0 a_k^t (a_k^f/a_0^f) Int L C_f C_t B_f B_t e^{-ik omega t} dJ dL
                                    / Int L C_f dL Int B_f dJ,

con A_f = exp(-sin^2(Q/2)/sp^2), B_f = J^2 exp(-J^2/sr^2), C_f = exp(-(L-l0)^2/sl^2),
A_t = exp(-sin^2(Q/2)/sq^2), B_t = J^2 exp(-(J-j)^2/sj^2), C_t = exp(-(L-lt)^2/slt^2),
a_k = (1/2pi) Int_0^2pi A e^{-ikQ} dQ, y las integrales sobre el soporte muestreado
[jminc,jmaxc] x [lminc,lmaxc].

Se calculan dos referencias:
  continuo  la integral, con Gauss-Legendre (--nj, --nl nodos);
  discreto  la misma suma sobre los nodos de punto medio que usa el código (Nrc en J,
            Npc en Q, Nlc en L), que aísla el error del integrador y del mapa.

Uso:  python3 tools/hk_exacto.py <directorio de la corrida> [--nj 4000] [--nl 400] [--fn 1]
"""
import argparse, os, sys
import numpy as np


def leer_par(fn):
    p = {}
    for linea in open(fn):
        linea = linea.split('#')[0].split('!')[0].strip()
        if '=' in linea:
            k, v = linea.split('=', 1)
            p[k.strip().lower()] = v.strip()
    return p


def a_k(sigma, kmax=4, n=4096):
    # Trapecio periódico: espectralmente exacto para este integrando analítico.
    q = 2*np.pi*np.arange(n)/n
    A = np.exp(-np.sin(q/2)**2/sigma**2)
    return np.array([np.mean(A*np.cos(k*q)) for k in range(kmax+1)])


def omega(J, L):
    return 1.0/(J + 0.5*(L + np.sqrt(L**2 + 4)))**3


def hk_continuo(t, P, fn, nj, nl, kmax=4):
    xj, wj = np.polynomial.legendre.leggauss(nj)
    xl, wl = np.polynomial.legendre.leggauss(nl)
    J = 0.5*(P['jmaxc']-P['jminc'])*xj + 0.5*(P['jmaxc']+P['jminc']); wJ = 0.5*(P['jmaxc']-P['jminc'])*wj
    L = 0.5*(P['lmaxc']-P['lminc'])*xl + 0.5*(P['lmaxc']+P['lminc']); wL = 0.5*(P['lmaxc']-P['lminc'])*wl
    return _hk(t, P, fn, J, wJ, L, wL, kmax, None)


def hk_discreto(t, P, fn, kmax=4):
    dJ = (P['jmaxc']-P['jminc'])/P['nrc']; dL = (P['lmaxc']-P['lminc'])/P['nlc']
    J = P['jminc'] + (np.arange(P['nrc'])+0.5)*dJ
    L = P['lminc'] + (np.arange(P['nlc'])+0.5)*dL
    return _hk(t, P, fn, J, np.full(J.size, dJ), L, np.full(L.size, dL), kmax, P['npc'])


def _hk(t, P, fn, J, wJ, L, wL, kmax, nq):
    s = str(fn)
    if nq is None:
        af, at = a_k(P['sp'], kmax), a_k(P['sq'+s], kmax)
    else:
        # Suma de trapecio en Q con los Npc nodos del código: a_k^f discreto; a_k^t es el
        # de Simpson 512 del código, igual al exacto a 1e-16.
        q = (np.arange(nq)+0.5)*2*np.pi/nq
        Af = np.exp(-np.sin(q/2)**2/P['sp']**2)
        af = np.array([np.mean(Af*np.exp(-1j*k*q)) for k in range(kmax+1)])
        at = a_k(P['sq'+s], kmax)
    Bf = J**2*np.exp(-J**2/P['sr']**2)
    Bt = J**2*np.exp(-(J-P['j'+s])**2/P['sj'+s]**2)
    Cf = np.exp(-(L-P['l0'])**2/P['sl']**2)
    Ct = np.exp(-(L-P['lt'+s])**2/P['slt'+s]**2)
    WJ = (Bf*Bt*wJ)[:, None]
    WL = (L*Cf*Ct*wL)[None, :]
    om = omega(J[:, None], L[None, :])
    den = np.sum(L*Cf*wL)*np.sum(Bf*wJ)*af[0]
    t = np.atleast_1d(t)
    h = np.empty((t.size, kmax+1), complex)
    for k in range(kmax+1):
        for i, ti in enumerate(t):
            h[i, k] = P['a0']*at[k]*af[k]*np.sum(WJ*WL*np.exp(-1j*k*om*ti))/den
    return h


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dir')
    ap.add_argument('--nj', type=int, default=4000)
    ap.add_argument('--nl', type=int, default=400)
    ap.add_argument('--fn', type=int, default=1, choices=[1, 2])
    a = ap.parse_args()

    raw = leer_par(os.path.join(a.dir, 'params_usados.par'))
    if raw['state'] != 'aa_quad' or raw['autointeraction'].lower() in ('.true.', 'true', 't'):
        sys.exit('solo para state=aa_quad sin autogravedad')
    if raw['bgtype'] != 'Isochrone':
        sys.exit('solo para BGtype=Isochrone')
    P = {}
    for k, v in raw.items():
        try:
            P[k] = float(v)
        except ValueError:
            pass
    for k in ('nrc', 'npc', 'nlc'):
        P[k] = int(P[k])
    if raw['cutoff'] not in ('0', '0.0'):
        print('aviso: cutoff distinto de cero; la referencia no incluye el corte')

    d = np.loadtxt(os.path.join(a.dir, f'hk{a.fn}_complex.tl'))
    t = d[:, 0]
    z = d[:, 1::2] + 1j*d[:, 2::2]
    hc = hk_continuo(t, P, a.fn, a.nj, a.nl)
    hd = hk_discreto(t, P, a.fn)

    np.savez(os.path.join(a.dir, f'hk{a.fn}_exacto.npz'), t=t, codigo=z, continuo=hc, discreto=hd)
    print(f'{a.dir}: N_J={P["nrc"]} N_Q={P["npc"]} N_L={P["nlc"]}, t hasta {t[-1]:g}')
    print('error relativo máximo por modo, |dh_k|/max_t|h_k(t)|:')
    esc = np.abs(hc).max(0)
    print('  código   - continuo :', ' '.join('%.2e' % x for x in (np.abs(z-hc)).max(0)/esc))
    print('  discreto - continuo :', ' '.join('%.2e' % x for x in (np.abs(hd-hc)).max(0)/esc))
    print('  código   - discreto :', ' '.join('%.2e' % x for x in (np.abs(z-hd)).max(0)/esc))


if __name__ == '__main__':
    main()
