"""h_k exacto para el estado aa_quad sin autogravedad en el isócrono, y comparación con la
salida del código.

Sin autogravedad J y L se conservan y Q = Q0 + omega(J,L) t, con omega = 1/(J+c)^3 y
c = (L + sqrt(L^2+4))/2. Con la distribución F0(Q,J,L) del estado aa_quad (dftype, la
misma de distribution.f90) y la función de prueba n del código (analysish.f90),

    h_k(t) = a0 a_k^t Int L C_t B_t Fk(J,L) e^{-ik omega t} dJ dL / Int L F0hat(J,L) dJ dL,
    Fk(J,L) = (1/2pi) Int_0^2pi F0(Q,J,L) e^{-ikQ} dQ,

con A_t = exp(-sin^2(Q/2)/sq^2), B_t = J^2 exp(-(J-j)^2/sj^2), C_t = exp(-(L-lt)^2/slt^2),
a_k = (1/2pi) Int A_t e^{-ikQ} dQ y las integrales sobre el soporte muestreado
[jminc,jmaxc] x [lminc,lmaxc].

Se calculan dos referencias:
  continuo  la integral: Gauss-Legendre en J (--nj) y L (--nl), trapecio periódico en Q
            (--nq); para king, J solo hasta el corte J_t, donde F0 tiene una esquina;
  discreto  la misma suma sobre los nodos de punto medio que usa el código (Nrc en J,
            Npc en Q, Nlc en L), que aísla el error del integrador y del mapa.

Uso:  python3 tools/hk_exacto.py <directorio de la corrida> [--nj 1000] [--nl 100] [--nq 512] [--fn 1]
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


# Constantes de distribution.f90
BIM = dict(b1=0.8, b2=0.3, Ja=0.10, sa=0.025, Jb=0.24, sb=0.035, wb=0.7)
SPI = dict(J0=0.15, s0=0.04, beta=50.0, sq=0.5)
KIN = dict(Jt=0.35, sE2=0.01, eps=0.5)


def F0(Q, J, L, P, dftype):
    CL = np.exp(-(L - P['l0'])**2/P['sl']**2)
    if dftype == 'gauss':
        return np.exp(-np.sin(0.5*Q)**2/P['sp']**2)*np.exp(-J**2/P['sr']**2)*J**2*CL
    if dftype == 'bimodal':
        b = BIM
        return ((1 + b['b1']*np.cos(Q) + b['b2']*np.cos(2*Q))
                *(np.exp(-(J-b['Ja'])**2/b['sa']**2) + b['wb']*np.exp(-(J-b['Jb'])**2/b['sb']**2))*CL)
    if dftype == 'spiral':
        s = SPI
        return np.exp(-(J-s['J0'])**2/s['s0']**2)*np.exp(-np.sin(0.5*(Q-s['beta']*J))**2/s['sq']**2)*CL
    if dftype == 'king':
        k = KIN
        E = lambda JJ: -1.0/(2*(JJ + 0.5*(L + np.sqrt(L**2 + 4)))**2)
        feq = np.where(J < k['Jt'], np.exp(-E(J)/k['sE2']) - np.exp(-E(k['Jt'])/k['sE2']), 0.0)
        return feq*(1 + k['eps']*np.cos(Q))*CL
    raise SystemExit(f'dftype desconocido: {dftype}')


def hk_referencia(t, P, fn, dftype, J, wJ, L, wL, q, kmax=4):
    """h_k(t) con reglas de cuadratura dadas en J y L, y nodos equiespaciados q en Q."""
    s = str(fn)
    at = a_k(P['sq'+s], kmax)
    Bt = J**2*np.exp(-(J-P['j'+s])**2/P['sj'+s]**2)
    num = np.zeros((np.atleast_1d(t).size, kmax+1), complex)
    den = 0.0
    t = np.atleast_1d(t)
    for Lk, wl in zip(L, wL):
        Ct = np.exp(-(Lk-P['lt'+s])**2/P['slt'+s]**2)
        Fq = F0(q[None, :], J[:, None], Lk, P, dftype)                     # (nJ, nQ)
        Fk = np.array([np.mean(Fq*np.exp(-1j*k*q[None, :]), axis=1) for k in range(kmax+1)])
        den += Lk*wl*np.sum(wJ*Fk[0].real)
        om = omega(J, Lk)
        for k in range(kmax+1):
            fase = np.exp(-1j*k*np.outer(t, om))                              # (nt, nJ)
            num[:, k] += Lk*wl*Ct*np.sum(wJ*Bt*Fk[k]*fase, axis=1)
    return P['a0']*at[None, :]*num/den


def hk_continuo(t, P, fn, dftype, nj, nl, nq):
    Jhi = min(P['jmaxc'], KIN['Jt']) if dftype == 'king' else P['jmaxc']
    xj, wj = np.polynomial.legendre.leggauss(nj)
    xl, wl = np.polynomial.legendre.leggauss(nl)
    J = 0.5*(Jhi-P['jminc'])*xj + 0.5*(Jhi+P['jminc']); wJ = 0.5*(Jhi-P['jminc'])*wj
    L = 0.5*(P['lmaxc']-P['lminc'])*xl + 0.5*(P['lmaxc']+P['lminc']); wL = 0.5*(P['lmaxc']-P['lminc'])*wl
    return hk_referencia(t, P, fn, dftype, J, wJ, L, wL, 2*np.pi*np.arange(nq)/nq)


def hk_discreto(t, P, fn, dftype):
    dJ = (P['jmaxc']-P['jminc'])/P['nrc']; dL = (P['lmaxc']-P['lminc'])/P['nlc']
    J = P['jminc'] + (np.arange(P['nrc'])+0.5)*dJ
    L = P['lminc'] + (np.arange(P['nlc'])+0.5)*dL
    q = (np.arange(P['npc'])+0.5)*2*np.pi/P['npc']
    return hk_referencia(t, P, fn, dftype, J, np.full(J.size, dJ), L, np.full(L.size, dL), q)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dir')
    ap.add_argument('--nj', type=int, default=1000)
    ap.add_argument('--nl', type=int, default=100)
    ap.add_argument('--nq', type=int, default=512)
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
    dftype = raw.get('dftype', 'gauss')
    hc = hk_continuo(t, P, a.fn, dftype, a.nj, a.nl, a.nq)
    hd = hk_discreto(t, P, a.fn, dftype)

    np.savez(os.path.join(a.dir, f'hk{a.fn}_exacto.npz'), t=t, codigo=z, continuo=hc, discreto=hd)
    print(f'{a.dir}: dftype={dftype} N_J={P["nrc"]} N_Q={P["npc"]} N_L={P["nlc"]}, t hasta {t[-1]:g}')
    # Modos que la distribución no tiene (h_k exacto = 0, p.ej. k=3,4 en bimodal) se miden
    # contra 1e-6 max|h_0| para que el cociente no divida por el redondeo.
    esc = np.maximum(np.abs(hc).max(0), 1e-6*np.abs(hc[:, 0]).max())
    print('error relativo máximo por modo, |dh_k|/max(max_t|h_k|, 1e-6 max_t|h_0|):')
    print('  código   - continuo :', ' '.join('%.2e' % x for x in (np.abs(z-hc)).max(0)/esc))
    print('  discreto - continuo :', ' '.join('%.2e' % x for x in (np.abs(hd-hc)).max(0)/esc))
    print('  código   - discreto :', ' '.join('%.2e' % x for x in (np.abs(z-hd)).max(0)/esc))


if __name__ == '__main__':
    main()
