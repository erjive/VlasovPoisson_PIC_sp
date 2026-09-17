"""Equilibrio autoconsistente F_eq(J,L) con L en rejilla, y condición inicial para state=checkpoint.

Un equilibrio es una función de las acciones VERDADERAS, las del potencial total
Phi = Phi_iso + Phi_self[F_eq]. Como la acción depende del potencial y el potencial de la
distribución, se itera

    Phi_self -> para cada L_k: tabla J(E) -> F_eq(J,L_k) -> rho(r) -> Poisson -> Phi_self

con F_eq = A J^2 exp(-J^2/sigma_J^2) exp(-(L-l0)^2/sl^2) y L en los puntos medios de la
rejilla del código (Nlc celdas en [lminc,lmaxc]), la misma que verá la simulación. Cada
L_k es un problema de un grado de libertad en el mismo potencial total:

    4 pi r^2 rho(r) = 8 pi^2 Sum_k L_k dL Int F_eq(J(E(r,p,L_k)), L_k) dp.

La masa no depende del potencial (dr dp_r = dQ dJ), así que
M = 8 pi^2 Sum_k L_k dL C(L_k) 2 pi Int A J^2 exp(-J^2/sigma_J^2) dJ fija A de entrada.

Con el equilibrio convergido se colocan nodos de cuadratura en una rejilla regular de las
(Q,J,L) verdaderas (puntos medios, Nrc en [0,J_max], Npc en Q, Nlc en L) y se invierte el
mapa. La perturbación es F = F_eq (1 + eps cos Q), que no cambia la masa. El archivo tiene
"r p_r L F" en el orden del código: indx = (k-1) Nrc Npc + (i-1) Npc + j.

Uso:
    python3 tools/equilibrio_L.py --a0 1e-2 --eps 0.1 --nrc 200 --npc 20 --nlc 8 \\
        --lminc 1.6 --lmaxc 2.4 --l0 2 --sl 0.2 --salida ic.dat
Generaliza reproducir/scripts/equilibrio.py de vlasov-poisson_PIC (L fija).
"""
import os, sys, argparse, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aa_numerico_L import MapaAA, phi_iso


class Equilibrio:
    def __init__(self, a0, sigma_j, J_max, L, dL, l0, sl, r_malla=np.arange(0.01, 25.0 + 1e-9, 0.01)):
        self.a0, self.sigma_j, self.J_max = a0, sigma_j, J_max
        self.L, self.dL = np.asarray(L, float), dL
        self.C = np.exp(-(self.L - l0)**2/sl**2)
        Jq = np.linspace(0, J_max, 200001)
        self.A = a0/(16*np.pi**3*np.sum(self.L*self.C*dL)*np.trapezoid(self.forma(Jq), Jq))
        self.r = r_malla
        self.phi_self = np.zeros_like(r_malla)

    def forma(self, J):
        return J**2*np.exp(-J**2/self.sigma_j**2)

    def mapa(self):
        return MapaAA(self.r, self.phi_self)

    def densidad(self, m, tablas, npm=1201):
        u = np.linspace(-1, 1, npm)
        acum = np.zeros_like(self.r)
        for Lk, Ck, (E_t, J_t) in zip(self.L, self.C, tablas):
            phief = m.phi_ef(self.r, Lk)
            pmax = np.sqrt(2*np.maximum(E_t[-1] - phief, 0.0))
            E = 0.5*(pmax[:, None]*u[None, :])**2 + phief[:, None]
            J = np.interp(E, E_t, J_t, right=10*self.J_max)
            F = self.A*self.forma(J)*Ck
            acum += Lk*self.dL*np.trapezoid(F, u, axis=1)*pmax
        return 8*np.pi**2*acum/(4*np.pi*self.r**2)

    def poisson(self, rho):
        r = self.r
        dM = 4*np.pi*r**2*rho
        M = np.concatenate([[0], np.cumsum(0.5*(dM[1:] + dM[:-1])*np.diff(r))])
        g = 4*np.pi*r*rho
        ext = np.concatenate([np.cumsum((0.5*(g[1:] + g[:-1])*np.diff(r))[::-1])[::-1], [0]])
        return -M/r - ext, M[-1]

    def iterar(self, tol=1e-12, maxit=60, verboso=True):
        for it in range(maxit):
            m = self.mapa()
            tablas = [m.tabla_J_de_E(Lk, 1.05*self.J_max) for Lk in self.L]
            rho = self.densidad(m, tablas)
            nuevo, M = self.poisson(rho)
            cambio = np.max(np.abs(nuevo - self.phi_self))
            self.phi_self = nuevo
            if verboso:
                print(f'  iteración {it:2d}: max|dPhi_self| = {cambio:.2e}   masa = {M:.10e}', flush=True)
            if cambio < tol:
                break
        self.tablas, self.rho, self.masa = tablas, rho, M
        return self


def condicion_inicial(eq, eps, nrc, npc):
    m = eq.mapa()
    Jn = (np.arange(nrc) + 0.5)*eq.J_max/nrc
    Qn = (np.arange(npc) + 0.5)*2*np.pi/npc
    LL, JJ, QQ = np.meshgrid(eq.L, Jn, Qn, indexing='ij')      # orden k, i, j del código
    CC = np.meshgrid(eq.C, Jn, Qn, indexing='ij')[0]
    LL, JJ, QQ, CC = LL.ravel(), JJ.ravel(), QQ.ravel(), CC.ravel()
    r, p = m.invertir(QQ, JJ, LL)
    F = eq.A*eq.forma(JJ)*CC*(1 + eps*np.cos(QQ))
    return r, p, LL, F, QQ, JJ


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--a0', type=float, default=1e-2)
    ap.add_argument('--eps', type=float, default=0.1)
    ap.add_argument('--sigma-j', type=float, default=0.1)
    ap.add_argument('--jmax', type=float, default=0.6)
    ap.add_argument('--nrc', type=int, default=200)
    ap.add_argument('--npc', type=int, default=20)
    ap.add_argument('--nlc', type=int, default=8)
    ap.add_argument('--lminc', type=float, default=1.6)
    ap.add_argument('--lmaxc', type=float, default=2.4)
    ap.add_argument('--l0', type=float, default=2.0)
    ap.add_argument('--sl', type=float, default=0.2)
    ap.add_argument('--salida', required=True)
    a = ap.parse_args()
    dL = (a.lmaxc - a.lminc)/a.nlc
    L = a.lminc + (np.arange(a.nlc) + 0.5)*dL
    print(f'equilibrio: a0={a.a0:g}, F_eq = A J^2 exp(-J^2/{a.sigma_j}^2) exp(-(L-{a.l0})^2/{a.sl}^2),'
          f' L en {a.nlc} puntos medios de [{a.lminc},{a.lmaxc}]')
    eq = Equilibrio(a.a0, a.sigma_j, a.jmax, L, dL, a.l0, a.sl).iterar()
    r, p, LL, F, Q, J = condicion_inicial(eq, a.eps, a.nrc, a.npc)
    Qc, Jc, _ = eq.mapa()(r, p, LL)
    peso = F > 1e-6*F.max()
    print(f'inversión: max|J - J_nodo| = {np.max(np.abs(Jc - J)[peso]):.1e},'
          f' max|Q - Q_nodo| = {np.abs(np.angle(np.exp(1j*(Qc - Q))))[peso].max():.1e} (nodos con peso)')
    os.makedirs(os.path.dirname(os.path.abspath(a.salida)), exist_ok=True)
    np.savetxt(a.salida, np.column_stack([r, p, LL, F]), fmt='%.17e')
    base = os.path.splitext(a.salida)[0]
    np.savez(base + '_equilibrio.npz', r=eq.r, phi_self=eq.phi_self, rho=eq.rho, A=eq.A, a0=a.a0,
             eps=a.eps, nrc=a.nrc, npc=a.npc, nlc=a.nlc, L=L, dL=dL, l0=a.l0, sl=a.sl,
             J_max=a.jmax, sigma_j=a.sigma_j, masa=eq.masa)
    print(f'escrito {a.salida} ({len(r)} partículas) y {base}_equilibrio.npz')
