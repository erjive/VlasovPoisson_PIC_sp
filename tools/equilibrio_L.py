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



# ---------------------------------------------------------------------------
# El esquema discreto del código, para que el equilibrio lo sea del sistema
# que se va a integrar y no del continuo (AUDITORIA_2026-09-20.md, punto 21).
# Depósito con W_n y V_i (density.f90), cuadratura de la masa encerrada e
# integración cerrada del potencial (poisson_rk.f90), e interpolación de vuelta
# con el mismo W_n, incluidos el espejo en el origen y la cola exterior.
# ---------------------------------------------------------------------------

def Wn(n, y):
    a = np.abs(y)
    if n == 1:
        return np.where(a < 1.0, 1.0 - a, 0.0)
    if n == 2:
        return np.where(a < 0.5, 0.75 - y*y,
               np.where(a < 1.5, 0.125*(3.0 - 2.0*a)**2, 0.0))
    if n == 3:
        return np.where(a < 1.0, 2.0/3.0 - y*y + a**3/2.0,
               np.where(a < 2.0, (2.0 - a)**3/6.0, 0.0))
    raise ValueError('bsplineorder debe ser 1, 2 o 3')


class PoissonCodigo:
    """Phi_self tal como lo calcula el código, a partir de las partículas."""

    def __init__(self, dr, rmax, n):
        self.dr, self.n = dr, n
        self.Nr = int(rmax/dr) + 1
        self.r = (np.arange(1, self.Nr + 1) - 0.5)*dr
        self.vol = 4*np.pi*dr*(self.r**2 + (n + 1)*dr**2/12.0)
        # pesos de la cuadratura de masa: M(r_i) = M(r_{i-1}) + A_i rho_{i-1} + B_i rho_i
        self.A = np.zeros(self.Nr)
        self.B = np.zeros(self.Nr)
        self.B[0] = 4*np.pi*self.r[0]**3/3.0
        a, b = self.r[:-1], self.r[1:]
        I0 = (b**3 - a**3)/3.0
        I1 = (b**4 - a**4)/4.0 - a*I0
        self.A[1:] = 4*np.pi*(I0 - I1/dr)
        self.B[1:] = 4*np.pi*I1/dr

    def densidad(self, r_part, m_part):
        dr, n, r = self.dr, self.n, self.r
        cut = 0.5*(n + 1)*dr
        rho = np.zeros(self.Nr)
        ic = np.clip(np.round((r_part - r[0])/dr).astype(int), -2, self.Nr + 1)
        w = (n + 2)//2
        for d in range(-w - 1, w + 2):
            i = ic + d
            ok = (i >= 0) & (i < self.Nr)
            if not ok.any():
                continue
            ii, rj, mj = i[ok], r_part[ok], m_part[ok]
            peso = np.where(np.abs(r[ii] - rj) < cut, Wn(n, (r[ii] - rj)/dr), 0.0) \
                 + np.where(np.abs(r[ii] + rj) < cut, Wn(n, (r[ii] + rj)/dr), 0.0)
            np.add.at(rho, ii, mj*peso)
        return rho/self.vol

    def resolver(self, rho):
        r, dr = self.r, self.dr
        M = np.zeros(self.Nr)
        pot = np.zeros(self.Nr)
        M[0] = self.B[0]*rho[0]
        pot[0] = 2.0/3.0*np.pi*rho[0]*r[0]**2
        for i in range(1, self.Nr):
            ra, rb = r[i-1], r[i]
            s = (rho[i] - rho[i-1])/dr
            c3 = 4*np.pi*(rho[i-1] - s*ra)/3.0
            c4 = np.pi*s
            c0 = M[i-1] - c3*ra**3 - c4*ra**4
            M[i] = M[i-1] + self.A[i]*rho[i-1] + self.B[i]*rho[i]
            pot[i] = pot[i-1] + (-c0/rb + c3*rb**2/2.0 + c4*rb**3/3.0) \
                              - (-c0/ra + c3*ra**2/2.0 + c4*ra**3/3.0)
        force = -M/r**2
        pot = pot - (pot[-1] - force[-1]*r[-1])
        return pot, M

    def interpolar(self, pot, r_eval):
        """El potencial que siente una partícula en r_eval: la misma suma de
        pesos del código, con el espejo para los nodos j <= 0 y la cola
        kepleriana más allá del último."""
        dr, n, r, Nr = self.dr, self.n, self.r, self.Nr
        cut = 0.5*(n + 1)*dr
        w = (n + 2)//2
        jc = np.round((r_eval - r[0])/dr).astype(int) + 1
        out = np.zeros_like(np.asarray(r_eval, float))
        for d in range(-w, w + 1):
            j = jc + d
            rj = r[0] + (j - 1)*dr
            peso = np.where(np.abs(r_eval - rj) < cut, Wn(n, (r_eval - rj)/dr), 0.0)
            m = np.where(j <= 0, 1 - j, j)                 # espejo: Phi es par
            rm = r[0] + (m - 1)*dr
            dentro = m <= Nr
            pj = np.where(dentro, pot[np.clip(m, 1, Nr) - 1], pot[-1]*r[-1]/np.maximum(rm, 1e-300))
            out += pj*peso
        return out

    def phi(self, r_part, m_part, r_eval):
        pot, M = self.resolver(self.densidad(r_part, m_part))
        return self.interpolar(pot, r_eval), pot, M[-1]


class Equilibrio:
    def __init__(self, a0, sigma_j, J_max, L, dL, l0, sl, nrc, npc, dr, rmax, bsplineorder,
                 r_malla=None):
        self.a0, self.sigma_j, self.J_max = a0, sigma_j, J_max
        self.L, self.dL = np.asarray(L, float), dL
        self.C = np.exp(-(self.L - l0)**2/sl**2)
        Jq = np.linspace(0, J_max, 200001)
        self.A = a0/(16*np.pi**3*np.sum(self.L*self.C*dL)*np.trapezoid(self.forma(Jq), Jq))
        self.nrc, self.npc = nrc, npc
        self.dJ, self.dQ = J_max/nrc, 2*np.pi/npc
        self.pc = PoissonCodigo(dr, rmax, bsplineorder)
        if r_malla is None:
            r_malla = np.arange(0.01, rmax + 1e-9, 0.01)
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

    def nodos(self, m, eps=0.0):
        """Los nodos de cuadratura que va a ver el código, en el orden del código."""
        Jn = (np.arange(self.nrc) + 0.5)*self.dJ
        Qn = (np.arange(self.npc) + 0.5)*self.dQ
        LL, JJ, QQ = np.meshgrid(self.L, Jn, Qn, indexing='ij')   # orden k, i, j
        CC = np.meshgrid(self.C, Jn, Qn, indexing='ij')[0]
        LL, JJ, QQ, CC = LL.ravel(), JJ.ravel(), QQ.ravel(), CC.ravel()
        r, p = m.invertir(QQ, JJ, LL)
        F = self.A*self.forma(JJ)*CC*(1 + eps*np.cos(QQ))
        return r, p, LL, F, QQ, JJ

    def masas(self, F, LL):
        """La masa por partícula, normalizada a a0 como hace initial_data.f90."""
        mp = 8*np.pi**2*self.dJ*self.dQ*self.dL*F*LL
        return mp*(self.a0/mp.sum())

    def iterar(self, tol=1e-12, maxit=60, verboso=True):
        """Punto fijo contra el solver DEL CÓDIGO: en cada paso se colocan los
        mismos nodos que se van a escribir, se depositan en la malla del código
        y se resuelve Poisson como allí. El punto fijo es entonces una
        distribución que es función de las acciones del potencial que el código
        va a calcular a partir de esas mismas partículas, que es lo que hace
        falta para que no evolucione. Iterarlo contra la Poisson del continuo
        dejaba un desajuste de ~1e-4 relativo en Phi_self
        (AUDITORIA_2026-09-20.md, punto 21)."""
        for it in range(maxit):
            m = self.mapa()
            r_part, _, LL, F, _, _ = self.nodos(m)
            nuevo, pot_malla, M = self.pc.phi(r_part, self.masas(F, LL), self.r)
            cambio = np.max(np.abs(nuevo - self.phi_self))
            self.phi_self = nuevo
            if verboso:
                print(f'  iteración {it:2d}: max|dPhi_self| = {cambio:.2e}   masa = {M:.10e}', flush=True)
            if cambio < tol:
                break
        self.pot_malla, self.masa = pot_malla, M
        m = self.mapa()
        self.tablas = [m.tabla_J_de_E(Lk, 1.05*self.J_max) for Lk in self.L]
        self.rho = self.densidad(m, self.tablas)
        return self


def condicion_inicial(eq, eps):
    return eq.nodos(eq.mapa(), eps)


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
    ap.add_argument('--dr', type=float, default=0.05, help='dr de la malla del código')
    ap.add_argument('--rmax', type=float, default=25.0, help='rmax de la malla del código')
    ap.add_argument('--bsplineorder', type=int, default=1)
    ap.add_argument('--salida', required=True)
    a = ap.parse_args()
    dL = (a.lmaxc - a.lminc)/a.nlc
    L = a.lminc + (np.arange(a.nlc) + 0.5)*dL
    print(f'equilibrio: a0={a.a0:g}, F_eq = A J^2 exp(-J^2/{a.sigma_j}^2) exp(-(L-{a.l0})^2/{a.sl}^2),'
          f' L en {a.nlc} puntos medios de [{a.lminc},{a.lmaxc}]')
    print(f'            Poisson discreta del código: dr={a.dr:g}, rmax={a.rmax:g},'
          f' bsplineorder={a.bsplineorder}')
    eq = Equilibrio(a.a0, a.sigma_j, a.jmax, L, dL, a.l0, a.sl,
                    a.nrc, a.npc, a.dr, a.rmax, a.bsplineorder).iterar()
    r, p, LL, F, Q, J = condicion_inicial(eq, a.eps)
    Qc, Jc, _ = eq.mapa()(r, p, LL)
    peso = F > 1e-6*F.max()
    print(f'inversión: max|J - J_nodo| = {np.max(np.abs(Jc - J)[peso]):.1e},'
          f' max|Q - Q_nodo| = {np.abs(np.angle(np.exp(1j*(Qc - Q))))[peso].max():.1e} (nodos con peso)')
    os.makedirs(os.path.dirname(os.path.abspath(a.salida)), exist_ok=True)
    np.savetxt(a.salida, np.column_stack([r, p, LL, F]), fmt='%.17e')
    base = os.path.splitext(a.salida)[0]
    np.savez(base + '_equilibrio.npz', r=eq.r, phi_self=eq.phi_self, rho=eq.rho, A=eq.A, a0=a.a0,
             eps=a.eps, nrc=a.nrc, npc=a.npc, nlc=a.nlc, L=L, dL=dL, l0=a.l0, sl=a.sl,
             J_max=a.jmax, sigma_j=a.sigma_j, masa=eq.masa,
             r_codigo=eq.pc.r, phi_codigo=eq.pot_malla, dr=a.dr, rmax=a.rmax,
             bsplineorder=a.bsplineorder)
    print(f'escrito {a.salida} ({len(r)} partículas) y {base}_equilibrio.npz')
