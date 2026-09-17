"""Variables ángulo-acción numéricas con L por partícula, en un potencial esférico cualquiera.

El potencial es el isócrono (masa y escala unitarias) más una corrección tabulada en una
malla radial, por ejemplo el potencial propio de una corrida autogravitante, interpolada
con un spline cúbico natural y continuada fuera de la tabla como -M/r. Para cada
partícula (r, p_r, L):

    E  = p_r^2/2 + Phi_ef(r),                 Phi_ef = Phi + L^2/(2 r^2)
    r_c : mínimo de Phi_ef (bisección sobre dPhi_ef/dr, con la derivada del spline)
    r_-, r_+ : raíces de Phi_ef(r) = E a cada lado de r_c (bisección)
    J  = (1/pi) Int_{r_-}^{r_+} sqrt(2(E - Phi_ef)) dr
    T  = 2 Int_{r_-}^{r_+} dr / sqrt(2(E - Phi_ef))
    Q  = (2 pi/T) Int_{r_-}^{r} dr'/sqrt(2(E - Phi_ef))   (p_r >= 0),  2 pi - eso  (p_r < 0)

Las integrales usan r = rm + ra sin(theta), que cancela la singularidad de raíz inversa en
los puntos de retorno, con Gauss-Legendre. Generaliza reproducir/scripts/aa_numerico.py de
vlasov-poisson_PIC (L fija).

Uso como módulo:
    m = MapaAA(r_malla, phi_self)      # o MapaAA() para el isócrono solo
    Q, J, E = m(r, p, L)
    r, p = m.invertir(Q, J, L)         # inversa, por tablas J(E) para cada L distinto
"""
import numpy as np

R_MIN, R_MAX = 1e-3, 80.0


def phi_iso(r):
    return -1.0/(1.0 + np.sqrt(1.0 + r**2))


def dphi_iso(r):
    s = np.sqrt(1.0 + r**2)
    return r/(s*(1.0 + s)**2)


class SplineCubico:
    """Spline cúbico natural en una malla uniforme (sin scipy), con su derivada."""

    def __init__(self, x, y):
        x = np.asarray(x, float); y = np.asarray(y, float)
        n = len(x); h = x[1] - x[0]
        assert np.allclose(np.diff(x), h), 'la malla debe ser uniforme'
        a = np.full(n-2, 1.0); b = np.full(n-2, 4.0); c = np.full(n-2, 1.0)
        d = 6.0*(y[2:] - 2*y[1:-1] + y[:-2])/h**2
        for i in range(1, n-2):                      # Thomas
            m = a[i]/b[i-1]
            b[i] -= m*c[i-1]; d[i] -= m*d[i-1]
        M = np.zeros(n)
        M[n-2] = d[-1]/b[-1]
        for i in range(n-4, -1, -1):
            M[i+1] = (d[i] - c[i]*M[i+2])/b[i]
        self.x0, self.h, self.y, self.M, self.n = x[0], h, y, M, n

    def _base(self, xq):
        xq = np.asarray(xq, float)
        i = np.clip(((xq - self.x0)//self.h).astype(int), 0, self.n-2)
        xi = self.x0 + i*self.h
        A = (xi + self.h - xq)/self.h
        return i, A, 1.0 - A

    def __call__(self, xq):
        i, A, B = self._base(xq)
        return (A*self.y[i] + B*self.y[i+1]
                + ((A**3 - A)*self.M[i] + (B**3 - B)*self.M[i+1])*self.h**2/6.0)

    def deriv(self, xq):
        i, A, B = self._base(xq)
        return ((self.y[i+1] - self.y[i])/self.h
                + (-(3*A**2 - 1)*self.M[i] + (3*B**2 - 1)*self.M[i+1])*self.h/6.0)


class MapaAA:
    def __init__(self, r_malla=None, phi_self=None, nodos=64):
        self.tabla = None
        if phi_self is not None:
            self.tabla = SplineCubico(r_malla, phi_self)
            self.r0_t, self.r_ult = float(r_malla[0]), float(r_malla[-1])
            # Fuera de la tabla el potencial propio es kepleriano: -M/r.
            self.M_self = -float(phi_self[-1])*self.r_ult
        self.x, self.w = np.polynomial.legendre.leggauss(nodos)

    # --- potencial -------------------------------------------------------------
    def phi(self, r):
        r = np.asarray(r, float)
        p = phi_iso(r)
        if self.tabla is None:
            return p
        dentro = self.tabla(np.clip(r, self.r0_t, self.r_ult))
        return p + np.where(r <= self.r_ult, dentro, -self.M_self/np.maximum(r, 1e-300))

    def dphi(self, r):
        r = np.asarray(r, float)
        d = dphi_iso(r)
        if self.tabla is None:
            return d
        dentro = self.tabla.deriv(np.clip(r, self.r0_t, self.r_ult))
        return d + np.where(r <= self.r_ult, dentro, self.M_self/np.maximum(r, 1e-300)**2)

    def phi_ef(self, r, L):
        return self.phi(r) + 0.5*L**2/r**2

    # --- órbitas ---------------------------------------------------------------
    def r_circular(self, L, it=100):
        """Mínimo de Phi_ef: dPhi/dr = L^2/r^3, por bisección."""
        L = np.asarray(L, float)
        a = np.full_like(L, R_MIN); b = np.full_like(L, R_MAX)
        for _ in range(it):
            m = 0.5*(a + b)
            neg = self.dphi(m) - L**2/m**3 < 0
            a = np.where(neg, m, a); b = np.where(neg, b, m)
        return 0.5*(a + b)

    def _raiz(self, E, L, a, b, it=80):
        """Bisección vectorizada de Phi_ef(r) = E en [a, b] (un cambio de signo)."""
        a = np.array(a, float, copy=True); b = np.array(b, float, copy=True)
        ga = self.phi_ef(a, L) - E
        for _ in range(it):
            m = 0.5*(a + b)
            gm = self.phi_ef(m, L) - E
            izq = np.sign(gm) == np.sign(ga)
            a = np.where(izq, m, a); ga = np.where(izq, gm, ga)
            b = np.where(izq, b, m)
        return 0.5*(a + b)

    def _integrales(self, E, L, rm, ra, th_hi):
        """Int_{-pi/2}^{th_hi} de sqrt(2(E-Phi_ef)) ra cos y de ra cos/sqrt(2(E-Phi_ef))."""
        x, w = self.x, self.w
        lo = -0.5*np.pi
        mid, half = 0.5*(th_hi + lo), 0.5*(th_hi - lo)
        th = mid[:, None] + half[:, None]*x[None, :]
        rr = rm[:, None] + ra[:, None]*np.sin(th)
        v = np.sqrt(np.maximum(2.0*(E[:, None] - self.phi_ef(rr, L[:, None])), 0.0))
        jac = ra[:, None]*np.cos(th)
        cociente = np.where(v > 0, jac/np.where(v > 0, v, 1.0), 0.0)
        return half*np.sum(w*v*jac, axis=1), half*np.sum(w*cociente, axis=1)

    def _extremos(self, E, L):
        rc = self.r_circular(L)
        r1 = self._raiz(E, L, np.full_like(E, R_MIN), rc)
        r2 = self._raiz(E, L, rc, np.full_like(E, R_MAX))
        return 0.5*(r1 + r2), 0.5*(r2 - r1)

    def __call__(self, r, p, L, lote=4000):
        r, p = np.asarray(r, float), np.asarray(p, float)
        L = np.broadcast_to(np.asarray(L, float), r.shape).copy()
        Q = np.empty_like(r); J = np.empty_like(r); E = np.empty_like(r)
        for s in range(0, r.size, lote):
            sl = slice(s, s+lote)
            Q[sl], J[sl], E[sl] = self._lote(r[sl], p[sl], L[sl])
        return Q, J, E

    def _lote(self, r, p, L):
        E = 0.5*p**2 + self.phi_ef(r, L)
        rm, ra = self._extremos(E, L)
        Ip, It = self._integrales(E, L, rm, ra, np.full_like(r, 0.5*np.pi))
        J = Ip/np.pi
        T = 2.0*It
        s = np.clip((r - rm)/np.where(ra > 0, ra, 1.0), -1.0, 1.0)
        _, t_r = self._integrales(E, L, rm, ra, np.arcsin(s))
        frac = np.where(T > 0, t_r/np.where(T > 0, T, 1.0), 0.0)     # órbita circular: Q indefinido
        Q = np.where(p >= 0, 2.0*np.pi*frac, 2.0*np.pi - 2.0*np.pi*frac)
        return np.mod(Q, 2.0*np.pi), J, E

    # --- inversa ---------------------------------------------------------------
    def tabla_J_de_E(self, L, J_max, n=4000):
        """J(E) a L fijo, sobre órbitas que pasan por r_c con distintos p_r, hasta cubrir J_max."""
        rc = float(self.r_circular(np.array([L]))[0])
        E0 = float(self.phi_ef(rc, L))
        c = 0.5*(L + np.sqrt(L**2 + 4))
        Emax = max(-0.5/(J_max + c)**2, E0 + 1e-12)      # isócrono como punto de partida
        for _ in range(60):
            _, Jt, _ = self(np.array([rc]), np.array([np.sqrt(2*(Emax - E0))]), L)
            if Jt[0] >= J_max:
                break
            Emax = 0.5*Emax if Emax < 0 else Emax
            Emax = Emax + 0.1*abs(E0)*(Emax >= 0)
        E = np.linspace(E0, Emax, n)
        _, J, _ = self(np.full_like(E, rc), np.sqrt(2*np.maximum(E - E0, 0)), L)
        J[0] = 0.0
        return E, J

    def invertir(self, Q, J, L, nb=60, lote=4000):
        """(Q,J,L) -> (r,p_r): E de la tabla J(E) de cada L distinto, y bisección en theta."""
        Q, J = np.asarray(Q, float), np.asarray(J, float)
        L = np.broadcast_to(np.asarray(L, float), Q.shape)
        r = np.empty_like(Q); p = np.empty_like(Q)
        for Lv in np.unique(L):
            sel = np.nonzero(L == Lv)[0]
            E_t, J_t = self.tabla_J_de_E(Lv, 1.05*J[sel].max())
            for s in range(0, sel.size, lote):
                ix = sel[s:s+lote]
                r[ix], p[ix] = self._invertir_lote(Q[ix], J[ix], np.interp(J[ix], J_t, E_t),
                                                   np.full(ix.size, Lv), nb)
        return r, p

    def _invertir_lote(self, Q, J, E, L, nb, newton=3):
        # La tabla da E(J) con interpolación lineal (~1e-8); se corrige con Newton,
        # dE/dJ = omega = 2 pi/T.
        for _ in range(newton):
            rm, ra = self._extremos(E, L)
            Ip, It = self._integrales(E, L, rm, ra, np.full_like(E, 0.5*np.pi))
            T = 2*It
            E = E + (J - Ip/np.pi)*np.where(T > 0, 2*np.pi/np.where(T > 0, T, 1.0), 0.0)
        rm, ra = self._extremos(E, L)
        _, It = self._integrales(E, L, rm, ra, np.full_like(E, 0.5*np.pi))
        T = 2*It
        Qm = np.mod(Q, 2*np.pi)
        ida = Qm <= np.pi
        objetivo = np.where(ida, Qm, 2*np.pi - Qm)*T/(2*np.pi)
        a = np.full_like(E, -0.5*np.pi); b = np.full_like(E, 0.5*np.pi)
        for _ in range(nb):
            c = 0.5*(a + b)
            bajo = self._integrales(E, L, rm, ra, c)[1] < objetivo
            a = np.where(bajo, c, a); b = np.where(bajo, b, c)
        r = rm + ra*np.sin(0.5*(a + b))
        p = np.sqrt(np.maximum(2*(E - self.phi_ef(r, L)), 0))
        return r, np.where(ida, p, -p)


def aa_isocrono(r, p, L):
    """Mapa analítico del isócrono, el mismo de analysish.f90 (para validar)."""
    E = phi_iso(r) + 0.5*L**2/r**2 + 0.5*p**2
    er1 = np.sqrt((1+E*(2+L**2)-np.sqrt(1+2*E*(2+2*E+L**2)))/(2*E**2))
    er2 = np.sqrt((1+E*(2+L**2)+np.sqrt(1+2*E*(2+2*E+L**2)))/(2*E**2))
    s1 = 1+np.sqrt(1+er1**2); s2 = 1+np.sqrt(1+er2**2); s = 1+np.sqrt(1+r**2)
    arg = np.clip((s1+s2-2*s)/(s2-s1), -1, 1)
    eta = np.where(p >= 0, np.arccos(arg), np.arccos(-arg)+np.pi)
    Q = eta - np.sqrt((-2*E)**3)*np.sqrt(np.maximum(-L**2-2*E-2-0.5/E, 0))/(-2*E)*np.sin(eta)
    J = 1/np.sqrt(-2*E) - 0.5*(L+np.sqrt(L**2+4))
    return np.mod(Q, 2*np.pi), J, E
