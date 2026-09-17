"""Figuras de docs/introduccion/vlasov_L_intro.tex, e impresión de las cifras que cita el texto.

Uso (desde cualquier directorio):
    python3 reproducir/figuras/generar_figuras.py            # todas
    python3 reproducir/figuras/generar_figuras.py espiral    # solo las que contienen el nombre

Lee las corridas de exe/rep/ (reproducir/correr.sh) y los archivos que escriben
tools/hk_exacto.py, tools/hk_numerico.py y tools/delta_phi.py. Las figuras 'frecuencias',
'envolvente' y 'convergencia' son analíticas o de cuadratura en Python.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(os.path.dirname(AQUI))
REP = os.path.join(RAIZ, 'exe', 'rep')
SAL = os.path.join(RAIZ, 'docs', 'introduccion', 'figuras')
sys.path.insert(0, os.path.join(RAIZ, 'tools'))
import hk_exacto as H                                   # noqa: E402

plt.rcParams.update({'font.size': 9, 'axes.titlesize': 9, 'legend.fontsize': 7.5,
                     'figure.dpi': 150, 'savefig.bbox': 'tight'})
AZUL, NAR, VER, GRIS, ROJO = '#1f77b4', '#ff7f0e', '#2ca02c', '0.55', '#d62728'


def guarda(fig, nombre):
    os.makedirs(SAL, exist_ok=True)
    fig.savefig(os.path.join(SAL, nombre + '.pdf'))
    plt.close(fig)
    print(f'  -> figuras/{nombre}.pdf')


def complejo(fn):
    a = np.loadtxt(fn)
    return a[:, 0], a[:, 1::2] + 1j*a[:, 2::2]


def params(d):
    raw = H.leer_par(os.path.join(d, 'params_usados.par'))
    P = {}
    for k, v in raw.items():
        try:
            P[k] = float(v)
        except ValueError:
            P[k] = v
    for k in ('nrc', 'npc', 'nlc'):
        P[k] = int(P[k])
    return P


def c_de(L):
    return 0.5*(L + np.sqrt(L**2 + 4))


# ---------------------------------------------------------------------------
def frecuencias():
    print('frecuencias')
    J = np.linspace(0, 0.6, 400)
    fig, ax = plt.subplots(1, 3, figsize=(10, 2.9))
    for L, col in [(1.6, AZUL), (2.0, 'k'), (2.4, NAR)]:
        c = c_de(L)
        ax[0].plot(J, (J + c)**-3, color=col, label=f'$L={L}$')
        ax[1].plot(J, -3*(J + c)**-4, color=col)
        ax[2].plot(J, -3*(J + c)**-4*0.5*(1 + L/np.sqrt(L**2 + 4)), color=col)
    ax[0].set(xlabel='$J$', ylabel=r'$\omega$', title=r'(a) $\omega(J,L)=(J+c)^{-3}$')
    ax[1].set(xlabel='$J$', ylabel=r'$\partial\omega/\partial J$', title='(b) gradiente en $J$')
    ax[2].set(xlabel='$J$', ylabel=r'$\partial\omega/\partial L$', title='(c) gradiente en $L$')
    ax[0].legend()
    for a in ax:
        a.axvspan(0, 0.35, color='0.93', zorder=0)
    guarda(fig, 'frecuencias')
    c = c_de(2.0)
    wJ = -3*(0.1 + c)**-4
    wL = wJ*0.5*(1 + 2/np.sqrt(8))
    print(f'  J=0.1, L=2: omega={(0.1+c)**-3:.4f} dw/dJ={wJ:.4f} dw/dL={wL:.4f} cociente={wL/wJ:.3f}')
    for L in (1.6, 2.4):
        print(f'  L={L}: omega en J=0, 0.35: {(c_de(L))**-3:.4f} {(0.35+c_de(L))**-3:.4f}')


# ---------------------------------------------------------------------------
def envolvente():
    """Exacto continuo para la gaussiana con distintos anchos en L; ajuste del coeficiente de t^2."""
    print('envolvente')
    P0 = dict(a0=1e-3, l0=2.0, sp=0.1, sr=0.1, jminc=0.0, jmaxc=0.6,
              sq1=0.1, sj1=0.1, j1=0.0, lt1=2.0)
    t = np.linspace(0, 400, 201)
    casos = [(0.002, 'k', r'$\sigma_L=0.002$ ($L$ casi fija)'), (0.05, AZUL, r'$\sigma_L=0.05$'),
             (0.1, VER, r'$\sigma_L=0.1$'), (0.2, NAR, r'$\sigma_L=0.2$')]
    fig, ax = plt.subplots(1, 2, figsize=(8.5, 3.0))
    coefs = {}
    for sl, col, lab in casos:
        P = dict(P0, sl=sl, slt1=sl, lminc=2 - 4*sl, lmaxc=2 + 4*sl)
        h = H.hk_continuo(t, P, 1, 'gauss', 600, 120, 512)
        y = np.abs(h[:, 1])/np.abs(h[0, 1])
        ax[0].semilogy(t, y, color=col, label=lab)
        m = t <= 150
        coefs[sl] = np.polyfit(t[m]**2, np.log(y[m]), 1)[0]
        ax[1].plot(t[m]**2, np.log(y[m]), color=col)
    ax[0].set(xlabel='$t$', ylabel=r'$|h_1(t)|/|h_1(0)|$', title='(a) solución exacta, $k=1$', ylim=(1e-2, 1.5))
    ax[0].legend()
    ax[1].set(xlabel='$t^2$', ylabel=r'$\ln|h_1/h_1(0)|$', title=r'(b) tramo inicial ($t\leq150$)')
    guarda(fig, 'envolvente')
    # Predicción lineal: coef = -(omega_J^2 s_J^2 + omega_L^2 s_L^2)/4, con los anchos efectivos
    # del integrando en J y L, que se miden como segundos momentos del peso.
    c = c_de(2.0)
    base = coefs[0.002]
    for sl in (0.05, 0.1, 0.2):
        L = np.linspace(2 - 4*sl, 2 + 4*sl, 4001)
        w = L*np.exp(-2*(L - 2)**2/sl**2)
        var = np.sum(w*(L - np.sum(w*L)/w.sum())**2)/w.sum()
        wL = -3*(0.1 + c)**-4*0.5*(1 + 2/np.sqrt(8))
        pred = -0.5*wL**2*var            # exp(-x^2/s^2) con s^2 = 2 var: -(wL^2 s^2)/4
        print(f'  sigma_L={sl}: coef t^2 = {coefs[sl]:.4e}; exceso sobre L fija = {coefs[sl]-base:.4e};'
              f' predicción lineal -omega_L^2 var/2 = {pred:.4e}')
    print(f'  L casi fija: coef t^2 = {base:.4e}')


# ---------------------------------------------------------------------------
def verificacion():
    print('verificacion')
    V = os.path.join(REP, '08_verificacion')
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.2))
    fig.subplots_adjust(wspace=0.38)
    e = np.load(os.path.join(V, 'gauss_y4', 'hk1_exacto.npz'))
    for k, col in zip(range(1, 5), [AZUL, NAR, VER, ROJO]):
        ax[0].semilogy(e['t'], np.abs(e['codigo'][:, k]), color=col, label=f'$k={k}$')
        ax[0].semilogy(e['t'][::10], np.abs(e['continuo'][::10, k]), 'o', ms=2.5, mfc='none', color='k')
    ax[0].set(xlabel='$t$', ylabel=r'$|h_k|$', title='(a) código (líneas) y exacto (círculos)')
    ax[0].legend(ncol=2)
    for nombre, col, lab in [('gauss_y4', AZUL, 'yoshida4, $\\Delta t=0.025$'), ('gauss_an', NAR, 'analytic')]:
        e = np.load(os.path.join(V, nombre, 'hk1_exacto.npz'))
        esc = np.abs(e['continuo'][:, 1]).max()
        ax[1].semilogy(e['t'][1:], np.abs(e['codigo'][1:, 1] - e['discreto'][1:, 1])/esc, color=col, label=lab)
        print(f'  {nombre}: max |codigo-discreto|/max|h_1| = {np.abs(e["codigo"][:,1]-e["discreto"][:,1]).max()/esc:.2e};'
              f' |discreto-continuo| = {np.abs(e["discreto"][:,1]-e["continuo"][:,1]).max()/esc:.2e}')
    ax[1].semilogy(e['t'][1:], np.abs(e['discreto'][1:, 1] - e['continuo'][1:, 1])/esc, color=GRIS,
                   label='suma discreta $-$ continuo')
    ax[1].set(xlabel='$t$', ylabel=r'error relativo de $h_1$', title='(b) separación de errores')
    ax[1].legend(fontsize=6.5)
    # Convergencia de cuadratura: gauss en N_L (J convergido), king en N_J (L exacto).
    P = params(os.path.join(V, 'gauss_y4'))
    t = np.array([0., 100., 200.])
    hc = H.hk_continuo(t, P, 1, 'gauss', 600, 120, 512)
    nls = np.array([4, 8, 16, 32, 64, 128])
    err = []
    for nl in nls:
        Q = dict(P, nlc=int(nl))
        hd = H.hk_discreto(t, Q, 1, 'gauss')
        err.append(np.abs(hd[:, 1] - hc[:, 1]).max()/np.abs(hc[:, 1]).max())
    ax[2].loglog(nls, err, 'o-', color=AZUL, label='gauss, en $N_L$')
    print('  gauss N_L:', ' '.join(f'{n}:{x:.2e}' for n, x in zip(nls, err)))
    Pk = dict(a0=1e-3, l0=2.0, sl=0.2, sp=0.1, sr=0.1, jminc=0.0, jmaxc=0.35, lminc=1.6, lmaxc=2.4,
              sq1=0.5, sj1=0.2, j1=0.1, lt1=2.0, slt1=0.2)
    hk = H.hk_continuo(t, Pk, 1, 'king', 1000, 200, 256)
    xl, wl = np.polynomial.legendre.leggauss(200)
    Lg, wL = 2 + 0.4*xl, 0.4*wl
    njs = np.array([25, 50, 100, 200, 400, 800])
    errk = []
    for nj in njs:
        dJ = 0.35/nj
        Jn = (np.arange(nj) + 0.5)*dJ
        hd = H.hk_referencia(t, Pk, 1, 'king', Jn, np.full(nj, dJ), Lg, wL, 2*np.pi*np.arange(256)/256)
        errk.append(np.abs(hd[:, 1] - hk[:, 1]).max()/np.abs(hk[:, 1]).max())
    ax[2].loglog(njs, errk, 's-', color=NAR, label='king, en $N_J$')
    print('  king N_J:', ' '.join(f'{n}:{x:.2e}' for n, x in zip(njs, errk)))
    n = np.array([4, 800])
    ax[2].loglog(n, 2e-2*(n/4.)**-2, ':', color=GRIS, label=r'$\propto N^{-2}$')
    ax[2].set(xlabel='nodos', ylabel='error de cuadratura', title='(c) convergencia de la regla')
    ax[2].legend(fontsize=6.5)
    guarda(fig, 'verificacion')
    for nombre in ('bim_an', 'bim_y4', 'king_an', 'king_y4'):
        tt, z = complejo(os.path.join(V, nombre, 'hk1_complex.tl'))
        e = np.load(os.path.join(V, nombre, 'hk1_exacto.npz'))
        esc = np.abs(e['continuo'][:, 1]).max()
        prohib = (3, 4) if nombre.startswith('bim') else (2, 3, 4)
        print(f'  {nombre}: modos prohibidos max|h_k|/max|h_1| =',
              ' '.join(f'k={k}:{np.abs(z[:,k]).max()/esc:.1e}' for k in prohib),
              f'; k=1 codigo-discreto {np.abs(e["codigo"][:,1]-e["discreto"][:,1]).max()/esc:.1e}')


# ---------------------------------------------------------------------------
def espiral():
    print('espiral')
    V = os.path.join(REP, '08_verificacion')
    fig, ax = plt.subplots(1, 2, figsize=(8.5, 3.0))
    e0 = np.load(os.path.join(V, 'spi_L0_an', 'hk1_exacto.npz'))
    t0 = e0['t']
    for k, col in zip((1, 2, 3), (AZUL, NAR, VER)):
        ax[0].semilogy(t0, np.abs(e0['codigo'][:, k])/np.abs(e0['codigo'][0, k]), color=col, label=f'$k={k}$')
        m = t0 > 50
        i = np.argmax(np.abs(e0['codigo'][m, k]))
        print(f'  L fija, k={k}: pico en t={t0[m][i]:.0f}, |h_k|/|h_k(0)|={np.abs(e0["codigo"][m,k][i])/np.abs(e0["codigo"][0,k]):.3g}')
    c = c_de(2.0)
    tstar = 50.0/(3*(0.15 + c)**-4)
    ax[0].axvline(tstar, color=GRIS, ls='--', lw=0.8)
    print(f'  t* lineal = {tstar:.0f}')
    ax[0].set(xlabel='$t$', ylabel=r'$|h_k(t)|/|h_k(0)|$', title='(a) espiral, $L$ casi fija (código)')
    ax[0].legend()
    e = np.load(os.path.join(V, 'spi_L16_an', 'hk1_exacto.npz'))
    t = e['t']
    ax[1].semilogy(t, np.abs(e['codigo'][:, 1])/np.abs(e['codigo'][0, 1]), color=NAR, label=r'código, $\sigma_L=0.2$, $N_L=16$')
    ax[1].semilogy(t, np.abs(e['continuo'][:, 1])/np.abs(e['continuo'][0, 1]), '--', color='k', label=r'exacto continuo, $\sigma_L=0.2$')
    P = params(os.path.join(V, 'spi_L16_an'))
    tt = np.linspace(0, 1500, 151)
    for sl, col in [(0.02, AZUL), (0.05, VER)]:
        Q = dict(P, sl=sl, slt1=sl, lminc=2 - 4*sl, lmaxc=2 + 4*sl)
        h = H.hk_continuo(tt, Q, 1, 'spiral', 400, 120, 256)
        ax[1].semilogy(tt, np.abs(h[:, 1])/np.abs(h[0, 1]), color=col, lw=0.9, label=rf'exacto, $\sigma_L={sl}$')
        i = np.argmax(np.abs(h[5:, 1])) + 5
        print(f'  sigma_L={sl}: max |h_1|/|h_1(0)|={np.abs(h[i,1])/np.abs(h[0,1]):.2f} en t={tt[i]:.0f}')
    for tv in (700, 1000, 1500):
        i = np.argmin(abs(t - tv))
        print(f'  sigma_L=0.2 t={t[i]:.0f}: codigo {np.abs(e["codigo"][i,1])/np.abs(e["codigo"][0,1]):.2e}'
              f' discreto {np.abs(e["discreto"][i,1])/np.abs(e["discreto"][0,1]):.2e}'
              f' continuo {np.abs(e["continuo"][i,1])/np.abs(e["continuo"][0,1]):.2e}')
    # Recurrencia por muestreo en L: nodos vecinos separados dL = 0.05 desfasan 2 pi cuando
    # t = 2 pi/(k |dw/dL| dL); el primero en hacerlo es donde |dw/dL| es mayor (L = 1.6).
    Ln = 1.6 + (np.arange(16) + 0.5)*0.05
    wn = (0.15 + c_de(Ln))**-3
    dw = np.abs(np.diff(wn))
    print(f'  N_L=16, J=0.15: Delta omega_L total={wn[0]-wn[-1]:.4f}; 2pi/max(dw)={2*np.pi/dw.max():.0f};'
          f' 2pi/min(dw)={2*np.pi/dw.min():.0f}; 2pi N_L/Delta omega={2*np.pi*16/(wn[0]-wn[-1]):.0f}')
    for tv in (1100, 1200, 1300, 1400):
        i = np.argmin(abs(t - tv))
        print(f'  t={t[i]:.0f}: codigo/continuo = {np.abs(e["codigo"][i,1])/np.abs(e["continuo"][i,1]):.2f}')
    ax[1].axvline(tstar, color=GRIS, ls='--', lw=0.8)
    ax[1].set(xlabel='$t$', ylabel=r'$|h_1(t)|/|h_1(0)|$', title='(b) con dispersión en $L$', ylim=(1e-6, 5))
    ax[1].legend(fontsize=6.5, loc='upper left', bbox_to_anchor=(1.02, 1.0))
    guarda(fig, 'espiral')


# ---------------------------------------------------------------------------
def equilibrio():
    print('equilibrio')
    E = os.path.join(REP, '11_equilibrio')
    fig, ax = plt.subplots(1, 2, figsize=(8.5, 3.0))
    etiquetas = {'eq_e0': ('equilibrio, $\\varepsilon=0$', 'k'), 'eq_e01': ('equilibrio, $\\varepsilon=0.1$', AZUL),
                 'iso_aq': ('acciones del isócrono', NAR)}
    for d, (lab, col) in etiquetas.items():
        dp = np.load(os.path.join(E, d, 'delta_phi.npz'))
        m = (dp['r'] >= 1) & (dp['r'] <= 15)
        y = np.abs(dp['dphi'][:, m]).max(1)
        ax[0].semilogy(dp['t'][1:], y[1:], color=col, label=lab)
        print(f'  {d}: max_t max_r |dPhi| = {y.max():.2e}')
    ax[0].set(xlabel='$t$', ylabel=r'$\max_r|\Phi(r,t)-\Phi(r,0)|$', title=r'(a) el potencial, $1\leq r\leq15$',
              ylim=(1e-9, 3e-3))
    ax[0].legend(loc='lower right')
    import h5py
    from aa_numerico_L import phi_iso
    eq = np.load(os.path.join(E, 'ic_e0_equilibrio.npz'))
    h = h5py.File(os.path.join(E, 'eq_e0', 'vlasov_output.h5'), 'r')
    r = h['grid/r'][()]
    m = (r > 1) & (r < 15)
    pot0 = h['step_0000000000']['potential'][()]
    py = np.interp(r, eq['r'], eq['phi_self']) + phi_iso(r)
    print(f'  t=0: max|Phi_codigo - Phi_python| (1<r<15) = {np.abs(pot0-py)[m].max():.2e}; max|Phi_self| = {np.abs(eq["phi_self"]).max():.2e}')
    for d, col in (('eq_e0', 'k'), ('eq_e01', AZUL)):
        tc, zc = complejo(os.path.join(E, d, 'hk1_complex.tl'))
        a = np.loadtxt(os.path.join(E, d, 'hk1_numerico_eq.tl'))
        tn, zn = a[:, 0], a[:, 1::2] + 1j*a[:, 2::2]
        ai = np.loadtxt(os.path.join(E, d, 'hk1_numerico.tl'))
        zi = ai[:, 1::2] + 1j*ai[:, 2::2]
        print(f'  {d}: marco instantáneo |h1|/h0 t=0,50,100: ' +
              ' '.join(f'{abs(zi[i,1])/abs(zi[i,0]):.3e}' for i in (0, len(zi)//2, len(zi)-1)))
        ax[1].semilogy(tc, np.abs(zc[:, 1])/np.abs(zc[:, 0]), color=col, ls='--', lw=0.9)
        ax[1].semilogy(tn, np.abs(zn[:, 1])/np.abs(zn[:, 0]), 'o-', ms=2.5, color=col, lw=0.9)
        idx = [np.argmin(abs(tc - x)) for x in tn]
        print(f'  {d}: |h1|/h0 isócrono t=0,50,100: ' +
              ' '.join(f'{abs(zc[i,1])/abs(zc[i,0]):.3e}' for i in (idx[0], idx[len(idx)//2], idx[-1])) +
              ' ; numérico (marco del equilibrio): ' + ' '.join(f'{abs(zn[i,1])/abs(zn[i,0]):.3e}' for i in (0, len(tn)//2, len(tn)-1)))
    ax[1].plot([], [], 'k--', lw=0.9, label='mapa del isócrono (código)')
    ax[1].plot([], [], 'ko-', ms=2.5, lw=0.9, label='mapa numérico, potencial de equilibrio')
    ax[1].set(xlabel='$t$', ylabel=r'$|h_1|/h_0$', title='(b) $h_1$ con los dos mapas')
    ax[1].legend(fontsize=6.5)
    guarda(fig, 'equilibrio')
    print(f'  a_1/a_0 * eps/2 (sigma_Q=0.1) = {H.a_k(0.1)[1]/H.a_k(0.1)[0]*0.05:.5f}')


# ---------------------------------------------------------------------------
def fondos():
    print('fondos (sin figura)')
    F = os.path.join(REP, '09_fondos')
    for bg in ('iso', 'isotrun', 'nfw', 'burkert'):
        errs = []
        for c in ('0.5', '0.25', '0.125'):
            E = np.loadtxt(os.path.join(F, f'{bg}_c{c}', 'vlasov_k_phi_e.tl'))[:, 3]
            errs.append(np.max(np.abs(E/E[0] - 1)))
        print(f'  {bg}: max|dE/E| ' + ' '.join(f'{x:.2e}' for x in errs) +
              f'  cocientes {errs[0]/errs[1]:.1f} {errs[1]/errs[2]:.1f}')
    G = os.path.join(REP, '10_energia')
    for d in ('aa_sg', 'gauss1_sg', 'sphere_sg'):
        E = np.loadtxt(os.path.join(G, d, 'vlasov_k_phi_e.tl'))
        print(f'  {d}: t hasta {E[-1,0]:g}, max|dE/E| = {np.max(np.abs(E[:,3]/E[0,3]-1)):.2e}')


def mezcla():
    """Solución exacta F(Q,J,L,t) = F0(Q - omega(J,L) t, J, L) para la gaussiana, en dos cortes."""
    print('mezcla')
    Q = np.linspace(0, 2*np.pi, 400)
    J = np.linspace(0.0, 0.35, 300)
    L = np.linspace(1.6, 2.4, 300)
    tiempos = (0, 400, 2000)
    fig, ax = plt.subplots(2, 3, figsize=(10, 5.2), sharex=True)
    for j, t in enumerate(tiempos):
        QQ, JJ = np.meshgrid(Q, J)
        w = (JJ + c_de(2.0))**-3
        F = np.exp(-np.sin(0.5*(QQ - w*t))**2/0.1**2)*np.exp(-JJ**2/0.1**2)*JJ**2
        ax[0, j].pcolormesh(Q, J, F, shading='auto', cmap='Blues', rasterized=True)
        ax[0, j].set(title=f'$t={t}$')
        QQ, LL = np.meshgrid(Q, L)
        w = (0.1 + c_de(LL))**-3
        F = np.exp(-np.sin(0.5*(QQ - w*t))**2/0.1**2)*np.exp(-(LL - 2)**2/0.2**2)
        ax[1, j].pcolormesh(Q, L, F, shading='auto', cmap='Oranges', rasterized=True)
        ax[1, j].set(xlabel='$Q$')
    ax[0, 0].set(ylabel='$J$  ($L=2$)')
    ax[1, 0].set(ylabel='$L$  ($J=0.1$)')
    guarda(fig, 'mezcla')
    for t in (400, 2000):
        dJ = 2*np.pi/(abs(-3*(0.1 + c_de(2.0))**-4)*t)
        dL = 2*np.pi/(abs(-3*(0.1 + c_de(2.0))**-4*0.5*(1 + 2/np.sqrt(8)))*t)
        print(f'  t={t}: separación de franjas k=1: dJ={dJ:.3f}, dL={dL:.3f}')


FIGURAS = [frecuencias, mezcla, envolvente, verificacion, espiral, equilibrio, fondos]

if __name__ == '__main__':
    filtro = sys.argv[1] if len(sys.argv) > 1 else ''
    for f in FIGURAS:
        if filtro in f.__name__:
            f()
