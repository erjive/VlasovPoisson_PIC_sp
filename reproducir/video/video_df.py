"""Video 3D de la evolución de la DF gaussiana en (r, p_r, L), (Q, J, L) y (x, y, z).

Lee las instantáneas HDF5 de la corrida reproducir/video/df_gauss.par (sin autogravedad,
integrador analytic) y dibuja, en cada instantánea, las partículas del código:

  (r, p_r, L)  tal como las guarda el código;
  (Q, J, L)    con el mapa ángulo-acción del isócrono (el mismo de analysish.f90), Q en
               [-pi, pi) para que la distribución inicial quede centrada;
  (x, y, z)    cada partícula representa una capa esférica: se dibuja con K estrellas en
               planos orbitales de orientación aleatoria (fija en el tiempo). El radio es el
               del código; el ángulo dentro del plano, psi, no lo calcula el código y se
               obtiene exacto de dpsi/dt = L/r^2: como r depende solo de Q,
                   dpsi/dQ = L / (omega r(Q)^2),
               que se integra una vez por cada par (J, L) en una rejilla de Q (trapecio
               periódico) y da psi(t) = psi0 + Psi(Q0 + omega t) - Psi(Q0).
               Se quita una cuña de 90 grados que mira a la cámara para ver el interior.

El color y la opacidad son el peso f de cada partícula (la DF), en escala logarítmica.

Uso:
    python3 reproducir/video/video_df.py --valida          # comprueba psi(t)
    python3 reproducir/video/video_df.py --cuadro 150      # un cuadro de prueba (PNG)
    python3 reproducir/video/video_df.py                   # todos los cuadros y el MP4
"""
import argparse, os, subprocess, sys, time
from multiprocessing import Pool
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import mpl_toolkits
# El mpl_toolkits del sistema puede no corresponder a la versión de matplotlib instalada.
mpl_toolkits.__path__[:] = [os.path.join(os.path.dirname(os.path.dirname(matplotlib.__file__)), 'mpl_toolkits')]
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401,E402
import matplotlib.pyplot as plt  # noqa: E402

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(os.path.dirname(AQUI))
RUN = os.path.join(RAIZ, 'exe', 'rep', 'video', 'df_gauss')
SAL = os.path.join(RAIZ, 'exe', 'rep', 'video')
sys.path.insert(0, os.path.join(RAIZ, 'tools'))
from aa_numerico_L import aa_isocrono  # noqa: E402

K = 6           # estrellas por capa en (x, y, z)
SEMILLA = 20260917


def c_de(L):
    return 0.5*(L + np.sqrt(L**2 + 4))


def r_de_QJL(Q, J, L):
    """Inversa del mapa del isócrono (misma fórmula que invert_QJ_to_rp en utils.f90)."""
    E = -1.0/(2.0*(J + c_de(L))**2)
    er1 = np.sqrt((1 + E*(2 + L**2) - np.sqrt(1 + 2*E*(2 + 2*E + L**2)))/(2*E**2))
    er2 = np.sqrt((1 + E*(2 + L**2) + np.sqrt(1 + 2*E*(2 + 2*E + L**2)))/(2*E**2))
    s1, s2 = 1 + np.sqrt(1 + er1**2), 1 + np.sqrt(1 + er2**2)
    ecc = np.sqrt((-2*E)**3)*np.sqrt(np.maximum(-L**2 - 2*E - 2 - 0.5/E, 0))/(-2*E)
    Qm = np.mod(Q, 2*np.pi)
    eta = Qm.copy()
    for _ in range(60):
        eta = eta - (eta - ecc*np.sin(eta) - Qm)/(1 - ecc*np.cos(eta))
    s = 0.5*(s1 + s2 - np.cos(eta)*(s2 - s1))
    return np.sqrt(np.maximum((s - 1)**2 - 1, 0))


class Azimut:
    """Psi(Qu) acumulado para cada par (J, L) distinto, con Qu el ángulo radial sin envolver."""

    def __init__(self, J, L, nq=1024):
        self.claves, self.inv = np.unique(np.round(np.column_stack([J, L]), 12), axis=0, return_inverse=True)
        self.inv = self.inv.ravel()
        Ju, Lu = self.claves[:, 0], self.claves[:, 1]
        q = 2*np.pi*np.arange(nq + 1)/nq
        om = (Ju + c_de(Lu))**-3
        r = r_de_QJL(q[None, :], Ju[:, None], Lu[:, None])
        g = Lu[:, None]/(om[:, None]*r**2)                         # dpsi/dQ
        acum = np.concatenate([np.zeros((Ju.size, 1)), np.cumsum(0.5*(g[:, 1:] + g[:, :-1])*(q[1] - q[0]), axis=1)], axis=1)
        self.q, self.tabla, self.periodo = q, acum, acum[:, -1]      # psi por periodo radial
        self.nq = nq

    def __call__(self, Qu, idx):
        """Psi en Qu para las partículas con claves idx."""
        n = np.floor(Qu/(2*np.pi))
        x = (Qu - 2*np.pi*n)/(2*np.pi)*self.nq
        i = np.clip(x.astype(int), 0, self.nq - 1)
        w = x - i
        t = self.tabla
        return n*self.periodo[idx] + (1 - w)*t[idx, i] + w*t[idx, i + 1]


def carga():
    h = h5py.File(os.path.join(RUN, 'vlasov_output.h5'), 'r')
    pasos = sorted([g for g in h if g.startswith('step_')], key=lambda g: int(g.split('_')[1]))
    return h, pasos


def preparar(h, pasos):
    g0 = h[pasos[0]]
    r0, p0, L = g0['r_part'][()], g0['p_part'][()], g0['l_part'][()]
    f = g0['fl'][()]/L
    Q0, J, _ = aa_isocrono(r0, p0, L)
    om = (J + c_de(L))**-3
    az = Azimut(J, L)
    rng = np.random.default_rng(SEMILLA)
    n = r0.size
    # Plano orbital aleatorio para cada una de las K copias: normal uniforme en la esfera,
    # y una dirección de referencia e1 en el plano; psi0 aleatorio.
    nrm = rng.normal(size=(n, K, 3)); nrm /= np.linalg.norm(nrm, axis=2, keepdims=True)
    a = rng.normal(size=(n, K, 3))
    e1 = a - np.sum(a*nrm, axis=2, keepdims=True)*nrm; e1 /= np.linalg.norm(e1, axis=2, keepdims=True)
    e2 = np.cross(nrm, e1)
    psi0 = rng.uniform(0, 2*np.pi, size=(n, K))
    return dict(L=L, f=f, Q0=Q0, J=J, om=om, az=az, e1=e1, e2=e2, psi0=psi0,
                Psi0=az(Q0, az.inv))


def valida():
    """psi(t) de las tablas frente a integrar L/r(t)^2 con r(t) del mapa inverso, paso fino."""
    h, pasos = carga()
    P = preparar(h, pasos)
    rng = np.random.default_rng(1)
    sel = rng.choice(P['L'].size, 200, replace=False)
    T, dt = 2400.0, 0.02
    t = np.arange(0, T + dt/2, dt)
    peor = 0.0
    for j in sel:
        Qu = P['Q0'][j] + P['om'][j]*t
        r = r_de_QJL(Qu, np.full_like(t, P['J'][j]), np.full_like(t, P['L'][j]))
        g = P['L'][j]/r**2
        num = np.concatenate([[0], np.cumsum(0.5*(g[1:] + g[:-1])*dt)])
        tab = P['az'](Qu, np.full(t.size, P['az'].inv[j])) - P['Psi0'][j]
        peor = max(peor, np.max(np.abs(tab - num)))
    print(f'200 partículas, t<=2400: max |psi(tabla) - psi(integración, dt={dt})| = {peor:.2e} rad')
    # Consistencia con el código: r(t) del mapa inverso frente a r_part de las instantáneas.
    dif = 0.0
    for k in (0, 150, 300, 600):
        g = h[pasos[k]]
        tk = g.attrs['time']
        rr = r_de_QJL(P['Q0'] + P['om']*tk, P['J'], P['L'])
        dif = max(dif, np.max(np.abs(rr - g['r_part'][()])))
    print(f'r(t) del mapa inverso frente al código en t=0,600,1200,2400: max|dr| = {dif:.2e}')


LIM = {}


def limites(h, pasos, P):
    rs = np.concatenate([h[pasos[k]]['r_part'][()] for k in range(0, len(pasos), 20)])
    ps = np.concatenate([h[pasos[k]]['p_part'][()] for k in range(0, len(pasos), 20)])
    rmax = np.percentile(rs, 99.9)
    return dict(r=(np.min(rs), rmax), p=(np.min(ps), np.max(ps)), L=(1.6, 2.4),
                J=(0, np.max(P['J'])), x=rmax)


ESTILO = dict(fondo='#0b0e14', texto='#d8dee9', rejilla='#2e3440')


def ejes(ax, titulo, etiquetas):
    ax.set_facecolor(ESTILO['fondo'])
    for eje in (ax.xaxis, ax.yaxis, ax.zaxis):
        eje.set_pane_color((0.07, 0.09, 0.13, 1.0))
        eje._axinfo['grid']['color'] = ESTILO['rejilla']
        eje.label.set_color(ESTILO['texto'])
    ax.tick_params(colors=ESTILO['texto'], labelsize=7, pad=0)
    ax.set_xlabel(etiquetas[0], labelpad=2, fontsize=10)
    ax.set_ylabel(etiquetas[1], labelpad=2, fontsize=10)
    ax.set_zlabel(etiquetas[2], labelpad=2, fontsize=10)
    ax.set_title(titulo, color=ESTILO['texto'], fontsize=15, pad=-6)


def colores(f):
    lf = np.log10(np.maximum(f/f.max(), 1e-3))
    u = (lf + 3)/3                                        # 0..1
    c = plt.get_cmap('inferno')(0.25 + 0.75*u)
    c[:, 3] = 0.05 + 0.55*u**1.5
    return c


def cuadro(k, h=None, pasos=None, P=None, lim=None, destino=None):
    g = h[pasos[k]]
    t = g.attrs['time']
    r, p, L, f = g['r_part'][()], g['p_part'][()], g['l_part'][()], P['f']
    Q, J, _ = aa_isocrono(r, p, L)
    Qc = np.mod(Q + np.pi, 2*np.pi) - np.pi
    col = colores(f)
    orden = np.argsort(f)                                  # los pesos grandes al frente

    fig = plt.figure(figsize=(19.2, 10.8), dpi=100, facecolor=ESTILO['fondo'])
    gira = 35 + 40*k/len(pasos)                            # rotación lenta de la cámara

    ax1 = fig.add_axes([-0.01, 0.06, 0.35, 0.84], projection='3d')
    ax1.scatter(r[orden], p[orden], L[orden], s=0.6, c=col[orden], linewidths=0, depthshade=False)
    ejes(ax1, r'$(r,\ p_r,\ L)$', (r'$r$', r'$p_r$', r'$L$'))
    ax1.set_xlim(*lim['r']); ax1.set_ylim(*lim['p']); ax1.set_zlim(*lim['L'])
    ax1.set_box_aspect(None, zoom=1.02)
    ax1.view_init(elev=22, azim=-60 + gira)

    ax2 = fig.add_axes([0.325, 0.06, 0.35, 0.84], projection='3d')
    ax2.scatter(Qc[orden], J[orden], L[orden], s=0.6, c=col[orden], linewidths=0, depthshade=False)
    ejes(ax2, r'$(Q,\ J,\ L)$', (r'$Q$', r'$J$', r'$L$'))
    ax2.set_xlim(-np.pi, np.pi); ax2.set_ylim(*lim['J']); ax2.set_zlim(*lim['L'])
    ax2.set_box_aspect(None, zoom=1.02)
    ax2.view_init(elev=22, azim=-60 + gira)

    # (x, y, z): K estrellas por capa
    Qu = P['Q0'] + P['om']*t
    psi = P['psi0'] + (P['az'](Qu, P['az'].inv) - P['Psi0'])[:, None]
    pos = r[:, None, None]*(np.cos(psi)[..., None]*P['e1'] + np.sin(psi)[..., None]*P['e2'])
    xyz = pos.reshape(-1, 3)
    cK = np.repeat(col, K, axis=0)
    fK = np.repeat(f*L, K)
    # Cuña de 90 grados centrada en la dirección de la cámara (azimut de la vista).
    a = np.deg2rad(-60 + gira)
    phi = np.arctan2(xyz[:, 1], xyz[:, 0])
    visible = np.abs(np.angle(np.exp(1j*(phi - a)))) > np.pi/4
    o = np.argsort(fK[visible])
    X, cc = xyz[visible][o], cK[visible][o]
    cc = cc.copy(); cc[:, 3] = np.minimum(1.0, cc[:, 3]*0.9)
    ax3 = fig.add_axes([0.66, 0.06, 0.35, 0.84], projection='3d')
    ax3.scatter(X[:, 0], X[:, 1], X[:, 2], s=0.5, c=cc, linewidths=0, depthshade=False)
    ejes(ax3, r'$(x,\ y,\ z)$', (r'$x$', r'$y$', r'$z$'))
    m = lim['x']
    ax3.set_xlim(-m, m); ax3.set_ylim(-m, m); ax3.set_zlim(-m, m)
    ax3.set_box_aspect((1, 1, 1), zoom=1.02)
    ax3.view_init(elev=22, azim=-60 + gira)

    T_r = 2*np.pi/np.median(P['om'])
    fig.text(0.5, 0.955, 'Phase mixing de la DF gaussiana en el isócrono (sin autogravedad)',
             ha='center', color=ESTILO['texto'], fontsize=18)
    fig.text(0.5, 0.915, rf'$t = {t:6.0f}$   ($\approx {t/T_r:4.1f}$ periodos radiales)',
             ha='center', color=ESTILO['texto'], fontsize=15)
    fig.text(0.5, 0.03,
             r'$F_0\propto \exp[-\sin^2(Q/2)/0.1^2]\,J^2\,\exp(-J^2/0.1^2)\,\exp[-(L-2)^2/0.2^2]$'
             '    ·    58 560 partículas de VP_PIC (integrador analytic)'
             '    ·    color y opacidad: peso $f$    ·    en $(x,y,z)$, 6 estrellas por capa, sin la cuña que mira a la cámara',
             ha='center', color='#8f9bb3', fontsize=10.5)
    fig.savefig(destino, facecolor=ESTILO['fondo'])
    plt.close(fig)


_G = {}


def _init():
    h, pasos = carga()
    P = preparar(h, pasos)
    _G.update(h=h, pasos=pasos, P=P, lim=limites(h, pasos, P))


def _trabajo(k):
    d = os.path.join(SAL, 'cuadros', f'c{k:04d}.png')
    if not os.path.exists(d):
        cuadro(k, destino=d, **_G)
    return k


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--valida', action='store_true')
    ap.add_argument('--cuadro', type=int)
    ap.add_argument('--procesos', type=int, default=6)
    ap.add_argument('--fps', type=int, default=30)
    a = ap.parse_args()
    if a.valida:
        valida()
        return
    os.makedirs(os.path.join(SAL, 'cuadros'), exist_ok=True)
    if a.cuadro is not None:
        _init()
        t0 = time.time()
        cuadro(a.cuadro, destino=os.path.join(SAL, f'prueba_{a.cuadro}.png'), **_G)
        print(f'cuadro {a.cuadro} en {time.time()-t0:.1f} s')
        return
    h, pasos = carga()
    n = len(pasos)
    h.close()
    t0 = time.time()
    with Pool(a.procesos, initializer=_init) as pool:
        for i, _ in enumerate(pool.imap_unordered(_trabajo, range(n))):
            if (i + 1) % 50 == 0:
                print(f'{i+1}/{n} cuadros, {time.time()-t0:.0f} s', flush=True)
    salida = os.path.join(SAL, 'df_gauss_3d.mp4')
    subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', str(a.fps),
                    '-i', os.path.join(SAL, 'cuadros', 'c%04d.png'),
                    '-c:v', 'libx264', '-preset', 'slow', '-crf', '17', '-pix_fmt', 'yuv420p',
                    '-movflags', '+faststart', salida], check=True)
    print('escrito', salida)


if __name__ == '__main__':
    main()
