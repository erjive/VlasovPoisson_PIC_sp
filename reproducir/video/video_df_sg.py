"""Video 3D de la DF gaussiana con autogravedad en (r, p_r, L), (Q, J, L) y (x, y, z).

Lee las instantáneas HDF5 de reproducir/video/df_gauss_sg.par (a0=1e-2, yoshida4,
dt=0.025, una instantánea cada 4 unidades) y dibuja los mismos tres paneles que
video_df.py, con dos diferencias impuestas por la autogravedad:

  (Q, J, L)  con el mapa ángulo-acción NUMÉRICO en el potencial de cada instantánea
             (isócrono + potencial propio guardado por el código, tools/aa_numerico_L.py).
             Con el mapa del isócrono aparecería el artefacto de coordenadas de las notas.
  (x, y, z)  el ángulo en el plano orbital ya no tiene forma cerrada: se integra
             dpsi/dt = L/r(t)^2 entre instantáneas, con r(t) interpolado por Hermite
             cúbico con los extremos r y dr/dt = p_r de cada una, y Gauss-Legendre de
             8 nodos por intervalo. El método se valida en la corrida sin autogravedad
             contra el psi exacto (--valida).

Uso:
    python3 reproducir/video/video_df_sg.py --valida
    python3 reproducir/video/video_df_sg.py --aa          # (Q,J) de todas las instantáneas
    python3 reproducir/video/video_df_sg.py --cuadro 300
    python3 reproducir/video/video_df_sg.py               # cuadros y MP4
"""
import argparse, os, subprocess, sys, time
from multiprocessing import Pool
import numpy as np
import h5py

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import video_df as V                         # noqa: E402  estilo, geometría y psi exacto
from video_df import plt                     # noqa: E402
from aa_numerico_L import MapaAA, phi_iso    # noqa: E402

RUN = os.path.join(V.RAIZ, 'exe', 'rep', 'video', 'df_gauss_sg')
SAL = os.path.join(V.RAIZ, 'exe', 'rep', 'video')
AA = os.path.join(SAL, 'aa_sg')


def carga(run=RUN):
    h = h5py.File(os.path.join(run, 'vlasov_output.h5'), 'r')
    pasos = sorted([g for g in h if g.startswith('step_')], key=lambda g: int(g.split('_')[1]))
    return h, pasos


# ---------------------------------------------------------------------------
xg, wg = np.polynomial.legendre.leggauss(8)


def integra_psi(h, pasos):
    """psi_k acumulado para todas las partículas: Hermite cúbico de r(t) entre instantáneas."""
    g = h[pasos[0]]
    L = g['l_part'][()]
    t0, r0, p0 = g.attrs['time'], g['r_part'][()], g['p_part'][()]
    psi = np.zeros((len(pasos), L.size))
    for k in range(1, len(pasos)):
        g = h[pasos[k]]
        t1, r1, p1 = g.attrs['time'], g['r_part'][()], g['p_part'][()]
        dt = t1 - t0
        suma = np.zeros_like(L)
        for x, w in zip(xg, wg):
            s = 0.5*(x + 1)
            h00, h10, h01, h11 = 2*s**3 - 3*s**2 + 1, s**3 - 2*s**2 + s, -2*s**3 + 3*s**2, s**3 - s**2
            r = h00*r0 + h10*dt*p0 + h01*r1 + h11*dt*p1
            suma += 0.5*w*L/r**2
        psi[k] = psi[k-1] + suma*dt
        t0, r0, p0 = t1, r1, p1
    return psi


def valida():
    h, pasos = V.carga()                     # corrida SIN autogravedad, psi exacto conocido
    P = V.preparar(h, pasos)
    psi = integra_psi(h, pasos)
    peor = 0.0
    for k in range(0, len(pasos), 50):
        t = h[pasos[k]].attrs['time']
        exacto = P['az'](P['Q0'] + P['om']*t, P['az'].inv) - P['Psi0']
        peor = max(peor, np.max(np.abs(psi[k] - exacto)))
    print(f'Hermite + Gauss 8, instantáneas cada 4: max |psi - psi exacto| hasta t=2400 = {peor:.2e} rad'
          f' ({V.K} estrellas por capa, {psi.shape[1]} capas)')


# ---------------------------------------------------------------------------
def _aa_trabajo(k):
    destino = os.path.join(AA, f'aa{k:04d}.npz')
    if os.path.exists(destino):
        return k
    h, pasos = carga()
    g = h[pasos[k]]
    rg = h['grid/r'][()]
    m = MapaAA(rg, g['potential'][()] - phi_iso(rg))
    r, p, L = g['r_part'][()], g['p_part'][()], g['l_part'][()]
    E = 0.5*p**2 + m.phi_ef(r, L)
    Q = np.full(r.size, np.nan); J = np.full(r.size, np.nan)
    lig = E < 0
    Q[lig], J[lig], _ = m(r[lig], p[lig], L[lig])
    np.savez(destino, Q=Q.astype(np.float32), J=J.astype(np.float32), no_ligadas=int((~lig).sum()))
    return k


def calcula_aa(procesos):
    os.makedirs(AA, exist_ok=True)
    h, pasos = carga()
    n = len(pasos)
    h.close()
    t0 = time.time()
    with Pool(procesos) as pool:
        for i, _ in enumerate(pool.imap_unordered(_aa_trabajo, range(n))):
            if (i + 1) % 50 == 0:
                print(f'(Q,J): {i+1}/{n} instantáneas, {time.time()-t0:.0f} s', flush=True)


# ---------------------------------------------------------------------------
def cuadro(k, h, pasos, G, destino):
    g = h[pasos[k]]
    t = g.attrs['time']
    r, p, L = g['r_part'][()], g['p_part'][()], g['l_part'][()]
    f = G['f']
    d = np.load(os.path.join(AA, f'aa{k:04d}.npz'))
    Q, J = d['Q'].astype(float), d['J'].astype(float)
    Qc = np.mod(Q + np.pi, 2*np.pi) - np.pi
    col = V.colores(f)
    orden = np.argsort(f)
    lim = G['lim']

    fig = plt.figure(figsize=(19.2, 10.8), dpi=100, facecolor=V.ESTILO['fondo'])
    gira = 35 + 40*k/len(pasos)

    ax1 = fig.add_axes([-0.01, 0.06, 0.35, 0.84], projection='3d')
    ax1.scatter(r[orden], p[orden], L[orden], s=0.6, c=col[orden], linewidths=0, depthshade=False)
    V.ejes(ax1, r'$(r,\ p_r,\ L)$', (r'$r$', r'$p_r$', r'$L$'))
    ax1.set_xlim(*lim['r']); ax1.set_ylim(*lim['p']); ax1.set_zlim(*lim['L'])
    ax1.set_box_aspect(None, zoom=1.02)
    ax1.view_init(elev=22, azim=-60 + gira)

    ax2 = fig.add_axes([0.325, 0.06, 0.35, 0.84], projection='3d')
    ok = np.isfinite(Q)
    o2 = orden[ok[orden]]
    ax2.scatter(Qc[o2], J[o2], L[o2], s=0.6, c=col[o2], linewidths=0, depthshade=False)
    V.ejes(ax2, r'$(Q,\ J,\ L)$  del potencial real', (r'$Q$', r'$J$', r'$L$'))
    ax2.set_xlim(-np.pi, np.pi); ax2.set_ylim(*lim['J']); ax2.set_zlim(*lim['L'])
    ax2.set_box_aspect(None, zoom=1.02)
    ax2.view_init(elev=22, azim=-60 + gira)

    psi = G['psi0'] + np.asarray(G['psi'][k])[:, None]
    pos = r[:, None, None]*(np.cos(psi)[..., None]*G['e1'] + np.sin(psi)[..., None]*G['e2'])
    xyz = pos.reshape(-1, 3)
    cK = np.repeat(col, V.K, axis=0)
    fK = np.repeat(f*L, V.K)
    a = np.deg2rad(-60 + gira)
    phi = np.arctan2(xyz[:, 1], xyz[:, 0])
    visible = np.abs(np.angle(np.exp(1j*(phi - a)))) > np.pi/4
    o = np.argsort(fK[visible])
    X, cc = xyz[visible][o], cK[visible][o].copy()
    cc[:, 3] = np.minimum(1.0, cc[:, 3]*0.9)
    ax3 = fig.add_axes([0.66, 0.06, 0.35, 0.84], projection='3d')
    ax3.scatter(X[:, 0], X[:, 1], X[:, 2], s=0.5, c=cc, linewidths=0, depthshade=False)
    V.ejes(ax3, r'$(x,\ y,\ z)$', (r'$x$', r'$y$', r'$z$'))
    mm = lim['x']
    ax3.set_xlim(-mm, mm); ax3.set_ylim(-mm, mm); ax3.set_zlim(-mm, mm)
    ax3.set_box_aspect((1, 1, 1), zoom=1.02)
    ax3.view_init(elev=22, azim=-60 + gira)

    fig.text(0.5, 0.955, r'Phase mixing de la DF gaussiana en el isócrono, con autogravedad ($a_0=10^{-2}$)',
             ha='center', color=V.ESTILO['texto'], fontsize=18)
    fig.text(0.5, 0.915, rf'$t = {t:6.0f}$   ($\approx {t/G["Tr"]:4.1f}$ periodos radiales)',
             ha='center', color=V.ESTILO['texto'], fontsize=15)
    fig.text(0.5, 0.042,
             r'$F_0\propto \exp[-\sin^2(Q/2)/0.1^2]\,J^2\,\exp(-J^2/0.1^2)\,\exp[-(L-2)^2/0.2^2]$'
             '    ·    58 560 partículas de VP_PIC (yoshida4, $\\Delta t=0.025$)'
             '    ·    color y opacidad: peso $f$',
             ha='center', color='#8f9bb3', fontsize=11)
    fig.text(0.5, 0.012,
             '$(Q,J)$: mapa ángulo-acción numérico en el potencial de cada instante'
             '    ·    en $(x,y,z)$, 6 estrellas por capa, sin la cuña que mira a la cámara',
             ha='center', color='#8f9bb3', fontsize=11)
    fig.savefig(destino, facecolor=V.ESTILO['fondo'])
    plt.close(fig)


_G = {}


def _init():
    h, pasos = carga()
    g0 = h[pasos[0]]
    L = g0['l_part'][()]
    f = g0['fl'][()]/L
    n = L.size
    rng = np.random.default_rng(V.SEMILLA)
    nrm = rng.normal(size=(n, V.K, 3)); nrm /= np.linalg.norm(nrm, axis=2, keepdims=True)
    a = rng.normal(size=(n, V.K, 3))
    e1 = a - np.sum(a*nrm, axis=2, keepdims=True)*nrm; e1 /= np.linalg.norm(e1, axis=2, keepdims=True)
    e2 = np.cross(nrm, e1)
    psi0 = rng.uniform(0, 2*np.pi, size=(n, V.K))
    # mmap: psi.npy son 281 MB y cada proceso lee solo la fila de su cuadro.
    psi = np.load(os.path.join(AA, 'psi.npy'), mmap_mode='r')
    J0 = np.load(os.path.join(AA, 'aa0000.npz'))['J']
    rs = np.concatenate([h[pasos[k]]['r_part'][()] for k in range(0, len(pasos), 20)])
    ps = np.concatenate([h[pasos[k]]['p_part'][()] for k in range(0, len(pasos), 20)])
    Js = np.concatenate([np.load(os.path.join(AA, f'aa{k:04d}.npz'))['J'] for k in range(0, len(pasos), 20)])
    rmax = np.percentile(rs, 99.9)
    lim = dict(r=(np.min(rs), rmax), p=(np.min(ps), np.max(ps)), L=(1.6, 2.4),
               J=(0, float(np.nanpercentile(Js, 99.9))), x=rmax)
    Tr = 2*np.pi/np.median((J0 + V.c_de(L))**-3)
    _G.update(h=h, pasos=pasos, G=dict(f=f, e1=e1, e2=e2, psi0=psi0, psi=psi, lim=lim, Tr=Tr))


def _trabajo(k):
    d = os.path.join(SAL, 'cuadros_sg', f'c{k:04d}.png')
    if not os.path.exists(d):
        cuadro(k, destino=d, **_G)
    return k


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--valida', action='store_true')
    ap.add_argument('--aa', action='store_true')
    ap.add_argument('--cuadro', type=int)
    ap.add_argument('--procesos', type=int, default=2)
    ap.add_argument('--fps', type=int, default=30)
    a = ap.parse_args()
    if a.valida:
        valida()
        return
    if a.aa:
        calcula_aa(a.procesos)
        h, pasos = carga()
        np.save(os.path.join(AA, 'psi.npy'), integra_psi(h, pasos))
        nl = [int(np.load(os.path.join(AA, f'aa{k:04d}.npz'))['no_ligadas']) for k in range(len(pasos))]
        print('psi integrado; partículas no ligadas por instantánea: máx', max(nl))
        return
    os.makedirs(os.path.join(SAL, 'cuadros_sg'), exist_ok=True)
    if a.cuadro is not None:
        _init()
        t0 = time.time()
        cuadro(a.cuadro, destino=os.path.join(SAL, f'prueba_sg_{a.cuadro}.png'), **_G)
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
    salida = os.path.join(SAL, 'df_gauss_sg_3d.mp4')
    subprocess.run(['ffmpeg', '-y', '-loglevel', 'error', '-framerate', str(a.fps),
                    '-i', os.path.join(SAL, 'cuadros_sg', 'c%04d.png'),
                    '-c:v', 'libx264', '-preset', 'slow', '-crf', '17', '-pix_fmt', 'yuv420p',
                    '-movflags', '+faststart', salida], check=True)
    print('escrito', salida)


if __name__ == '__main__':
    main()
