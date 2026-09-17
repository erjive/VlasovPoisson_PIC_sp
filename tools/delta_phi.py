"""Perturbación del potencial, delta Phi(r,t), desde las instantáneas HDF5 de una corrida.

delta Phi(r,t) = Phi(r,t) - Phi(r,0) no depende de ninguna elección de coordenadas en el
espacio fase, a diferencia de h_k, así que es el primer control de si "algo se mueve".
Con --referencia se resta además la misma cantidad de una corrida sin perturbación en la
misma rejilla (el equilibrio, eps=0), lo que cancela el desajuste determinista entre el
equilibrio calculado y el Poisson discreto del código:

    delta Phi_pert(r,t) = [Phi(r,t) - Phi(r,0)] - [Phi_ref(r,t) - Phi_ref(r,0)].

Uso:  python3 tools/delta_phi.py <corrida> [--referencia <corrida eps=0>] [--rmin 1 --rmax 15]
Imprime max_r |delta Phi| por instantánea y escribe <corrida>/delta_phi.npz (t, r, dphi).
"""
import argparse, os
import numpy as np, h5py


def serie(d):
    h = h5py.File(os.path.join(d, 'vlasov_output.h5'), 'r')
    pasos = sorted([g for g in h if g.startswith('step_')], key=lambda g: int(g.split('_')[1]))
    if 'potential' not in h[pasos[0]]:
        raise SystemExit(f'{d}: sin potencial en la malla (solo se guarda con autointeraction)')
    t = np.array([h[g].attrs['time'] for g in pasos])
    P = np.array([h[g]['potential'][()] for g in pasos])
    return t, h['grid/r'][()], P


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dir')
    ap.add_argument('--referencia')
    ap.add_argument('--rmin', type=float, default=0.0)
    ap.add_argument('--rmax', type=float, default=np.inf)
    a = ap.parse_args()

    t, r, P = serie(a.dir)
    dphi = P - P[0]
    if a.referencia:
        tr, rr, Pr = serie(a.referencia)
        if not (np.array_equal(r, rr) and np.allclose(t, tr)):
            raise SystemExit('la referencia debe tener la misma malla y los mismos tiempos')
        dphi = dphi - (Pr - Pr[0])
    m = (r >= a.rmin) & (r <= a.rmax)
    for ti, fila in zip(t, dphi):
        print(f'{ti:10.3f}  max_r|dPhi| = {np.max(np.abs(fila[m])):.3e}')
    np.savez(os.path.join(a.dir, 'delta_phi.npz'), t=t, r=r, dphi=dphi, referencia=a.referencia or '')


if __name__ == '__main__':
    main()
