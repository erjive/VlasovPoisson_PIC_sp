"""h_k en las variables ángulo-acción verdaderas, a partir de las instantáneas HDF5 de una corrida.

analysish.f90 usa el mapa analítico del isócrono aun con autogravedad. Con masa propia eso
deja en h_k una meseta estática proporcional a la masa, que no es un modo (vlasov-poisson_PIC,
notas, secciones 9 y 10). Aquí, para cada instantánea, el potencial propio es la tabla
`potential` menos el isócrono, y el mapa numérico de aa_numerico_L da (Q,J) para cada
partícula con su L. h_k usa la función de prueba n de la corrida (params_usados.par):

    h_k = 8 pi^2 drc dpc dlc Sum_j f_j L_j B(J_j) C(L_j) a_k exp(-ik Q_j),

igual que el código; drc dpc dlc se reconstruye de los parámetros. Sin autogravedad
coincide con hk{n}_complex.tl hasta la precisión del mapa.

Uso:  python3 tools/hk_numerico.py <directorio> [--fn 1] [--cada 1]
Escribe <directorio>/hk{n}_numerico.tl con t, Re h_0, Im h_0, ..., Re h_4, Im h_4.
"""
import argparse, os, sys
import numpy as np, h5py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aa_numerico_L import MapaAA, phi_iso
from hk_exacto import leer_par, a_k


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dir')
    ap.add_argument('--fn', type=int, default=1, choices=[1, 2])
    ap.add_argument('--cada', type=int, default=1, help='usar una de cada n instantáneas')
    a = ap.parse_args()

    raw = leer_par(os.path.join(a.dir, 'params_usados.par'))
    P = {}
    for k, v in raw.items():
        try:
            P[k] = float(v)
        except ValueError:
            pass
    if raw['bgtype'] != 'Isochrone':
        sys.exit('solo para BGtype=Isochrone (el mapa suma el isócrono a la tabla)')
    sg = raw['autointeraction'].lower() in ('.true.', 'true', 't')
    s = str(a.fn)
    drc = (P['rmaxc'] - P['rminc'])/P['nrc']
    dpc = (P['pmaxc'] - P['pminc'])/P['npc']
    dlc = (P['lmaxc'] - P['lminc'])/P['nlc']
    ak = a_k(P['sq'+s])

    h = h5py.File(os.path.join(a.dir, 'vlasov_output.h5'), 'r')
    r_malla = h['grid/r'][()]
    pasos = sorted([g for g in h if g.startswith('step_')], key=lambda g: int(g.split('_')[1]))[::a.cada]
    filas = []
    for g in pasos:
        G = h[g]
        r, p, L, fl = G['r_part'][()], G['p_part'][()], G['l_part'][()], G['fl'][()]
        if sg:
            m = MapaAA(r_malla, G['potential'][()] - phi_iso(r_malla))
        else:
            m = MapaAA()
        E = 0.5*p**2 + m.phi_ef(r, L)
        lig = E < 0
        Q, J, _ = m(r[lig], p[lig], L[lig])
        w = fl[lig]*J**2*np.exp(-(J - P['j'+s])**2/P['sj'+s]**2)*np.exp(-(L[lig] - P['lt'+s])**2/P['slt'+s]**2)
        hk = 8*np.pi**2*drc*dpc*dlc*np.array([ak[k]*np.sum(w*np.exp(-1j*k*Q)) for k in range(5)])
        filas.append([G.attrs['time']] + [x for z in hk for x in (z.real, z.imag)])
        print(f"t={G.attrs['time']:9.3f}  |h_0|={abs(hk[0]):.6e}  |h_1|={abs(hk[1]):.6e}  no ligadas={np.sum(~lig)}",
              flush=True)
    salida = os.path.join(a.dir, f'hk{a.fn}_numerico.tl')
    np.savetxt(salida, np.array(filas), fmt='%24.16e')
    print('escrito', salida)


if __name__ == '__main__':
    main()
