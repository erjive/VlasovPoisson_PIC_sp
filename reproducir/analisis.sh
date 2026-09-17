#!/bin/bash
# Análisis que alimentan las figuras (después de reproducir/correr.sh con los cuatro grupos).
# Desde la raíz del repositorio: reproducir/analisis.sh; luego reproducir/figuras/generar_figuras.py
set -eu
RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$RAIZ/exe/rep"
T="$RAIZ/tools"
for d in gauss_y4 gauss_an bim_an bim_y4 king_an king_y4; do python3 $T/hk_exacto.py 08_verificacion/$d; done
python3 $T/hk_exacto.py 08_verificacion/spi_L0_an --nl 4
python3 $T/hk_exacto.py 08_verificacion/spi_L16_an --nj 400 --nl 200 --nq 256
E=11_equilibrio
for d in eq_e0 eq_e01 iso_aq; do
  python3 $T/hk_numerico.py $E/$d > $E/$d.hkn.log
  python3 $T/delta_phi.py $E/$d --rmin 1 --rmax 15 > $E/$d.dphi.log
done
for d in eq_e0 eq_e01; do
  python3 $T/hk_numerico.py $E/$d --equilibrio $E/ic_e0_equilibrio.npz > $E/$d.hkn_eq.log
done
