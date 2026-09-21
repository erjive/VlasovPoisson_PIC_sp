#!/bin/bash
# Compara dos ejecutables bit a bit en tres casos cortos (sin autogravedad,
# con autogravedad y leapfrog, con autogravedad y yoshida4). Sirve para
# demostrar que un cambio que se declara neutro (refactor, optimizacion,
# limpieza) no cambia ningun bit de la salida.
#
#   verificacion/regresion.sh <VP_PIC_referencia> <VP_PIC_nuevo>
#
# Mismo numero de hilos para ambos (las sumas de energia y h_k son
# deterministas para un numero de hilos dado, no entre numeros distintos).
set -u
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
A="$(realpath "$1")"; B="$(realpath "$2")"
W="${VP_VERIF:-/tmp/vp_verificacion}/regresion"; rm -rf "$W"; mkdir -p "$W"; cd "$W"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
base="dr = 0.1
Nrc = 30
Npc = 30
Nlc = 6
courant = 0.5
dt_switch = fix
Nt = 20
rmin = 0.0
rmax = 20.0
rminc = 0.5
rmaxc = 10.0
pminc = -0.8
pmaxc = 0.8
lminc = 0.0
lmaxc = 2.4
pmax = 2.0
time_output = 100
spatial_output = 5
output_format = hdf5
a0 = 0.1
l0 = 1.0
sr = 0.1
sp = 0.3
sl = 0.8
state = aa
cutoff = 0.001
bsplineorder = 2
BGtype = Isochrone"
distintos=0
for caso in "fijo leapfrog .false." "sg_lf leapfrog .true." "sg_y4 yoshida4 .true."; do
  set -- $caso
  for v in A B; do
    exe=$A; [ $v = B ] && exe=$B
    printf '%s\nintegrator = %s\nautointeraction = %s\ndirectory = %s_%s\n' "$base" "$2" "$3" "$1" "$v" > "$1_$v.par"
    "$exe" "$1_$v.par" > "$1_$v.log" 2>&1 || { echo "  $1: el ejecutable $v fallo"; tail -3 "$1_$v.log"; }
  done
  printf '  %-6s' "$1"
  if python3 "$AQUI/py/h5igual.py" "$1_A/vlasov_output.h5" "$1_B/vlasov_output.h5"; then :; else distintos=$((distintos+1)); fi
  for f in hk1_complex.tl hk2_complex.tl; do cmp -s "$1_A/$f" "$1_B/$f" || { echo "         $f distinto"; distintos=$((distintos+1)); }; done
done
[ $distintos = 0 ] && echo "IDENTICOS bit a bit" || echo "DIFIEREN ($distintos)"
exit $distintos
