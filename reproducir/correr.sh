#!/bin/bash
# Repite las simulaciones de las notas (docs/introduccion/vlasov_L_intro.tex).
#
#   reproducir/correr.sh <grupo> [patrón]      corre, una detrás de otra
#   reproducir/correr.sh <grupo> [patrón] -n   solo muestra los comandos
#
# <grupo>: una carpeta de reproducir/corridas (08_verificacion, 09_fondos, 10_energia,
# 11_equilibrio). [patrón]: opcional, solo los .par cuyo nombre lo contiene.
#
# Las corridas se ejecutan desde exe/ y escriben en exe/rep/<grupo>/<nombre>. Las de
# 11_equilibrio leen su estado inicial de exe/rep/11_equilibrio/ic_*.dat, que se genera
# aquí con tools/equilibrio_L.py si no existe. Nunca corre dos simulaciones a la vez.

set -u
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAIZ="$(dirname "$AQUI")"
grupo="${1:?uso: correr.sh <grupo> [patrón] [-n]}"
patron=""; seco=0
for a in "${@:2}"; do
  if [ "$a" = "-n" ]; then seco=1; else patron="$a"; fi
done
dir="$AQUI/corridas/$grupo"
[ -d "$dir" ] || { echo "no existe el grupo $grupo"; exit 1; }

cd "$RAIZ/exe" || exit 1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4} OMP_PLACES=cores OMP_PROC_BIND=close

ejecuta () { if [ $seco = 1 ]; then echo "  $*"; else "$@"; fi; }

if [ "$grupo" = "11_equilibrio" ]; then
  mkdir -p rep/11_equilibrio
  for e in 0 0.1; do
    nombre=ic_e$(echo $e | tr -d .)
    if [ ! -f "rep/11_equilibrio/$nombre.dat" ] || [ $seco = 1 ]; then
      ejecuta python3 "$RAIZ/tools/equilibrio_L.py" --a0 1e-2 --eps $e --nrc 200 --npc 20 --nlc 8 \
              --lminc 1.6 --lmaxc 2.4 --l0 2 --sl 0.2 --salida "rep/11_equilibrio/$nombre.dat"
    fi
  done
fi

for par in "$dir"/*"$patron"*.par; do
  [ -f "$par" ] || continue
  nombre=$(basename "$par" .par)
  destino=$(awk -F= '/^directory/ {gsub(/ /,"",$2); print $2}' "$par")
  if [ $seco = 1 ]; then
    echo "  ./VP_PIC $par    # -> exe/$destino"
    continue
  fi
  mkdir -p "$(dirname "$destino")"
  rm -rf "$destino"
  t0=$SECONDS
  if ./VP_PIC "$par" > "$destino.log" 2>&1; then
    echo "OK    $nombre ($((SECONDS-t0)) s)"
  else
    echo "FALLO $nombre"; tail -3 "$destino.log"
  fi
done
