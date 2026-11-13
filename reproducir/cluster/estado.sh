#!/bin/bash
# Estado de las corridas de la nube en el cluster. Desde la raíz del repositorio:
#
#   bash reproducir/cluster/estado.sh [<tanda o corrida> ...] > estado.txt 2>&1
#
# Imprime los trabajos del usuario en la cola y, de cada corrida, hasta dónde llegó (corrida,
# analizada, medida), cuánto tardó y su tamaño en disco. Al final, la tabla de medidas.

source "$(dirname "${BASH_SOURCE[0]}")/comun.sh"
cd "$RAIZ" || exit 1
echo "fecha: $(date)   máquina: $(hostname)   HEAD: $(git rev-parse --short HEAD 2>/dev/null)   vlasov-poisson_PIC: $(cd "$FIJO" && git rev-parse --short HEAD 2>/dev/null)"
if command -v squeue > /dev/null 2>&1; then
  echo "== trabajos en la cola"
  squeue -u "$USER" -o "%.10i %.12j %.8T %.10M %.6C %.20S %R" 2>/dev/null
fi
corridas=$(python3 reproducir/scripts/nube_L.py lista "$@") || { echo "$corridas"; exit 1; }
echo "== corridas"
for n in $corridas; do
  e="sin correr"
  [ -f "$NUBE/$n.ok" ] && e="corrida"
  [ -f "$NUBE/$n/bulto.npz" ] && e="analizada"
  [ -f "$NUBE/$n/resumen.npz" ] && e="medida"
  s=$(grep -h "^OK *$n " "$NUBE"/slurm_${n}_*.out 2>/dev/null | tail -1 | sed 's/^OK *//')
  d=$([ -d "$NUBE/$n" ] && du -sh "$NUBE/$n" | cut -f1)
  err=$(cat "$NUBE"/slurm_${n}_*.err 2>/dev/null | grep -v -i warn | tail -1)
  printf "%-12s %-10s %-28s %6s  %s\n" "$n" "$e" "${s:--}" "${d:--}" "$err"
done
echo "== modelo reducido"
logs=$(ls -t "$FIJO"/log.reducido*.out 2>/dev/null | head -6)
[ -n "$logs" ] || echo "sin salida todavía"
for ult in $logs; do
  echo "${ult#$FIJO/}:"
  grep -E "^(inicio|fin|== serie|OK |FALTA|medidas de)" "$ult" | tail -n 12
  tail -n 1 "$ult"
  err="${ult%.out}.err"
  [ -s "$err" ] && { echo "${err#$FIJO/} (últimas líneas):"; tail -n 5 "$err"; }
done
echo "== tabla"
python3 reproducir/scripts/nube_L.py tabla 2>&1 | tail -n +6
