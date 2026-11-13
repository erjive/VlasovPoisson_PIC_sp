#!/bin/bash
# Envía a Slurm corridas de la nube con L distribuida, un trabajo por corrida. Desde la raíz del
# repositorio:
#
#   reproducir/cluster/enviar.sh [-c cpus] [-t hh:mm:ss] [-p partición] [-m memoria] [-n] <tanda o corrida> ...
#
# Las tandas son primera, resto y convergencia (TANDAS de nube_L.py); sin nombres, todas. Hacen
# falta el .par de cada corrida (versionado) y su dato inicial en exe/rep/12_nube/ic/. -n muestra
# los comandos sin enviar. Una corrida ya medida, o con un trabajo del mismo nombre en la cola,
# se salta. Por omisión -c 8, que en "olin" son cuatro núcleos: las corridas no ganan con más, y
# así caben seis a la vez en el nodo. La partición por omisión es VP_PARTICION, de entorno.sh.

source "$(dirname "${BASH_SOURCE[0]}")/comun.sh"
cpus=8; tiempo=24:00:00; part="$VP_PARTICION"; mem=8G; seco=0
while getopts "c:t:p:m:n" o; do
  case $o in
    c) cpus=$OPTARG ;; t) tiempo=$OPTARG ;; p) part=$OPTARG ;; m) mem=$OPTARG ;; n) seco=1 ;;
    *) sed -n '2,13p' "${BASH_SOURCE[0]}"; exit 1 ;;
  esac
done
shift $((OPTIND-1))

cd "$RAIZ" || exit 1
corridas=$(python3 reproducir/scripts/nube_L.py lista "$@") || { echo "$corridas"; exit 1; }
mkdir -p "$NUBE"
for n in $corridas; do
  [ -f "reproducir/corridas/12_nube/$n.par" ] || { echo "$n: falta su .par"; continue; }
  [ -f "$NUBE/ic/$n.dat" ] || { echo "$n: falta el dato inicial exe/rep/12_nube/ic/$n.dat"; continue; }
  [ -f "$NUBE/$n/resumen.npz" ] && { echo "$n: ya está medida"; continue; }
  if command -v squeue > /dev/null 2>&1 && squeue -h -u "$USER" -n "$n" 2>/dev/null | grep -q .; then
    echo "$n: ya hay un trabajo con ese nombre en la cola"; continue
  fi
  cmd=(sbatch -J "$n" -c "$cpus" -t "$tiempo" -p "$part" --mem="$mem"
       -o "exe/rep/12_nube/slurm_%x_%j.out" -e "exe/rep/12_nube/slurm_%x_%j.err"
       --export="ALL,VP_RAIZ=$RAIZ,VP_PIC_FIJO=$FIJO" reproducir/cluster/nube.slurm "$n")
  if [ $seco = 1 ]; then echo "${cmd[*]}"; else "${cmd[@]}"; fi
done
