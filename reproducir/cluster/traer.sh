#!/bin/bash
# Trae del cluster al portátil las medidas de las corridas de la nube hechas allá. En el portátil,
# desde la raíz del repositorio:
#
#   reproducir/cluster/traer.sh [-o origen] [-g] [-f] [-n] [<tanda o corrida> ...]
#
# origen   este repositorio en el cluster, usuario@máquina:/ruta (o la variable VP_CLUSTER_SP).
# -g       trae también bulto.npz, los mapas (Q, J) de todas las instantáneas (16 kB por partícula:
#          0.7 GB con 42000). Sin él vienen solo los archivos pequeños, que bastan para
#          "nube_L.py tabla": medidas.npz, capas.npz, campo.npz, resumen.npz, nucleo.npz,
#          params_usados.par, y el .ok, el .meta, el .log y la salida de Slurm de cada corrida.
# -f       pisa una corrida que ya exista en el portátil (por omisión se niega).
# -n       muestra lo que traería, sin copiar.
#
# Quedan en exe/rep/12_nube/, donde nube_L.py las espera. Una sola conexión para todas.

AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAIZ="$(cd "$AQUI/../.." && pwd)"
origen="${VP_CLUSTER_SP:-}"; grande=0; forzar=0; seco=""
while getopts "o:gfn" o; do
  case $o in
    o) origen=$OPTARG ;; g) grande=1 ;; f) forzar=1 ;; n) seco="--dry-run" ;;
    *) sed -n '2,16p' "${BASH_SOURCE[0]}"; exit 1 ;;
  esac
done
shift $((OPTIND-1))
[ -n "$origen" ] || { echo "falta el origen: -o usuario@máquina:/ruta/de/VlasovPoisson_PIC_sp, o la variable VP_CLUSTER_SP"; exit 1; }

cd "$RAIZ" || exit 1
corridas=$(python3 reproducir/scripts/nube_L.py lista "$@") || { echo "$corridas"; exit 1; }
destino="$RAIZ/exe/rep/12_nube"
mkdir -p "$destino"
filtros=()
for n in $corridas; do
  if [ -e "$destino/$n.ok" ] && [ $forzar = 0 ] && [ -z "$seco" ]; then
    echo "$n ya existe en el portátil: no se pisa (use -f)"; continue
  fi
  filtros+=(--include="/$n.ok" --include="/$n.meta" --include="/$n.log" --include="/slurm_${n}_*" --include="/$n/"
            --include="/$n/medidas.npz" --include="/$n/capas.npz" --include="/$n/campo.npz"
            --include="/$n/resumen.npz" --include="/$n/nucleo.npz" --include="/$n/params_usados.par")
  [ $grande = 1 ] && filtros+=(--include="/$n/bulto.npz")
done
[ ${#filtros[@]} -gt 0 ] || { echo "nada que traer"; exit 0; }
rsync -av $seco --prune-empty-dirs "${filtros[@]}" --exclude='*' "$origen/exe/rep/12_nube/" "$destino/"
rc=$?
echo "destino: $destino"
exit $rc
