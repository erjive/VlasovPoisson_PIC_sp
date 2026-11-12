#!/bin/bash
# Trae del cluster al portátil las medidas de las corridas de la nube hechas allá. En el portátil,
# desde la raíz del repositorio:
#
#   reproducir/cluster/traer.sh [-o origen] [-r] [-g] [-f] [-n] [<tanda o corrida> ...]
#
# origen   este repositorio en el cluster, usuario@máquina:/ruta (o la variable VP_CLUSTER_SP).
# -r       trae también las medidas del modelo reducido de las series Lb y control
#          (<corrida>_diag.npz y <corrida>_capas.npz), de la copia de vlasov-poisson_PIC que está
#          junto a esta en el cluster (o en la variable VP_CLUSTER_FIJO) a la del portátil. Nunca
#          pisa un archivo del portátil.
# -g       trae también los archivos grandes: bulto.npz de cada corrida del PIC, los mapas (Q, J)
#          de todas las instantáneas (16 kB por partícula: 0.7 GB con 42000), y con -r las corridas
#          del modelo (8 kB por partícula). Sin él vienen solo los archivos pequeños, que bastan
#          para "nube_L.py tabla": medidas.npz, capas.npz, campo.npz, resumen.npz, nucleo.npz,
#          params_usados.par, y el .ok, el .meta, el .log y la salida de Slurm de cada corrida.
# -f       pisa una corrida del PIC que ya exista en el portátil (por omisión se niega).
# -n       muestra lo que traería, sin copiar.
#
# Las corridas del PIC quedan en exe/rep/12_nube/, donde nube_L.py las espera.

AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAIZ="$(cd "$AQUI/../.." && pwd)"
origen="${VP_CLUSTER_SP:-}"; modelo=0; grande=0; forzar=0; seco=""
while getopts "o:rgfn" o; do
  case $o in
    o) origen=$OPTARG ;; r) modelo=1 ;; g) grande=1 ;; f) forzar=1 ;; n) seco="--dry-run" ;;
    *) sed -n '2,22p' "${BASH_SOURCE[0]}"; exit 1 ;;
  esac
done
shift $((OPTIND-1))
[ -n "$origen" ] || { echo "falta el origen: -o usuario@máquina:/ruta/de/VlasovPoisson_PIC_sp, o la variable VP_CLUSTER_SP"; exit 1; }
origen="${origen%/}"

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
rc=0
if [ ${#filtros[@]} -gt 0 ]; then
  rsync -av $seco --prune-empty-dirs "${filtros[@]}" --exclude='*' "$origen/exe/rep/12_nube/" "$destino/" || rc=1
  echo "destino: $destino"
else
  echo "del PIC no hay nada que traer"
fi

if [ $modelo = 1 ]; then
  fijo="${VP_PIC_FIJO:-$(dirname "$RAIZ")/vlasov-poisson_PIC}"
  origen_fijo="${VP_CLUSTER_FIJO:-${origen%/*}/vlasov-poisson_PIC}"
  nombres=$(cd "$fijo/reproducir/scripts" && python3 -c "
import reducido as R
print(' '.join(n for s in ('Lb', 'control') for n, *_ in R.SERIES[s]))") || { echo "no pude leer las series de $fijo"; exit 1; }
  fm=()
  for n in $nombres; do
    fm+=(--include="/${n}_diag.npz" --include="/${n}_capas.npz")
    [ $grande = 1 ] && fm+=(--include="/$n.npz")
  done
  mkdir -p "$fijo/exe/demo_eta/reducido"
  rsync -av $seco --ignore-existing "${fm[@]}" --exclude='*' "$origen_fijo/exe/demo_eta/reducido/" "$fijo/exe/demo_eta/reducido/" || rc=1
  echo "destino del modelo: $fijo/exe/demo_eta/reducido"
fi
exit $rc
