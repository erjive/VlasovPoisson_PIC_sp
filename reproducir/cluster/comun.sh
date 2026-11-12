# Lo que comparten los guiones de esta carpeta: dónde están los dos repositorios y el entorno del
# cluster. El entorno (módulos, compilador, núcleos físicos) es el de vlasov-poisson_PIC,
# reproducir/cluster/entorno.sh: es el único archivo que hay que editar al cambiar de máquina.
# vlasov-poisson_PIC debe estar junto a este repositorio, o su ruta en la variable VP_PIC_FIJO.

RAIZ="${VP_RAIZ:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
FIJO="${VP_PIC_FIJO:-$(dirname "$RAIZ")/vlasov-poisson_PIC}"
if [ ! -f "$FIJO/reproducir/cluster/entorno.sh" ]; then
  echo "no encuentro vlasov-poisson_PIC en $FIJO: clónelo junto a este repositorio o fije VP_PIC_FIJO"
  exit 1
fi
export VP_PIC_FIJO="$FIJO"
source "$FIJO/reproducir/cluster/entorno.sh"
NUBE="$RAIZ/exe/rep/12_nube"
