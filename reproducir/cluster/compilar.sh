#!/bin/bash
# Compila VP_PIC (el código con L distribuida) en el cluster y lo prueba (probar.sh).
#
#   bash reproducir/cluster/compilar.sh
#
# Usa los módulos de vlasov-poisson_PIC/reproducir/cluster/entorno.sh. El Makefile compila a
# través de h5fc si está en el PATH y envuelve al compilador elegido; si no, usa las rutas
# VP_HDF5_INC y VP_HDF5_LIBS de ese archivo.

source "$(dirname "${BASH_SOURCE[0]}")/comun.sh"
set -e
cd "$RAIZ"

args=(FC="$VP_FC")
if [ -n "${VP_HDF5_INC:-}" ]; then
  args+=(HDF5_WRAPPER= "HDF5_INC=$VP_HDF5_INC" "HDF5_LIBS=$VP_HDF5_LIBS")
fi
echo "compilador: $(command -v "$VP_FC" || echo "$VP_FC no está en el PATH")"
echo "h5fc:       $(command -v h5fc || echo "no está; se usan VP_HDF5_INC y VP_HDF5_LIBS")"
echo "make ${args[*]}"
make clean
make "${args[@]}"
echo "bibliotecas de exe/VP_PIC que no se encuentran:"
ldd exe/VP_PIC | grep "not found" || echo "  ninguna"
bash "$RAIZ/reproducir/cluster/probar.sh"
