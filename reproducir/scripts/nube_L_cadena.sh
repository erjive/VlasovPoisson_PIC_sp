#!/bin/bash
# Corre y analiza tandas de nube_L.py, una detrás de otra y nunca dos a la vez:
#
#   reproducir/scripts/nube_L_cadena.sh primera resto convergencia
#
# Cada corrida se analiza al terminar, así que lo ya hecho sirve aunque la cadena se interrumpa, y
# al volver a lanzarla solo se hace lo que falta. El avance queda en exe/rep/12_nube/cadena.log,
# que termina con una línea FIN.
set -u
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAL="$AQUI/../../exe/rep/12_nube"
mkdir -p "$SAL"
for t in "$@"; do
  echo "== tanda $t: $(date '+%F %T')" >> "$SAL/cadena.log"
  python3 "$AQUI/nube_L.py" tanda "$t" >> "$SAL/cadena.log" 2>&1 || echo "FALLO en la tanda $t" >> "$SAL/cadena.log"
done
echo "FIN $(date '+%F %T')" >> "$SAL/cadena.log"
