#!/bin/bash
# Toda la serie de la nube con L distribuida en el cluster, con un solo comando. En el nodo de
# entrada, desde la raíz de este repositorio, con vlasov-poisson_PIC en la carpeta hermana:
#
#   bash reproducir/cluster/todo.sh [-n] > todo.txt 2>&1
#
#  1. compila VP_PIC y lo prueba (compilar.sh); si la prueba falla, se detiene. Si hay trabajos de
#     esta serie en la cola no recompila: "make clean" borraría el ejecutable que van a usar;
#  2. escribe los datos iniciales que falten (nube_L.py preparar, un minuto);
#  3. envía el modelo reducido de las series Lb y control (reducido.slurm de vlasov-poisson_PIC,
#     16 CPUs, unas 2.5 h). Va primero porque es el trabajo más largo y la cola arranca los
#     trabajos en el orden de llegada;
#  4. envía las catorce corridas (enviar.sh primera resto convergencia): un trabajo de 8 CPUs y
#     8 GB por corrida, que la corre, la analiza y la mide.
#
# -n, o una máquina sin sbatch, hace 1 y 2 y solo muestra los sbatch de 3 y 4. Lo hecho o ya en la
# cola se salta, así que se puede volver a correr. Después: bash reproducir/cluster/estado.sh.

source "$(dirname "${BASH_SOURCE[0]}")/comun.sh"
seco=0; [ "${1:-}" = "-n" ] && seco=1
if ! command -v sbatch > /dev/null 2>&1; then
  [ $seco = 0 ] && echo "no hay sbatch en esta máquina: solo se muestran los comandos"
  seco=1
fi
opc=(); [ $seco = 1 ] && opc=(-n)
cd "$RAIZ" || exit 1
corridas=$(python3 reproducir/scripts/nube_L.py lista) || { echo "ALTO: nube_L.py no carga"; echo "$corridas"; exit 1; }
echo "inicio: $(date)   máquina: $(hostname)   HEAD: $(git rev-parse --short HEAD)   vlasov-poisson_PIC: $(cd "$FIJO" && git rev-parse --short HEAD)"

echo; echo "== 1. compilación y prueba"
encola=""
if command -v squeue > /dev/null 2>&1; then
  encola=$(squeue -h -u "$USER" -o %j 2>/dev/null | grep -x -F -f <(printf '%s\n' $corridas reducido))
fi
if [ -n "$encola" ]; then
  echo "hay trabajos de esta serie en la cola ($(echo $encola)): no se recompila"
else
  bash reproducir/cluster/compilar.sh || { echo "ALTO: falló la compilación o la prueba"; exit 1; }
fi

echo; echo "== 2. datos iniciales"
faltan=$(for n in $corridas; do [ -f "$NUBE/ic/$n.dat" ] || echo $n; done)
if [ -n "$faltan" ]; then
  python3 reproducir/scripts/nube_L.py preparar || { echo "ALTO: falló nube_L.py preparar"; exit 1; }
fi
for n in $corridas; do [ -f "$NUBE/ic/$n.dat" ] || { echo "ALTO: falta el dato inicial de $n"; exit 1; }; done
echo "están los datos iniciales de las $(echo $corridas | wc -w) corridas"

echo; echo "== 3. modelo reducido (series Lb y control)"
if (cd "$FIJO/reproducir/scripts" && python3 -c "
import os, sys, reducido as R
sys.exit(0 if any(not os.path.exists(os.path.join(R.DATOS, n + '_capas.npz')) for s in ('Lb', 'control') for n, *_ in R.SERIES[s]) else 1)"); then
  if [ $seco = 0 ] && squeue -h -u "$USER" -n reducido 2>/dev/null | grep -q .; then
    echo "ya hay un trabajo reducido en la cola"
  elif [ $seco = 1 ]; then
    echo "(cd $FIJO && sbatch reproducir/cluster/reducido.slurm Lb control)"
  else
    (cd "$FIJO" && sbatch reproducir/cluster/reducido.slurm Lb control)
  fi
else
  echo "ya están hechas"
fi

echo; echo "== 4. corridas de la nube"
reproducir/cluster/enviar.sh "${opc[@]}" primera resto convergencia
echo; echo "fin: $(date)"
