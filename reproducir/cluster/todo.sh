#!/bin/bash
# Toda la serie de la nube con L distribuida en el cluster, con un solo comando. En el nodo de
# entrada, desde la raíz de este repositorio, con vlasov-poisson_PIC en la carpeta hermana:
#
#   bash reproducir/cluster/todo.sh [-n | -t] > todo.txt 2>&1
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
# -n, o una máquina sin sbatch, hace 1 y 2 y solo muestra los sbatch de 3 y 4. -t hace 1 y 2 y, en
# vez de la serie, envía la prueba de diez minutos (prueba.slurm). Los trabajos van a la partición
# VP_PARTICION de entorno.sh. Lo hecho o ya en la cola se salta, así que se puede volver a correr.
# Después: bash reproducir/cluster/estado.sh.

source "$(dirname "${BASH_SOURCE[0]}")/comun.sh"
seco=0; solo_prueba=0
case "${1:-}" in -n) seco=1 ;; -t) solo_prueba=1 ;; "") ;; *) sed -n '2,21p' "${BASH_SOURCE[0]}"; exit 1 ;; esac
if ! command -v sbatch > /dev/null 2>&1; then
  [ $seco = 0 ] && echo "no hay sbatch en esta máquina: solo se muestran los comandos"
  seco=1
fi
opc=(); [ $seco = 1 ] && opc=(-n)
cd "$RAIZ" || exit 1
# Si Python no carga nube_L.py se compila de todos modos, para saber también si eso funciona.
err_py=$(mktemp)
if corridas=$(python3 reproducir/scripts/nube_L.py lista 2> "$err_py"); then py_bien=1; else py_bien=0; corridas=""; fi
echo "inicio: $(date)   máquina: $(hostname)   partición: $VP_PARTICION   HEAD: $(git rev-parse --short HEAD)   vlasov-poisson_PIC: $(cd "$FIJO" && git rev-parse --short HEAD)"

echo; echo "== 1. compilación y prueba"
encola=""
if command -v squeue > /dev/null 2>&1; then
  encola=$(squeue -h -u "$USER" -o %j 2>/dev/null | grep -x -F -f <(printf '%s\n' $corridas reducido prueba))
fi
if [ -n "$encola" ]; then
  echo "hay trabajos de esta serie en la cola ($(echo $encola)): no se recompila"
else
  bash reproducir/cluster/compilar.sh || comp_mal=1
fi
if [ $py_bien = 0 ]; then
  echo; echo "ALTO: $(command -v python3 || echo 'no hay python3') ($(python3 --version 2>&1)) no carga nube_L.py:"
  tail -n 6 "$err_py"
  echo "Hacen falta numpy y scipy: ver la sección Python de diagnostico.txt y los módulos de entorno.sh."
  rm -f "$err_py"; exit 1
fi
rm -f "$err_py"
[ -z "${comp_mal:-}" ] || { echo "ALTO: falló la compilación o la prueba"; exit 1; }

echo; echo "== 2. datos iniciales"
faltan=$(for n in $corridas; do [ -f "$NUBE/ic/$n.dat" ] || echo $n; done)
if [ -n "$faltan" ]; then
  python3 reproducir/scripts/nube_L.py preparar || { echo "ALTO: falló nube_L.py preparar"; exit 1; }
fi
for n in $corridas; do [ -f "$NUBE/ic/$n.dat" ] || { echo "ALTO: falta el dato inicial de $n"; exit 1; }; done
echo "están los datos iniciales de las $(echo $corridas | wc -w) corridas"

if [ $solo_prueba = 1 ]; then
  echo; echo "== prueba de diez minutos en un nodo de cálculo"
  if [ $seco = 1 ]; then
    echo "sbatch -p $VP_PARTICION reproducir/cluster/prueba.slurm"
  else
    sbatch -p "$VP_PARTICION" reproducir/cluster/prueba.slurm && echo "cuando termine, el resultado está en log.prueba.<trabajo>.out"
  fi
  echo; echo "fin: $(date)"
  exit 0
fi

echo; echo "== 3. modelo reducido (series Lb y control)"
if (cd "$FIJO/reproducir/scripts" && python3 -c "
import os, sys, reducido as R
sys.exit(0 if any(not os.path.exists(os.path.join(R.DATOS, n + '_capas.npz')) for s in ('Lb', 'control') for n, *_ in R.SERIES[s]) else 1)"); then
  if [ $seco = 0 ] && squeue -h -u "$USER" -n reducido 2>/dev/null | grep -q .; then
    echo "ya hay un trabajo reducido en la cola"
  elif [ $seco = 1 ]; then
    echo "(cd $FIJO && sbatch -p $VP_PARTICION reproducir/cluster/reducido.slurm Lb control)"
  else
    (cd "$FIJO" && sbatch -p "$VP_PARTICION" reproducir/cluster/reducido.slurm Lb control)
  fi
else
  echo "ya están hechas"
fi

echo; echo "== 4. corridas de la nube"
reproducir/cluster/enviar.sh "${opc[@]}" primera resto convergencia
echo; echo "fin: $(date)"
