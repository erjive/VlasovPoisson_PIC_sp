#!/bin/bash
#SBATCH --job-name=prueba
#SBATCH --output=log.%x.%j.out
#SBATCH --error=log.%x.%j.err
#SBATCH --partition=olin
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G

# Prueba de unos 10 minutos en un nodo de cálculo, con los mismos recursos que un trabajo de la
# serie, antes de enviar la serie. La envía "todo.sh -t" (con la partición de entorno.sh), después
# de compilar y escribir los datos iniciales. A mano, desde la raíz de este repositorio:
#
#   sbatch [-p partición] reproducir/cluster/prueba.slurm
#
#  1. el nodo: sistema, procesador, memoria, CPUs asignadas, núcleos físicos, módulos y versiones;
#  2. el ejecutable frente al portátil, y el Python del nodo (probar.sh);
#  3. la cadena completa con una corrida corta de 12600 partículas y seis valores de L: correr,
#     analizar y medir (nube_L.py prueba);
#  4. el modelo reducido: compila su biblioteca en C en el nodo y pasa sus pruebas (reducido_pruebas.py);
#  5. el costo: 2000 pasos de cada una de las catorce corridas (nube_L.py tiempo).
#
# Termina con un resumen PASA/FALLA de cada paso. Escribe en log.prueba.<trabajo>.out. No deja
# archivos de corridas: los de las pruebas se borran.

VP_RAIZ="${VP_RAIZ:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "$VP_RAIZ/reproducir/cluster/comun.sh"
hilos=$(vp_nucleos)
[ -n "${SLURM_CPUS_PER_TASK:-}" ] && [ "$hilos" -gt "$SLURM_CPUS_PER_TASK" ] && hilos=$SLURM_CPUS_PER_TASK
export VP_HILOS="${VP_HILOS:-$hilos}" VP_PROCESOS="${VP_PROCESOS:-$hilos}"
export OMP_NUM_THREADS="$VP_HILOS" OMP_PLACES=cores OMP_PROC_BIND=close    # para el modelo reducido
ulimit -s unlimited 2>/dev/null
export OMP_STACKSIZE=256M
cd "$RAIZ" || exit 1

resumen=()
paso () {               # paso <título> <comando ...>
  local titulo="$1" t0=$SECONDS; shift
  echo; echo "== $titulo"
  if "$@"; then resumen+=("PASA  $titulo ($((SECONDS-t0)) s)"); else resumen+=("FALLA $titulo ($((SECONDS-t0)) s)"); fi
}

echo "== 1. el nodo"
echo "máquina: $(hostname)   trabajo: ${SLURM_JOB_ID:-sin Slurm}   partición: ${SLURM_JOB_PARTITION:-?}   $(date)"
grep -h PRETTY_NAME /etc/os-release 2>/dev/null
LC_ALL=C lscpu | grep -E "^Model name|^CPU\(s\)|^Thread|^Core|^Socket"
free -g | head -2
echo "CPUs asignadas: $(taskset -cp $$ 2>/dev/null | sed 's/.*: *//')   núcleos físicos: $(vp_nucleos)   hilos de las corridas: $VP_HILOS   procesos del análisis: $VP_PROCESOS"
type module > /dev/null 2>&1 && module list 2>&1 | head -12
echo "$(gcc --version 2>&1 | head -1)  |  $(gfortran --version 2>&1 | head -1)  |  $(python3 --version 2>&1)"
echo "repositorios: $(git rev-parse --short HEAD) y vlasov-poisson_PIC $(cd "$FIJO" && git rev-parse --short HEAD)"

paso "2. ejecutable y Python (probar.sh)" bash reproducir/cluster/probar.sh
paso "3. cadena corta: correr, analizar y medir (nube_L.py prueba)" python3 reproducir/scripts/nube_L.py prueba
paso "4. modelo reducido (reducido_pruebas.py)" bash -c "cd '$FIJO/reproducir/scripts' && python3 reducido_pruebas.py"
paso "5. costo de las catorce corridas (nube_L.py tiempo)" python3 reproducir/scripts/nube_L.py tiempo

echo; echo "===== Resumen ($(hostname), $(date))"
printf '%s\n' "${resumen[@]}"
! printf '%s\n' "${resumen[@]}" | grep -q "^FALLA"
