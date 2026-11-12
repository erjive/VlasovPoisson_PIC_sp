#!/bin/bash
#SBATCH --job-name=nube
#SBATCH --output=log.%x.%j.out
#SBATCH --error=log.%x.%j.err
#SBATCH --partition=olin
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G

# Corridas de la nube con L distribuida (reproducir/scripts/nube_L.py): cada una se corre, se
# analiza y se mide, una detrás de otra. Desde la raíz del repositorio:
#
#   sbatch -J <nombre> reproducir/cluster/nube.slurm <corrida> [<corrida> ...]
#
# enviar.sh manda un trabajo por corrida. La memoria se pide porque la partición da 2 GB por
# omisión y el análisis de una corrida de 63000 partículas guarda en memoria unos 3 GB. Lo ya
# hecho se salta, así que un trabajo interrumpido se puede volver a enviar.
# Fuera de Slurm corre igual: bash reproducir/cluster/nube.slurm <corrida>.

VP_RAIZ="${VP_RAIZ:-${SLURM_SUBMIT_DIR:-$PWD}}"
source "$VP_RAIZ/reproducir/cluster/comun.sh"
[ -x "$RAIZ/exe/VP_PIC" ] || { echo "no hay $RAIZ/exe/VP_PIC: compile con reproducir/cluster/compilar.sh"; exit 1; }

# Un hilo por núcleo físico asignado (vp_nucleos, del entorno), sin pasar de las CPUs pedidas.
# El análisis reparte las instantáneas entre el mismo número de procesos.
hilos=$(vp_nucleos)
[ -n "${SLURM_CPUS_PER_TASK:-}" ] && [ "$hilos" -gt "$SLURM_CPUS_PER_TASK" ] && hilos=$SLURM_CPUS_PER_TASK
export VP_HILOS="${VP_HILOS:-$hilos}" VP_PROCESOS="${VP_PROCESOS:-$hilos}"
ulimit -s unlimited 2>/dev/null
export OMP_STACKSIZE=256M

cd "$RAIZ" || exit 1
echo "máquina: $(hostname)   trabajo: ${SLURM_JOB_ID:-sin Slurm}   hilos: $VP_HILOS   procesos del análisis: $VP_PROCESOS"
echo "inicio: $(date)"
rc=0
for n in "$@"; do
  python3 reproducir/scripts/nube_L.py correr --solo "$n" &&
  python3 reproducir/scripts/nube_L.py analizar --solo "$n" &&
  python3 reproducir/scripts/nube_L.py medir --solo "$n" || rc=1
  [ -f "$NUBE/$n/resumen.npz" ] || { echo "FALTA $n"; rc=1; }
done
echo "fin: $(date)   código de salida $rc"
exit $rc
