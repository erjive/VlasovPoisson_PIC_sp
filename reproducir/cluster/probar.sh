#!/bin/bash
# Prueba corta, de segundos, de que el cluster puede hacer la serie de la nube con L distribuida:
#
#  1. el ejecutable: 10240 nodos en (Q, J, L) con autogravedad, 400 pasos, 2 hilos; |h_0| ... |h_4|
#     del final frente a los del portátil (gfortran 13.3). Otro compilador no da los mismos bits;
#     se pide acuerdo relativo a 1e-10.
#  2. Python: numpy y scipy, el repositorio vlasov-poisson_PIC y la lista de corridas.
#  3. los datos iniciales de las corridas (exe/rep/12_nube/ic/), que se copian del portátil.
#
# La corre compilar.sh; se puede repetir sola: bash reproducir/cluster/probar.sh

source "$(dirname "${BASH_SOURCE[0]}")/comun.sh"
cd "$RAIZ/exe" || exit 1
ulimit -s unlimited 2>/dev/null

# Última línea de hk1.tl (t y |h_0| ... |h_4|) de esta misma corrida en el portátil.
REF="1.0000000000000075E+01  1.0704919835348240E-06  1.0462378148497553E-06  9.7672029653178660E-07  8.7096062981538485E-07  7.4183123634273383E-07"

fallos=0
rm -rf rep/cluster_prueba
OMP_NUM_THREADS=2 OMP_PLACES=cores OMP_PROC_BIND=close \
  ./VP_PIC ../reproducir/video/df_gauss_sg.par Nrc=40 Npc=32 Nlc=8 Nt=400 spatial_output=100 field_output=400 \
           time_output=1000000 output_format=raw directory=rep/cluster_prueba > rep_cluster_prueba.log 2>&1
if [ -s rep/cluster_prueba/hk1.tl ] && [ -s rep/cluster_prueba/vlasov_output.raw ]; then
  tail -1 rep/cluster_prueba/hk1.tl | awk -v ref="$REF" '
    { n = split(ref, r, " "); peor = 0
      for (i = 1; i <= n; i++) { d = ($i - r[i])/(r[i] == 0 ? 1 : r[i]); if (d < 0) d = -d; if (d > peor) peor = d }
      printf "  ejecutable: diferencia relativa máxima con el portátil %.1e\n", peor
      exit (peor < 1e-10 ? 0 : 1) }' || fallos=$((fallos+1))
  grep -h "Time step fixed" rep_cluster_prueba.log
else
  echo "  ejecutable: la corrida falló"; tail -5 rep_cluster_prueba.log; fallos=$((fallos+1))
fi
rm -rf rep/cluster_prueba rep_cluster_prueba.log

cd "$RAIZ" || exit 1
if corridas=$(python3 reproducir/scripts/nube_L.py lista 2>&1); then
  echo "  python: $(python3 -c 'import numpy, scipy, sys; print("python", sys.version.split()[0], "numpy", numpy.__version__, "scipy", scipy.__version__)')"
  faltan=""
  for n in $corridas; do [ -f "$NUBE/ic/$n.dat" ] || faltan="$faltan $n"; done
  if [ -z "$faltan" ]; then
    echo "  datos iniciales: están los de las $(echo $corridas | wc -w) corridas"
  else
    echo "  datos iniciales: faltan los de$faltan"
    echo "     cópielos del portátil (rsync de exe/rep/12_nube/ic/) o genérelos aquí con"
    echo "     python3 reproducir/scripts/nube_L.py preparar"
  fi
else
  echo "  python: falla la carga de nube_L.py:"; echo "$corridas" | tail -3; fallos=$((fallos+1))
fi
[ $fallos = 0 ] && echo "PASA la prueba" || echo "FALLA la prueba ($fallos)"
exit $fallos
