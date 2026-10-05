#!/bin/bash
# Pruebas de verificacion de VP_PIC (AUDITORIA_CIENTIFICA_2026-09-21.md).
#
#   verificacion/correr.sh [--rapido] [VP_PIC]
#
# Compila las pruebas unitarias contra objs/ (hace falta haber corrido make)
# y las corre una detras de otra, nunca dos a la vez. --rapido omite las de
# integracion (I1 a I6), que corren el ejecutable completo (~20 s).
# Escribe en ${VP_VERIF:-/tmp/vp_verificacion}. Codigo de salida: numero de
# pruebas que fallan. U7, U8, I3 e I4 son de AUDITORIA_FISICA_2026-09-21.md.
set -u
AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAIZ="$(dirname "$AQUI")"
rapido=0; exe="$RAIZ/exe/VP_PIC"
for a in "$@"; do
  if [ "$a" = "--rapido" ]; then rapido=1; else exe="$(realpath "$a")"; fi
done
W="${VP_VERIF:-/tmp/vp_verificacion}"; mkdir -p "$W"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
FC=$(command -v h5fc || command -v gfortran)
FL="-O3 -ffree-form -ffree-line-length-none -fopenmp"
OBJS=$(ls "$RAIZ"/objs/*.o | grep -v '/main.o$')
fallos=0
resumen=()
anota () { resumen+=("$1"); case "$1" in FALLA*) fallos=$((fallos+1));; esac; }

for t in t_forma t_origen t_poisson_nodal t_poisson_orden t_fondos t_aa t_r0 t_masa; do
  $FC $FL -I"$RAIZ/objs" "$AQUI/src/$t.f90" $OBJS -o "$W/$t" || { anota "FALLA $t (no compila)"; continue; }
done

cd "$W"
echo "== U1";  ./t_forma | tee u1.txt;          anota "$(grep -E '^(PASA|FALLA)' u1.txt)"
echo "== U2";  ./t_origen | grep -v -E 'Number of|Memory|^ *$' | tee u2.txt;  anota "$(grep -E '^(PASA|FALLA)' u2.txt)"
echo "== U3"
for caso in "2.0 0.05" "20.0 0.01"; do
  set -- $caso; ./t_poisson_nodal $1 $2 > /dev/null
  n=$(python3 -c "print(int($1/$2+0.5))"); tol=$([ $n -gt 500 ] && echo 1e-13 || echo 1e-14)
  python3 "$AQUI/py/ref_nodal.py" nodal_$n.dat $tol | tee u3_$n.txt; anota "$(grep -E '^(PASA|FALLA)' u3_$n.txt)"
done
echo "== U4";  ./t_poisson_orden | grep -v -E "Number of|Memory|^ *$" | tee u4.txt; anota "$(grep -E '^(PASA|FALLA)' u4.txt)"
echo "== U5";  ./t_fondos | grep -v -E 'Number of|Memory|^ *$' | tee u5.txt;  anota "$(grep -E '^(PASA|FALLA)' u5.txt)"
echo "== U6";  ./t_aa | grep -v -E 'Number of|Memory|^ *$|integrator|particles are|angle-action|Aborting' | tee u6.txt; anota "$(grep -E '^(PASA|FALLA)' u6.txt)"

echo "== U7";  ./t_r0 | grep -v -E 'Number of|Memory|^ *$' | tee u7.txt; anota "$(grep -E '^(PASA|FALLA)' u7.txt)"
echo "== U8";  ./t_masa | grep -v -E 'Number of|Memory|^ *$' | tee u8.txt; anota "$(grep -E '^(PASA|FALLA)' u8.txt)"

if [ $rapido = 0 ]; then
  echo "== I1";  python3 "$AQUI/py/integradores.py" "$exe" "$W/i1" 2>&1 | grep -v -i warn | tee i1.txt; anota "$(grep -E '^(PASA|FALLA)' i1.txt)"
  echo "== I2";  python3 "$AQUI/py/colapso.py" "$exe" "$W/i2" 2>&1 | grep -v -i warn | tee i2.txt;     anota "$(grep -E '^(PASA|FALLA)' i2.txt)"
  echo "== I3";  python3 "$AQUI/py/medida.py" "$exe" "$W/i3" 2>&1 | grep -v -i warn | tee i3.txt;      anota "$(grep -E '^(PASA|FALLA)' i3.txt)"
  echo "== I4";  python3 "$AQUI/py/limite_L0.py" "$exe" "$W/i4" 2>&1 | grep -v -i warn | tee i4.txt;   for x in $(grep -n -E '^(PASA|FALLA)' i4.txt | cut -d: -f1); do anota "$(sed -n ${x}p i4.txt)"; done
  echo "== I5";  python3 "$AQUI/py/salida_raw.py" "$exe" "$W/i5" 2>&1 | grep -v -i warn | tee i5.txt;   anota "$(grep -E '^(PASA|FALLA)' i5.txt)"
  echo "== I6";  python3 "$AQUI/py/muestreo.py" "$exe" "$W/i6" 2>&1 | grep -v -i warn | tee i6.txt;     anota "$(grep -E '^(PASA|FALLA)' i6.txt)"
fi

echo; echo "===== Resumen ($exe)"
printf "%s\n" "${resumen[@]}" | grep -v "^ *$"
exit $fallos
