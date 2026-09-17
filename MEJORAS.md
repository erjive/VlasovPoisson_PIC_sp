# Mejoras para `VlasovPoisson_PIC_sp` a partir de lo aprendido en `vlasov-poisson_PIC`

Inventario de cambios a hacer en este código (distribución en momento angular L,
`l_part`), sacado de comparar su estado actual (rama `fix/bugs-mejoras`, commit
`74234b1`) con lo que se corrigió, midió o descubrió en `vlasov-poisson_PIC`
(L fijo, rama `clean/comentarios`). Cada punto dice qué hay hoy aquí, qué se aprendió
allá y qué cambiar.

## Estado (rama `mejoras/portar`)

Cada punto se implementó en su propio commit y se verificó antes de seguir; las cifras
de cada verificación están en `BUGS_TODO.md`.

| punto | estado | commit |
|---|---|---|
| A1 energía ½ | hecho | `abb8bc9` |
| A2 mapa ángulo-acción con L | hecho como herramienta (`tools/aa_numerico_L.py`, `tools/hk_numerico.py`); `analysish` sigue con el isócrono | `ed6e9cb` |
| A3 cuadratura en (Q,J,L) | hecho (`state=aa_quad`) | `065519a` |
| A4 cutoff | medido; ya era 0 por omisión | `0fa4e10` |
| A5 nquad | hecho | `5e11235` |
| A6 funciones de prueba | hecho | `8203a63` |
| A7 eps | decidido: eps=0; medido | `0fa4e10` |
| A8 paso fijo | hecho (`dt_switch`) | `828b79b` |
| A9 errores dormidos | hechos (estados, forcetype, fondos); nota Wn/drc revisada | `e128127`, `8f079fd` |
| B1 parámetros por nombre | hecho | `9aad262` |
| B2 código de salida | hecho | `a5993e4` |
| B3 directorio | hecho | `e9cb7a9` |
| B4 semilla | no aplica: no hay estados aleatorios | — |
| B5 Makefile | hecho (`-w` fuera, `h5fc`) | `506a3df`, `6d8b997` |
| B6 legado | hecho | `fd39130` |
| B7 salida de h_k | hecho, con `field_output` | `8203a63` |
| B8 comentarios | hecho | `23a30ac` |
| C1 a_k una vez | hecho | `6200fe2` |
| C2 rendimiento | hecho: 2.16× con autogravedad, 2.39× sin ella; `build_cell_list` sigue serial | `5219f6f`..`75c5779` |
| C3 integrador analytic | hecho | `7bc1f93` |
| D1 checkpoint | hecho (y `l_part` en HDF5) | `894a120` |
| D2 equilibrio F(J,L) | hecho (`tools/equilibrio_L.py`) | `ed6e9cb` |
| D3 exacto y pruebas agnósticas | hecho (`dftype`, `tools/hk_exacto.py`) | `065519a`, `94e697d` |
| D4 controles | `tools/delta_phi.py`; el solver lineal con L queda pendiente | `0f446e4` |

Además, encontrados al verificar: h_k NaN con partículas no ligadas (`732250f`) y sumas
OpenMP no deterministas en energía y h_k (`e6d2011`).

Referencias: `vlasov-poisson_PIC/BUGS_TODO.md`, `vlasov-poisson_PIC/PREGUNTAS_ABIERTAS.md`
y `vlasov-poisson_PIC/docs/introduccion/vlasov_intro.tex`.

---

## A. Errores que afectan resultados

### A1. La energía cuenta dos veces la autoenergía gravitatoria
**Hoy:** `energy.f90` suma `(p^2/2 + pot_part) * f * l_part`, con `pot_part` que incluye
el potencial propio completo.
**Aprendido:** la energía de interacción lleva un factor 1/2 (cada par se cuenta desde
ambos lados). Sin él, la energía "cambia" en ½ΔW_self cuando la masa se redistribuye:
un error de -4.96e-4 con a0=1e-3, inmune a dr, orden del spline y N, y lineal en a0.
Corregido, el error real fue 1.7e-7 (commit `b0f3915`).
**Cambio:** guardar el potencial propio por partícula (`potself_part`) justo después de
`poisson_rk` y usar `pot_part - potself_part/2` en la energía.

### A2. Con autogravedad, h_k se calcula con las variables ángulo-acción equivocadas
**Hoy:** `analysish.f90` usa el mapa analítico del isócrono aunque haya autogravedad.
**Aprendido:** eso produce una meseta estática en h_1 proporcional a la masa
(1.77e-3 h_0 con a0=1e-3, 1.23e-2 con 1e-2), sin fase, que imita un modo que no decae
e invalida cualquier criterio de "mezclado" por debajo de ~2 a0 h_0. Con el mapa del
potencial real desaparece (baja 2200 veces).
**Cambio:** mapa ángulo-acción numérico con L por partícula. Con L variable, J(E,L)
es una tabla en dos variables, o un mapa unidimensional por cada valor de L de la
rejilla. Base: `vlasov-poisson_PIC/reproducir/scripts/aa_numerico.py`. Afecta directamente la tabla de
"mixed yes/no" del borrador `Vlasov_Poisson_evolutions/main.md`.

### A3. El estado `aa` desplaza las partículas dentro de su celda
**Hoy:** el estado `aa` es una rejilla regular en (r,p,L) donde se evalúa F(Q,J,L).
(Corrección: el desplazamiento con secuencia de Weyl solo existe en la rama
`experimento/jitter`, no en la principal.)
**Aprendido:** la recurrencia es el límite de Nyquist en J. Cualquier desplazamiento
destruye la convergencia de cuadratura y deja el error en 1/sqrt(N) (medido y revertido,
commit `c4df269`). Lo que funciona es una rejilla regular en las variables
ángulo-acción, sin desplazamientos (`aa_quad`, commit `7b9dff4`): error de h_k al nivel
del integrador (1.2e-9 de h_0) con N=5000.
**Cambio:** estado de cuadratura en (Q, J, L) con nodos en puntos medios de celda, e
inversión del mapa para obtener (r, p_r). Requisito de resolución, ahora en dos
direcciones:
- N_J ≳ 4 k Δω_J t_max / 2π
- N_L ≳ 4 k Δω_L t_max / 2π, con ∂ω/∂L = -3 (J+c)^-4 · ½(1 + L/sqrt(L²+4)).
  En J=0.1, L=2: ∂ω/∂J = -0.0751 y ∂ω/∂L = -0.0641, así que una dispersión en L mezcla
  casi tanto como una en J.

### A4. El corte `cutoff` trunca la distribución
**Hoy:** `aa` y `gaussian1` eliminan las partículas con f ≤ `cutoff`·f_max. Las
entradas del artículo usan `cutoff=0.01`.
**Aprendido:** en el código de L fijo, un corte análogo (`r0=0.01`) dejaba el error
clavado en ~1e-2 sin importar N, porque descartaba la mayoría de los nodos y truncaba
la cuadratura (commit `e547021`).
**Cambio:** `cutoff=0` por omisión; si se quiere ahorrar partículas, medir antes el
sesgo que introduce.

### A5. La cuadratura angular de la función de prueba es imprecisa
**Hoy:** `analysish.f90` usa Simpson con `nquad=20`.
**Aprendido (medido):** con σ_Q=0.1, el valor de las entradas del artículo, Simpson de
20 intervalos tiene un error relativo de 1.3e-2 (k=1) y 3.1e-2 (k=4) en el coeficiente
a_k. Con 512 intervalos es ≤2.4e-15. Ese error reescala h_k y sesga la comparación con
la solución exacta.
**Cambio:** `nquad=512`, calculado una sola vez (ver C1).

### A6. La función de prueba está atada a los parámetros de la condición inicial
**Hoy:** Φ usa `sp`, `sr`, `sl`, `l0` (los mismos anchos que la distribución) y J_0=0.
**Aprendido:** la envolvente de mezcla depende del ancho de la ventana de la función de
prueba; conviene poder variarla sin cambiar la física. El código de L fijo usa
parámetros propios (`j1`, `sj1`, `sq1`) y dos funciones de prueba.
**Cambio:** parámetros independientes para Φ, incluida su parte en L.

### A7. Sin suavizado centrífugo, y L puede acercarse a 0
**Hoy:** `eps=0` fijo (la regresión documentada en `BUGS_TODO.md`); con dispersión en L,
las partículas de L pequeño llegan cerca de r=0, donde se documentó ~350% de error de
energía.
**Aprendido:** `eps≠0` es inconsistente con el mapa ángulo-acción analítico, que no
incluye el suavizado, y creó un piso en h_k de ~5e-12 (commits `8954066`, `23e03de`).
**Cambio:** mantener `eps=0` y acotar el soporte en L lejos de cero (`lminc>0` con
margen), o, si hace falta suavizar, que el mapa numérico use el mismo potencial
suavizado. Decidir antes de producir resultados.

### A8. El paso de tiempo se adapta durante las corridas con autogravedad
**Hoy:** `set_timestep` se vuelve a llamar durante la evolución si
`autointeraction=.true.` (commit `64c47b3`).
**Aprendido:** un paso variable rompe la conservación simpléctica del integrador. Con
paso fijo y yoshida4, los resultados son independientes de Δt entre 0.05 y 0.2: la parte
estática de h_1 cambia 1e-4 relativo y h_1 difiere ≤1.6e-6 relativo hasta t=800.
**Cambio:** paso fijo; elegirlo con una prueba de Δt.

### A9. Errores pendientes de este código (de su propio `BUGS_TODO.md`)
- `compact`, `compact2` y `Plummer` dejan `l_part=0`, y los diagnósticos, ponderados
  por `l_part`, dan cero.
- `state="gaussian"` (el valor por omisión) no corresponde a ninguna rama: f queda en
  cero sin aviso.
- `forcetype="self"` no hace nada.
- Los fondos `iso`, `isotrun`, `nfw` y `burkert` no actualizan `pot_part`/`force_part`.
- `Wn` usa `dr` como ancho de la función de forma y `rho` usa `drc`. **Revisar también
  en el código de L fijo:** `vlasov-poisson_PIC` tiene la misma forma
  (`Wn(bsplineorder,(r-r_part)/dr)` en `density.f90` y `poisson_rk.f90`). Allí la
  meseta no cambió con dr entre 0.2 y 0.025, pero la consistencia no se ha verificado.

---

## B. Reproducibilidad y uso

### B1. Parámetros por posición desde la entrada estándar
**Hoy:** `read(*,*)` de líneas en orden fijo; un parámetro fuera de lugar cambia otro
sin aviso.
**Cambio:** archivo `clave = valor` con overrides en la línea de comandos y volcado de la
configuración resuelta (`params_usados.par`) en cada corrida (`paramfile.f90`, commits
`95820d0`, `7787b51`).
**Aprendido además:** el volcado de `vlasov-poisson_PIC` omitía `dftype`, y las corridas
con una distribución distinta de la de omisión no se podían repetir desde su archivo.
Al portar, comprobar que se escribe cada parámetro de la lista `pname` y repetir una
corrida desde su `params_usados.par`, comparando h_k.

### B2. Los errores terminan con código de salida 0
**Hoy:** `stop` sin código en todas las rutas de error.
**Aprendido:** eso escondió 15 corridas fallidas seguidas en un script.
**Cambio:** `stop 1` en errores (commit `45d2a58`).

### B3. El directorio de salida se trunca a 20 caracteres
**Hoy:** `character(20) :: directory`.
**Cambio:** `character(100)` y argumentos de longitud asumida en las subrutinas de
escritura.

### B4. Semilla reproducible
**Hoy:** no hay estados aleatorios.
**Cambio:** si se agrega un estado Monte Carlo, parámetro `seed` escrito en
`params_usados.par` (gfortran siembra desde el sistema operativo por omisión).

### B5. Makefile
**Hoy:** `-w -Wall`; el `-w` silencia todas las advertencias.
**Aprendido:** `-w` escondía un truncamiento real de una cadena de 512 a 256 caracteres.
**Cambio:** quitar `-w`; detectar `h5fc` y comprobar las banderas (commits `225eabe`,
`8b3f706`).

### B6. Archivos legado en `src/`
`analysish`, `analysish_old.f90`, `poisson`, `poisson_ps`, `reduce_arrays`: no se
compilan o están desincronizados. Mover a `legacy/` o borrar.

### B7. Salida de h_k
**Hoy:** solo |h_k| en `hk.tl`, con una función de prueba.
**Aprendido:** la fase de h_k fue decisiva: distinguió un residuo estático de un modo, y
permitió separar la parte que gira y ajustar ω y γ.
**Cambio:** escribir Re y Im de h_k (`hk1_complex.tl`), dos funciones de prueba, y
cadencia de instantáneas (`field_output`) separada de la de h_k (commits `8b511e3`,
`aaaf27d`).

### B8. Comentarios
Limpiar a solo método y física, sin bitácora de cambios en el código (commit `fb75fc2`).

---

## C. Rendimiento

### C1. La cuadratura angular de la función de prueba se repite por partícula
**Hoy:** `analysish` evalúa la integral en Q (21 nodos × 5 modos) para cada partícula.
**Aprendido:** la función de prueba es separable, así que a_k es una constante de toda
la corrida. Sacarla del ciclo hizo la corrida completa 10.4× más rápida (commit
`285daa2`). Aquí también separa: A(Q)·B(J)·C(L).
**Cambio:** calcular a_k una vez y por partícula solo B(J)·C(L)·e^{-ikQ}.

### C2. Pasadas seriales y lista de celdas
En el código de L fijo, eliminar las pasadas seriales por partícula dio 23×
(commit `1d4c9b1`) y paralelizar `build_cell_list`, 3.08× (commit `bb57363`). Ninguno de
los dos commits está en este repositorio; medir antes y después con ABBA.

### C3. Integradores
**Hoy:** euler, leapfrog, yoshida4, rk4.
**Cambio:** agregar el integrador `analytic` (avance exacto en variables ángulo-acción,
con J y L por partícula) como referencia de validación. yoshida6 no hizo falta.

---

## D. Métodos y diagnósticos nuevos

### D1. Estado `checkpoint`
Leer partículas "r p_r L f" de un archivo, para arrancar desde estados que el código no
construye (commit `601c6d2`).

### D2. Equilibrio autoconsistente con perturbación separada
Generalizar `vlasov-poisson_PIC/reproducir/scripts/equilibrio.py` a F_eq(J, L): con L en rejilla, cada valor
de L es un problema de un grado de libertad en el mismo potencial total, e iterar
Poisson sobre todos. Es el montaje necesario para medir amortiguamiento de Landau
(pregunta abierta 2 del otro código).

### D3. Validación con solución exacta y distribuciones de prueba
- h_k(t) = a_0 a_k ∫L dL ∫dJ B(J,L) c_k(J,L) e^{-ikω(J,L)t} / ∫L dL ∫dJ c_0.
- Repetir las pruebas agnósticas del otro código: modos prohibidos, espiral que se
  desenrolla en t* y convergencia h² de King.

### D4. Diagnósticos y controles que resultaron decisivos
- δΦ(r,t) respecto del equilibrio: no depende de coordenadas.
- Restar la corrida sin perturbación (misma rejilla) para cancelar el desajuste
  determinista del Poisson discreto.
- Pruebas de resolución (N_J, N_Q, N_L, Δr, Δt) antes de interpretar cualquier señal a
  tiempos largos: después de t≈2000, la componente que gira de h_1 resultó ruido de
  discreción.
- Referencia de phase mixing libre en el potencial de equilibrio.
- Solver lineal sin partículas (en validación en el otro código).
