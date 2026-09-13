# Lista de bugs y mejoras — seguimiento

Encontrados en la revisión de código de la rama `fix/bugs-mejoras`. Un
commit por ítem resuelto; esta lista se va tachando a medida que se
avanza.

## Resueltos

- [x] **Makefile: `FLAGS` vacío para gfortran.** Todas las líneas de
  `FLAGS` del bloque `gfortran` estaban comentadas, rompiendo el build
  (faltaba `-Jobjs`, se perdía `-fopenmp`). — commit `fix(build):
  activate gfortran FLAGS in Makefile`
- [x] **`analysish.f90`: línea de 133 columnas.** Excedía el límite de
  132 columnas del formato libre de Fortran, error de compilación real
  detectado al probar con gfortran instalado. — commit `fix(build):
  wrap overlong line in analysish.f90 and rebuild binary`
- [x] **`set_timestep` (`utils.f90`): `dt = dtr` anulaba `min(dtr,dtp)`.**
  El límite de paso de tiempo por fuerza nunca se aplicaba. De paso se
  reemplazó `dtp = courant*dpc/Fmax` por el criterio de aceleración
  estilo GADGET-2 `dtp = courant*sqrt(2*drc/Fmax)`, menos restrictivo y
  físicamente más correcto para leapfrog (dt·ω ≲ 2). — commit
  `fix(timestep): restore force-based dt bound with a less restrictive,
  physically motivated criterion`
- [x] **Makefile: `-fallow-argument-mismatch` innecesario.** Se había
  añadido preventivamente al arreglar `FLAGS`; verificado que no hace
  falta (gfortran no chequea el tamaño de arreglos de forma explícita
  cuando la dimensión depende de un argumento en tiempo de ejecución).
  — commit `fix(build): drop unneeded -fallow-argument-mismatch`
- [x] **`density.f90`: relleno de ghost zones incorrecto** en `density`
  y `avg_density` (`rho(i-1)=rho(i)` en vez de `rho(1-i)=rho(i)`),
  corrompía el primer punto físico y nunca llenaba el punto fantasma
  externo. Afecta directamente `poisson_rk.f90` (siembra la integración
  RK justo en el origen). — commit `fix(density): correct ghost-zone
  mirroring in density/avg_density`

## Resueltos (cont.)

- [x] **`save1Ddata` recibía arreglos con ghost cells en un argumento
  mudo de forma explícita más chico (`dimension(1:Nr)`).** Por
  asociación de secuencia de Fortran, el arreglo recibido quedaba
  corrido `ghost` posiciones: los `.rl` de salida (`vlasov_density`,
  `vlasov_avg_density`, `vlasov_curr`, y si `autointeraction`,
  `vlasov_force`/`vlasov_potential`) incluían al principio los puntos
  fantasma (r negativo) y perdían los últimos `ghost` puntos físicos
  reales cerca de `r=rmax`. Arreglado limitando cada llamada al rango
  físico `(1:Nr)` en `utils.f90`. — commit `fix(io): correct index
  shift when saving ghost-augmented grid arrays`

## Resueltos (cont. 2)

- [x] **`utils.f90` `reduce_arrays`: operadores de comparación
  inconsistentes** entre el conteo (`r_part(i)<=rmax`) y la copia
  (`r_aux(i)<rmax`) — una partícula justo en `r=rmax` dejaba una
  entrada sin inicializar en los arreglos reasignados. Arreglado
  usando `<=rmax` en ambos lados. — commit `fix(reduce_arrays): use
  consistent <=rmax in count and copy loops`

## Resueltos (cont. 3)

- [x] **`density.f90` / `poisson_rk.f90`: `collapse(2)` sin protección
  (condición de carrera OpenMP).** Los bucles más caros del código
  (`O(Nr×Npart)`) acumulan en `rho(i)`/`curr(i)`/`avg_rho(i)` y
  `pot_part(i)`/`force_part(i)`, indexados solo por el índice externo,
  pero paralelizaban con `collapse(2)` sobre externo+interno sin
  `atomic`/`reduction` — dos hilos podían caer en el mismo `i` (con
  distinto `j`) a la vez. Estaba dormido mientras `-fopenmp` estaba
  roto en el Makefile; quedó activo al arreglarlo. Arreglado
  paralelizando solo en el índice externo (mismo patrón que ya usa
  correctamente `avg_density` en el mismo archivo). Verificado:
  3 corridas con 8 hilos dan salida idéntica byte a byte. — commit
  `fix(openmp): remove collapse(2) data race in density() and
  poisson_rk()`

## Resueltos (cont. 4)

- [x] **`utils.f90` `deallocate_mem` era código muerto y roto**:
  intentaba desasignar `p_part_hp` (nunca asignado — typo por
  `p_part_h`), `res` (nunca asignado, con o sin `conv_test`),
  desasignaba `force`/`pot`/`dev_pot` sin comprobar `autointeraction`,
  y desasignaba `pot`/`dev_pot` dos veces si `autointeraction=.true.`.
  Además le faltaba `l_part` (fuga). No se llamaba desde ningún lado.
  Arreglado para reflejar exactamente `alloc_mem_set0`, y se agregó
  la llamada real al final de `main.f90` para poder validarlo en
  ejecución. Probado en ambas configuraciones (`autointeraction`
  `.true.`/`.false.`): terminan limpio con "Memory deallocated". —
  commit `fix(memory): repair deallocate_mem and actually call it`

## Resueltos (cont. 5)

- [x] **`functions.f90` `Sn`/`Wn` no rechazaban ningún orden inválido
  salvo `n` mayor al máximo** (`else if (n>4)`/`else if (n>3)`, no
  `else`): con `n<1` ninguna rama coincidía y se usaba el valor de
  retorno sin inicializar. Cambiado a un `else` genérico que atrapa
  cualquier orden inválido, y corregido de paso el mensaje de error
  de `Wn` (decía "greater than 4" cuando la condición real era
  `n>3`). Verificado con `bsplineorder=5` (sí alcanzable en la
  práctica): aborta limpio con el mensaje correcto. Nota: `n<1`
  resultó no ser alcanzable por ningún sitio de llamada actual (los
  6 están todos protegidos por `abs(distancia)<=cutoff`, y con
  `bsplineorder<=0` ese cutoff ya es `<=0`), así que no pude
  reproducirlo como fallo real — el arreglo queda como
  endurecimiento de robustez para cualquier llamada futura directa.
  — commit `fix(functions): Sn/Wn now reject any invalid order, not
  just n>4/n>3`

## Pendientes

- [ ] **`grav_force.f90`: condición `r_part(i)<1.d0` en el fondo
  `sphere`** sin `abs()` — cualquier partícula con `r_part` negativo
  transitorio (antes de la reflexión de simetría al final de cada
  paso) usa la fórmula interior sin importar su magnitud real.
- [ ] **`grav_force.f90`: fondos `iso`, `isotrun`, `nfw`, `burkert` no
  actualizan `pot_part`/`force_part`**, solo los arreglos de malla
  `pot`/`force` — las partículas nunca sienten esa fuerza de fondo
  (solo el término centrífugo). Además esos arreglos solo se reservan
  si `autointeraction=.true.`, así que con `autointeraction=.false.`
  (el caso normal para un fondo fijo) se escribe en memoria no
  reservada.
- [ ] **`initial_data.f90` / `parameters.f90`: `state="gaussian"` (el
  valor por defecto) no coincide con ninguna rama** (`initial_data.f90`
  solo reconoce `"gaussian1"`, no `"gaussian"`), y no hay `else`/`stop`
  de captura — con la configuración de fábrica, `f` queda en cero en
  silencio.
- [ ] **`grav_force.f90` / `parameters.f90`: `forcetype="self"` es una
  opción documentada que no hace nada.** `main.f90` bifurca sobre
  `forcetype=="self"`, pero `grav_force.f90` solo mira `forcetype=="bg"`
  y el booleano independiente `autointeraction`. Con
  `forcetype="self"` y `autointeraction=.false.` las partículas no
  sienten fuerza radial (documentación/código desincronizados).
- [ ] **Carpeta `src/`: archivos legado sin extensión `.f90`**
  (`analysish`, `poisson`, `poisson_ps`, `reduce_arrays`) no se
  compilan, están desincronizados de sus homónimos activos, y
  `poisson_ps` referencia un módulo `chebyshev` inexistente. Mover a
  `legacy/` o borrar.

## Mejoras de rendimiento

- [x] **`density()`/`avg_density()`/`poisson_rk()`: búsqueda de vecinos
  por fuerza bruta `O(Nr×Npart)`.** Reemplazada por una lista de
  celdas (`build_cell_list` en `utils.f90`, `O(Nr+Npart)`). Validado
  bit a bit contra la versión sin optimizar en 3 escenarios (incluido
  autogravitante). Benchmark autogravitante completo
  (`input_parameters`, 17084 partículas, 100000 pasos,
  `autointeraction=.true.`): main 4 hilos = 2467.0s vs fix 1 hilo =
  849.5s (**2.90× más rápido con ¼ de los hilos**). Energía y
  diagnósticos finales coinciden a varias cifras significativas;
  diferencias por partícula (~10⁻⁶ mediana) consistentes con caos
  numérico, no con un bug. — commit `perf(density): replace
  O(Nr*Npart) brute-force deposit/interpolation with a cell list`.
  Nota: con la optimización, más hilos ya no ayuda claramente a
  Npart moderado (ver el commit) — posible trabajo futuro: paralelizar
  `build_cell_list` o reusar sus buffers entre llamadas.
- [ ] **`avg_density()` corre cada paso en modo autogravitante**
  (vía `grav_force()→poisson_rk()`), no solo cada `spatial_output`
  como `density()`. No es un bug — es inherente al método de campo
  autoconsistente — pero sigue siendo el costo dominante en corridas
  autogravitantes largas. Ya mitigado en gran parte por el punto
  anterior; posible mejora futura si se necesita más velocidad.

## Mejoras de diseño (no son bugs)

- [ ] `test_consistency` (`utils.f90`) nunca se llama desde `main.f90`
  y referencia `lmin`/`lmax` que no existen (son `lminc`/`lmaxc`) — no
  compilaría si se descomenta tal cual.
- [ ] Integrador `rk4` declarado como opción válida pero no
  implementado (aborta con mensaje).
- [ ] Limpiar el código muerto/comentado en el bucle principal de
  `main.f90` (ecuación de continuidad, paso de tiempo adaptativo para
  autogravedad, etc.).
