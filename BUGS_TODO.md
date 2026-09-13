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
- [ ] **`eps` (longitud de suavizado del término centrífugo) está fija
  en `0.0` siempre** (`set_grid_size` en `utils.f90`, línea con
  `eps = 0.0D0`, con el cálculo real comentado justo arriba). Esto
  significa que el potencial centrífugo $L^2/(2r^2)$ es genuinamente
  singular en $r=0$ sin ningún suavizado — encontrado al validar el
  `set_timestep` adaptativo: un encuentro cercano con $r\approx 0$ y
  momento angular chico sigue dando ~350% de deriva de energía
  incluso con `dt` adaptativo, porque no hay forma de "suavizar" la
  fuerza cerca de la singularidad. Activar `eps` (o exponerlo como
  parámetro configurable) podría mejorar sustancialmente la
  estabilidad de corridas con partículas de bajo momento angular.

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
  Npart moderado (ver el commit).
- [x] **Intentado y revertido: reusar los buffers de `build_cell_list`**
  en vez de reservarlos/liberarlos en cada llamada (para evitar el
  `allocate`/`deallocate` de arreglos `Npart` en cada paso, en modo
  autogravitante). Implementado moviendo `cell_start`/`particle_order`
  y los arreglos auxiliares a nivel de módulo (`arrays.f90`),
  persistentes entre llamadas. **Empeoró el rendimiento en vez de
  mejorarlo.** Con instrumentación directa (500 llamadas): el overhead
  real de `allocate`/`deallocate` que se quería eliminar era chico
  (0.227s de 2.13s totales, ~10%), y el propio bucle de depósito
  — algorítmicamente idéntico en ambas versiones — pasó de 1.902s a
  3.063s al leer de arreglos a nivel de módulo en vez de argumentos
  mudos locales con `intent()` explícito (hipótesis: gfortran pierde
  margen para asumir ausencia de aliasing y vectoriza peor). Revertido
  por completo (`arrays.f90`/`density.f90` de vuelta al commit
  anterior); no vale la pena perseguir esta variante de la idea.
- [x] **Bug nuevo encontrado de paso: `r(0)` se escribía fuera de los
  límites del arreglo cuando `rmin>0`** (`construct_grid` llena
  `r(i)` para `i=0,Nr` en ese caso, pero `alloc_mem_set0` solo
  reservaba `r(1:Nr)`, ya que `ghost=0` cuando `rmin>0`). Corrompía
  el heap en silencio — nunca se detectaba porque nada validaba la
  integridad del heap hasta que algo lo liberaba, y `deallocate_mem`
  era código muerto hasta hace pocos commits. Encontrado por casualidad
  al validar el intento de reuso de buffers de arriba, con un caso de
  prueba `rmin>0` — la primera vez en esta sesión que se ejercitó
  `rmin>0` junto con `deallocate_mem` realmente ejecutándose. — commit
  `fix(grid): allocate r(0:Nr) for rmin>0, not r(1:Nr)`
- [ ] **`avg_density()` corre cada paso en modo autogravitante**
  (vía `grav_force()→poisson_rk()`), no solo cada `spatial_output`
  como `density()`. No es un bug — es inherente al método de campo
  autoconsistente — pero sigue siendo el costo dominante en corridas
  autogravitantes largas. Ya mitigado en gran parte por el punto
  anterior; posible mejora futura si se necesita más velocidad.
- [x] **Output en ASCII de texto plano, pesado y lento de escribir**
  (`vlasov_fdist.2D` llegó a 838 MB en una corrida de 100 000 pasos).
  Agregado output HDF5 opcional (`output_format="hdf5"`), un archivo
  `.h5` por corrida, un grupo por snapshot, comprimido con gzip.
  Validado bit a bit contra ASCII en 2 escalas. A escala realista
  (50 000 partículas, 21 snapshots): **34% más rápido, 62% más chico**.
  Rama `feature/hdf5-output`, commit `feat(io): add optional HDF5
  output, selected via output_format parameter`. `hk.tl` y
  `vlasov_rhomix.tl` se quedaron en ASCII (fuera de alcance, chicos).
- [x] **Integrador simpléctico de orden superior** (Yoshida 1990,
  4º orden, 3 sub-pasos tipo leapfrog). Agregado como opción nueva
  `integrator="yoshida4"` (no reemplaza `leapfrog`, hay que elegirlo
  explícitamente). Validado con cuidado porque las primeras pruebas
  autogravitantes *explotaron* (tanto `leapfrog` como `yoshida4`) —
  no era un bug del integrador nuevo, era `dt` insuficiente cerca de
  la barrera centrífuga (`L²/r²`, ver el ítem de `set_timestep`
  abajo). Con potencial de fondo suave y momento angular alejado de
  `r=0`: `yoshida4` da **6600× menos deriva de energía** que
  `leapfrog` al mismo `dt` (9.6×10⁻⁹ vs 6.3×10⁻⁵), a ~2× el costo por
  paso (no 4×, la evaluación de fuerza no es todo el costo por paso).
  Con autogravedad (partículas lejos de `r=0` para que no explote):
  **sin diferencia medible** entre ambos (1.85% de deriva los dos) —
  el ruido de discretización del propio método PIC (densidad estimada
  con partículas finitas) domina sobre el error de integración
  temporal, así que un orden temporal mayor no ayuda ahí a menos que
  suba la resolución (más partículas). Recomendación: usar
  `yoshida4` para corridas con `forcetype="bg"` (fondo fijo, sin
  autointeracción); para `autointeraction=.true.` el beneficio no
  está probado y cuesta ~2× igual. — commit `feat(integrator): add
  optional 4th-order symplectic integrator (yoshida4)`
- [x] **`set_timestep` solo se llamaba una vez, antes del bucle
  principal** — en modo autogravitante `Fmax` puede crecer si el
  sistema colapsa. Se activó la reevaluación periódica (cada paso,
  cuando `autointeraction=.true.`), reemplazando el `if
  (forcetype=="self")` comentado que nunca se hubiera disparado (ver
  el ítem de `forcetype="self"` sin efecto, arriba) por `if
  (autointeraction)`, que es el flag que de verdad controla si la
  fuerza cambia con el tiempo.

  De paso, al validar esto contra el caso que hizo explotar a
  `yoshida4` (momento angular chico, partículas cerca de `r=0`)
  apareció un **segundo bug independiente** que anulaba por completo
  la reevaluación: `set_timestep()` solo calculaba `Fmax` si
  `BGtype/="null"` — con `BGtype="null"` y autogravedad pura
  (`autointeraction=.true.`), `Fmax` nunca se tocaba y se quedaba
  fijo en `0.0` toda la corrida, sin importar qué tan grande se
  pusiera la fuerza real. Confirmado con instrumentación directa.
  Arreglado a `BGtype/="null" .or. autointeraction`.

  Con ambos arreglos: el caso que antes divergía a ~10¹⁴ ya no
  explota (queda acotado, ~350% de deriva) — `dt` se ve encogiendo
  en los primeros pasos (0.025→0.0091→0.0176→de vuelta a 0.025)
  en respuesta al encuentro cercano inicial. El ~350% de deriva que
  queda es un caso genuinamente extremo (encuentro cercano casi
  singular, sin suavizado — `eps` está fijo en `0.0` en
  `set_grid_size`, posible mejora futura), no una promesa de
  conservación perfecta de energía. Verificado sin regresión en el
  caso suave (sin autogravedad: deriva idéntica bit a bit) y en un
  caso autogravitante bien resuelto (partículas lejos de `r=0`: sin
  cambio, 1.85% en ambos). — commit `feat(timestep): adapt dt
  periodically during self-gravitating runs`
- [ ] Paralelizar la generación de partículas iniciales en
  `initial_data.f90` (estados `aa`/`gaussian1`) — el `!$OMP` ya está
  escrito pero deshabilitado por la dependencia secuencial del índice
  `indx`; se resuelve calculándolo con una fórmula cerrada. Medido en
  ~1-2s para 2M candidatos hoy, así que no es urgente.
- [ ] Cachear la evaluación de `Sn(...)` en `density()` — se calcula
  dos veces con los mismos argumentos (una para `rho`, otra para
  `curr`) en el bucle más caliente del código.
- [x] **Intentado y revertido: `-march=native` / `-flto`.** Probado y
  medido con cuidado (binario de referencia sin flags, comparación
  intercalada, aislando cada flag por separado en dos builds
  distintos para no confundirlos entre sí — mismo rigor que con el
  intento de reuso de buffers). Resultado en este CPU (i7 Haswell,
  AVX2/FMA) y este código: `-march=native` solo **no dio ninguna
  ganancia medible** (idéntico a la base en 50k y 200k partículas, a
  veces hasta un poco peor dentro del ruido). `-flto` solo fue una
  **regresión clara de ~30%** (0.65s vs 0.50s, repetible). Combinados,
  el mismo ~30% de regresión (dominado por `-flto`). Revertido el
  `Makefile` por completo. No vale la pena perseguir esta idea en este
  código/máquina tal como está. `-ffast-math` ni se probó, dado que
  las otras dos ya no rindieron.
- [ ] Evitar copiar arreglos completos cada paso (`r_part_p=r_part`,
  `p_part_p=p_part` en `main.f90`) con un esquema ping-pong de índices
  en vez de copia.
- [ ] Vectorizar/agrupar la cuadratura de `phik` en `analysish.f90`
  (hoy hace Simpson de 21 puntos por partícula y por modo, con
  llamadas a función una por una).

## Mejoras de diseño (no son bugs)

- [ ] **Forma más amigable de pasar los parámetros de la simulación.**
  Hoy `read_initial_param` (`utils.f90`) lee `input_parameters` de
  forma puramente posicional (`read(*,*) x` sin nombres) — agregar,
  quitar o reordenar un parámetro rompe silenciosamente cualquier
  archivo de entrada existente que no se actualice a la par (ya nos
  pasó al agregar `output_format`: cualquier `input_parameters` viejo
  deja de funcionar hasta agregarle la línea nueva al final). Un
  formato con claves (namelist de Fortran, TOML, YAML, o JSON) sería
  mucho más robusto y auto-documentado, y permitiría detectar
  parámetros faltantes/mal escritos con un mensaje claro en vez de
  un `read` que falla de forma críptica o lee el valor equivocado.

- [ ] `test_consistency` (`utils.f90`) nunca se llama desde `main.f90`
  y referencia `lmin`/`lmax` que no existen (son `lminc`/`lmaxc`) — no
  compilaría si se descomenta tal cual.
- [ ] Integrador `rk4` declarado como opción válida pero no
  implementado (aborta con mensaje).
- [ ] Limpiar el código muerto/comentado en el bucle principal de
  `main.f90` (ecuación de continuidad, paso de tiempo adaptativo para
  autogravedad, etc.).
