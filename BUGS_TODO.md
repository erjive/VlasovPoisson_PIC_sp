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

## Pendientes

- [ ] **`save1Ddata` recibe arreglos con ghost cells en un argumento
  mudo de forma explícita más chico (`dimension(1:Nr)`).** Por
  asociación de secuencia de Fortran, el arreglo recibido queda corrido
  `ghost` posiciones: los `.rl` de salida (`vlasov_density`,
  `vlasov_avg_density`, `vlasov_curr`, y si `autointeraction`,
  `vlasov_force`/`vlasov_potential`) incluyen al principio los puntos
  fantasma (r negativo) y pierden los últimos `ghost` puntos físicos
  reales cerca de `r=rmax`. No afecta la física de la simulación, sí
  los datos de diagnóstico/gráficas. Afecta las llamadas en
  `utils.f90` dentro de `save_data`.
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
- [ ] **`utils.f90` `reduce_arrays`: operadores de comparación
  inconsistentes** entre el conteo (`r_part(i)<=rmax`) y la copia
  (`r_aux(i)<rmax`) — una partícula justo en `r=rmax` deja una entrada
  sin inicializar en los arreglos reasignados.
- [ ] **`arrays.f90` `deallocate_mem` es código muerto y roto**: intenta
  desasignar `p_part_hp` (nunca asignado), desasigna `force`/`pot`/
  `dev_pot` sin comprobar `autointeraction`, y los desasigna dos veces
  si `autointeraction=.true.`. Hoy no se llama desde ningún lado, pero
  crashea en cuanto alguien lo use.
- [ ] **`grav_force.f90` / `parameters.f90`: `forcetype="self"` es una
  opción documentada que no hace nada.** `main.f90` bifurca sobre
  `forcetype=="self"`, pero `grav_force.f90` solo mira `forcetype=="bg"`
  y el booleano independiente `autointeraction`. Con
  `forcetype="self"` y `autointeraction=.false.` las partículas no
  sienten fuerza radial (documentación/código desincronizados).
- [ ] **`functions.f90` `Sn`/`Wn` no validan `n<1`** (p. ej.
  `bsplineorder=0` o negativo): ninguna rama coincide y se usa la
  variable de retorno sin inicializar.
- [ ] **Carpeta `src/`: archivos legado sin extensión `.f90`**
  (`analysish`, `poisson`, `poisson_ps`, `reduce_arrays`) no se
  compilan, están desincronizados de sus homónimos activos, y
  `poisson_ps` referencia un módulo `chebyshev` inexistente. Mover a
  `legacy/` o borrar.

## Mejoras de diseño (no son bugs)

- [ ] `test_consistency` (`utils.f90`) nunca se llama desde `main.f90`
  y referencia `lmin`/`lmax` que no existen (son `lminc`/`lmaxc`) — no
  compilaría si se descomenta tal cual.
- [ ] Integrador `rk4` declarado como opción válida pero no
  implementado (aborta con mensaje).
- [ ] Limpiar el código muerto/comentado en el bucle principal de
  `main.f90` (ecuación de continuidad, paso de tiempo adaptativo para
  autogravedad, etc.).
