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

## Resueltos (cont. 6)

- [x] **La energía contaba dos veces la autoenergía gravitatoria**
  (MEJORAS A1). `energy.f90` sumaba el potencial propio completo por
  partícula; la energía de interacción lleva ½. Ahora `grav_force`
  guarda `potself_part` justo después de `poisson_rk` y la energía usa
  `pot_part - potself_part/2`. Medido con los casos de referencia
  (Isócrono + autogravedad, a0=1e-3, yoshida4, dt=0.025):
  t2 (`aa`, t=100) máx |E/E0-1| de 5.87e-4 a 2.52e-7;
  t3 (`gaussian1`, t=50) de 4.91e-4 a 3.46e-7; t1 (sin autogravedad)
  sin cambio. Trayectorias, densidades y h_k idénticos bit a bit.
  Caveat: con `BGtype=sphere` y autogravedad, `grav_force` sobrescribe
  `pot_part` (y `force_part`) con el fondo y descarta el potencial
  propio (bug de A9, sin corregir aquí), así que ahí la energía resta
  ½·potself sin haberlo sumado.

- [x] **El paso de tiempo se readaptaba en cada paso con autogravedad**
  (MEJORAS A8). Un paso variable rompe el carácter simpléctico de
  leapfrog/yoshida4. Ahora el paso se fija con la fuerza inicial
  (`dt_switch=fix`, valor por omisión) y solo se recalcula con
  `dt_switch=var`. Verificado: t1/t2/t3 idénticos bit a bit (ahí dt
  nunca cambiaba: manda el límite de Courant); con `pmax=0.01 Nrc=400`
  (manda el criterio de fuerza), `var` reproduce bit a bit el binario
  anterior (dt entre 0.570 y 0.982 por bloque) y `fix` mantiene
  dt=0.53768.

- [x] **h_k salía NaN toda la corrida si había una partícula no ligada.**
  Encontrado al verificar C1: el caso de referencia t3 (`gaussian1`,
  r0=5, L en [1.8,2.2]) tiene 10 de 782 partículas con E>0, para las
  que J y Q no existen; `analysish` sumaba NaN (ya pasaba en la versión
  original). Ahora se omiten (Φ=0), que es la extensión continua de
  B(J)=J² exp(-J²/sr²) cuando E→0⁻. Verificado: t1 y t2 idénticos bit a
  bit; h_k(t=0) de t3 coincide a 1e-15 con un cálculo independiente en
  Python sobre las 772 partículas ligadas.

- [x] **Cuadratura angular de h_k imprecisa** (MEJORAS A5, con C1).
  Simpson con 20 intervalos y sp=0.1 daba a_k con errores relativos de
  1.24e-2 (k=0) a 3.15e-2 (k=4) frente a la forma cerrada
  a_k = e^{-x} I_k(x), x=1/(2 sp²); con 512, ≤2.2e-16. Verificado: en
  t1/t2/t3, h_k(512)/h_k(20) es en cada fila la constante predicha
  a_k(512)/a_k(20) (1.01253 … 1.03249) hasta los 12 decimales de la
  predicción; las demás salidas, idénticas bit a bit.

- [x] **Función de prueba atada a la condición inicial y h_k solo en
  magnitud** (MEJORAS A6, B7). Dos funciones de prueba con parámetros
  propios (`j1 sj1 sq1 lt1 slt1`, ídem con 2; sin darlos toman sp, sr,
  0, l0, sl), salida `hk1.tl`, `hk2.tl`, `hk1_complex.tl`,
  `hk2_complex.tl` con 17 cifras (reemplaza `hk.tl`), y `field_output`
  para la cadencia de instantáneas. Verificado: con valores por omisión
  hk1 = hk anterior a ≤1.7e-15 y hk2 = hk1; con parámetros propios y
  p0=0.05, h_k(t=0) de ambas funciones coincide a ≤2.6e-15 (fases a 12
  decimales) con Python usando a_k exacto (Bessel); `field_output=80`
  deja las instantáneas pares idénticas y `60` termina con código 1.

- [x] **Energía y h_k cambiaban a 1e-15 entre corridas idénticas.**
  Las reducciones OpenMP de `energy` y `analysish` combinan las sumas
  parciales en el orden en que terminan los hilos (lo muestra un
  programa mínimo: la misma suma con 4 hilos dio dos valores en 6
  repeticiones, también con `schedule(static)`); con 8 cifras de salida
  no se veía. Ahora cada hilo acumula en local y las partes se suman en
  orden de hilo. Verificado: tres corridas de t2 (HDF5) y de t3 (ascii)
  idénticas bit a bit con 4 hilos; `field_output=80` frente a 40,
  idénticas en las instantáneas comunes y en h_k; tiempo en ABBA
  30.6 s → 30.0 s. Con otro número de hilos el redondeo sigue siendo
  distinto (≤1e-15), como es de esperar.

- [x] **Revisado: `Wn` con `dr` frente a `Sn` con `drc` en la densidad**
  (nota de A9). Poisson usa `avg_rho` (núcleo `Wn` de ancho `dr`), que
  conserva la masa: ∫4πr²ρ dr = a0 a ≤4.2e-5 en t1/t2/t3 a lo largo de
  la corrida. `rho` (núcleo `Sn` de ancho `drc` muestreado cada `dr`) solo
  se escribe como diagnóstico y no conserva la masa cuando drc≠dr: en t2
  (drc=0.225, dr=0.1) varía de -3.9% a +4.8%. La dinámica es consistente;
  `rho` no debe usarse para cantidades integradas.

- [x] **Sesgo de `cutoff`** (MEJORAS A4). El valor por omisión ya era 0;
  las entradas del artículo usan 0.01. Medido frente a `cutoff=0`:
  t1 (`aa`, sin autogravedad) 0.01 → |Δh_k|/|h_k| de 3.7e-3 (k=0) a
  7.5e-3 (k=4), 0.001 → 5e-4 a 1.3e-3; t2 (autogravedad) 0.01 → 1.7e-2
  a 2.5e-2, 0.001 → 1.3e-3 a 2.2e-3. A cambio, `cutoff=0` usa 24–34 veces
  más partículas (t1: 66/91/1575; t2: 160/272/5434) y tarda 10–18 veces
  más. `exe/input_parameters` conserva 0.01 (configuración del artículo).

- [x] **Estado `aa_quad`: cuadratura en (Q,J,L)** (MEJORAS A3). Rejilla de
  puntos medios en J (`Nrc`, [`jminc`,`jmaxc`], por omisión [0, 6 sr]),
  Q (`Npc`) y L (`Nlc`), invertida a (r,p_r) con `invert_QJ_to_rp`
  (`utils.f90`). Verificado: (1) ida y vuelta del mapa en 480 mil puntos
  (L=0.05, 0.6, 2.4; J hasta 3): Q a ≤4.6e-12, J a ≤6.7e-15; (2) sin
  autogravedad, yoshida4, dt=0.025, t≤200, N_J=200, N_Q=40, N_L=16:
  h_k del código frente a la suma discreta exacta sobre los mismos
  nodos, ≤2.6e-11 (integrador y mapa); suma discreta frente a la
  integral continua, 2.0e-4, dominado por N_L con orden 2 (5.1e-5 con
  32, 1.3e-5 con 64; N_J ya convergido en 200) porque C_f está truncada
  en lminc/lmaxc. Herramienta: `tools/hk_exacto.py`. t1/t2/t3 idénticos
  bit a bit.

- [x] **Estado `checkpoint`** (MEJORAS D1). Lee "r p_r L f" por partícula
  (exactamente Nrc·Npc·Nlc líneas) y normaliza la masa a a0 como
  `aa_quad`. Las instantáneas HDF5 ahora incluyen `l_part`, así que
  cualquier instantánea sirve de punto de partida. Verificado con
  `aa_quad` (16000 partículas, 2000 pasos, HDF5): sin autogravedad,
  recargar t=0 y reanudar desde el paso 1000 dan r y p idénticos bit a
  bit; f difiere en un factor 7.4e-15 (la renormalización de masa) y
  energía y h_k ≤8.5e-15. Con autogravedad ese factor entra a la
  fuerza: r, p ≤1.2e-14, h_k ≤4.2e-14. Archivos con menos o más líneas
  terminan con código 1.

- [x] **Integrador `analytic`** (MEJORAS C3). Avance exacto
  Q = Q0 + ω(J,L) t con J y L por partícula e inversión
  `invert_QJ_to_rp`; solo sin autogravedad, isócrono y eps=0, sin
  `reduceparticles`, y aborta si hay partículas no ligadas. Verificado
  con `aa_quad` (N_J=200, N_Q=40, N_L=16, t≤200): h_k frente a la suma
  discreta exacta ≤2.9e-12 (3.5e-14 en k=0); yoshida4 (dt=0.025) frente
  a analytic 1.2e-11–2.6e-11, que es el error del integrador.

- [x] **yoshida4 evaluaba 4 fuerzas por paso; la trayectoria usa 3**
  (MEJORAS C2, paso 1). La fuerza tras el último drift solo la usan los
  diagnósticos de cada `spatial_output` y `dt_switch=var`; ahora se
  calcula solo entonces (igual en `analytic`). Perfil previo (128 mil
  partículas `aa_quad`, 1000 pasos, 4 hilos): con autogravedad 92 s, de
  ellos `avg_density` 54 s, interpolación a partículas 20 s, fondo 8 s,
  centrífugo 4.6 s; sin autogravedad 15.6 s, fondo 7.2 s y centrífugo
  4.2 s. Verificado idéntico bit a bit: t1/t2/t3, `dt_switch=var`,
  leapfrog, analytic y todos los datos HDF5. ABBA (128 mil partículas):
  con autogravedad 25.04 → 19.50 s (1.28×), 300 pasos; sin autogravedad
  10.38 → 9.24 s (1.12×), 1000 pasos.

- [x] **Fondo y término centrífugo en pasadas seriales** (MEJORAS C2,
  paso 2). Para `Isochrone` y `Central` (y el centrífugo de los demás
  fondos) ahora es una sola pasada paralela por partícula, con
  sqrt(1+r²) y r²+eps² formados una vez y las mismas expresiones.
  Verificado idéntico bit a bit en t1/t2/t3 y en 14 casos más
  (Isochrone, Central, sphere, null, iso y nfw, con y sin autogravedad;
  leapfrog; analytic). ABBA (128 mil partículas): con autogravedad
  19.19 → 17.90 s, sin autogravedad 8.51 → 5.67 s (1.50×).

- [x] **Ventanas del núcleo más anchas que su soporte** (MEJORAS C2,
  paso 3). `avg_density` recorría 7 celdas y la interpolación de
  `poisson_rk` 7 nodos por partícula; `Wn` de orden n solo es distinto
  de cero en floor((n+2)/2) celdas a cada lado (1, 2, 2). Lo demás sumaba
  ceros exactos. Además el denominador de la cáscara se forma una vez por
  punto y `Wn` se evalúa una vez para potencial y fuerza. Verificado
  idéntico bit a bit: t1/t2/t3 y, con autogravedad, `bsplineorder` 1, 2 y
  3 en t2, t3 y `aa_quad`, más nfw y null. ABBA (128 mil partículas,
  autogravedad, 300 pasos): 17.49 → 12.42 s (1.41×).

- [x] **Avance de yoshida4 y reflexión en el origen seriales** (MEJORAS
  C2, paso 4). Seis pasadas de arreglo más dos copias por paso pasan a
  cuatro ciclos paralelos (cada kick fusionado con el drift siguiente,
  sin copias), y la reflexión en r=0 es paralela y evalúa `rmin==0` una
  vez. Verificado idéntico bit a bit: t1/t2/t3, bsplineorder 2,
  leapfrog, euler, `dt_switch=var`, `aa_quad` con autogravedad, null, y
  un caso con L∈[0,0.02] que cruza el origen (caótico: cualquier
  diferencia se amplificaría). ABBA (128 mil partículas): con
  autogravedad 12.40 → 11.58 s, sin autogravedad 5.59 → 4.34 s.
  Total de C2: con autogravedad 25.04 → 11.58 s (2.16×), sin autogravedad
  10.38 → 4.34 s (2.39×). Nota para A7: ese caso de L∈[0,0.02] da un
  error de energía de 1.4e14 (órbitas casi radiales por el origen).

- [x] **Mapa ángulo-acción numérico con L y equilibrio autoconsistente**
  (MEJORAS A2, D2). `tools/aa_numerico_L.py` (mapa directo e inverso en un
  potencial isócrono + tabla, L por partícula), `tools/equilibrio_L.py`
  (F_eq(J,L) con L en la rejilla del código, iterado con Poisson, y
  condición inicial "r p_r L F" para `checkpoint`) y
  `tools/hk_numerico.py` (h_k en variables verdaderas desde instantáneas
  HDF5). Verificado:
  - mapa sin tabla frente al isócrono analítico, 14556 partículas con
    L∈[0.1,3]: J a 8.9e-15, Q a ≤2.7e-8 (peor caso órbitas casi
    circulares, J≈1e-4); inversa con Newton: J a 1.2e-15, Q a 4.8e-9;
  - equilibrio a0=1e-2, L en 8 nodos de [1.6,2.4]: converge a 1.6e-13 en
    6 iteraciones; nodos recuperados a 2.8e-16 en J y 1.7e-9 en Q;
  - con el código (32000 partículas, autogravedad, t=100): Φ del código
    en t=0 frente al de Python, 1.3e-7 (|Φ_self|≤1.6e-3);
    max_r|Φ(t)-Φ(0)| 4.1e-7 en equilibrio, 4.3e-5 con eps=0.1 y 9.6e-4
    con el mismo F en acciones del isócrono (`aa_quad`);
  - h_k numérico sin autogravedad = h_k del código a ≤4.3e-11; con
    autogravedad en equilibrio, |h_1|/h_0 = 4.96e-3 constante con el mapa
    del isócrono (el artefacto de coordenadas) y 2.9e-7 → 1.6e-5 con el
    numérico; con eps=0.1, 0.04936 en t=0 frente a (a_1/a_0)·eps/2=0.0494
    esperado (isócrono: 0.0545).

- [x] **Rutas de HDF5 fijas a Debian/Ubuntu** (MEJORAS B5). El Makefile
  compila a través de `h5fc` si está en el camino y envuelve el
  compilador elegido; `make HDF5_WRAPPER=` usa las rutas explícitas.
  Verificado: sin wrapper el binario es idéntico byte a byte al anterior;
  con wrapper (aquí enlaza HDF5 estático, 3.8 MB) t1/t2/t3 idénticos y
  los 338 datasets y atributos de una salida HDF5 idénticos.

- [x] **Distribuciones de prueba y pruebas agnósticas** (MEJORAS D3).
  `distribution.f90` (portado del código de L fija) con `dftype` = gauss,
  bimodal, spiral, king, todas por C(L)=exp(-(L-l0)²/sl²) y king con
  E(J,L); `jmaxc` por omisión es la ventana de cada una.
  `tools/hk_exacto.py` generalizado a F0(Q,J,L) no separable. Sin
  autogravedad, `aa_quad`, L en [1.6,2.4] salvo donde se indica:
  - gauss: neutral bit a bit (t1/t2/t3 y `aa_quad` con y sin
    autogravedad); la herramienta repite 2.02e-4 / 2.6e-11 / 2.9e-12;
  - bimodal (N_J=200, N_Q=40, N_L=16, t≤200): |h_3|,|h_4| ≤ 5.2e-13
    max|h_1| con analytic y ≤ 7.6e-12 con yoshida4; k≤2 = suma discreta
    a ≤8.1e-13 (analytic) y ≤2.5e-11 (yoshida4);
  - spiral con L≈2 fijo (N_J=400): |h_k| llega al máximo en t=710 para
    k=1..4 (2.2, 25, 1.5e3 y 4.7e5 veces h_k(0)); la estimación lineal
    t*=-β/ω_J(J0)=721; código = discreto a ≤5e-11;
  - spiral con sl=0.2 (N_J=200, N_L=16, analytic): la dispersión en L
    borra el desenrollado, |h_1(700)|/|h_1(0)| = 2.18e-3 (con sl=0.02 el
    continuo da 2.0; con 0.05, 1.6); en t=1500 el código da 1.87e-2,
    igual a la suma discreta, frente a 5.0e-5 del continuo: recurrencia
    por muestreo en L (Nyquist en N_L);
  - king (N_J=100, N_L=16, t≤100): código = discreto a 1.3e-13 (analytic)
    y 5.2e-11 (yoshida4) en k=0,1; el error de cuadratura en J converge
    como h² (cociente 4.00 de N_J=50 a 1600, con L integrado exacto).

## Pendientes

- [ ] **`grav_force.f90`: condición `r_part(i)<1.d0` en el fondo
  `sphere`** sin `abs()` — cualquier partícula con `r_part` negativo
  transitorio (antes de la reflexión de simetría al final de cada
  paso) usa la fórmula interior sin importar su magnitud real.
- [x] **`grav_force.f90`: fondos `iso`, `isotrun`, `nfw`, `burkert` no
  actualizan `pot_part`/`force_part`**, solo los arreglos de malla
  `pot`/`force` — las partículas nunca sienten esa fuerza de fondo
  (solo el término centrífugo). Además esos arreglos solo se reservan
  si `autointeraction=.true.`, así que con `autointeraction=.false.`
  (el caso normal para un fondo fijo) se escribe en memoria no
  reservada.
  *Hecho (A9):* `add_background` en `grav_force.f90` aplica `sphere`,
  `iso`, `isotrun`, `nfw` y `burkert` a las partículas, los suma a la
  autogravedad (antes `sphere` la reemplazaba) y toca la malla solo si
  existe. Las parejas potencial/fuerza se verificaron con mpmath
  (F=-dΦ/dr a ≤1.3e-28). t1/t2/t3 (Isochrone) y `sphere`/`null` sin
  autogravedad, idénticos bit a bit. Con t1: antes las partículas
  escapaban (r_min≈405, error de energía 2.4e6); ahora quedan ligadas y
  el error máximo de energía es 3.2e-5 (iso), 7.7e-7 (isotrun), 4.6e-5
  (nfw), 8.1e-6 (burkert), con convergencia de cuarto orden al reducir
  dt (cocientes 15–16×; 10–26× en iso y burkert). `sphere` con
  autogravedad (t2): 1.43e-4 → 2.0e-7.
- [x] **`initial_data.f90` / `parameters.f90`: `state="gaussian"` (el
  valor por defecto) no coincide con ninguna rama** (`initial_data.f90`
  solo reconoce `"gaussian1"`, no `"gaussian"`), y no hay `else`/`stop`
  de captura — con la configuración de fábrica, `f` queda en cero en
  silencio.
  *Hecho (A9):* el valor por omisión es `gaussian1`; `read_parameters`
  rechaza estados desconocidos y `initial_data` tiene un `else` que
  aborta.
- [x] **`grav_force.f90` / `parameters.f90`: `forcetype="self"` es una
  opción documentada que no hace nada.** `main.f90` bifurca sobre
  `forcetype=="self"`, pero `grav_force.f90` solo mira `forcetype=="bg"`
  y el booleano independiente `autointeraction`. Con
  `forcetype="self"` y `autointeraction=.false.` las partículas no
  sienten fuerza radial (documentación/código desincronizados).
  *Hecho (A9):* `forcetype` solo acepta `bg`; la autogravedad es
  `autointeraction` y sin fondo es `BGtype=null`. La rama de `main.f90`
  queda en calcular densidad y energía cada `spatial_output` (neutral:
  era la rama que tomaban todas las corridas válidas).
- [ ] **Carpeta `src/`: archivos legado sin extensión `.f90`**
  (`analysish`, `poisson`, `poisson_ps`, `reduce_arrays`) no se
  compilan, están desincronizados de sus homónimos activos, y
  `poisson_ps` referencia un módulo `chebyshev` inexistente. Mover a
  `legacy/` o borrar.
- [x] **`initial_data.f90`: los estados `compact`, `compact2` y
  `Plummer` dejan `l_part` en cero para todas las partículas (nunca lo
  asignan), y `density()`/`energy()`/`analysish()` pesan todas sus
  sumas por `l_part(j)` — con `l_part≡0` esas tres rutinas devuelven
  **idénticamente cero** (`rho`, `avg_rho`, energía cinética/potencial/
  total, y los cinco $h_k$), aunque `f` sea distinto de cero. Mientras
  tanto, la normalización de masa de `compact`/`compact2`
  (`f = f*drc*dpc*8π²`, sin `dlc` ni `l_part`) sí asume un sistema sin
  dependencia en $L$ — inconsistente con lo que exigen las rutinas de
  diagnóstico. Encontrado comparando contra el código histórico
  (`Lfix` escalar), que manejaba exactamente este caso con una rama
  explícita `if (Lfix==0.0d0) then factor=1.0 else ... endif`
  ("Zero Angular Momentum" vs "Include Angular Momentum") en
  `density.f90`/`energy.f90` — esa rama sigue **comentada** en el
  código actual (nunca se restauró al migrar de `Lfix` escalar a
  `l_part` arreglo), y `compact`/`compact2`/`Plummer` nunca se
  actualizaron para poblar `l_part` con algo sensato. Bono menor:
  `Plummer` (línea ~112) además indexa `l_part(i)` con `i` el índice
  de malla en `r` (1..Nrc), no el índice real de partícula
  `(i-1)*Npc+j` — con `l_part` en cero en todos lados da lo mismo
  numéricamente, pero es la forma equivocada de indexar. Bug real pero
  dormido: ningún run de esta sesión usó estos tres estados (todo fue
  `aa`/`gaussian1`, que sí asignan `l_part` correctamente).
  *Hecho (A9):* retirados del código activo a
  `legacy/initial_data_estados_2D.f90`; `read_parameters` los rechaza.
  Recuperarlos exige decidir su dependencia en L (y `Plummer` usaba la
  energía del isócrono).
- [x] **`eps` (longitud de suavizado del término centrífugo) está fija
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
  **Confirmado como regresión real** comparando contra código
  histórico rescatado por el usuario
  (`/home/erik/Documentos/old_VlasovPoisson_PIC_sp/`, versión previa
  al refactor 2D→3D del estado `aa`): ahí `eps` sí se calculaba y se
  usaba activamente,
  ```fortran
  eps = Lfix/(10.0D0*pmax)                                    ! utils.f90
  pot_part   = pot_part + 0.5d0*Lfix**2/(r_part**2 + eps*eps)  ! grav_force.f90
  force_part = force_part + Lfix**2*r_part/(r_part**2+eps*eps)**2
  ```
  y también se usaba en `initial_data.f90` para calcular la energía
  inicial. El suavizado se perdió en algún punto del refactor de
  `Lfix` (escalar) a `l_part` (arreglo distribuido en la malla
  `(r,p,l)`) — nadie escribió el equivalente con `l_part`. Pendiente
  de decidir la forma funcional correcta de `eps` en términos de
  `l_part`/`pmax` para reintroducirlo.
  *Decidido (MEJORAS A7):* se mantiene `eps=0`. Con `eps≠0` el
  potencial deja de ser el del mapa ángulo-acción (en el código de L fija
  creó un piso de h_k de ~5e-12). Con los estados activos L ≥ lminc+dlc
  > 0 (los estados con `l_part=0` se retiraron). Medido con `aa`, sin
  autogravedad, yoshida4, Nlc=8, cutoff=0: L en [0,0.4] (r_min=0.061)
  da max|dE/E| 6.35e-8 con dt=0.025 y 4.23e-9 con dt/2 (15×); con L en
  [0.4,0.8], [1,1.4] y [1.8,2.2], 0 a 8 cifras. El caso del ~350% no se
  puede reproducir (sus parámetros no quedaron registrados). Si se
  estudian L muy pequeños, verificar con una prueba de dt. El valor por
  omisión `lminc=-2` (inválido) pasó a 0.
- [x] **Investigado y descartado: ¿el corte de la ventana de depósito
  en `density()` (`bsplineorder*drc`/`bsplineorder*dr`) recorta la
  cola del B-spline respecto al código histórico
  (`(bsplineorder+1)*drc`/`dr`)?** No. `functions.f90` (donde viven
  `Sn`/`Wn`) es idéntico entre el código viejo y `main`/`fix` — el
  soporte compacto real es $|y|<n/2$ para $S_n$ y $|y|<(n+1)/2$ para
  $W_n=S_{n+1}$ (leído directamente del código). El corte actual
  (`n·dr_c` para $S_n$, `n·dr` para $W_n$) sigue siendo **el doble**
  del radio de soporte real en los 4 órdenes soportados — nunca
  recorta nada; el `+1` del código viejo solo era margen extra sin
  necesidad matemática. Verificado con la tabla soporte-vs-corte para
  $n=1..4$ ($S_n$) y $n=1..3$ ($W_n$). Única nota menor (sin efecto
  numérico): dentro de `density.f90` el corte de `avg_rho` usa
  `bsplineorder*dr` mientras que la subrutina separada
  `avg_density()` todavía usa `(bsplineorder+1)*dr` — inconsistencia
  de estilo entre dos rutinas que calculan lo mismo, ambos cortes son
  seguros igual.
- [ ] **`functions.f90`/`density.f90`/`poisson_rk.f90`: `Wn` asume
  implícitamente `drc==dr`, pero ningún `input_parameters` real
  cumple eso.** El paper define
  $W_m(R_k-r_j):=\int S_m(r-r_j)\,b_0\!\left(\frac{r-R_k}{\Delta
  r}\right)dr$ — la convolución genuina entre el shape function de la
  partícula (ancho $\Delta$, que en el código es `drc`, ver
  `Sn(bsplineorder,(r(i)-r_part(j))/drc,drc)` en `density()`) y la
  caja de la malla espacial (ancho $\Delta r$=`dr`). El comentario en
  `functions.f90` ("the weight function W_n are just the S_(n+1)
  function multiplied by dx") sólo es válido cuando `drc==dr` — ahí
  la convolución colapsa a la recursión estándar de B-splines. Pero
  `avg_density()`/`poisson_rk()` llaman `Wn(bsplineorder,(r-r_part)/dr)`,
  es decir usan como ancho del shape function únicamente `dr`,
  **ignorando por completo `drc`**. Confirmado con los
  `input_parameters` reales: `input_parameters` tiene `dr=0.03125`
  vs `drc=(rmaxc-rminc)/Nrc=(5-1)/400=0.01`; `input_article_N1e3`
  tiene `dr=0.1` vs `drc=(10-1)/240=0.0375` — en ambos casos `drc` y
  `dr` difieren por un factor ~2-3×, no son iguales en ningún caso de
  uso real inspeccionado. El resultado es que la misma partícula, con
  el mismo shape function nominal de ancho `drc`, se trata con dos
  anchos distintos según qué cantidad se calcule: `rho` (diagnóstico
  puntual, `vlasov_density`) usa correctamente `drc`, pero `avg_rho`
  — la fuente del lado derecho de la ecuación de Poisson, y también
  el kernel usado para interpolar `pot`/`force` de vuelta a las
  partículas en `poisson_rk.f90` — usa `dr` en su lugar. No es
  catastrófico (la integral/masa total de la convolución no depende
  de qué ancho se use, sólo la forma), pero sí afecta la resolución
  espacial real del depósito y de la fuerza autogravitante sentida
  por las partículas: con `drc` bastante más chico que `dr` (caso
  típico visto arriba), el kernel efectivamente usado es más
  ancho/suave que el que debería representar a la partícula.
  Pendiente de decidir la forma correcta de `Wn` cuando `drc≠dr`
  (convolución explícita de dos anchos distintos, no la recursión de
  un solo ancho) y cuantificar el efecto numérico antes de corregir.
  Encontrado en una revisión de las ecuaciones del código contra
  `Vlasov_Poisson_evolutions/main.md` (el paper de referencia).

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
- [x] **Paralelizar la generación de partículas iniciales** en
  `initial_data.f90` (estados `aa`/`gaussian1`) — el índice `indx` se
  incrementaba a mano (dependencia secuencial que bloqueaba el
  `!!$OMP` ya escrito pero deshabilitado); reemplazado por la fórmula
  cerrada `indx=(k-1)*Nrc*Npc+(i-1)*Npc+j`, que hace que cada hilo
  escriba solo en su propio índice (a diferencia del bug de
  `collapse(2)` en `density()`, acá no hay accumulation compartida,
  es un scatter puro — seguro de colapsar del todo). Validado:
  datos idénticos bit a bit contra la versión serial (gaussian1 y aa),
  y determinístico entre 1/4/8 hilos. A escala de producción (2M
  candidatos, estado `aa`): **~1.7× más rápido** (0.53s→0.31-0.32s,
  con retornos decrecientes más allá de los 4 núcleos físicos, como
  se esperaba). — commit `perf(initial_data): parallelize gaussian1/aa
  candidate generation`
- [x] **Intentado y revertido: cachear `Sn(...)` en `density()` +
  evitar copiar `r_part_p`/`p_part_p` cada paso (ping-pong en
  `main.f90`).** Ambos implementados, validados bit a bit contra la
  versión anterior en 4 escenarios (euler, leapfrog, yoshida4,
  autogravitante) — correctos. Al medir tiempo: **cada uno por
  separado, sin regresión medible** (10 corridas intercaladas cada
  uno contra la referencia, dentro del ruido). **Combinados: ~20% más
  lento, de forma consistente y reproducible** (10 corridas
  intercaladas, cero superposición entre las dos distribuciones —
  23-26s vs 29.5-31s en una corrida de 3000 pasos autogravitante).
  No encontré una explicación clara — probablemente alguna interacción
  del optimizador de gfortran entre las dos subrutinas modificadas a
  la vez (parecido en espíritu a la pesimización de arreglos a nivel
  de módulo que ya vimos con el intento de reuso de buffers, pero acá
  ninguno de los dos cambios toca arreglos de módulo). Revertidos
  ambos por completo — ninguno de los dos demostró una ganancia
  medible por sí solo como para justificar el riesgo de la
  interacción al combinarlos.
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
- [x] **Vectorizar/agrupar la cuadratura de `phik` en `analysish.f90`.**
  `phik` integraba $g(J,l,Q)\cdot e^{-i\cdot\text{modo}\cdot Q}$ por
  Simpson (21 puntos), llamada una vez por (partícula, modo) — 5 veces
  por partícula — recalculando $g(J,l,Q)$ (la parte que NO depende del
  modo) desde cero cada vez. Reestructurado a partícula-externo,
  calculando $g$ una sola vez por partícula y reusándola en los 5
  modos. De paso: el resultado de `phik` siempre fue real (el código
  original se quedaba con `real(...)`, descartando la parte imaginaria
  en silencio) — reemplazado `exp(-i·modo·Q)` por `cos(modo·Q)`
  directamente, evitando calcular el `sin` que de todos modos se
  tiraba. También paralelizado sobre partículas con
  `REDUCTION(+:hk)` (acumulador de solo 5 elementos, sin atomics).
  Encontré y corregí un bug propio antes de integrarlo (el extremo
  $Q=\pi$ necesita el factor $\cos(\text{modo}\cdot\pi)=(-1)^{modo}$,
  no asumir que es 1 como en $Q=0$) probando la cuadratura aislada
  contra una referencia en Python del algoritmo *original* antes de
  tocar el archivo real.
  Medido (200k partículas, sin autogravedad, 10 llamadas a
  `analysish`, 1 hilo): mediana 14.2s→9.5s, **~1.5× más rápido** solo
  por reusar $g(Q)$; más hilos no dieron ganancia adicional en este
  benchmark pero tampoco regresión. — commit `perf(analysish): reuse
  the mode-independent quadrature weight across modes`
  **CORRECCIÓN (encontrada comparando main vs fix en una corrida
  completa):** la afirmación de "hk.tl idéntico byte a byte contra la
  versión anterior" era **incorrecta** — de hecho no lo era, y esto
  reveló que la reescritura corrigió, sin que me diera cuenta en su
  momento, un **bug real preexistente** (presente en `main` y en todo
  el historial de `fix` hasta este commit): el código original llamaba
  `phik(Jr(j),l_part(j),l0,mode,sp,sr,sl)` dentro de
  `do i=0,mode ... end do`, pasando la variable **`mode`** (constante,
  siempre 4, el número total de modos) en vez de **`i`** (el modo
  actual de la iteración). Por lo tanto `phik` se evaluaba con
  $\cos(4\cdot Q)$ (o $e^{-i4Q}$) para **todos** los modos 0..3, no con
  $\cos(i\cdot Q)$ — solo el modo 4 (el último) daba el resultado
  correcto por coincidencia. Confirmado reconstruyendo el binario justo
  antes de este commit (`2f345bc^`) y comparando `hk.tl` contra el
  binario actual con una corrida corta idéntica: h0..h3 difieren
  (p.ej. h0: 4.06344320E-06 → 4.27107028E-06, ~5.1%), h4 es idéntico
  (como se espera, ya que ahí `mode`==`i`==4 en el código viejo). Esto
  también explica exactamente la diferencia de ~5% en h0 vista en la
  comparación main-vs-fix de una corrida larga (100k pasos, sin
  autointeracción): main nunca tuvo esta corrección. No es un bug
  nuevo en la lista — es una reclasificación: era un bug real en
  `phik`, no solo una optimización, y ya está corregido en `fix`.

## Mejoras de diseño (no son bugs)

- [x] **Forma más amigable de pasar los parámetros de la simulación.**
  *Hecho (MEJORAS B1):* `paramfile.f90` lee `nombre = valor` con overrides en la
  línea de comandos, valida opciones y escribe `params_usados.par`;
  `tools/posicional_a_par.py` convierte los archivos viejos.
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

- [x] `test_consistency` (`utils.f90`) nunca se llama desde `main.f90`
  y referencia `lmin`/`lmax` que no existen (son `lminc`/`lmaxc`) — no
  compilaría si se descomenta tal cual.
  *Hecho:* se corrigió y activó en `506a3df`; con B1 sus comprobaciones pasaron a
  `validate` en `paramfile.f90` y la subrutina se eliminó.
- [ ] Integrador `rk4` declarado como opción válida pero no
  implementado (aborta con mensaje).
- [ ] Limpiar el código muerto/comentado en el bucle principal de
  `main.f90` (ecuación de continuidad, paso de tiempo adaptativo para
  autogravedad, etc.).
