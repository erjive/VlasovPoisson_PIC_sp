# Auditoría agnóstica del código — 2026-09-20

Revisión completa de `src/` buscando implementaciones no probadas, no justificadas o
que rompan la física que se quiere describir. Las afirmaciones marcadas como *medido*
se comprobaron replicando el esquema (depósito, RK2 de Poisson, interpolación) en un
programa independiente, no leyendo el código.

Estado de cada punto. **Reevaluado el 2026-09-20 tras corregir los seis primeros**: dos
enunciados no sobrevivieron a su propia prueba de aceptación y uno se confirmó midiendo.
Los puntos se escribieron leyendo el código; solo la medición decide.

| # | punto | estado |
|---|---|---|
| 1 | autofuerza | **corregido**, medido |
| 2 | masa perdida cerca del origen | **corregido**, medido |
| 3 | `BGtype=null` acumulaba fuerza | **corregido**, medido |
| 4 | paridad de los fondos en r | **corregido**, medido |
| 5 | radicandos sin proteger | **corregido**, medido |
| 6 | arreglos de tamaño `Npart` | **corregido**, pero el riesgo enunciado no existía |
| 7 | `courant > 1` sin verificar | pendiente; **peor de lo dicho**: 6 archivos, no 1 |
| 8 | paso de tiempo sin frecuencia orbital | pendiente; **demostrado** al corregir el 4 |
| 9 | mapa ángulo–acción triplicado | pendiente, sigue en tres archivos |
| 10 | pruebas de verificación degeneradas | **parcialmente atendido** |
| 11 | "densidad uniforme exacta" con n=1 | pendiente, la afirmación sigue en el fuente |
| 12 | `cutoff` no reporta la masa descartada | pendiente |
| 13 | constantes de forma compiladas | pendiente |
| 14 | nombres sobrecargados | pendiente; la mitad de `eps` decae con el punto 19 |
| 15 | dos ramas del isócrono | **corregido**, medido |
| 16 | fondo en los puntos fantasma | pendiente |
| 17 | Newton sin aviso | pendiente |
| 18 | lista de celdas y factor de Courant | pendiente |
| 19 | `eps` | **enunciado erróneo**, severidad rebajada |
| 20 | convenio del depósito | medido, decisión pendiente |
| 21 | el equilibrio no lo es del sistema discreto | **confirmado**, medido por tres vías |


---

## SEVERIDAD ALTA

### 1. El esquema tenía autofuerza, y el código afirmaba lo contrario
**Estado: CORREGIDO** (2026-09-20)

`functions.f90:7` y la cabecera de `density.f90` afirman que usar el mismo `W_n` para
depositar e interpolar evita la autofuerza. *Medido*: una cáscara aislada de masa 1
siente sobre sí misma

| n | r_j | F_interp(r_j) | −m/(2r_j²) | cociente |
|---|---|---|---|---|
| 1 | 0.5 | −1.980e0 | −2.000e0 | 0.990 |
| 1 | 1.0 | −4.975e−1 | −5.000e−1 | 0.995 |
| 1 | 2.0 | −1.2469e−1 | −1.250e−1 | 0.998 |
| 2 | 1.0 | −4.975e−1 | −5.000e−1 | 0.995 |
| 3 | 1.0 | −4.971e−1 | −5.000e−1 | 0.994 |

Es exactamente la autogravedad clásica de una cáscara. El argumento de cancelación por
paridad que vale en PIC cartesiano no vale en simetría esférica: el campo de una cáscara
en su propio radio es la mitad del exterior, no cero.

Es un sesgo **sistemático hacia adentro**, no ruido; escala como 1/N y no mejora
refinando la malla. La corrección estándar (restar el campo propio de la partícula) no
está aplicada ni mencionada.

Magnitud en `eq_e01.par` (a0=0.01, Npart=32000): ~2e−5 de la fuerza de autogravedad.
No explica la meseta ya medida, pero contamina los estudios de convergencia en N,
porque introduce un término 1/N sistemático donde se busca ruido de muestreo.

**Por qué aparece.** Con w_i los pesos del depósito y K(i',i) la fuerza en el punto i'
por masa unidad en i, la autofuerza es la forma cuadrática m_j Σ w_i' K(i',i) w_i. En
malla cartesiana K ∝ −sgn(i'−i) es antisimétrica y la forma cuadrática se anula: ese es
el teorema que justifica usar el mismo W para depositar e interpolar. En simetría
esférica el teorema de las capas da K(i',i) = −H(i'−i)/r_i'² con H el escalón, y
H = ½ + ½ sgn: la parte sgn se cancela igual que antes, pero la parte ½ sobrevive y deja

    F_auto = −(m_j/r_j²)·½·(Σ w_i)² = −m_j/(2 r_j²)

porque los pesos suman uno. Eso explica que el valor medido no dependa ni del orden del
B-spline ni de la posición dentro de la celda.

**Por qué no basta restar la forma cerrada.** Conviene recordar primero que la autofuerza
no aproxima nada: en Vlasov vale cero exactamente, y −m/(2r²) es la autogravedad de una
cáscara de masa finita, o sea el mismo artefacto idealizado. Lo que sí debe converger al
orden del esquema es el campo de una densidad suave, y converge (ver el punto 2).

La autofuerza discreta se aparta de −m/(2r²) en **primer orden** en dr/r, porque es una
integral sobre el ancho propio de la partícula, y a lo largo de ese ancho varían en
O(dr/r) tanto la medida 4πr²dr (la mitad exterior del soporte pesa más) como el 1/r² con
que se interpola de vuelta. Las dos reducen |F_auto|. Medido con la partícula en r=1 y
misma fase en la celda, el cociente (1 − F/F_ref)/(dr/r) para dr = 0.2, 0.1, 0.05, 0.025,
0.0125 da 0.709, 0.765, 0.798, 0.815, 0.824 con n=1, es decir tiende a 5/6 ≈ 0.833 (para
n=3 tiende a ≈ 8/9). El coeficiente depende además de la fase dentro de la celda.

Restar la forma cerrada dejaría entonces un residuo de (5/6)(dr/r) veces el artefacto
— con dr=0.05 y r≈1, un 4 % — que sigue escalando como 1/N. Se resta **exacta**: la densidad propia de la partícula se
construye sobre los puntos de su soporte con los pesos del depósito, su masa encerrada
con los mismos `mcoefA`/`mcoefB` que usa el solver, y el resultado se interpola de vuelta
con el mismo W_n. La imagen en −r_j pertenece a la misma partícula y entra también.

El autopotencial se resta igual, porque si no la energía deja de corresponder a la
dinámica: lo que queda es la suma sobre pares con j ≠ k.

Prueba de aceptación en el binario — una partícula con L=1, r₀=1, p₀=0, sin fondo, con
autogravedad, cuya solución exacta es r(t) = √(1+t²). A t = 10, r exacto = 10.04987562112:

| | r(10) | error relativo |
|---|---|---|
| con autofuerza | 10.010672361 | 3.9e−3 |
| sin autofuerza | 10.049875621 | **8e−12** |

Efecto sobre `10_energia/aa_sg.par`:

| | max\|ΔE/E₀\| | \|h₁\| final |
|---|---|---|
| RK2 original | 1.9371e−07 | 8.43542392e−08 |
| + punto 2 | 2.4851e−07 | 8.43542124e−08 |
| + punto 1, solo la fuerza | 1.0535e−05 | 8.43527526e−08 |
| + punto 1, fuerza y potencial | **2.0773e−07** | 8.43527526e−08 |

La fila tercera es el diagnóstico de que la resta del potencial hace falta: quitar la
fuerza sin quitar la autoenergía empeora la conservación 40 veces. Con las dos, la
conservación queda mejor que en el código original, y eso explica de paso el
empeoramiento que se había anotado en el punto 2: parte de él era la misma
inconsistencia, en menor grado.

|h₁| cambia 1.8e−5 en esa corrida. El costo añadido no se distingue de la variación entre
corridas. La salida sigue siendo idéntica bit a bit entre corridas con el mismo número
de hilos.

### 2. Poisson perdía la masa de las partículas cercanas al origen
**Estado: CORREGIDO** (2026-09-20)

*Medido*: masa que ve el campo, M = −F(Nr)·r(Nr)², para una partícula de masa 1.

| n | r_j/dr | M vista | error |
|---|---|---|---|
| 1 | 0.5 | 0.000000 | −100 % |
| 1 | 1.0 | 0.583 | −42 % |
| 1 | 2.0 | 1.096 | +9.6 % |
| 1 | 5.0 | 1.00399 | +4.0e−3 |
| 1 | 20.0 | 1.000202 | +2.0e−4 |
| 2 | 20.0 | 0.999994 | −6.5e−6 |

Una partícula sobre el primer punto de malla desaparece por completo del campo con CIC.
Verificado también a mano: con ρ solo en la celda 1, el arranque
`dev_pot(1) = (4/3)πρ₀r(1)` (`poisson_rk.f90:84`) y el primer paso RK2 dan
`dev_pot(2) = 0` exactamente, y sin más fuente el campo es nulo en toda la malla.

Lejos del origen el error es O((dr/r)²), la convergencia de segundo orden esperada; el
problema se limita a r ≲ 5dr. No se dispara en producción porque L ∈ [1.6, 2.4] mantiene
los pericentros en r ≳ 20dr, pero invalida cualquier corrida con autogravedad y órbitas
de L pequeño o radiales.

**Corrección.** `poisson_rk.f90` ahora integra la masa encerrada en lugar de dPhi/dr:

    dM/dr = 4 pi r**2 rho,      dPhi/dr = M/r**2

con rho reconstruida como la recta entre valores de malla, de modo que M es una cuártica
y ambas integrales se hacen en forma cerrada sobre cada intervalo. Desaparece el término
2/r y con él la rigidez en el origen. Para el peso lineal (bsplineorder = 1) la recta
reproduce exactamente el perfil depositado, así que cada partícula aporta toda su masa
sea cual sea su radio.

Medido en el binario, no en la réplica (una partícula de masa 1, `dr` = 0.01):

| n | r_j/dr | antes | después |
|---|---|---|---|
| 1 | 0.5 | −0.00000000 | **1.00000000** |
| 1 | 1.0 | 0.58341142 | **1.00000000** |
| 1 | 2.0 | 1.09610630 | **1.00000000** |
| 1 | 20.0 | 1.00020220 | **1.00000000** |
| 2 | 0.5 | 0.14099109 | 0.85000000 |
| 2 | 20.0 | 0.99999356 | 0.99979141 |
| 3 | 0.5 | 0.18192399 | 0.75115207 |
| 3 | 20.0 | 0.99978488 | 0.99958264 |

Para n ≥ 2 la reconstrucción lineal no reproduce el B-spline, así que queda un error
O((dr/r)²) del mismo orden que antes, pero sin el fallo catastrófico cerca del origen.

Esfera uniforme de masa 1 y radio 5 (2000 cáscaras, `dr` = 0.05), error relativo de la
fuerza:

| | antes | después |
|---|---|---|
| dentro (r < 4.95) | 5.0006e−4 | 5.0006e−4 |
| fuera (r > 5.25) | 8.8595e−5 | **4.3144e−9** |
| masa total vista | 0.9999937120 | **1.0000000003** |

El error interno es idéntico porque lo domina el sesgo de muestreo de CIC (punto 11),
no el solver.

Efecto sobre la física, corrida `10_energia/aa_sg.par` (a0=1e−3, autogravedad, 4000 pasos):

| | antes | después |
|---|---|---|
| E(0) | −7.846893890e−05 | −7.846891980e−05 |
| max\|ΔE/E₀\| | 1.9371e−07 | 2.4851e−07 |
| \|h₁\| final | 8.43542392e−08 | 8.43542124e−08 |

Los resultados cambian ~3e−7 en esa corrida. La conservación de energía queda
ligeramente peor (mismo orden); no se investigó el motivo, y conviene vigilarlo.

Las corridas sin autogravedad no cambian en absoluto: `poisson_rk` solo se llama desde
`grav_force` bajo `if (autointeraction)`.

### 3. `BGtype="null"` sin autogravedad acumulaba fuerza sin límite
**Estado: CORREGIDO** (2026-09-20)

En `grav_force.f90:119-141`, con `BGtype="null"` y `autointeraction=.false.` no se
ejecuta `add_background` ni ninguna otra asignación, y el bucle siguiente hace
`pot_part(i) = pot_part(i) + centrífugo` sobre valores nunca reinicializados. Con
yoshida4 son 3 o 4 acumulaciones por paso. `null` es una opción aceptada
(`paramfile.f90:387`) y corresponde a la prueba exacta más elemental del integrador:
partículas libres con momento angular.

**Corrección.** La rama añade el caso `null` sin autogravedad y pone ambos arreglos a
cero antes de sumar el término centrífugo. Con autogravedad no hace falta, porque
`poisson_rk` acaba de asignarlos.

Prueba de aceptación: una partícula con L=1, r₀=1, p₀=0, `BGtype=null`,
`autointeraction=.false.`, solución exacta r(t) = √(1+t²). A t = 10, r exacto =
10.04987562112:

| | r(10) |
|---|---|
| antes | 6972.33 |
| después | **10.049875621** |

---

## SEVERIDAD MEDIA

### 4. Los fondos NFW y Burkert no tenían la paridad correcta en r
**Estado: CORREGIDO** (2026-09-20)

El código deja que r sea negativo dentro de un paso (la reflexión está en `main.f90:239`,
después del integrador) y la interpolación de malla lo maneja con el espejo; los fondos
analíticos no. *Medido*, fuerza en r = −0.05 contra la extensión impar correcta:

| fondo | F(−0.05) | −F(+0.05) | |
|---|---|---|---|
| sphere, isotrun | | | OK |
| iso | | | fuerza OK; el **potencial** `3 log(x)` da NaN para r<0 |
| nfw | −8.565 | +7.495 | signo invertido |
| burkert | 0.2306 | 0.2139 | discrepa 8 % |

Solo se dispara con L≈0, que es el régimen en que se usaría una cúspide.

**Corrección.** `bgpot` y `bgforce` evalúan ahora las formas cerradas en |r| y la fuerza
lleva el signo de r. El fondo `Central`, que vive fuera de esas dos funciones, se arregló
igual (su potencial era impar y su fuerza par, las dos al revés de lo que deben ser).

Prueba de aceptación: una partícula con L ≈ 0 soltada en reposo en r = 1 en el fondo
`nfw`, sin autogravedad, que oscila cruzando el origen unas treinta veces. Error relativo
máximo de la energía:

| | max\|ΔE/E₀\| |
|---|---|
| antes | 6.32e−01 |
| después | **4.00e−02** |

El 4 % que queda no es paridad: es el paso de tiempo, que no conoce la frecuencia orbital
de una órbita que se hunde en una cúspide (punto 8).

### 5. Raíz sin proteger en `analysish`: un NaN de órbita casi circular envenenaba los h_k
**Estado: CORREGIDO** (2026-09-20)

`analysish.f90:81-82` calcula `sqrt(-L²-2E-2-0.5/E)` sin protección; ese radicando vale
cero exactamente en una órbita circular y el redondeo puede volverlo negativo. La misma
expresión sí lleva `max(...,0)` en `utils.f90:733` y `utils.f90:781`. El filtro del bucle
solo descarta `energy(j) >= 0` (`analysish.f90:147`), así que una partícula ligada con NaN
contamina el acumulador completo. `initial_data.f90:122` tiene el mismo hueco, enmascarado
por un `f /= f` que atribuye el NaN a nodos no ligados cuando puede ser una órbita
circular válida a la que se le pone masa cero en silencio.

**No era un riesgo teórico.** Evaluando los radicandos para 500 órbitas exactamente
circulares del isócrono, con radios repartidos en [0.5, 5]:

| radicando | negativos por redondeo | mínimo |
|---|---|---|
| discriminante 1+2E(2+2E+L²) | **87 de 500** | −4.4e−16 |
| excentricidad −L²−2E−2−0.5/E | **120 de 500** | −1.8e−15 |

Es decir, una de cada cinco órbitas circulares producía un NaN.

**Corrección.** Todos los radicandos van con `max(...,0)`, la división por s₂−s₁ se
protege (vale cero exactamente en una órbita circular) y las partículas no ligadas se
saltan antes de construir el mapa en vez de después. El caso degenerado es inocuo por sí
mismo: en una órbita circular la acción radial se anula y la función de prueba lleva un
factor J², así que esa partícula no aporta nada sea cual sea el ángulo que se le dé.

En `initial_data.f90` (estado `aa`) se separan los dos casos que antes tragaba el mismo
`f /= f`: los nodos no ligados se cuentan y se reportan, y los casi circulares dejan de
perder su masa en silencio.

Prueba de aceptación: las 500 órbitas circulares de arriba, en una corrida real.

| | h₀ | h₁ … h₄ |
|---|---|---|
| antes | 5.04e−39 | **NaN** |
| después | 5.04e−39 | 2.42e−40, 1.52e−39, 3.24e−41, 1.12e−40 |

h₀ sobrevivía porque no usa el ángulo; a partir de k=1 el factor exp(−ikQ) propagaba el
NaN a todo el acumulador.

### 6. Arreglos automáticos de tamaño `Npart` (el riesgo que enuncié no existe)
**Estado: CORREGIDO, pero por otra razón que la que dije** (2026-09-20)

`analysish.f90` declaraba nueve arreglos automáticos de longitud `Npart`. Escribí que
eran 72 MB de pila con 10⁶ partículas y un segfault dependiente del `ulimit -s`.
**Eso no se reproduce**: con 10⁶ partículas y el límite de pila por omisión (8192 kB) la
versión con arreglos corre sin problema, porque gfortran con estas banderas no los pone
en la pila.

Lo que sí es real es el consumo. Con 10⁶ partículas, pico de memoria del proceso:

| | pico | |
|---|---|---|
| con arreglos | 152 420 kB | |
| sin ellos | **90 044 kB** | −41 % |

**Corrección.** Nada se comparte entre partículas, así que el mapa ángulo–acción se
calcula con escalares dentro del lazo paralelo y los nueve arreglos desaparecen. La
salida es **idéntica bit a bit** (hk1.tl, hk1\_complex.tl, hk2.tl y la energía, corrida de
equilibrio con 8000 partículas y 8000 pasos).

El precio es tiempo: con 10⁵ partículas y 200 llamadas a `analysish`, el mejor de tres
pasa de 4.74 s a 5.10 s, un 7 % más. Las expresiones sobre arreglos completos
vectorizaban mejor que el lazo escalar con sus ramas. En producción `analysish` corre
cada `spatial_output` pasos (40 en las corridas de referencia), así que ese 7 % es una
fracción pequeña del total, y a cambio el programa ocupa un 41 % menos en el caso que
importa para escalar N.

### 7. `courant = 80.0` en un archivo de parámetros distribuido
**Estado: PENDIENTE**

`reproducir/video/df_gauss.par` fija `courant = 80.0` para un parámetro documentado como
fracción de Courant (`parameters.f90:45`). Solo tiene sentido con `integrator = analytic`,
donde dt únicamente espacia las instantáneas. Nada lo verifica: cambiar el integrador en
ese archivo da dt = 4 sin aviso.

**Reevaluación:** es más extendido de lo que escribí. Contando los `.par` distribuidos,
`courant` vale 20.0 en cinco de ellos y 80.0 en uno; solo catorce usan 0.5 y ocho valores
menores. Seis archivos, no uno.

### 8. El criterio de paso de tiempo desconoce la frecuencia orbital
**Estado: PENDIENTE (documentado en el fuente)**

`utils.f90:251-255` lo documenta, incluido el caso nfw que divergió. La prueba de
convergencia queda a cargo del usuario y nada la exige.

**Reevaluación:** dejó de ser una objeción teórica. La prueba de aceptación del punto 4
— una órbita radial en el fondo `nfw` — conserva la energía al 4 % con el paso que el
código elige solo, y ese 4 % es enteramente el paso de tiempo. Es el mismo número que
mediría cualquiera que use una cúspide.

### 9. El mapa ángulo–acción está triplicado
**Estado: PENDIENTE**

La misma transformación directa en `initial_data.f90:107-124`, `utils.f90:763-785` y
`analysish.f90:63-82`, con tres implementaciones que ya divergieron (la protección `max`).

---

## SEVERIDAD BAJA, Y DUDAS

### 10. Las pruebas de verificación son degeneradas respecto de lo que deben probar
**Estado: PARCIALMENTE ATENDIDO**

La esfera uniforme es exactamente el caso que la fórmula
V_i = 4πdr(r_i² + (n+1)dr²/12) hace exacto por construcción, y también el caso en que el
arranque del RK2 (densidad constante dentro de r(1)) es exacto. Un esquema con los
defectos 1 y 2 pasa esa prueba sin marcas. Hacen falta casos que el diseño no privilegie:
cáscara única cerca del origen, órbitas radiales, convergencia en N a dr fijo.

**Reevaluación:** el diagnóstico está ya escrito en §7.2 de
`docs/introduccion/vlasov_L_intro.tex`, y las correcciones 1 a 6 produjeron justamente
las pruebas que faltaban: cáscara aislada a distintos radios, movimiento libre
r(t)=√(1+t²) con y sin autogravedad, órbita radial que cruza el origen en un fondo con
cúspide, 500 órbitas exactamente circulares, esfera uniforme, Plummer, y un millón de
partículas. Lo que falta es **moverlas al repositorio**: hoy viven en el directorio de
trabajo de la sesión, no en `reproducir/`, así que nadie más las puede repetir.

### 11. "Una densidad uniforme se reproduce exactamente" es falso para el orden por omisión
**Estado: PENDIENTE**

*Medido*, densidad uniforme ρ₀=1 muestreada con espaciamiento D, máx |ρ_i/ρ₀ − 1|:

| n | D=dr | D=dr/2 | D=dr/4 |
|---|---|---|---|
| 1 | 1.79e−4 | 2.24e−5 | 5.60e−6 |
| 2 | 4.9e−15 | 4.7e−15 | 3.3e−15 |
| 3 | 4.7e−15 | 4.7e−15 | 2.9e−15 |

La afirmación de `density.f90:15` vale a precisión de máquina para n ≥ 2 y falla para
n = 1, que es el valor por omisión (`parameters.f90:81`) y el de producción. El segundo
momento *discreto* del B-spline lineal oscila entre 0 y 1/4 según la posición dentro de
la celda, en vez de valer 1/6; para n ≥ 2 sí es constante. El error es una modulación con
la periodicidad de la malla, coherente y no aleatoria, de tamaño 2e−4.

### 12. El corte `cutoff` altera la DF y nunca se reporta cuánto
**Estado: PENDIENTE**

`initial_data.f90:202` descarta nodos y renormaliza la masa a a0 sin imprimir la fracción
descartada. La truncación recorta las colas en J y en L, o sea cambia var_J y var_L, que
son las cantidades que predicen la envolvente de mezcla. Verificar que
`tools/hk_exacto.py` aplique el mismo corte.

### 13. Constantes de forma compiladas
**Estado: PENDIENTE**

`spi_beta = 50`, `bim_Ja = 0.10`, etc. en `distribution.f90:30-56` son `parameter`.
Cambiarlas exige recompilar y `params_usados.par` no las registra: la reproducibilidad
depende del binario, no del archivo de parámetros.

### 14. Nombres sobrecargados
**Estado: PENDIENTE**

`sp` y `sr` son anchos en p y r en `gaussian1`, pero en Q y J en `aa` y en `df0`
(`distribution.f90:121`). `eps` en el código es el suavizado centrífugo, mientras que en
las notas ε es la amplitud de la perturbación. Ninguna colisión está advertida.

**Reevaluación:** la mitad de `eps` pierde fuerza con el punto 19: como no se puede poner
en un `.par`, nadie tropieza con la colisión de nombres en la práctica. Queda la de
`sp`/`sr`, que sí son parámetros de entrada y sí cambian de significado según el estado.

### 15. Las dos ramas del isócrono no eran idénticas bit a bit
**Estado: CORREGIDO** (2026-09-20)

`grav_force.f90:64` escribe la fuerza como `-r/(sq*(1+sq)**2)` y `grav_force.f90:80` como
`-r/sq*pot**2`: iguales en aritmética exacta, no en punto flotante. El comentario de las
líneas 52-53 afirma identidad bit a bit, pero esa afirmación cubre el refactor de arreglos
a bucles, no la comparación entre ramas.

**Medido.** Las dos expresiones difieren en el último bit en 222 113 de 400 002 radios
probados, con diferencia relativa de hasta 5.7e−16. En el código: dos corridas idénticas
salvo `autointeraction`, con a0 = 1e−300 para que la autogravedad no influya, dan
`hk1_complex.tl` **distintos**. (El archivo de partículas sale idéntico, pero eso es un
espejismo: `save2Ddata_particles` escribe solo ocho cifras.)

Vale la pena arreglarlo porque la propiedad perdida es útil: con masa despreciable, una
corrida con autogravedad debería reproducir exactamente la corrida sin ella, y eso es
justo el control que hemos usado varias veces en esta auditoría.

**Corrección.** Las dos ramas escriben ahora la fuerza igual,
`-r/(sq*(1+sq)**2)`. Prueba de aceptación: las mismas dos corridas de arriba dan
`hk1_complex.tl`, `hk2_complex.tl`, la energía y el archivo de partículas **idénticos bit
a bit**.

### 16. El fondo se suma a los puntos fantasma en unas ramas y no en otras
**Estado: PENDIENTE**

Líneas 70-71 y 102-103 operan sobre el arreglo completo; `add_background` solo sobre
`1:Nr`. Sin efecto porque la interpolación no lee los fantasmas
(`poisson_rk.f90:177-183`), pero sobreviviría a un cambio del interpolador.

### 17. Newton sin aviso de no convergencia
**Estado: PENDIENTE**

`utils.f90:725-731` itera 50 veces y sale si |g|<1e−14; si no converge devuelve el último
valor en silencio.

### 18. La lista de celdas acopla el factor de Courant con la corrección del depósito
**Estado: PENDIENTE**

`utils.f90:186-189` afirma que el depósito aplica el soporte exacto del peso para las
partículas archivadas en la celda del extremo. Vale para el término directo, no para el
término imagen: una partícula en r < −W_cell·dr tiene su imagen sobre puntos que no
escanean la celda 1 y su masa desaparecería. Con courant ≤ 0.5 nunca ocurre, pero nada
lo verifica.

### 19. `eps`: código muerto (el enunciado original de este punto era erróneo)
**Estado: ENUNCIADO CORREGIDO, severidad rebajada** (2026-09-20)

Lo que escribí primero — que `paramfile.f90` acepta `eps` del archivo de parámetros y
`set_grid_size` lo descarta en silencio — **es falso**. `eps` no está en la lista de
nombres leíbles de `paramfile.f90`: un `.par` que lo ponga es rechazado con
"Did you mean Npc?". No hay ningún valor del usuario que se descarte.

Lo que sí queda, y es mucho menor: `eps` vale cero siempre (su valor por omisión en
`parameters.f90`, reafirmado en `set_grid_size`), de modo que el suavizado del término
centrífugo `den = r**2 + eps*eps` de `grav_force.f90` y `energy.f90` es código muerto, y
la comprobación `eps /= 0` de `main.f90:54` no puede fallar nunca. La declaración en
`parameters.f90:42` ya lo documenta correctamente.

Opciones, ninguna urgente: dejarlo como está (una funcionalidad latente), quitarlo del
todo (sería neutral bit a bit, porque `r**2 + 0.0*0.0` es exactamente `r**2`), o
convertirlo en parámetro de verdad, lo que exigiría rehacer los mapas ángulo–acción con
el potencial suavizado.

### 20. El código y el artículo de referencia usan convenios de depósito distintos
**Estado: MEDIDO, decisión pendiente** (hallado al documentar el integrador)

`Vlasov_Poisson_evolutions/main.md` define el volumen de celda como el geométrico exacto,
ΔV_k = (4/3)π[(R_{k+1/2})³ − (R_{k−1/2})³] = 4πΔr(R_k² + Δr²/12), y la función de peso
como W_m = S_m * b₀, la forma de la partícula convolucionada con el promedio de celda. El
código usa V_i = 4πΔr(r_i² + (n+1)Δr²/12) y aplica `Wn` directamente. Como
b_m * b₀ = b_{m+1}, el `Wn(n)` del código **sí** es el W_m del artículo con m = n−1; lo que
no coincide es el volumen.

Son dos convenios coherentes cada uno consigo mismo y ambos de segundo orden:
ρ_k es el promedio exacto en la celda (artículo) o un estimador puntual insesgado de
ρ(R_k) (código). El primero conserva la masa exactamente y sesga el perfil de densidad;
el segundo hace lo contrario. `BUGS_TODO.md` registra el cambio a (n+1)/12 como una
corrección de error, cuando en realidad alejó el código del artículo.

Medido en la réplica, con Δr = 0.01:

*Masa total vista de una partícula aislada* — artículo exacto (1.0000000000) para todo n
y todo radio; código exacto solo para n=1 (para n=3 va de 0.751 a 0.99958).

*Esfera uniforme, error relativo de la fuerza dentro:*

| n | código | artículo |
|---|---|---|
| 1 | 1.28e−4 | 8.32e−3 |
| 2 | 4.9e−16 | 1.64e−2 |
| 3 | 9.6e−16 | 2.46e−2 |

*Plummer (masa concentrada), error relativo de la fuerza:*

| n | región | código | artículo |
|---|---|---|---|
| 1 | r < 0.2 | 3.50e−4 | 3.88e−2 |
| 2 | r < 0.2 | 3.49e−4 | 7.65e−2 |
| 3 | r < 0.2 | 4.51e−4 | 1.15e−1 |
| 1 | 0.2 < r < 2 | 2.22e−4 | 3.48e−4 |
| 1 | 2 < r < 3.9 | 9.94e−6 | 8.70e−6 |

El sesgo del convenio del artículo es n·Δr²/(12r²): crece hacia el centro y con el orden
del B-spline, hasta el 11 % dentro de r = 0.2 con n = 3. El del código lo deconvoluciona
y gana por factores de 100 a 250 donde está la masa; el del artículo gana por factores de
1.5 a 5 en la región exterior de baja densidad.

**Recomendación:** quedarse con el convenio del código y corregir en el artículo tanto
ΔV_k como el párrafo que describe el Runge–Kutta (§7.2 de las notas explica por qué ese
algoritmo perdía masa). Queda por comprobar si las corridas ya publicadas tenían masa
dentro de ~5Δr del origen; con L₀ = 2 la barrera centrífuga probablemente las protege,
pero no está verificado.

**Mejora posible:** para n ≥ 2, reconstruir ρ entre nodos con el B-spline del mismo orden
en vez de con una recta recuperaría también la conservación exacta de masa sin perder la
densidad insesgada. Para n = 1 no hace falta: la recta ya es el sombrero.

**No medido:** energía y h₁ en una corrida con autogravedad bajo el convenio del
artículo, que exigiría implementarlo en el código.

### 21. El equilibrio de `tools/equilibrio_L.py` no es un equilibrio del sistema discreto
**Estado: CONFIRMADO** (2026-09-20)

Al medir el ruido con un equilibrio autoconsistente y `eps = 0` — una configuración que
no debería evolucionar — la señal |h₁(t)−h₁(0)| a t=100 resultó ser:

| N | n=1 | n=2 |
|---|---|---|
| 8 000 | 2.118e−11 | 2.084e−11 |
| 16 000 | 1.857e−11 | 1.808e−11 |
| 32 000 | 1.740e−11 | 1.690e−11 |
| 64 000 | 1.687e−11 | 1.636e−11 |

Dos lecturas. La primera, buscada: **n=2 no compra nada**, un 3 %, lo que cierra la
decisión sobre el orden del B-spline y deja la reconstrucción de orden coincidente como
mejora disponible pero sin uso.

La segunda, no buscada y más importante: la señal **apenas baja con N** (factor 1.26 al
multiplicar N por 8, con la razón entre pasos consecutivos tendiendo a 0.97). El ruido de
discreción debería caer con N; esto satura. La explicación más plausible es que el
equilibrio que construye el generador no es equilibrio del sistema **discreto**: su
Φ_self, calculado con el mapa ángulo–acción en Python, difiere del que calcula la malla
del código, y esa diferencia produce una evolución sistemática independiente de N.

**Confirmado por tres vías.**

*Uno: el desajuste existe y se mide.* Comparando el Φ_self que el generador deja en
`ic_*_equilibrio.npz` con el que calcula la malla del código para esa misma condición
inicial (corrida con Nt=0, restando el isócrono del potencial de salida):

| | |
|---|---|
| profundidad del pozo Φ_self | 1.597e−03 |
| max\|Φ_código − Φ_generador\| | 1.439e−07 |
| relativo a la profundidad | **9.0e−05** |

*Dos: no es el integrador temporal.* Repitiendo la corrida con `courant` = 0.5, 0.25 y
0.125 (es decir dt, dt/2 y dt/4, con Nt escalado para llegar al mismo t=100):

| courant | \|h₁(fin)−h₁(0)\| |
|---|---|
| 0.5 | 2.118272622610e−11 |
| 0.25 | 2.118269538849e−11 |
| 0.125 | 2.118268579241e−11 |

Coinciden a **seis cifras** y convergen al refinar: la señal es una evolución real del
sistema discreto que el integrador resuelve bien, no un error suyo.

*Tres: el tamaño cuadra.* δΦ = 1.44e−07 desplaza la acción en δJ ≈ δΦ/ω ≈ 2.9e−07, o sea
δω/ω ≈ 3·δJ/(J+c) ≈ 3.5e−07. Sobre t=100 eso es un desfase δQ ≈ 1.4e−05 rad, y con
|h₀| = 1.06e−06 da |h₁| ≈ 1.5e−11, del mismo orden que los 2.1e−11 medidos.

**Consecuencias.** (a) Esa prueba **no mide el piso de discreción** y no se puede citar
como tal; hace falta otro diseño para medirlo. (b) Toda corrida que entre por
`state=checkpoint` con un equilibrio arrastra ese desajuste, que no mejora al subir N
porque es un efecto de la malla, no del muestreo.

**Arreglo propuesto.** Iterar el equilibrio contra el solver **del código** en vez del de
Python. No hace falta llamar al binario: la réplica del esquema (depósito con V_i, la
cuadratura de masa de §7.2 y la interpolación) reproduce al código a ocho cifras, así que
basta con sustituir la Poisson continua de `tools/equilibrio_L.py` por esa. Mientras no se
haga, conviene decir en las notas que los equilibrios son aproximados al nivel de 1e−4
relativo.

---

## A favor del código

- La advertencia de masa fuera de la malla (`density.f90:46-58`) es más severa de lo
  necesario: por el teorema de las capas una partícula exterior no ejerce fuerza sobre
  las interiores, así que la dinámica interior es correcta. El error real se limita a las
  fuerzas mutuas entre partículas que ya salieron.
- Como el depósito y la interpolación usan el mismo W, se cumple exactamente
  Σ_j m_j Φ(r_j) = Σ_i ρ_i Φ_i V_i, de modo que la energía potencial de `energy.f90:66`
  es exactamente la energía del campo discreto. El factor ½ está bien puesto.
- Los B-splines de `functions.f90` tienen los segundos momentos continuos correctos
  ((n+1)/12 para n=1,2,3), que es lo que justifica la fórmula de V_i.
- El manejo de espejo y exterior en la interpolación (`poisson_rk.f90:177-191`) es
  correcto y robusto para radios negativos y más allá de la malla.
