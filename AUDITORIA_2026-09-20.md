# Auditoría agnóstica del código — 2026-09-20

Revisión completa de `src/` buscando implementaciones no probadas, no justificadas o
que rompan la física que se quiere describir. Las afirmaciones marcadas como *medido*
se comprobaron replicando el esquema (depósito, RK2 de Poisson, interpolación) en un
programa independiente, no leyendo el código.

Estado de cada punto: PENDIENTE / EN CURSO / CORREGIDO (con el commit).

---

## SEVERIDAD ALTA

### 1. El esquema tiene autofuerza, y el código afirma lo contrario
**Estado: PENDIENTE**

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

### 2. Poisson pierde la masa de las partículas cercanas al origen
**Estado: PENDIENTE**

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

### 3. `BGtype="null"` sin autogravedad acumula fuerza sin límite
**Estado: PENDIENTE**

En `grav_force.f90:119-141`, con `BGtype="null"` y `autointeraction=.false.` no se
ejecuta `add_background` ni ninguna otra asignación, y el bucle siguiente hace
`pot_part(i) = pot_part(i) + centrífugo` sobre valores nunca reinicializados. Con
yoshida4 son 3 o 4 acumulaciones por paso. `null` es una opción aceptada
(`paramfile.f90:387`) y corresponde a la prueba exacta más elemental del integrador:
partículas libres con momento angular.

---

## SEVERIDAD MEDIA

### 4. Los fondos NFW y Burkert no tienen la paridad correcta en r
**Estado: PENDIENTE**

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

### 5. Raíz sin proteger en `analysish`: un NaN de órbita casi circular envenena los h_k
**Estado: PENDIENTE**

`analysish.f90:81-82` calcula `sqrt(-L²-2E-2-0.5/E)` sin protección; ese radicando vale
cero exactamente en una órbita circular y el redondeo puede volverlo negativo. La misma
expresión sí lleva `max(...,0)` en `utils.f90:733` y `utils.f90:781`. El filtro del bucle
solo descarta `energy(j) >= 0` (`analysish.f90:147`), así que una partícula ligada con NaN
contamina el acumulador completo. `initial_data.f90:122` tiene el mismo hueco, enmascarado
por un `f /= f` que atribuye el NaN a nodos no ligados cuando puede ser una órbita
circular válida a la que se le pone masa cero en silencio.

### 6. Arreglos automáticos de tamaño `Npart` en la pila
**Estado: PENDIENTE**

`analysish.f90:31-32` declara nueve arreglos automáticos de longitud `Npart`: 72 MB de
pila con 10⁶ partículas. Segfault dependiente del `ulimit -s`, justo al escalar N.

### 7. `courant = 80.0` en un archivo de parámetros distribuido
**Estado: PENDIENTE**

`reproducir/video/df_gauss.par` fija `courant = 80.0` para un parámetro documentado como
fracción de Courant (`parameters.f90:45`). Solo tiene sentido con `integrator = analytic`,
donde dt únicamente espacia las instantáneas. Nada lo verifica: cambiar el integrador en
ese archivo da dt = 4 sin aviso.

### 8. El criterio de paso de tiempo desconoce la frecuencia orbital
**Estado: PENDIENTE (documentado en el fuente)**

`utils.f90:251-255` lo documenta, incluido el caso nfw que divergió. La prueba de
convergencia queda a cargo del usuario y nada la exige.

### 9. El mapa ángulo–acción está triplicado
**Estado: PENDIENTE**

La misma transformación directa en `initial_data.f90:107-124`, `utils.f90:763-785` y
`analysish.f90:63-82`, con tres implementaciones que ya divergieron (la protección `max`).

---

## SEVERIDAD BAJA, Y DUDAS

### 10. Las pruebas de verificación son degeneradas respecto de lo que deben probar
**Estado: PENDIENTE**

La esfera uniforme es exactamente el caso que la fórmula
V_i = 4πdr(r_i² + (n+1)dr²/12) hace exacto por construcción, y también el caso en que el
arranque del RK2 (densidad constante dentro de r(1)) es exacto. Un esquema con los
defectos 1 y 2 pasa esa prueba sin marcas. Hacen falta casos que el diseño no privilegie:
cáscara única cerca del origen, órbitas radiales, convergencia en N a dr fijo.

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

### 15. Las dos ramas del isócrono no son idénticas bit a bit
**Estado: PENDIENTE**

`grav_force.f90:64` escribe la fuerza como `-r/(sq*(1+sq)**2)` y `grav_force.f90:80` como
`-r/sq*pot**2`: iguales en aritmética exacta, no en punto flotante. El comentario de las
líneas 52-53 afirma identidad bit a bit, pero esa afirmación cubre el refactor de arreglos
a bucles, no la comparación entre ramas.

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
