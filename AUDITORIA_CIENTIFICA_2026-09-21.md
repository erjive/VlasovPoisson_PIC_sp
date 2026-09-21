# Auditoría científica y numérica de VlasovPoisson_PIC_sp

**Fecha:** 2026-09-21. **Versión auditada:** rama `mejoras/portar`, commit `66bc5bb` (HEAD).
**Referencias:** versión original `84fc3b2` (2024-07-05), 96 commits intermedios, y el
manuscrito `Vlasov_Poisson_evolutions/main.md`.

**Regla de esta auditoría.** No se modificó `src/`. Toda afirmación de la forma "X es correcto"
o "X está mal" lleva al lado la prueba que la sostiene y su número. Donde no hubo prueba, se dice.
Las pruebas se corrieron una por una, nunca dos simulaciones a la vez. Las comparaciones de
tiempo se hicieron en orden ABBA. Las copias modificadas del código que se usaron como evidencia
(instrumentación, variante propuesta) se compilaron aparte y no tocan el repositorio.

**Qué se agregó al repositorio:** este documento y el directorio `verificacion/` (pruebas
automatizadas, sección E).

**Estado (actualizado el 2026-09-21, después de la auditoría).** D1 está **revertido**:
`poisson_rk.f90`, `arrays.f90` y `utils.f90` volvieron a su estado de `5c0e07f`. El
ejecutable nuevo reproduce a `5c0e07f` bit a bit en el colapso frío (82 datasets, 36
atributos), I2 pasa (orden 1.95 en 1/N) y U2 ya no depende del número de hilos (la lectura
de memoria sin inicializar estaba en el bloque revertido). La capa aislada de la antigua
prueba de aceptación converge con orden 1 en Δr a la solución de una cáscara delgada,
r̈ = L²/r³ − m/(2r²) (error 1.9·10⁻⁵ con Δr = 0.01). Las notas (`vlasov_L_intro.tex`,
"La autofuerza") y `AUDITORIA_2026-09-20.md` (punto 1) están corregidas. Siguen abiertos
D2 (solo la parte del depósito), D3, D4 y el resto de la tabla B. Las referencias de
línea y a `mcoefA/B` de este documento son de `66bc5bb`: la reversión devolvió el avance de
la masa a la forma c₀ + c₃r³ + c₄r⁴ de `5c0e07f`, igual de exacta (U3: 6.6·10⁻¹⁶ con 41
nodos, 1.1·10⁻¹⁴ con 2001), y con ella D12 cambia de expresión pero no de tamaño.

---

## A. Resumen ejecutivo

### Lo que quedó demostrado

| Afirmación | Prueba | Resultado |
|---|---|---|
| Las funciones de forma W₁, W₂, W₃ son las B-splines correctas: partición de la unidad, área 1, segundo momento (n+1)/12, simetría exacta. | U1 | error ≤ 3·10⁻¹⁶ |
| El depósito conserva la masa y reproduce una densidad uniforme exactamente (r > −Δr). | U2, T3c | \|ΔM\|/m ≤ 2·10⁻¹⁶ |
| El integrador radial de Poisson es **exacto** para su interpolante lineal por tramos. | U3 | 4·10⁻¹⁶ (41 nodos), 10⁻¹⁴ (2001 nodos) |
| El campo en la malla converge con **orden 2.00** en Δr para n = 1, 2, 3. | U4 | razón 4.0 exacta en 3 niveles |
| Los siete fondos cumplen F = −dΦ/dr y la paridad Φ(−r) = Φ(r), F(−r) = −F(r) **bit a bit**. | U5 | 10⁻¹² a 10⁻¹³; paridad 0 |
| Las fórmulas del isócrono (E(J,L), Ω, J, Q) coinciden con cuadratura independiente. | T6b | \|ΔJ\| ≤ 10⁻³⁰, \|ΔQ\| ≤ 3·10⁻¹⁶ |
| Euler, leapfrog y yoshida4 convergen con **orden 1, 2 y 4** frente a la solución exacta. | I1 | 1.02, 2.00, 4.00 |
| Sin autogravedad, el código actual reproduce al original de 2024 **en las 8 cifras impresas** durante 2000 pasos. | R1 | 8192/8192 partículas idénticas |
| Siete commits de optimización o limpieza son **neutros bit a bit**. | R2 | 46 datasets + 20 atributos idénticos por par |

### Lo que quedó demostrado como defecto

| Id | Severidad | Hallazgo | Evidencia |
|---|---|---|---|
| **D1** | **HIGH** | La "corrección" del punto 1 de la auditoría anterior (commit `bfd9b29`) quita una fuerza **física**: la autogravedad −m/(2r²) de cada cáscara. Con ella, el error del colapso frío frente al continuo es O(1/N); sin ella, O(1/N²). En N = 800 el código actual es **2700 veces peor** que la versión anterior (`5c0e07f`) y **30 veces peor que el original de 2024**. Además duplica el costo (×2.09, ABBA). | I2, tabla §8.2 |
| **D2** | **HIGH** | Partícula en r ≤ −Δr durante un paso: la resta de la autofuerza **lee memoria no inicializada** (resultado depende del número de hilos) y el depósito **pierde masa** (hasta 100 %). Alcanzable sin aviso con opciones admitidas (`BGtype=Central`, L pequeño): ≥ 3363 evaluaciones afectadas en 500 pasos. | U2, §8.4 |
| **D3** | **MEDIUM** | La iteración de Newton de la ecuación de Kepler en `invert_QJ_to_rp` **diverge** (residuo 10²⁵) para e ≳ 0.98 en 1–3 % de los ángulos, y entrega un r plausible pero falso, sin aviso. Afecta `aa_quad` y el integrador `analytic`. | U6, T6c |
| **D4** | LOW | `invert_QJ_to_rp` e `init_action_angle` no acotan los radicandos de los puntos de retorno: con L = 0 el redondeo da NaN y r = 0 para todo Q (650/1600 casos). Sin masa (f·L = 0), pero trayectorias falsas. | U6 |

(La lista completa, con 17 hallazgos, está en la sección B.)

### Veredicto por categoría

- **Matemáticamente consistente:** sí en cada pieza verificada por separado (formas, depósito,
  Poisson, fondos, mapa ángulo-acción, integradores), **salvo** D2 (simetría en el origen rota
  para r ≤ −Δr), D3 y D4 (inversión del mapa fuera de su dominio seguro).
- **Numéricamente consistente:** sí. Órdenes medidos iguales a los esperados en Δr y Δt. El error
  de energía con autogravedad no depende de Δt (demostrado) y baja con N (medido); es una
  propiedad del esquema PIC, no un defecto.
- **Físicamente consistente:** **no del todo, desde `bfd9b29` (2026-09-20)**. D1 cambia la física
  de cada elemento de masa. Sin autogravedad (el régimen de las corridas del artículo) la física
  está intacta y demostrada contra la solución exacta.
- **Reproducible:** sí para un número de hilos dado (sumas deterministas), **excepto** cuando D2 se
  activa (lectura de memoria no inicializada). Entre números de hilos distintos, las energías y h_k
  difieren en el último bit (documentado en el código).

### Lo que no se pudo verificar (sección F)

Evolución de largo plazo con autogravedad frente a una teoría (no hay solución analítica en el
régimen del artículo); h_k frente al semianalítico en esta pasada (lo hicieron sesiones previas);
`reduceparticles`, `dt_switch=var`, estado `checkpoint` más allá de su uso en las pruebas;
herramientas de `tools/`.

### Recomendación principal

Revertir `bfd9b29` (volver a `5c0e07f` en `poisson_rk.f90`), corregir D2 evaluando la fuerza en
|r| y D3 con un arranque robusto de Newton, y **después** volver a correr cualquier resultado con
autogravedad producido desde el 2026-09-20. Las corridas del artículo (sin autogravedad) no
dependen de ninguno de estos cuatro defectos.

---

## 1. Modelo matemático reconstruido

Reconstruido del código (no de los nombres) y contrastado con `main.md`.

**Ecuaciones.** Vlasov–Poisson con simetría esférica y distribución en el momento angular. Con
F(t, r, p_r, L) la función de distribución reducida:

    ∂F/∂t + p_r ∂F/∂r + (−∂Φ_tot/∂r) ∂F/∂p_r = 0,
    Φ_tot = Φ_ext(r) + L²/(2r²) + Φ_gas(r),
    (1/r²) d/dr (r² dΦ_gas/dr) = 4π ρ_gas,
    ρ_gas(r) = (2π/r²) ∫∫ F L dp_r dL,     d³x d³v = 8π² L dr dp_r dL.

L es constante de movimiento y rotula a cada partícula; no hay ecuación para L.

**Unidades.** G = M_iso = b = 1 (`main.md` §2.2). La masa del gas es `a0` en unidades de M_iso;
el cociente m/M_iso se absorbe en F.

**Signos.** Fuerza radial = −∂Φ/∂r (atractiva negativa). En el código `force = −dev_pot` y
`force_part` es la fuerza total por unidad de masa, centrífugo incluido (+L²/r³).

**Variables.** `r_part`, `p_part`, `l_part`: (r, p_r, L) de cada partícula. `f`: valor de F en el
nodo. Masa de la partícula m_j = 8π² Δr_c Δp_c ΔL_c f_j L_j (regla del punto medio de la medida).

**Densidad, potencial, campo.** `avg_rho(i)` = ρ en el nodo r_i. `pot`, `force` = Φ_gas y −dΦ_gas/dr
en la malla (más el fondo, sumado después). `pot_part`, `force_part` = totales en la partícula.

**Condiciones iniciales.** Nodos en los puntos medios de una malla regular en (r, p, L)
(`gaussian1`, `aa`) o en (J, Q, L) (`aa_quad`), o leídos de archivo (`checkpoint`). Se descartan
los nodos con f ≤ cutoff·f_max y se renormaliza la masa a `a0`.

**Condiciones de frontera.** En el origen, la simetría F(r, p) = F(−r, −p): puntos fantasma con
ρ par y dΦ/dr impar, imágenes de cada partícula en el depósito, y reflexión (r, p, F) → (−r, −p, −F)
al final de cada paso (`main.f90:239-249`). En el borde exterior, Φ + r dΦ/dr = 0 (solución
exterior −M/r) y, para partículas más allá de la malla, la solución exterior de la masa en la
malla. No hay periodicidad.

**Poisson.** Integración hacia afuera de dM/dr = 4πr²ρ y dΦ/dr = M/r² con ρ lineal entre nodos
(constante en [0, r₁]) y ambas integrales en forma cerrada (`poisson_rk.f90:108-136`, pesos
`mcoefA`, `mcoefB` de `utils.f90:193-205`). Ya no es el RK2 que describe `main.md`.

**Avance temporal.** `euler` (orden 1), `leapfrog` KDK (orden 2), `yoshida4` (composición de
Yoshida, orden 4), `analytic` (solución exacta en ángulo-acción, solo isócrono sin autogravedad).
Paso fijo por omisión.

**Partícula → malla.** ρ_i = Σ_j m_j [W_n((r_i−r_j)/Δr) + W_n((r_i+r_j)/Δr)] / V_i con
V_i = 4πΔr (r_i² + (n+1)Δr²/12), el volumen que cubre W_n (`density.f90:136-160`).

**Malla → partícula.** Φ(r_j) = Σ_i W_n((r_j−r_i)/Δr) Φ_i con los puntos r_i ≤ 0 tomados del
espejo (Φ par, F impar) y los r_i > r_Nr de la solución exterior (`poisson_rk.f90:201-234`).
Desde `bfd9b29` se resta además la contribución de la propia partícula (§10, D1).

**Magnitudes conservadas esperadas.** Masa total (exacta salvo recortes explícitos); L de cada
partícula (exacta, es un rótulo); peso f de cada partícula (Liouville, exacto por construcción);
energía total (exacta en el continuo; en el discreto, acotada con integradores simplécticos sin
autogravedad, y no exacta con autogravedad, ver §6). No hay momento lineal que conservar
(simetría esférica).

## 2. Mapa ecuación → código

| # | Ecuación | Discretización | Implementación | Supuestos | Discrepancias |
|---|---|---|---|---|---|
| E1 | ṙ = p_r | Drift de cada integrador | `main.f90:170,180,195-219` | Paso fijo | Ninguna. `main.md` (Leapfrog2) escribe Δt/2 en el drift: errata del manuscrito, el código usa Δt. |
| E2 | ṗ_r = −∂Φ_tot/∂r | Kicks | `main.f90:171,179,184,202-218` | F evaluada en posiciones nuevas | `main.md` (Leapfrog3) escribe F(t_n, r_{n+1}); el código usa t_{n+1}, que es lo correcto. |
| E3 | Φ_ext isócrono = −1/(1+√(1+r²)) | Analítica | `grav_force.f90:64-95` | Par en r | Verificado U5. |
| E4 | Otros fondos | Analítica, en \|r\| | `grav_force.f90:189-249` | Par/impar | Verificado U5. |
| E5 | Φ_L = L²/(2r²), fuerza L²/r³ | Analítica | `grav_force.f90:74-75,149-155` | eps = 0 forzado | Verificado U5. |
| E6 | ρ = (2π/r²)∫∫F L dp dL | Suma de partículas × W_n / V_i, con imágenes | `density.f90:95-175` | V_i con (n+1)/12 | `main.md` usa W_m = S_m∗b₀ y ΔV_k geométrico (punto 20 de la auditoría anterior, sin decidir). Pierde masa si r_j ≤ −1.5Δr (D2). |
| E7 | dM/dr = 4πr²ρ | ρ lineal entre nodos, integral cerrada | `poisson_rk.f90:117-131`, `utils.f90:193-205` | Densidad par en el origen | Exacto para el interpolante (U3). `main.md` describe RK2 sobre (Φ, Φ′). |
| E8 | dΦ/dr = M/r² | Integral cerrada de M/r² | `poisson_rk.f90:122-136` | — | Exacto (U3). |
| E9 | Φ + r Φ′ = 0 en r_max | Desplazamiento de Φ | `poisson_rk.f90:162` | Toda la masa dentro de la malla | Aviso una vez si sale masa (`density.f90:46-58`). |
| E10 | Φ(r_j), F(r_j) | Interpolación W_n con espejo y exterior | `poisson_rk.f90:201-234` | Mismo W_n que el depósito | Verificado U4 (orden 2). |
| E11 | (sin ecuación) resta de la autointeracción | Densidad, masa y campo propios por la malla | `poisson_rk.f90:252-338` | "La autofuerza es un artefacto" | **Falso en simetría esférica (D1).** Ventana mal calculada para r ≤ −Δr (D2). |
| E12 | F(r,p) = F(−r,−p) | Reflexión tras el paso | `main.f90:239-249` | Fuerza impar exacta | Exacta para fondos y centrífugo (U5); rota para la autofuerza con r ≤ −Δr (D2). |
| E13 | E = ∫(p²/2 + L²/2r² + Φ_ext + Φ_gas/2) F d³xd³v | Suma de partículas | `energy.f90:60-67` | Factor 1/2 solo en Φ_gas | Correcto. El reparto K/W se redefinió en `abb8bc9` (centrífugo en K). |
| E14 | h_k = 8π²∫F Φ̂_k(J)* e^{−ikQ} dr dp L dL | a_k por Simpson (512 intervalos), suma de partículas | `analysish.f90` | Mapa del isócrono aun con autogravedad | Coincide con `main.md` (Eq. h00k). El comentario de cabecera (línea 4) dice Φ en vez de Φ̂_k. |
| E15 | Mapa (r,p,L) ↔ (Q,J,L) del isócrono | Fórmulas cerradas + Newton de Kepler | `utils.f90:728-825`, `analysish.f90`, `initial_data.f90` | e < 1 | Tres copias (punto 9 anterior): solo dos tienen los radicandos acotados (D4). Newton sin salvaguarda (D3). |

## 3. Auditoría de las transformaciones y optimizaciones

**Método.** Para cada commit que se declaró optimización, limpieza o refactor: (i) diff del código
sin comentarios (script que conserva las directivas `!$OMP`), (ii) prueba algebraica de la
equivalencia, (iii) compilación del padre y del hijo con las mismas banderas y corrida de dos casos
con autogravedad (leapfrog y yoshida4, 5400 partículas, 20 pasos, salida HDF5 cada 5 pasos),
(iv) comparación bit a bit de los 46 datasets y 20 atributos (script `h5igual.py`). Las banderas
del Makefile son `-O3 -funroll-loops -fopenmp` sin `-ffast-math`: el compilador no reasocia y
x86-64 base no contrae a FMA, así que misma secuencia de operaciones ⇒ mismos bits.

| Commit | Qué hace | A. ¿Equivalente? | B. ¿Exacta? | C. ¿Siempre? | D. ¿Cambia bits? | E. ¿Estabilidad/conservación? | F. ¿Cambio semántico oculto? | Evidencia |
|---|---|---|---|---|---|---|---|---|
| `f762437` | Quita `collapse(2)` | Sí (elimina una carrera) | — | Sí | Sí, quita no determinismo | No | No; corrige | Diff |
| `505a32d` | Lista de celdas en el depósito | Sí | Sí en ℝ | Sí para r ≥ 0 | **Sí**: cambia el orden de suma | No | No en su momento. Con las imágenes de `a81f4cf` su ventana se volvió insuficiente para r ≤ −1.5Δr (D2) | Mensaje del commit: 10⁻⁶ en posiciones tras 10⁵ pasos caóticos |
| `9860975` | Datos iniciales en paralelo | Sí (índice directo) | Sí | Sí | No | No | No; ciclo por elemento, sin compactación | Diff |
| `2f345bc` | "Reusar el peso de cuadratura" | **No** | — | — | Sí | — | **Sí**: corrigió de paso φ̂_k, que usaba siempre k = 4; se documentó después en `e5847a0` | Diff |
| `6200fe2` | a_k una sola vez | Sí: Σ w_k A_k B C cos = B C Σ w_k A_k cos | Sí en ℝ | Sí | Sí, último bit de h_k | No | No | h_k idéntico en 8 cifras; datos de partículas idénticos; las 12 diferencias en energía son ruido de la reducción OpenMP del padre (dos corridas del padre ya difieren) |
| `e6d2011` | Sumas por hilo en orden fijo | Sí | Sí en ℝ | Sí | Sí, ≤ 6.4·10⁻¹⁶ relativo | No | No; hace el resultado determinista para un número de hilos | R2 |
| `5219f6f` | Omitir la fuerza final de yoshida4 si no hay diagnóstico | Sí: esa fuerza nunca se lee en pasos sin diagnóstico | Sí | Sí (incluye `dt_switch=var`) | No | No | No | **Idéntico bit a bit** |
| `4f5c86b` | Fondo + centrífugo en una pasada; √(1+r²) una vez | Sí, misma secuencia por elemento | Sí | Sí (`forcetype` solo admitía `bg`) | No | No | No | **Idéntico bit a bit** |
| `568bee7` | Ventanas del núcleo a (n+2)/2 | Sí: las celdas quitadas solo sumaban ceros exactos | Sí | Sí para r ≥ 0 en su momento | No | No | No | **Idéntico bit a bit** |
| `75c5779` | Kick+drift de yoshida4 fundidos | Sí, misma secuencia por elemento | Sí | Sí | No | No | No | **Idéntico bit a bit** |
| `23a30ac` | "Solo comentarios" | Sí | Sí | Sí | No | No | Quitó variables sin uso y una rama `rk4` ya rechazada por `paramfile` | **Idéntico bit a bit** |
| `ea928cc` | Código muerto | Sí | Sí (`m0` = 1.0 no era legible; ×1.0 es exacto) | Sí | No | No | No | **Idéntico bit a bit** |
| `a646629` | Mapa AA con escalares | Sí, mismas expresiones | Sí | Sí | No | No | No | **Idéntico bit a bit** |
| `86b3c80` | Una sola expresión de la fuerza del isócrono | −r/√(1+r²)·Φ² = −r/(s(1+s)²) | Sí en ℝ | Sí | Sí, último bit (222113 de 400002 radios, medido antes) | No | No | Auditoría anterior |
| `5c0e07f` | Integrar M en vez de (Φ, Φ′) | **No es optimización**: cambia el método | — | — | Sí | Mejora (§8.2) | Intencional y documentado | U3, I2 |
| `bfd9b29` | Restar autofuerza y autoenergía | **No es optimización**: cambia la física | — | — | Sí | **Empeora** (D1) | Intencional, pero la justificación es falsa | I2 |

En el código actual quedan además estas reescrituras algebraicas respecto del original, todas
revisadas (tabla C): el volumen del depósito, los pesos `mcoefA/B`, los coeficientes c₀, c₃, c₄, la
potencia compleja `expv**i`, y la separación K/W de la energía.

## 4. Discretización

**Órdenes esperados y medidos.**

| Componente | Esperado | Medido | Prueba |
|---|---|---|---|
| Depósito de un perfil suave | O(Δr²) (sesgo de suavizado ∝ (n+1)Δr²/12) | incluido en U4 | U4 |
| Integrador radial | Exacto para el interpolante; O(Δr²) frente al continuo | 4·10⁻¹⁶; orden 2.00 | U3, U4 |
| Interpolación a partículas | O(Δr²) | orden 2 en la malla; en partículas el error dominante es O(1/N) por D1 | T3 |
| Borde discontinuo (esfera uniforme) | O(Δr) local | orden 1 en el borde (2.3·10⁻² → 2.4·10⁻³ para Δr 0.04 → 0.005) | I2 |
| Tiempo | 1, 2, 4 | 1.02, 2.00, 4.00 | I1 |
| Cuadratura de partículas (punto medio) | O(h²) para F suave | masa del Plummer 1.0089 → 1.0002 al refinar ×2 | §6 |

**Lista de errores clásicos de PIC revisados.**

| Error clásico | Resultado |
|---|---|
| Off-by-one en la malla escalonada | No. r_i = r_min + (i−½)Δr; fantasmas r_{1−k} = −r_k verificados (U2 para \|x\| < Δr). |
| Índices desplazados | No en el código actual. El original escribía la densidad corrida 2 nodos (§9.1), corregido en `41facf0`. |
| Factores Δx, Δt | Correctos. Masa = 8π²ΔrcΔpcΔLc f L en depósito, energía, h_k y `rhomix`. |
| Normalizaciones | Consistentes (U2, U4). |
| Doble conteo | No. Imagen y directa suman W(r_i−r_j) + W(r_i+r_j), que no se solapan para r_i, r_j > 0 más allá del soporte. |
| Pesos | Correctos (U1). |
| Inconsistencia depósito/interpolación | Mismo W_n y misma ventana para r ≥ 0. **Para r ≤ −Δr la ventana de la autofuerza no cubre el soporte de la imagen (D2).** |
| Signo de la fuerza | Correcto (U5, I1, I2). |
| Fronteras | Origen: correcto salvo D2. Exterior: Φ = −M/r verificado (U4). |
| Paso de tiempo | Criterio de Courant con `pmax` del usuario y criterio de aceleración con la fuerza de t = 0; no se revisa después (D10). |
| Tamaño de malla | `Nr = int((rmax−rmin)/dr)+1` trunca por redondeo (D11). |

## 5. Pruebas con solución conocida

Todas con errores cuantificados; la tabla completa está en D.

- **Campo exactamente cero.** Una partícula aislada: con la versión actual, fuerza propia 0 a
  redondeo para r > 0 (esto es lo que D1 impone, y es físicamente incorrecto para una cáscara).
  Con r ≤ −Δr, fuerza espuria de hasta 0.157 para m = 0.01 (D2).
- **Densidad uniforme.** Interior de la esfera uniforme: el depósito es exacto y F = −r también;
  el error restante en partículas es 2πρ₀Δr/K (exclusión de la propia masa, D1), que coincide con
  lo medido (9.1·10⁻⁴ con Δr = 0.04, K = 64).
- **Poisson analítico.** ρ = ρ₀(1−r²)³: orden 2.00 en Φ y F para n = 1, 2, 3 (U4).
- **Interpolante lineal.** Exactitud a redondeo (U3).
- **Simetría conocida.** Paridad de fondos exacta (U5); simetría del depósito exacta para
  \|r\| < Δr y rota más allá (U2).
- **Equilibrio conocido.** Plummer isótropo autogravitante (§6).
- **Solución exacta no lineal.** Colapso frío de una esfera uniforme antes del cruce de cáscaras
  (§8.2).
- **Perturbación sinusoidal / solución linealizada.** No aplicable directamente: el código es
  radial y no tiene un problema lineal con solución cerrada fuera del mezclado de fase, que
  sesiones anteriores ya contrastaron con `tools/hk_exacto.py`. No se repitió aquí.

## 6. Propiedades físicas y conservación

| Magnitud | ¿Debe conservarse? | Resultado | Tipo de error |
|---|---|---|---|
| Masa total | Sí, salvo recortes explícitos | Depósito exacto a 2·10⁻¹⁶ para r > −Δr; pérdida de hasta 100 % de una partícula en r ≤ −1.5Δr (D2); la masa más allá de r_Nr no entra en Poisson (aviso único) | algoritmo (D2) |
| L de cada partícula | Sí, exacta | Nunca se modifica; L mín/máx idénticos al inicio y al final en todas las corridas | — |
| Peso f | Sí, exacta | Constante por construcción | — |
| Energía, fondo fijo | Acotada con simplécticos | leapfrog 1.3·10⁻⁴ → 1.3·10⁻⁷ (Δt 0.2 → 0.00625), yoshida4 1.3·10⁻⁵ → 3.2·10⁻¹² | discretización temporal |
| Energía, autogravedad | No exacta en el semidiscreto | Plummer N = 19356: 4.6·10⁻⁴ **igual para leapfrog y yoshida4 con Δt de 0.0167 a 0.0021**; N = 154270: 8.4·10⁻⁶ (t = 50). Con Δr sin tendencia limpia (1.3·10⁻⁴ a 4.6·10⁻⁴) | espacial/muestreo, ver abajo |
| Virial 2K/\|W\| | ≈ 1 en equilibrio | Plummer N = 154270: 0.992 → 0.997; K(0) = 0.14757 frente a 0.14726 analítico | muestreo inicial |
| Momento lineal | No aplica (simetría esférica) | — | — |

**Por qué la energía no se conserva con autogravedad.** Lo demostrado es que el error no baja con
Δt, así que el sistema semidiscreto (tiempo continuo) no conserva la energía que calcula
`energy.f90`. La explicación estándar es que la fuerza interpolada Σ W F_i no es el gradiente de la
energía discreta ½ Σ m Φ(r_j) (esquema "que conserva momento", no "que conserva energía"). No
separé esa contribución de la del ruido de muestreo; para hacerlo haría falta implementar la
variante variacional (fuerza = −d/dr_j de Σ W Φ_i) y comparar.

**Distinción de errores en este código.** Físico: D1 (se quitó una fuerza real). Discretización:
O(Δr²) en la malla, O(Δt^k) en el tiempo, O(Δr) en bordes discontinuos. Estadístico/muestreo:
O(1/N) con D1, O(1/N²) sin él (cuadratura en malla regular, no Monte Carlo). Redondeo: 10⁻¹⁴ en
el integrador radial con 2000 nodos (D12), último bit por orden de suma. Algoritmo: D2, D3, D4.

## 7. Convergencia

| Variable | Experimento | Resultado |
|---|---|---|
| Δr | U4 (Poisson analítico), 3 niveles | orden 2.00 (n = 1, 2, 3) |
| Δr | Colapso frío, Δr 0.04 → 0.005, dt fijo | interior independiente de Δr (densidad uniforme: exacta); borde orden 1 |
| Δt | I1, 6 niveles | 1.02 / 2.00 / 4.00 |
| Δt | Plummer autogravitante | error de energía independiente de Δt (piso espacial) |
| N | Colapso frío, N = 50 … 800 | **HEAD: orden 1.00 en 1/N; `5c0e07f`: orden ≈ 2** |
| Tamaño del dominio | — | **No se varió.** La condición exterior se verificó (U4) pero no la sensibilidad a r_max con masa saliendo. |
| Resolución de la malla de partículas | Plummer 19356 → 154270 | masa de cuadratura 1.0089 → 1.0002; energía 5·10⁻⁴ → 8·10⁻⁶ (dos N, sin orden) |

## 8. Física de las simulaciones

### 8.1 Sin autogravedad (régimen del artículo)

Contra la solución exacta del isócrono (integrador `analytic`), 60 órbitas con L ∈ {0.3, 1, 2}:
leapfrog y yoshida4 convergen con su orden, sin deriva de energía (I1). No hay crecimiento,
amortiguamiento ni calentamiento artificiales atribuibles al código en este régimen. Las
corridas `exe/input_article_*` (isócrono, sin autogravedad, L ≈ 2, leapfrog) no alcanzan D1, D2,
D3 ni D4.

### 8.2 Autogravedad: colapso frío y qué commit cambió la física

Esfera uniforme R = 1, M = 1, en reposo, L = 10⁻⁴, n = 1, Δr = 0.01, yoshida4 con dt = 0.0025,
t = 0.8 (r/r₀ = 0.634). Antes del cruce de cáscaras cada cáscara sigue la cicloide del continuo
r = r₀cos²θ, t = (θ + sen θ cos θ)/√2. Error relativo mediano en r para 0.1 < r₀ < 0.9:

| Versión | Qué cambia | N = 200 | N = 800 | Orden en 1/N |
|---|---|---|---|---|
| `894a120`, `ea928cc` | RK2 y acoplamiento original | 9.6·10⁻⁵ | 8.0·10⁻⁵ | ~0 (piso de la malla) |
| `a81f4cf` | acoplamiento en el origen corregido | 1.9·10⁻⁵ | 2.5·10⁻⁶ | 1.5 |
| `5c0e07f` | integración de la masa | 1.7·10⁻⁵ | **9.1·10⁻⁷** | 2.1 |
| `66bc5bb` (HEAD) | **resta de la autofuerza** | 9.9·10⁻³ | **2.5·10⁻³** | 1.00 |
| variante de prueba (dt = 0.005) | resta por malla + −m/(2r²) analítico | 1.9·10⁻⁴ | 4.3·10⁻⁵ | 1.1 |

Con dt fijo y Δr de 0.04 a 0.005 (N = 800), el error de HEAD no cambia (2.4·10⁻³ a 2.5·10⁻³) y el
de `5c0e07f` tampoco (0.9·10⁻⁶ a 1.1·10⁻⁶): **el error de HEAD es un sesgo O(1/N), no un error de
malla.** La conservación de energía es parecida en todas (2·10⁻⁵ a 8·10⁻⁵).

**Por qué.** Una cáscara delgada de masa m siente el promedio del campo a ambos lados de sí misma,
−m/(2r²): si la cáscara tiene espesor δ y masa interior m(s) creciente de 0 a m,
⟨F⟩ = −(1/m)∫ m(s)/s² dm(s) → −m/(2r²). Es física newtoniana, no un artefacto: una cáscara de polvo
en reposo colapsa por su propio peso. Cada partícula de este código es una cáscara (un elemento de
cuadratura de una F suave), y la mitad de su masa está "dentro" de su radio. El campo de la malla
(depósito + Poisson + interpolación) evalúa el campo continuo de la densidad reconstruida,
incluida la propia; esa es la estimación consistente. En la malla cartesiana la autofuerza se
anula sola porque el campo de una nube simétrica es cero en su centro; en simetría esférica el
mismo principio da −m/(2r²), no cero. La variante que resta la parte de malla y suma el término
analítico mejora a HEAD pero sigue en orden 1: lo que importa es tratar la propia masa igual que
a las vecinas que el núcleo cuenta parcialmente.

La prueba de aceptación que motivó `bfd9b29` (una partícula aislada debe moverse libre,
r(t) = √(1+t²)) supone que una cáscara con masa no se atrae a sí misma, que es falso.

### 8.3 Autogravedad: equilibrio de Plummer

Plummer isótropo f(E) = C(−E)^{7/2}, cuadratura regular en (r, p, L) hasta r = 20, `BGtype=null`,
leapfrog, t = 50 (unos 20 tiempos dinámicos en el radio de media masa). Con N = 154270 ambas
versiones mantienen el equilibrio dentro de la resolución del muestreo inicial: radios del 10, 50
y 90 % de la masa 0.55/1.25/3.65 → 0.51/1.29/3.61 (los valores iniciales están cuantizados en los
nodos, Δr_c = 0.1), energía 8.4·10⁻⁶ (HEAD) y 1.5·10⁻⁵ (`5c0e07f`). A este N el efecto de D1 (~1/N)
no se distingue. Con N = 19356 el estado inicial no está en equilibrio a mejor que 2.5 % (virial
0.976) y no sirve para comparar versiones.

### 8.4 ¿Se alcanza D2 en la práctica?

Copia instrumentada (cuenta partículas con r ≤ −Δr en cada evaluación de fuerza), 6400 nodos
`gaussian1` con L ∈ [0, 0.2], Δr = 0.05, leapfrog, autogravedad `a0 = 0.01`:

- `BGtype=Isochrone`: 0 casos (la velocidad en el isócrono no pasa de 1 < `pmax`).
- `BGtype=Central`: **≥ 3363 evaluaciones partícula-fuerza en r ≤ −Δr en las primeras 501
  evaluaciones**, sin ningún aviso (el aviso de `|p| > pmax` solo se revisa en t = 0).

## 9. Regresión

### 9.1 Original (`84fc3b2`, 2024) frente al actual, sin autogravedad

Mismo estado inicial construido a mano para que coincida bit a bit (los nodos del original están
en el borde derecho de la celda y los actuales en el punto medio; se desplazó el soporte medio
paso con valores diádicos), isócrono, leapfrog, 8192 partículas, 2000 pasos.

| Salida | Resultado | Clasificación |
|---|---|---|
| r, p, f·L de cada partícula, t = 0 … 50 | **idénticas en las 8 cifras impresas** | — |
| Energía total | idéntica | — |
| Cinética y potencial por separado | difieren 7.7·10⁻³ | cambio legítimo: el centrífugo pasó de W a K (`abb8bc9`) |
| Densidad en la malla | difiere hasta 40 %; alineada 2 nodos, 7·10⁻⁷ | **error del original**: salida corrida 2 nodos por los fantasmas (corregido en `41facf0`); el resto es el volumen (n+1)/12 |
| h_k | NaN en el original | **error del original**: partículas no ligadas (corregido en `732250f`); además φ̂_k usaba siempre k = 4 (corregido en `2f345bc`) |

La reflexión en el origen del original no invertía la fuerza, y leapfrog usaba en el medio paso
siguiente una fuerza de signo contrario para las partículas que cruzan r = 0 (corregido en
`a81f4cf`). No aparece en esta prueba porque con L ≥ 0.5 nadie cruza.

### 9.2 Con autogravedad

Tabla de §8.2. Las diferencias se deben a: (1) error del original (acoplamiento en el origen),
(2) mejora numérica (integración de la masa), (3) **cambio involuntario de la física** (`bfd9b29`).

### 9.3 Commits de optimización

Tabla de §3: siete idénticos bit a bit (`5219f6f`, `4f5c86b`, `568bee7`, `75c5779`, `a646629`,
`23a30ac`, `ea928cc`), dos al nivel del redondeo por orden de suma (`6200fe2`, `e6d2011`), uno que
cambió resultados por corregir un error (`2f345bc`), uno con orden de suma distinto (`505a32d`) y
uno neutro por construcción no corrido (`9860975`).

## 10. Problemas, en detalle

### D1 · HIGH · Se quitó la autogravedad física de cada cáscara

**Estado: revertido el 2026-09-21.**


- **Archivo/función:** `src/poisson_rk.f90`, `poisson_rk`, líneas 236-338 (bloque completo), en
  particular 337-338.
- **Expresión matemática:** fuerza sobre la cáscara j = −[M_int(r_j) + m_j/2]/r_j² (+ fondo); en
  el PIC, la interpolación del campo de la densidad reconstruida completa.
- **Expresión implementada:** `force_part(i) = force_part(i) - fself` con `fself` el campo de la
  propia partícula por la malla; ídem `pot_part - pself`.
- **Problema:** quita −m/(2r²)(1 − (5/6)Δr/r + …), que es la autogravedad de una cáscara, no un
  artefacto. Pasa de O(1/N²) a O(1/N) frente al continuo; duplica el costo de la fuerza.
- **Evidencia:** §8.2 (N = 50 … 800, Δr 0.04 … 0.005); ABBA: 7.66 s frente a 3.67 s (medianas de 4).
  La auditoría anterior (`AUDITORIA_2026-09-20.md`, punto 1) y las notas (`vlasov_L_intro.tex`
  §7.4) presentan esto como corrección: ambos textos quedan invalidados en ese punto.
- **Impacto:** grande en régimen dominado por autogravedad con N moderado (2.5·10⁻³ en r con
  N = 800); pequeño en el régimen del artículo si se activa la autogravedad (fracción ~1/N de una
  fuerza que ya es `a0` veces el fondo); nulo sin autogravedad.
- **Corrección:** revertir `bfd9b29` en `poisson_rk.f90` (y `mcoefA/B` pueden quedarse, los usa
  el integrador). Sustituir la prueba de la partícula libre por el colapso frío (I2).

### D2 · HIGH · Partículas en r ≤ −Δr durante un paso

- **Archivo/función:** `src/poisson_rk.f90:252-253,303-322` (ventana `mlo:mhi` calculada con `jc`
  de r, no de \|r\|); `src/utils.f90:242-243` (`build_cell_list` archiva toda partícula con r < 0
  en la celda 1); `src/density.f90:131` (ventana (n+2)/2).
- **Expresión matemática:** por la simetría F(r,p) = F(−r,−p), el depósito y la fuerza de una
  partícula en −x deben ser los de una en +x con la fuerza cambiada de signo.
- **Implementado:** para r_j ≤ −Δr, `mhi = jc+Wgrid+1` queda por debajo de los nodos que toca la
  imagen en +x; `Mself(m-mlo)` se lee fuera de lo calculado en esta iteración (valor de otra
  partícula del mismo hilo o sin inicializar). En el depósito, los nodos que toca la imagen quedan
  fuera de la ventana de la celda 1 cuando \|r_j\| > 1.5Δr (n = 1, 3) o > 2Δr (n = 2).
- **Evidencia:** U2: fuerza espuria de hasta 0.157 para m = 0.01 (n = 1, x = 1.2Δr); masa perdida
  10 %, 40 %, 90 %, 100 % en x = 1.6, 1.9, 2.4, 3.7 Δr (n = 1); el número de casos fallidos cambia
  con el número de hilos (17 con 1 hilo, 12 a 17 con 4). §8.4: alcanzable con `BGtype=Central`.
- **Cuándo ocurre:** cuando \|p\|Δt > Δr cerca del origen, es decir si `pmax` subestima la
  velocidad real (se revisa solo en t = 0) o `courant` > 1. No ocurre en el isócrono con
  `pmax ≥ 1`.
- **Corrección:** evaluar la autogravedad en \|r\| y aplicar el signo (F(r) = sign(r)·F(\|r\|));
  archivar cada partícula en la celda de \|r\|. Alternativa mínima: `jc` a partir de \|r_part\|
  en `poisson_rk` y `ic(j)` a partir de \|r_part(j)\| en `build_cell_list`. Con D1 revertido,
  la parte de la autofuerza desaparece pero la del depósito sigue.

### D3 · MEDIUM · Newton de Kepler diverge para e ≳ 0.98

- **Archivo/función:** `src/utils.f90:751-758`, `invert_QJ_to_rp`.
- **Matemática:** Q = η − e sen η, e ∈ [0, 1).
- **Implementado:** Newton desde η₀ = Q, 50 iteraciones, sin salvaguarda ni aviso.
- **Evidencia:** U6/T6c: fallas desde e = 0.98 (4/2000 ángulos) hasta 61/2000 en e = 1−10⁻⁸;
  residuo máximo 7.7·10²⁵; ida y vuelta con L = 0.01, J = 7.44 (e = 0.986): \|ΔQ\| = 0.36.
- **Cuándo ocurre:** e = √((1+2E)² + 2EL²) ≥ 0.98 exige E ≳ −0.01 (órbitas muy poco ligadas).
  No en las distribuciones de `distribution.f90` con sus ventanas por omisión.
- **Corrección:** arranque η₀ = Q + 0.85 e·sign(sen Q) (Danby) o η₀ = π para e > 0.8, con
  bisección de respaldo en [0, 2π] y aviso si no converge.

### D4 · LOW · Radicandos sin acotar en dos de las tres copias del mapa AA

- **Archivo/función:** `src/utils.f90:742-743` (`invert_QJ_to_rp`) y `:800-801`
  (`init_action_angle`).
- **Problema:** con L = 0 el radicando de r₁ es exactamente 0 y el redondeo lo hace negativo:
  NaN, y `max(NaN,0)` devuelve 0, así que r = 0 para todo Q. Es consecuencia directa del punto 9
  de la auditoría anterior (tres copias): las protecciones del punto 5 se agregaron en
  `analysish.f90` e `initial_data.f90` pero no aquí.
- **Evidencia:** U6: 650/1600 casos con L = 0 fallan; con L ≥ 1, \|ΔQ\| ≤ 1.3·10⁻¹¹.
- **Impacto:** esas partículas no tienen masa (m ∝ f·L), pero su trayectoria es falsa.
- **Corrección:** `max(...,0)` como en las otras copias; mejor, una sola rutina (punto 9 anterior).

### D5 · INFO · La energía no se conserva con autogravedad aun con Δt → 0

Ver §6. Esperado en un PIC con fuerza interpolada; documentarlo en las notas con los números.

### D6 · INFO · Costo de la resta de autofuerza

×2.09 (ABBA, 4+4 corridas, N = 19356). Desaparece con la corrección de D1.

### D7 · MEDIUM · El manuscrito ya no describe el código

`main.md` describe RK2 sobre (Φ, Φ′), el depósito con W_m = S_m∗b₀ y ΔV_k geométrico, L fija, y
tiene dos erratas en el leapfrog (Δt/2 en el drift; F(t_n, r_{n+1})). El código integra la masa,
usa W_n con V_i = 4πΔr(r_i² + (n+1)Δr²/12), tiene distribución en L y hace lo correcto en el
leapfrog. Corrección: actualizar §4 del manuscrito o documentar la diferencia (punto 20 anterior
pendiente de decisión).

### D8 · LOW · El mapa ángulo-acción del isócrono se usa aun con autogravedad

`analysish.f90` calcula h_k con (Q, J) del isócrono aunque el potencial incluya Φ_gas. Es una
elección de marco fijo, válida si se declara; el aviso de `main.f90:121` solo cubre
`BGtype ≠ Isochrone`. Corrección: aviso también con `autointeraction` y nota en el texto.

### D9 · INFO · Comentario de h_k impreciso

`analysish.f90:4` dice ∫FΦe^{−ikQ}; lo implementado es ∫F Φ̂_k(J)* e^{−ikQ} con Φ̂_k = a_k B C,
que es lo que define `main.md`. Corregir el comentario.

### D10 · LOW · El paso de tiempo no se revisa después de t = 0

`set_timestep` usa `pmax` del usuario y la fuerza de t = 0; el aviso `|p| > pmax` es solo en t = 0
(`main.f90:115`). Es la puerta de entrada de D2. Corrección: revisar max\|p\| en cada diagnóstico
y avisar (sin cambiar dt).

### D11 · LOW · `Nr` trunca por redondeo

`utils.f90:27`: `int((rmax-rmin)/dr)+1`. Con rmax = 0.3, dr = 0.1 el cociente es
2.9999999999999996 y la malla termina en 0.25 en vez de 0.35 (también rmax = 0.7). Corrección:
`nint` con tolerancia, o `ceiling(x − 1e-9)`.

### D12 · INFO · Cancelación en los pesos de la cuadratura radial

`utils.f90:200` forma I₁ como diferencia de cuartas potencias: pierde ~log₁₀(r/Δr) cifras. Medido:
10⁻¹⁴ relativo con 2001 nodos (U3). Si se quiere: I₁ = Δr²(r_a²/2 + 2r_aΔr/3 + Δr²/4).

### D13 · INFO · Salidas ASCII con 8 cifras

`save0Ddata`, `save_energy`, `save1Ddata`: `ES16.8`. No permiten medir conservación mejor que
10⁻⁸; las series de energía de HDF5 y `hk*_complex.tl` tienen precisión completa.

### D14 · LOW · Recorte inicial y renormalización

Los nodos con f ≤ cutoff·f_max se quitan y la masa restante se reescala a `a0`: la forma de F
cambia en la fracción recortada (punto 12 anterior, pendiente). Informar la masa recortada.

### D15 · INFO · Reproducibilidad entre números de hilos

Energías y h_k son deterministas para un número de hilos dado, no entre números distintos
(diferencias ≤ 6.4·10⁻¹⁶). Documentado en el código.

### D16 · INFO · Binario versionado

`exe/VP_PIC` está en git y aparece modificado; un binario versionado no corresponde a ningún
commit en particular. Recomendación: sacarlo del control de versiones.

### D17 · LOW · Pendientes de la auditoría anterior que siguen abiertos

7 (`courant > 1`), 9 (mapa AA triplicado; causa de D4), 12 (masa recortada, D14), 16 (fondo en
fantasmas: comprobado inocuo, la interpolación no lee los fantasmas), 17 (aviso de Newton; ahora
D3), 20 (convenio del depósito, D7).

---

## B. Tabla de problemas

| Id | Severidad | Archivo | Función | Problema | Evidencia | Corrección |
|---|---|---|---|---|---|---|
| D1 | HIGH | `poisson_rk.f90:236-338` | `poisson_rk` | Quita la autogravedad −m/2r² de cada cáscara: error O(1/N), costo ×2 | I2: 2.5·10⁻³ frente a 9·10⁻⁷ (N = 800); ABBA ×2.09 | **Revertido** |
| D2 | HIGH | `poisson_rk.f90:252-322`, `utils.f90:242`, `density.f90:131` | `poisson_rk`, `build_cell_list` | r ≤ −Δr: memoria no inicializada en la autofuerza, pérdida de masa en el depósito | U2; §8.4 (≥ 3363 casos con `Central`) | Evaluar en \|r\| y aplicar el signo |
| D3 | MEDIUM | `utils.f90:751-758` | `invert_QJ_to_rp` | Newton diverge para e ≳ 0.98 | U6: residuo 7.7·10²⁵; ΔQ = 0.36 | Arranque robusto + bisección + aviso |
| D7 | MEDIUM | `main.md` §4 | — | El manuscrito describe otro solver y otro depósito; erratas del leapfrog | Lectura | Actualizar el manuscrito |
| D4 | LOW | `utils.f90:742-743,800-801` | `invert_QJ_to_rp`, `init_action_angle` | Radicandos sin acotar: L = 0 da r = 0 | U6: 650/1600 | `max(...,0)`; unificar el mapa |
| D8 | LOW | `analysish.f90`, `main.f90:121` | `analysish` | Marco del isócrono con autogravedad sin aviso | Lectura | Aviso y nota |
| D10 | LOW | `main.f90:115`, `utils.f90:268` | `set_timestep` | \|p\| > pmax solo se revisa en t = 0 | §8.4 | Revisar en cada diagnóstico |
| D11 | LOW | `utils.f90:27` | `set_grid_size` | `int()` trunca Nr | Cálculo directo | Redondeo con tolerancia |
| D14 | LOW | `initial_data.f90` | `initial_data` | Recorte + renormalización sin informar | Lectura | Informar fracción recortada |
| D17 | LOW | varios | — | Pendientes 7, 9, 12, 16, 17, 20 de la auditoría anterior | — | Ver AUDITORIA_2026-09-20.md |
| D5 | INFO | `energy.f90`, esquema | — | Energía no conservada con autogravedad aun con Δt → 0 | Plummer: 4.6·10⁻⁴ independiente de Δt | Documentar; variante variacional si hace falta |
| D6 | INFO | `poisson_rk.f90` | — | Costo ×2 de la autofuerza | ABBA | **Resuelto** con D1 |
| D9 | INFO | `analysish.f90:4` | — | Comentario Φ en vez de Φ̂_k | Lectura | Corregir comentario |
| D12 | INFO | `utils.f90:200` | `construct_grid` | Cancelación en I₁ | U3: 10⁻¹⁴ | Forma expandida |
| D13 | INFO | `utils.f90` | `save*` | ASCII con 8 cifras | — | Usar HDF5 para conservación |
| D15 | INFO | `energy.f90`, `analysish.f90` | — | Bits dependen del número de hilos | R2 | Documentado |
| D16 | INFO | `exe/VP_PIC` | — | Binario versionado | `git status` | Quitar del repositorio |

## C. Tabla de equivalencias

| Expresión original | Expresión actual | ¿Equivalentes? | Condiciones | Impacto numérico |
|---|---|---|---|---|
| `-r/sqrt(1+r²)*pot²`, pot = −1/(1+s) | `-r/(sq*(1+sq)**2)` | Sí, exacta en ℝ | Todo r | Último bit en 55 % de los radios; hace idénticas las dos ramas |
| `-1/r`, `-1/r²` (Central) | `-1/abs(r)`, `-sign(1,r)/r²` | Sí para r > 0 (bit a bit); extiende par/impar a r < 0 | r ≠ 0 | Ninguno para r > 0; corrige el signo para r < 0 |
| Fondos evaluados en r | Evaluados en \|r\|, fuerza con signo de r | Sí para r > 0 | r ≠ 0 | Corrige r < 0 (`nfw` invertía el signo) |
| `f/(dr+dr³/(12r²))/r²` | `f/(r²dr + (n+1)dr³/12)` | **No**: cambia 1/12 por (n+1)/12 | — | Relativo (n)Δr²/(12r²); legítimo (volumen correcto del núcleo) |
| RK2 sobre (Φ, Φ′) | M con pesos `mcoefA/B`, Φ con c₀, c₃, c₄ | **No**: otro método | — | Exacto para el interpolante; mejora cerca del origen |
| M(r_b) por c₀ + c₃r_b³ + c₄r_b⁴ | M(r_b) por `mcoefA`ρ_{i−1} + `mcoefB`ρ_i | Sí en ℝ | — | ≤ 10⁻¹⁴ con 2000 nodos (U3) |
| Σ_k w_k A_k B C cos(kQ_k) por partícula | a_k = Σ_k w_k A_k cos(kQ_k) una vez, × B C | Sí en ℝ | A independiente de la partícula | Último bit |
| `exp(-i k Q)` | `expv**k`, expv = exp(−iQ) | Sí en ℝ | k ≤ 4 | ~k·ε |
| Reducción OpenMP | Sumas por hilo en orden fijo | Sí en ℝ | — | ≤ 6.4·10⁻¹⁶; determinista por número de hilos |
| Suma j = 1…N por nodo | Suma por celdas, luego por índice | Sí en ℝ | r ≥ −1.5Δr | Último bit; **no** equivalente para r < −1.5Δr (D2) |
| `r = r_p + c·dt·p_p` (arreglos) | `r(i) = r(i) + c*dt*p(i)` (ciclo) | Sí, misma secuencia | — | Ninguno (bit a bit) |
| `sqrt(1+r²)` tres veces | `sq` una vez | Sí | — | Ninguno (bit a bit) |
| E = Σ(p²/2 + pot)fL, sin ½ | E = Σ(p²/2 + pot − ½pot_self)fL | **No**: corrige la autoenergía | — | Legítimo (`abb8bc9`) |
| K = Σ p²/2, W = Σ pot | K incluye L²/2r², W no | **No**: redefine el reparto | — | Total idéntico; K y W cambian (R1: 7.7·10⁻³) |
| Fuerza de la malla interpolada | La misma menos la de la propia partícula | **No** | — | **Cambio de física (D1)** |

## D. Tabla de validación

| Prueba | Resultado esperado | Resultado obtenido | Error | Estado |
|---|---|---|---|---|
| U1 Funciones de forma | Σ W = 1, área 1, m₂ = (n+1)/12, simetría | todas | ≤ 3·10⁻¹⁶ (m₂ de n = 1: 1.7·10⁻⁹ por la cuadratura del quiebre) | PASA |
| U2 Simetría en el origen | depósito y fuerza de −x = los de +x | rota para x ≥ Δr | fuerza hasta 0.157; masa hasta 100 % | **FALLA (D2)** |
| U3 Poisson exacto para el interpolante | error de redondeo | 41 nodos: 3.9·10⁻¹⁶ (F), 1.2·10⁻¹⁵ (Φ); 2001 nodos: 9.1·10⁻¹⁵, 1.0·10⁻¹⁴ | redondeo | PASA |
| U4 Poisson analítico ρ₀(1−r²)³ | orden 2 | 2.00 / 2.00 / 2.00 (n = 1, 2, 3) | L∞(F) = 1.4·10⁻⁴ (n = 1, Δr = 0.005) | PASA |
| T3c Esfera uniforme, interior | F = −r exacta en la malla | exacta; en partículas 2πρ₀Δr/K | 9.1·10⁻⁴ → 1.2·10⁻⁴ (sesgo D1) | PASA en malla, D1 en partículas |
| U5 Fondos | F = −dΦ/dr, paridad exacta | 7 fondos | ≤ 6.7·10⁻¹²; paridad 0 | PASA |
| T6b Fórmulas del isócrono | = cuadratura | J, Ω, Q | ≤ 3·10⁻¹⁶ | PASA |
| U6 Ida y vuelta del mapa AA | \|ΔQ\| < 10⁻⁸ | L ≥ 1: 1.3·10⁻¹¹; L = 0: 3.06; e > 0.98: 0.36 | — | **FALLA (D3, D4)** |
| I1 Orden temporal | 1, 2, 4 | 1.02, 2.00, 4.00 | — | PASA |
| I1 Energía, fondo fijo | acotada, ∝ Δt^k | leapfrog ×4, yoshida4 ×16 por mitad de Δt | 3.2·10⁻¹² (yoshida4, Δt = 0.00625) | PASA |
| I2 Colapso frío | error ≪ 1/N | HEAD: 2.5·10⁻³ (N = 800), orden 1.00 | — | **FALLA (D1)** |
| I2 Colapso frío con `5c0e07f` | error ≪ 1/N | 9.1·10⁻⁷ (N = 800), orden 2.1 | — | PASA |
| Plummer N = 154270 | estacionario, E conservada | radios estables dentro del muestreo; \|ΔE/E\| = 8.4·10⁻⁶ en t = 50 | — | PASA (con el límite de D5) |
| Plummer, barrido en Δt | error de E ∝ Δt² si fuera temporal | 4.6·10⁻⁴ constante | — | Demuestra D5 |
| R1 Original frente a actual, sin autogravedad | trayectorias iguales | 8 cifras iguales, 2000 pasos | — | PASA |
| R2 Pares de commits de optimización | bit a bit | 7 idénticos, 2 a redondeo, 1 cambio semántico documentado | — | PASA con observaciones |

## E. Pruebas automatizadas

En `verificacion/` (nuevo, sin confirmar en git):

    make                                  # compila objs/ y exe/VP_PIC
    verificacion/correr.sh [--rapido]     # U1-U6 (segundos) + I1, I2 (~1 min)
    verificacion/regresion.sh A B         # dos ejecutables, tres casos, bit a bit

| Prueba | Archivo | Qué detecta | Estado en HEAD |
|---|---|---|---|
| U1 | `src/t_forma.f90` | Pesos W_n | PASA |
| U2 | `src/t_origen.f90` | D2 | FALLA (tras revertir D1, solo los 10 casos del depósito) |
| U3 | `src/t_poisson_nodal.f90` + `py/ref_nodal.py` | Exactitud del integrador radial (mpmath, 50 dígitos) | PASA |
| U4 | `src/t_poisson_orden.f90` | Orden 2 de Poisson frente a solución analítica | PASA |
| U5 | `src/t_fondos.f90` | F = −dΦ/dr y paridad de los fondos | PASA |
| U6 | `src/t_aa.f90` | D3, D4 | FALLA |
| I1 | `py/integradores.py` | Orden de euler, leapfrog, yoshida4 frente a `analytic` | PASA |
| I2 | `py/colapso.py` | D1 (física de la autogravedad frente al continuo) | FALLA antes de revertir D1; PASA después |
| R | `regresion.sh` + `py/h5igual.py` | Neutralidad bit a bit de un cambio | — |

Las pruebas unitarias enlazan los objetos de `objs/` (todo menos `main.o`) y llaman a las
rutinas del código, así que prueban exactamente lo que se compila. Cada una imprime `PASA` o
`FALLA` con sus números; `correr.sh` devuelve el número de fallas. Con D1, D2, D3 y D4 corregidos,
las cuatro fallas deberían pasar sin tocar las pruebas.

## F. Lo que no se pudo verificar y qué haría falta

| Qué | Por qué no | Prueba necesaria |
|---|---|---|
| Evolución larga con autogravedad en el régimen del artículo | No hay solución analítica | Convergencia en N, Δr, Δt de h_k(t) con autogravedad y una referencia de mayor resolución |
| h_k frente al semianalítico | Lo hicieron sesiones anteriores (`reproducir/corridas/08_verificacion`); no se repitió | Volver a correr `08_verificacion` y comparar con `tools/hk_exacto.py` |
| Origen del error de energía con autogravedad | No se implementó la variante variacional | Fuerza = −∂/∂r_j Σ W Φ_i y comparar el piso de energía |
| Sensibilidad al tamaño del dominio con masa saliendo | No se varió r_max | Barrido en r_max con una distribución que se expande |
| `reduceparticles`, `dt_switch=var` | No se probaron | Casos dedicados |
| n = 2, 3 en dinámica | Solo en Poisson, depósito y regresión | Colapso frío con n = 2, 3 |
| `tools/*.py` (generador de equilibrio, h_k exacto) | Fuera del alcance de esta pasada | Auditoría propia |
| Regresión del original con autogravedad paso a paso | El original tiene RK2, depósito, energía y acoplamiento distintos | Se hizo por commits intermedios (§8.2), no con el original directo |

## G. Relación con AUDITORIA_2026-09-20.md

- **Punto 1 (autofuerza): invalidado.** La premisa ("artefacto") es falsa en simetría esférica;
  la corrección aplicada es D1 de esta auditoría. La prueba de aceptación de la partícula libre
  tiene la referencia equivocada.
- **Punto 2 (masa perdida en Poisson): confirmado y bien corregido** (`5c0e07f`; U3 e I2).
- **Puntos 3, 4, 15 (fondos `null`, paridad, dos ramas del isócrono): confirmados** (U5, R2).
- **Punto 5 (radicandos): incompleto.** Se protegieron dos copias del mapa y quedó la tercera (D4).
- **Punto 6 (memoria en `analysish`): neutro bit a bit** (R2 para `a646629`).
- **Punto 17 (aviso de Newton): subestimado.** No es solo falta de aviso: Newton diverge (D3).
- **Punto 20 (convenio del depósito): sigue pendiente** (D7).
- **`vlasov_L_intro.tex` §7.4:** documenta D1 como corrección; hay que reescribirlo.
