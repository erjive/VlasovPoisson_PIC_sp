# Auditoría del modelo, la física y la numérica de VP_PIC en (r, p_r, L)

**Fecha:** 2026-09-21. **Versión auditada:** rama `mejoras/portar`, commit `612f281`
(después de revertir la resta de la autofuerza). **Referencias:** `66bc5bb` (antes de la
reversión), `5c0e07f`, el original `84fc3b2` (2024) y el manuscrito
`Vlasov_Poisson_evolutions/main.md`.

**Relación con los documentos anteriores.** `AUDITORIA_2026-09-20.md` fue una revisión por
lectura; `AUDITORIA_CIENTIFICA_2026-09-21.md` auditó las optimizaciones y encontró que la resta
de la autofuerza era un error (ya revertido). Esta auditoría es la tercera: reconstruye el modelo
desde la medida del espacio de fases, se concentra en r = 0 y L = 0, vuelve a medir todo lo que
depende de la autogravedad con el binario actual y corre las simulaciones de `reproducir/`. Las
cifras que se reutilizan de la auditoría anterior se marcan y se justifica por qué siguen
valiendo.

**Reglas.** No se modificó `src/`. Una simulación a la vez. Cada afirmación lleva su prueba y su
número; donde no hubo prueba, se dice. Las copias instrumentadas se compilaron aparte.

Las notas del código que acompañan este informe están en `docs/codigo/notas_codigo.tex` (PDF al
lado).

---

## Resumen

| Categoría | Qué entra |
|---|---|
| **Matemáticamente demostrado** | Medida d³x d³v = 8π² L dr dp_r dL y ρ = (2π/r²)∫∫F L dp dL; características ṙ = p_r, ṗ_r = L²/r³ − Φ′, L̇ = 0; invariancia de la medida (divergencia nula y L constante); equivalencia de las formas de Poisson; límite regular en r → 0 del isócrono y de la esfera; exactitud del integrador radial para su interpolante |
| **Validado numéricamente** | Densidad de Plummer desde su F(E) (orden ≈ 2 en Δr, Δp, ΔL); campo analítico (orden 2.00) y fuerza sobre las partículas (orden 2.00); órbitas con L > 0 y L = 0 frente a RK4 independiente (orden 2.00 y 4.00); L de cada partícula y masa total constantes bit a bit; energía acotada con fondo fijo; colapso frío O(1/N²); equilibrio de Plummer autogravitante; `reproducir/` frente a sus soluciones exactas |
| **Consistente, no probado** | Evolución larga con autogravedad en el régimen del artículo; convergencia de h_k con autogravedad |
| **Sospechoso** | Energía de las partículas de L pequeño con paso fijo (medido: hasta 10⁶ relativo; el código no lo detecta) |
| **Incorrecto** | r = 0 exacto con L = 0 da NaN (debería ser el límite regular); `nfw`, `burkert` e `isotrun` pierden todas las cifras cerca de r = 0; depósito pierde masa para r ≤ −1.5Δr; Newton de Kepler diverge para e ≳ 0.98; mapa AA con L = 0 en `utils.f90` |

Ninguno de los defectos nuevos alcanza a las corridas del artículo (isócrono, sin autogravedad,
L ≈ 2).

---

## A. Modelo matemático reconstruido

### A.1 Qué es f en el código

Leído de `initial_data.f90`, `density.f90`, `energy.f90`, `analysish.f90`:

- `f(j)` es el **valor de la función de distribución tridimensional** en el nodo j:
  f(x, v) = F(r, p_r, L). No es una densidad respecto de dr dp dL: el factor L y el 8π² se
  aplican aparte en cada suma.
- La masa de la partícula j es **m_j = 8π² Δr_c Δp_c ΔL_c f_j L_j** (en `density.f90`
  aparece como 4π·`factor`·f·L con `factor` = 2π Δr_c Δp_c ΔL_c; en `energy.f90` como
  `factor` = 8π² Δr_c Δp_c ΔL_c por f·L).
- La normalización se fija en `initial_data`: f ← a0·f/(8π² Δr_c Δp_c ΔL_c Σ f L), de modo que
  Σ m_j = `a0`. En `aa_quad` y `checkpoint` el producto Δr_c Δp_c ΔL_c es una constante
  convencional que se cancela (la cuadratura real es en (J, Q, L) o la del archivo).
- r ≥ 0 es el radio; p_r = v_r es el momento radial por unidad de masa; L = |x × v| ≥ 0 es el
  módulo del momento angular por unidad de masa (m ≡ 1 en las ecuaciones de movimiento).

### A.2 Unidades y acoplamiento

G = M_iso = b = 1 (`main.md` §2.2): Φ_iso = −1/(1 + √(1+r²)), tiempo t_c = √(b³/GM). La
ecuación de Poisson es ∇²Φ = 4πρ (el 4π está en c₃, c₄ de `poisson_rk.f90` y en el
(4/3)π r₁³ del primer nodo). La masa del gas es `a0` en unidades de M_iso.

### A.3 La medida del espacio de fases (demostración)

En coordenadas esféricas, d³x = r² sen θ dr dθ dφ. En el espacio de velocidades, con v_r radial
y la componente tangencial en coordenadas polares (v_t, ψ), d³v = dv_r v_t dv_t dψ. Como
L = r v_t a r fijo, v_t dv_t = L dL / r². Entonces

    d³x d³v = sen θ dθ dφ dψ · dr dp_r L dL  →  (4π)(2π) · L dr dp_r dL = 8π² L dr dp_r dL.

Con simetría esférica f no depende de (θ, φ, ψ), y

    M   = 8π² ∫∫∫ F L dr dp_r dL,
    ρ(r) = ∫ F d³v = (2π/r²) ∫∫ F L dp_r dL,
    4πr² ρ(r) dr = 8π² (∫∫ F L dp_r dL) dr   (consistente con la anterior).

**Respecto de qué medida está definida f:** de la medida 8π² L dr dp_r dL. El factor r² de d³x
se cancela exactamente con el 1/r² de d²v_t; por eso **no** aparece r² en la masa de las
partículas, y sí aparece 1/r² en la densidad. El código hace exactamente esto: m_j ∝ f_j L_j,
sin r² (`density.f90:154`), y ρ_i = Σ m_j W/V_i con V_i ≈ 4πr_i²Δr, que es la forma discreta de
dM = 4πr²ρ dr.

**Invariancia.** El flujo (ṙ, ṗ_r, L̇) = (p_r, L²/r³ − Φ′, 0) tiene divergencia
∂ṙ/∂r + ∂ṗ_r/∂p_r + ∂L̇/∂L = 0 + 0 + 0, así que dr dp_r dL se conserva; como L es constante a lo
largo de cada trayectoria, también L dr dp_r dL. Por eso el peso f_j y la masa m_j de cada
partícula son constantes: no hace falta ningún término extra en la ecuación reducida.

**Verificación numérica (prueba M).** Se muestreó la F(E) de Plummer isótropo,
F = (24√2/7π³)(−E)^{7/2}, en mallas regulares de (r, p, L) y se corrió el código sin avanzar
(Nt = 0). La densidad depositada frente a ρ = (3/4π)(1+r²)^{−5/2} para 0.3 < r < 7:

| (N_r, N_p, N_L), Δr | máx rel | rms rel | K (3π/64 = 0.14726) | W (−3π/32 = −0.29452) |
|---|---|---|---|---|
| (100, 40, 20), 0.4 | 1.17·10⁻¹ | 3.23·10⁻² | 0.14822 | −0.29730 |
| (200, 80, 40), 0.2 | 3.90·10⁻² | 1.09·10⁻² | 0.14757 | −0.29597 |
| (400, 160, 80), 0.1 | 1.12·10⁻² | 2.67·10⁻³ | 0.14732 | −0.29489 |

Orden 1.58 y 1.80 en el máximo (tiende a 2). Un factor L, 2π o r² equivocado daría un error O(1)
que no baja al refinar. La diferencia de potencial Φ(r) − Φ(r₀) converge igual
(4.2·10⁻² → 9.1·10⁻³ → 1.7·10⁻³). **La medida está bien en todo el camino: condiciones
iniciales → masa → depósito → densidad → Poisson → energía.**

### A.4 Ecuaciones características

Hamiltoniano por unidad de masa, a L fijo: H = p_r²/2 + L²/(2r²) + Φ(r, t), con
Φ = Φ_ext + Φ_gas. Ecuaciones de Hamilton:

    ṙ = ∂H/∂p_r = p_r,
    ṗ_r = −∂H/∂r = L²/r³ − ∂Φ/∂r,
    L̇ = 0  (Φ esférico: no hay torque).

Vlasov reducida: ∂_t F + p_r ∂_r F + (L²/r³ − ∂_rΦ) ∂_{p_r} F = 0.

**En el código** (`grav_force.f90`): `force_part` = fondo + autogravedad + `l_part²·r_part/den²`
con den = r² + eps² y eps = 0 forzado (`utils.f90:35`), es decir +L²/r³. `pot_part` suma
`0.5·l_part²/den` = L²/(2r²). Signos, factores y dependencias correctos. Dimensiones: [L²/r³] =
[v²/r] = aceleración. L nunca se modifica en ningún integrador (`main.f90`).

**Verificación numérica (prueba C).** Tres órbitas en el isócrono (L = 0.1, 0.3, 1; r₀ = 1.5,
p₀ = 0) hasta t = 50 frente a una integración RK4 independiente de las mismas ecuaciones con
h = 5·10⁻⁵ (autoconsistencia 10⁻¹³):

| L | leapfrog Δt = 0.02, 0.01, 0.005 | orden | yoshida4 | orden |
|---|---|---|---|---|
| 0.1 | 3.3·10⁻⁴, 8.3·10⁻⁵, 2.1·10⁻⁵ | 2.00 | 5.5·10⁻⁶, 3.4·10⁻⁷, 2.1·10⁻⁸ | 4.00 |
| 0.3 | 8.9·10⁻⁵, 2.2·10⁻⁵, 5.6·10⁻⁶ | 2.00 | 1.0·10⁻⁷, 6.3·10⁻⁹, 3.9·10⁻¹⁰ | 4.00 |
| 1.0 | 1.1·10⁻⁴, 2.8·10⁻⁵, 7.0·10⁻⁶ | 2.00 | 1.9·10⁻⁸, 1.2·10⁻⁹, 7.6·10⁻¹¹ | 3.99 |

Esta prueba no depende del mapa ángulo-acción del código.

### A.5 Poisson esférica

Formas equivalentes (para ρ integrable y Φ regular en el origen):

    (1/r²) d/dr (r² dΦ/dr) = 4πρ
    ⇔ r² Φ′(r) = ∫₀^r 4πs²ρ ds ≡ M(r)       (integrando una vez, con r²Φ′ → 0 en r = 0)
    ⇔ dM/dr = 4πr²ρ,  dΦ/dr = M/r²,  g(r) = −dΦ/dr = −M(r)/r².

El código usa la segunda forma (`poisson_rk.f90:108-136`): ρ lineal entre nodos, M(x) =
c₀ + c₃x³ + c₄x⁴ con c₃ = (4π/3)(ρ_{i−1} − s r_{i−1}), c₄ = πs, y ∫M/x² dx = −c₀/x + c₃x²/2 +
c₄x³/3 en forma cerrada. `dev_pot` guarda M y al final M/r²; `force = −dev_pot` = −M/r² = g.
Signo correcto (atractivo).

**Condiciones de frontera.** En r = 0: M(0) = 0 y densidad constante en [0, r₁] (la densidad es
par y el nodo fantasma r₀ = −r₁ lleva el mismo valor), de donde M(r₁) = (4/3)πρ₁r₁³ y
Φ(r₁) − Φ(0) = (2/3)πρ₁r₁², la solución regular Φ = Φ₀ + (2π/3)ρ₀r² + O(r⁴) de `main.md`. En
el borde exterior: Φ(r_N) = −M(r_N)/r_N, impuesto desplazando Φ (`poisson_rk.f90:162`), y más
allá la solución exterior Φ = Φ(r_N)r_N/r, g = g(r_N)(r_N/r)².

**Verificación.** (U3) exacto para su interpolante: 5.5·10⁻¹⁶ con 41 nodos y 1.1·10⁻¹⁴ con 2001
(referencia en 50 dígitos). (U4) ρ = ρ₀(1−r²)³: orden 2.00 en g y Φ para n = 1, 2, 3.

### A.6 El origen, r → 0

Todas las divisiones por r en el código:

| Expresión | Dónde | ¿Puede r → 0? | Límite correcto | Implementado |
|---|---|---|---|---|
| V_i = 4πΔr(r_i² + (n+1)Δr²/12) | `density.f90:142` | No: r_i ≥ Δr/2 | — | — |
| M/r² | `poisson_rk.f90:140` | No: nodos | — | — |
| −c₀/x | `poisson_rk.f90:133` | No: x ≥ r₁ | — | — |
| L²/(2r²), L²r/r⁴ | `grav_force.f90`, `energy.f90` | **Sí** (partículas) | L = 0: 0; L > 0: barrera, inalcanzable | **0/0 = NaN con L = 0 y r = 0 exacto** |
| Isócrono −r/(s(1+s)²) | `grav_force.f90` | Sí | 0, regular | Exacto hasta r = 10⁻¹⁶ |
| `Central` −1/\|r\|, −sgn/r² | `grav_force.f90` | Sí | singular (físico) | Inf, correcto |
| `sphere` | `grav_force.f90` | Sí | regular | Exacto |
| `iso` 3 ln r | `grav_force.f90` | Sí | singular (físico) | Correcto |
| `isotrun` (a − atan a)/a² | `grav_force.f90:236` | Sí | −(10/9)a | **0 para a ≤ 10⁻⁸** (cancelación) |
| `nfw` (ln(1+a) − a/(1+a))/a² | `grav_force.f90:238` | Sí | −8 | **+1.72 en a = 10⁻⁸ (signo invertido), −1.3·10⁴ en 10⁻¹⁰** |
| `burkert` | `grav_force.f90:240` | Sí | −(40/9)a → 0 | **+10.7 en a = 10⁻⁸, −5.5·10³ en 10⁻¹⁰**; error relativo 5·10⁻⁵ ya en a = 10⁻⁴ |
| 0.5 L²/r² en el mapa AA | `analysish.f90:364`, `utils.f90:793` | Sí | 0 con L = 0 | NaN con r = 0 |

El primer nodo de la malla y el tratamiento de Poisson en [0, r₁] tienen el límite correcto
(A.5). Los problemas están en expresiones que se evalúan **en la partícula**. Con r = 0 exacto y
L = 0, `pot_part` y `force_part` salen NaN para **los siete fondos**, incluidos el isócrono y la
esfera, cuyo límite es regular (Φ = −½, g = 0). La partícula tiene masa 0 (m ∝ L), así que el NaN
no entra al depósito, pero sí a la energía y a h_k: NaN·0 = NaN (prueba R0, §D).

**Muestreo cerca del origen.** A radio r solo son accesibles L ≤ r·v_esc(r). Con una malla
uniforme en L de paso ΔL, los primeros radios (r ≲ ΔL/v_esc) tienen uno o ningún nodo de L
dentro del intervalo permitido, y la cuadratura en L no converge uniformemente allí: en la
prueba M el error del primer nodo fue 1.5·10⁻¹, 1.3·10⁻², 6.6·10⁻² en las tres resoluciones
(no monótono). No es un defecto del solver sino del muestreo regular en L.

### A.7 L = 0 y L → 0

- **Ecuaciones:** con L = 0 el término centrífugo desaparece y la órbita pasa por el centro. En
  (r, p_r), pasar por el centro es exactamente la reflexión (r, p_r) → (−r, −p_r) que aplica
  `main.f90:239-249` (junto con g → −g, porque g es impar). **Validado:** órbita L = 0 con 5
  cruces por el centro frente a RK4 en x (1D, x = ±r): leapfrog 1.50·10⁻⁵, 3.74·10⁻⁶,
  9.35·10⁻⁷ (orden 2.00), yoshida4 3.1·10⁻¹⁰, 1.9·10⁻¹¹, 1.2·10⁻¹² (orden 4.00); energía
  7.5·10⁻¹² con yoshida4, Δt = 0.01.
- **Medida:** el peso es ∝ L, así que L = 0 tiene masa nula; los nodos en punto medio evitan L = 0
  en los estados `gaussian1`, `aa`, `aa_quad`, no en `checkpoint`.
- **Interpolaciones:** no dependen de L.
- **Singularidades artificiales:** r = 0 exacto con L = 0 (A.6); el mapa AA de `utils.f90` con
  L = 0 (D4: el radicando de r₁ es exactamente 0 y el redondeo lo hace negativo).
- **L → 0 con paso fijo (el hallazgo principal de esta sección).** La órbita con L pequeño pasa
  por un pericentro r_p ~ L/v en un tiempo ~ L/v², donde la fuerza centrífuga L²/r³ es enorme. Un
  integrador de paso fijo tiene que resolver ese tiempo; los criterios de `set_timestep`
  (Courant con `pmax`, aceleración en t = 0) no lo ven. Mismo r₀ = 1.5, p₀ = 0, isócrono, t = 50,
  máx |ΔE/E₀|:

| L | 0 | 10⁻⁴ | 10⁻³ | 10⁻² | 3·10⁻² | 0.1 | 0.3 | 1 |
|---|---|---|---|---|---|---|---|---|
| leapfrog, Δt = 0.01 | 1.8·10⁻⁶ | 5.6·10⁶ | 3.5 | 1.7·10⁻³ | 2.7·10⁻⁴ | 2.5·10⁻⁵ | 1.2·10⁻⁶ | 2.1·10⁻⁶ |
| yoshida4, Δt = 0.01 | 7.5·10⁻¹² | 38 | 200 | 5.4·10⁻⁴ | 9.7·10⁻⁶ | 7.5·10⁻⁸ | 3.1·10⁻¹⁰ | 7.7·10⁻¹¹ |
| yoshida4, Δt = 0.005 | 4.7·10⁻¹³ | 69 | 200 | 3.5·10⁻⁵ | 6.0·10⁻⁷ | 4.7·10⁻⁹ | 2.0·10⁻¹¹ | 4.8·10⁻¹² |

  El límite L → 0 **no es uniforme** a Δt fijo: L = 0 es exacto y regular, L ≥ 0.1 converge con
  su orden, y en medio la energía de la órbita se destruye. Con `lminc = 0` (valor por omisión)
  los nodos más bajos de L caen en esa franja; llevan poca masa (∝ L), pero su dinámica está mal
  y el código no lo avisa. El criterio que falta es Δt ≲ η r_p²/L (tiempo de paso por el
  pericentro).

---

## B. Correspondencia matemática → código

| Ecuación | Discretización | Archivo:línea | Verificación |
|---|---|---|---|
| m_j = 8π² Δr_c Δp_c ΔL_c f_j L_j | punto medio de la medida | `initial_data.f90:81,178,239,285` | Prueba M |
| ρ = (2π/r²)∫∫F L dp dL | Σ m_j [W(r_i−r_j)+W(r_i+r_j)]/V_i | `density.f90:136-172` | M, U2, U4 |
| V_i = ∫W 4πr² dr | 4πΔr(r_i² + (n+1)Δr²/12) | `density.f90:142` | U1 (segundo momento), T-A |
| dM/dr = 4πr²ρ, Φ′ = M/r² | ρ lineal, integral cerrada | `poisson_rk.f90:108-140` | U3, U4 |
| M(r₁) = (4/3)πρ₁r₁³ | densidad par | `poisson_rk.f90:117-118` | U3 |
| Φ(r_N) = −M/r_N | desplazamiento | `poisson_rk.f90:162` | U4 |
| g(r_j) | Σ W g_i con espejo y exterior | `poisson_rk.f90:173-213` | U4, T-A, I2 |
| ṙ = p_r | drift | `main.f90:170,180,195-219` | I1, C |
| ṗ_r = L²/r³ − Φ′ | kick | `main.f90:171,179,184,202-218`; `grav_force.f90` | I1, C, U5 |
| L̇ = 0 | L no se toca | todo `main.f90` | Prueba D (bit a bit) |
| F(r,p) = F(−r,−p) | reflexión | `main.f90:239-249` | C (L = 0) |
| E = Σ m(p²/2 + L²/2r² + Φ_ext + ½Φ_self) | suma | `energy.f90:60-67` | E, §F |
| h_k = 8π²∫F Φ̂_k* e^{−ikQ} L dr dp dL | a_k por Simpson; mapa del isócrono | `analysish.f90` | 08_verificacion (§13) |

---

## C. Auditoría de optimizaciones y diff científico

Las transformaciones algebraicas se auditaron en `AUDITORIA_CIENTIFICA_2026-09-21.md` (§3, tabla
C) compilando padre e hijo y comparando bit a bit. Ninguno de esos commits tocó el camino sin
autogravedad después de `86b3c80`, y la reversión solo afectó `poisson_rk.f90`,
`arrays.f90`, `utils.f90` (pesos `mcoef`); por eso aquellas pruebas siguen valiendo, salvo lo que
se marca.

| Original | Actual | Equivalentes | Condiciones | Riesgo numérico |
|---|---|---|---|---|
| `-r/sqrt(1+r²)*pot²` | `-r/(sq*(1+sq)**2)` | Sí (ℝ) | todo r | último bit; ambos regulares en r = 0 |
| `-1/r`, `-1/r²` | `-1/abs(r)`, `-sign(1,r)/r²` | Sí para r > 0 (bit a bit) | r ≠ 0 | extiende la paridad a r < 0 |
| fondos en r | fondos en \|r\|, g con signo de r | Sí para r > 0 | r ≠ 0 | `nfw`, `burkert`, `isotrun` mal condicionados en r → 0 (ya lo estaban en r) |
| `L²/(2r²)`, `L²r/(r²+eps²)²` | igual (eps = 0 forzado) | — | r ≠ 0 | 0/0 con L = 0, r = 0 (también en el original) |
| `f/(dr+dr³/(12r²))/r²` | `f/(r²dr+(n+1)dr³/12)` | **No** (1/12 → (n+1)/12) | — | corrige el volumen del núcleo; legítimo |
| RK2 en (Φ, Φ′) | M con c₀, c₃, c₄ cerrados | **No** (otro método) | — | exacto para el interpolante; sin 2/r |
| (pesos `mcoefA/B`, `66bc5bb`) | c₀ + c₃r³ + c₄r⁴ (`5c0e07f`, actual) | Sí (ℝ) | — | ambos 10⁻¹⁴ con 2001 nodos (U3) |
| g propio restado (`66bc5bb`) | g propio incluido (actual) | **No** | — | la versión actual es la correcta (colapso) |
| a_k·Σ B C e^{−ikQ} por partícula en cuadratura | a_k una vez | Sí (ℝ) | A independiente de la partícula | último bit |
| `exp(-ikQ)` | `expv**k` | Sí (ℝ) | k ≤ 4 | ~kε |
| reducción OpenMP | sumas por hilo en orden fijo | Sí (ℝ) | — | 6·10⁻¹⁶; determinista por número de hilos |
| bucles de arreglos | bucles fundidos (yoshida4, fondos) | Sí, misma secuencia | — | idéntico bit a bit |
| E sin ½ | E con ½Φ_self | **No** | — | corrige la autoenergía |

**Diff científico frente a versiones anteriores**, medido con el binario actual:

| Cambio | Ecuación efectiva | ¿Cambió el resultado físico? | Evidencia |
|---|---|---|---|
| Nodos en el borde → punto medio de celda | misma F, otra cuadratura | Sí, O(Δ) en las condiciones iniciales | R1 de la auditoría anterior |
| Reflexión con g → −g (`a81f4cf`) | misma | Sí, corrige el medio kick tras cruzar r = 0 | Prueba C con L = 0 |
| Imágenes en el depósito, espejo en la interpolación (`a81f4cf`) | misma | Sí, corrige ρ y g cerca de r = 0 | colapso: 8·10⁻⁵ → 2.5·10⁻⁶ |
| Volumen (n+1)/12 | misma | Sí, ρ uniforme exacta | T-A |
| Poisson por masa (`5c0e07f`) | misma, sin 2/r | Sí, masa completa cerca del origen | U3; colapso 2.5·10⁻⁶ → 9.1·10⁻⁷ |
| Autofuerza restada y revertida | **cambiaba la física**; revertido | — | colapso 2.5·10⁻³ frente a 9.1·10⁻⁷ |
| ½ en la autoenergía | energía correcta | diagnóstico | `abb8bc9` |
| φ̂_k con el modo correcto, 512 intervalos | h_k correcto | diagnóstico | `2f345bc`, `5e11235` |
| Fuerza del isócrono en una expresión | misma | no (último bit) | 08: ≤ 3.5·10⁻¹² en h_k |

---

## D. Errores encontrados

**Estado de las correcciones** (rama `mejoras/portar`, desde el 2026-10-04). Son las de la
auditoría de `vlasov-poisson_PIC` (`AUDITORIA_L0_2026-09-21.md`, hallazgos E) que faltaban
aquí. Cada una lleva su medición. "Neutralidad" quiere decir: las 26 configuraciones de
`reproducir/corridas` recortadas a 200 pasos (283 archivos de salida) y los tres casos de
`verificacion/regresion.sh`, con el binario anterior al cambio y con el nuevo, a 4 hilos.

| Id | Estado | Medición |
|---|---|---|
| D3, D4 | corregido | Una función `kepler_eta(Q, e, tol)` en `utils` (la de `vlasov-poisson_PIC`, E12) resuelve Q = η − e sen η para `invert_QJ_to_rp`: primero el Newton de siempre, con las mismas operaciones; se acepta si salió por tolerancia con un último paso \|g/g′\| ≤ √tol y η ∈ [Q − e, Q + e], y si no, se resuelve con Newton acotado y bisección. Radicandos con `max(·, 0)` en `invert_QJ_to_rp` e `init_action_angle`, y fase de la órbita circular (s₁ = s₂) fijada en `init_action_angle`, como ya hacían `analysish` y el estado `aa` (E13). (a) U6, ecuación de Kepler: de 190 de 10000 ángulos sin converger (residuo hasta 7.7·10²⁵) a ninguno (residuo ≤ 8.9·10⁻¹⁶). La parte (b) de U6 probaba una copia del Newton viejo y ahora llama a `kepler_eta`. (b) U6, ida y vuelta (Q, J, L) → (r, p) → (Q, J): fases erradas en 650, 2 y 2 de 1600 órbitas con L = 0, 0.01 y 0.5 → en ninguna; max\|ΔQ\| ≤ 1.3·10⁻¹¹ en los seis valores de L. (c) Órbita circular exacta (J = 0) con L = 0.25: antes r = 0 en los 40 ángulos probados; ahora r_c = 0.80410872. (d) Neutralidad: idénticos bit a bit (9 de las 26 configuraciones usan `aa_quad` y 5 el integrador `analytic`, que invierte cada partícula en cada paso). (e) Costo del integrador `analytic` (`gauss_an`, 128 000 partículas, 2000 pasos, orden ABBA dos veces): 31.9, 31.5 y 31.4 s antes; 31.1, 31.1, 31.4 y 31.3 s después (la primera corrida de la serie, 27.1 s con el binario anterior, se descarta: el procesador aún no estaba caliente) |

Nuevos (N) y abiertos de la auditoría anterior (D).

| Id | Severidad | Archivo | Función | Problema | Evidencia | Corrección |
|---|---|---|---|---|---|---|
| N1 | **MEDIUM** | `utils.f90:268-330`, `main.f90` | `set_timestep` | El paso fijo no resuelve el pericentro de las órbitas de L pequeño: energía destruida para 10⁻⁴ ≤ L ≤ 10⁻² (hasta 10⁶ relativo), sin aviso | A.7, figura `limite_L0.pdf` | Cota Δt ≲ η·min_j(r_p,j²/L_j) o aviso con el error de energía por partícula; o `lminc` > 0 elegido con esa cota |
| D2′ | MEDIUM | `utils.f90:242-243`, `density.f90:131` | `build_cell_list`, `deposit` | Partícula en r ≤ −1.5Δr (n = 1, 3) o −2Δr (n = 2) durante un paso: su masa no llega a la malla (hasta 100 %) | U2: 10 casos, determinista | Archivar por \|r\| |
| D3 | MEDIUM | `utils.f90:751-758` | `invert_QJ_to_rp` | Newton de Kepler diverge para e ≳ 0.98 | U6 | Arranque robusto + bisección + aviso |
| N2 | MEDIUM | `grav_force.f90:74-75,90-91,106-107,122-123,151-153`; `energy.f90:64-66`; `analysish.f90:364`; `utils.f90:793` | varias | r = 0 exacto con L = 0: L²/(2r²) y L²r/r⁴ dan 0/0 = NaN en los siete fondos; en una corrida real vuelve NaN la energía total y todos los h_k de toda la corrida, sin aviso (código de salida 0). Probabilidad baja (exige r = 0.0 exacto, típicamente desde `checkpoint`), efecto total | R0, U7 | `if (l_part(i) /= 0)` para el término centrífugo |
| N3 | LOW | `grav_force.f90:207-212,236-240` | `bgpot`, `bgforce` | `nfw`, `burkert`, `isotrun` pierden todas las cifras cerca de r = 0 (fuerza de signo invertido en r = 10⁻⁸) | A.6 (mpmath, 60 dígitos) | Serie de Taylor para a < 10⁻³ |
| D4 | LOW | `utils.f90:742-743,800-801` | `invert_QJ_to_rp`, `init_action_angle` | Radicandos sin acotar: L = 0 da r = 0 | U6 | `max(...,0)`; una sola rutina |
| N4 | INFO | `initial_data.f90` | estados de malla | Malla uniforme en L: la cuadratura no converge uniformemente en r ≲ ΔL/v_esc | Prueba M, primeros nodos | Muestrear en v_t = L/r o refinar L cerca de 0 |
| N5 | INFO | `utils.f90:377-381`, `hdf5_io.f90:163-165` | `save_data` | `vlasov_density`, `avg_rho`, `rho`, `curr` guardan r²ρ, no ρ | lectura | Documentarlo (las notas lo hacen) |
| D7 | MEDIUM | `main.md` §4 | — | El manuscrito describe RK2, otro depósito, L fijo; erratas del leapfrog | lectura | Actualizar |
| D8, D10, D11, D12, D13, D14, D15, D16 | LOW/INFO | — | — | Ver `AUDITORIA_CIENTIFICA_2026-09-21.md` | — | — |

**Resueltos desde la auditoría anterior:** D1 (autofuerza: revertido), D6 (costo),
la parte de D2 que leía memoria sin inicializar.

### Detalle de N1 (paso de tiempo y L → 0)

- **Matemática:** en el pericentro de una órbita con L pequeño, la curvatura del potencial
  efectivo es U″ = 3L²/r_p⁴ + Φ″ ≈ 3L²/r_p⁴, con r_p ≈ L/v. Un integrador explícito necesita
  Δt·√U″ ≲ 1, es decir Δt ≲ r_p²/(√3 L) ≈ L/v².
- **Implementado:** Δt = min(courant·Δr/pmax, courant·√(2Δr/F_max(t=0))).
- **Diferencia:** ninguno de los dos criterios depende de L_min ni del pericentro; el segundo se
  evalúa en t = 0, donde una partícula que empieza en su apocentro no muestra su fuerza máxima.
- **Evidencia:** tabla de A.7.
- **Impacto:** las partículas de L pequeño llevan poca masa; en las corridas del artículo
  (L ≈ 2) no hay ninguna. En una distribución con `lminc = 0` sí.

### Detalle de N2 (r = 0 exacto)

`force_part(i) = force_part(i) + l_part(i)**2*r_part(i)/den**2` con den = r² = 0 y L = 0 da
0·0/0 = NaN. El límite correcto para L = 0 es 0 (el término no existe). Para L > 0 el punto es
inalcanzable en la dinámica exacta. Prueba R0 en §E.

---

## E. Validación analítica

| Test | Predicción analítica | Resultado numérico | Error | Estado |
|---|---|---|---|---|
| A. Esfera uniforme: densidad | ρ₀ en el interior | exacta para una densidad continua (V_i = ∫W 4πr²dr); con partículas discretas queda el error de su colocación | — | PASA |
| A. Esfera uniforme: campo en la malla (r < 0.9), 64 cáscaras por Δr | g = −r | n = 1: 2.8·10⁻⁶ → 7.0·10⁻⁷ (Δr 0.04 → 0.01); n = 3: 1.1·10⁻⁶ → 2.7·10⁻⁷ | orden 1: cada subcáscara va en su punto medio, no en su centro de masa (desplazamiento O(h²/r), mayor cerca del origen) | PASA |
| A. Esfera uniforme: potencial (r < 0.9) | Φ = −(3−r²)/2 | 8.0·10⁻⁴ → 5.0·10⁻⁵ | orden 2 (lo fija el borde) | PASA |
| A. Esfera uniforme: fuerza en las partículas (r < 0.9) | g = −M(<r)/r² | 2.8·10⁻⁶ → 3.5·10⁻⁷ (Δr 0.04 → 0.005) | orden 1, sobre la corrección 5/6·Δr/r de la capa | PASA |
| B. ρ₀(1−r²)³: fuerza en las partículas | g = −M(<r)/r² | 1.3·10⁻² → 2.0·10⁻⁴ | orden 2.00 | PASA |
| B. ρ₀(1−r²)³ (U4) | Φ, g cerrados | orden 2.00 (n = 1, 2, 3) | 1.4·10⁻⁴ (n = 1, Δr = 0.005) | PASA |
| B. Plummer desde F(E) (M) | ρ, g, K, W de Plummer | orden → 2 | 1.1·10⁻² máx, 2.7·10⁻³ rms | PASA |
| U3 interpolante lineal | integral exacta (50 dígitos) | 5.5·10⁻¹⁶ / 1.1·10⁻¹⁴ | redondeo | PASA |
| C. Órbitas L > 0 frente a RK4 | trayectoria | orden 2.00 / 4.00 | 7.6·10⁻¹¹ (yoshida4, L = 1) | PASA |
| C. Órbita L = 0 por el centro frente a RK4 | trayectoria | orden 2.00 / 4.00 | 1.2·10⁻¹² | PASA |
| C. Integradores frente a `analytic` (I1) | trayectoria | 1.02 / 2.00 / 4.00 | — | PASA |
| D. L constante | L(t) = L(0) | `l_part` idéntico bit a bit en las 11 instantáneas del Plummer autogravitante (N = 154270, t = 50) | 0 | PASA |
| D. Distribución de L | invariante | idéntica (consecuencia de lo anterior) | 0 | PASA |
| Masa | Σ f L constante | `fl` idéntico bit a bit en todas las instantáneas | 0 | PASA |
| E. Energía con fondo fijo | acotada, ∝ Δt^k | leapfrog ×4, yoshida4 ×16 por mitad de Δt | 3.2·10⁻¹² | PASA |
| E. Energía, L → 0 a Δt fijo | acotada | hasta 10⁶ para 10⁻⁴ ≤ L ≤ 10⁻² | — | **FALLA (N1)** |
| F. Convergencia | ver §G | | | |
| Colapso frío (I2) | cicloide del continuo | 4.5·10⁻⁶ (N = 400), orden 1.95 | — | PASA |
| Capa aislada con masa | r̈ = L²/r³ − m/2r² | 1.9·10⁻⁵ (Δr = 0.01), orden 1 en Δr | — | PASA |
| R0. r = 0 exacto, L = 0 | Φ_iso = −½, g = 0 | NaN | — | **FALLA (N2)** |
| U2 simetría en el origen | depósito par | 10 casos pierden masa | hasta 100 % | **FALLA (D2′)** |
| U6 mapa AA | ida y vuelta | e > 0.98 y L = 0 fallan | — | **FALLA (D3, D4)** |

---

## F. Conservación

Qué debe conservarse, en qué sentido, y qué se midió (figura `docs/codigo/figuras/conservacion.pdf`):

| Magnitud | Esperado | Medido | Tipo |
|---|---|---|---|
| Masa total Σ m_j | exacta (los pesos no cambian; solo `cutoff` en t = 0 y `reduceparticles` quitan partículas) | `fl` idéntico bit a bit en todas las instantáneas | exacta |
| Masa vista por Poisson | exacta mientras toda la masa esté en r < r_N | aviso único si sale (`density.f90:46-58`) | — |
| L de cada partícula | exacta | `l_part` idéntico bit a bit | exacta |
| Distribución de L | exacta | idéntica | exacta |
| Energía, fondo fijo | acotada, ∝ Δt^k | 60 órbitas del isócrono, t = 200, Δt = 0.0125: leapfrog ≤ 2.1·10⁻⁶ (1.8·10⁻⁶ en la primera mitad, 2.1·10⁻⁶ en la segunda: sin deriva secular visible), yoshida4 ≤ 9.4·10⁻⁹ (5.4·10⁻⁹ y 9.4·10⁻⁹: compatible con acotada, este registro no lo demuestra) | discretización temporal |
| Energía por órbita, L pequeño | acotada | hasta 10⁶ relativo (N1) | paso de tiempo insuficiente |
| Energía, autogravedad | no exacta en el semidiscreto | Plummer N = 154270, t = 50: máx 1.9·10⁻⁵ (muestreo cada 0.25); N = 19356: 5.0·10⁻⁴ idéntico con 7 pasos de tiempo y 2 integradores | espacial / muestreo, no temporal |
| Energía, equilibrio de `11_equilibrio` | acotada | eq_e0 5.1·10⁻⁹, eq_e01 9.4·10⁻⁹, iso_aq 4.9·10⁻⁷ (t = 100) | — |
| Virial 2K/\|W\| | 1 en equilibrio | Plummer 0.992 → 0.997 | muestreo inicial |

**Por qué la energía no se conserva exactamente con autogravedad.** La fuerza sobre la partícula
es la interpolación del campo de la malla, no el gradiente respecto de r_j de la energía discreta
½Σ m_jΦ(r_j). El sistema semidiscreto no es hamiltoniano en esas variables; lo demostrado es que
el error no baja con Δt. La parte debida al ruido de muestreo (que baja con N: 5·10⁻⁴ → 2·10⁻⁵
entre N = 19356 y 154270, con estados iniciales de calidad distinta) no se separó de la del
esquema.

## G. Convergencia

| Variable | Experimento | Orden observado | Esperado |
|---|---|---|---|
| Δr (Poisson) | U4, 3 niveles | 2.00 | 2 |
| Δr_c, Δp_c, ΔL_c (muestreo) | Prueba M, refinamiento conjunto ×2 | 1.58, 1.80 (máx) | 2 (F(E) es C³ en E = 0) |
| Δt | I1, 6 niveles; C | 1.02 / 2.00 / 4.00 | 1 / 2 / 4 |
| N | colapso frío, N = 50 … 800 | ≈ 2 | ≥ 1 |
| Δr (capa aislada) | 4 niveles | 1.00 | 1 (corrección 5/6·Δr/r) |
| Δt con autogravedad | Plummer N = 19356, 7 corridas | 0 (piso 5.0·10⁻⁴ independiente de Δt) | piso espacial |
| Δr con autogravedad (energía) | Plummer N = 19356, Δr 0.1 … 0.0125 | sin tendencia (1.4·10⁻⁴, 5.0·10⁻⁴, 3.4·10⁻⁴, 1.7·10⁻⁴) | — |
| Δr en la fuerza sobre partículas | perfil (1−r²)³ | 2.00 | 2 |

## Simulaciones del repositorio (`reproducir/`)

Se corrieron los cuatro grupos con el binario actual (8 + 12 + 3 + 3 corridas, más
`analisis.sh`) y se compararon con los resultados del 17 de septiembre (respaldados antes):

| Grupo | Qué mide | Resultado actual | Frente al anterior |
|---|---|---|---|
| 08_verificacion | h_k frente al exacto (sin autogravedad) | código − suma discreta exacta: 10⁻¹⁴ a 10⁻¹¹ (`analytic`), hasta 3·10⁻⁶ en modos que deben anularse (yoshida4); discreta − continuo: 2·10⁻⁴ (cuadratura de los nodos) | `analytic` idéntico bit a bit; yoshida4 ≤ 3.5·10⁻¹² (último bit de `86b3c80`) |
| 09_fondos | energía en 4 fondos × 3 pasos | ×16 por mitad de Δt (orden 4) en los cuatro fondos | idéntico en las cifras impresas |
| 10_energia | energía con autogravedad | 7·10⁻⁸ a 2.5·10⁻⁷ | 1.7·10⁻⁷ a 2.7·10⁻⁷ (RK2 original) |
| 11_equilibrio | deriva del equilibrio F(J,L) | \|h₀\| 2·10⁻⁵ en t = 100; δΦ 1.2·10⁻⁷; energía 5·10⁻⁹ | δΦ 1.3·10⁻⁷, energía 5.6·10⁻⁹ |

`spi_L16_an` muestra la recurrencia por Nyquist en L (diferencia con el continuo hasta 10³ en
k = 4), documentada en las notas; no es un defecto del código. En `11_equilibrio`, el
\|h₁\|/\|h₀\| de 5·10⁻³ que da el diagnóstico en línea es el artefacto del mapa del isócrono en un
potencial autogravitante; con el mapa verdadero (`tools/hk_numerico.py`) es 1.6·10⁻⁵.

No se buscó ni se encontró calentamiento, amortiguamiento o crecimiento artificiales en estas
corridas más allá de lo que explican el error de cuadratura y el paso de tiempo. No hay
violación de simetría (la dinámica es radial por construcción y la paridad se verificó bit a
bit). Los artefactos cerca de r = 0 son los de A.6 y A.7.

## Pruebas automatizadas

`verificacion/correr.sh` suma tres pruebas nuevas:

| Prueba | Archivo | Detecta | Estado |
|---|---|---|---|
| U7 | `src/t_r0.f90` | N2 (r = 0, L = 0) | FALLA |
| I3 | `py/medida.py` | factor de la medida (Plummer) | PASA (umbral fijado tras la primera corrida, justificado en el script) |
| I4a | `py/limite_L0.py` | órbita L = 0 por el centro | PASA |
| I4b | `py/limite_L0.py` | N1 (energía con L = 10⁻³) | FALLA |

Estado completo con el binario actual: pasan U1, U3 (×2), U4, U5, I1, I2, I3, I4a; fallan U2
(D2′), U6 (D3, D4), U7 (N2), I4b (N1).

## H. Conclusión

**Matemáticamente demostrado.** El modelo que resuelve el código es Vlasov–Poisson esférico con
F(r, p_r, L) respecto de la medida 8π² L dr dp_r dL; las características ṙ = p_r,
ṗ_r = L²/r³ − ∂Φ/∂r, L̇ = 0 están implementadas con el signo, el factor y la dependencia correctos;
la medida es invariante y por eso los pesos son constantes; la forma de Poisson por masa
encerrada es equivalente a la de segundo orden con regularidad en el origen; el primer nodo
reproduce el límite regular Φ = Φ₀ + (2π/3)ρ₀r²; la reflexión en el origen es la simetría exacta
de la órbita radial.

**Validado numéricamente.** La cadena F → partículas → densidad → Poisson → fuerza → energía,
con un sistema autogravitante de solución conocida (Plummer: densidad, campo, K y W convergen con
orden ≈ 2; colapso frío con orden 2 en 1/N). La dinámica de prueba con L > 0 y L = 0 frente a
integraciones independientes (orden 2 y 4). La conservación exacta de L y de la masa. La energía
acotada con fondo fijo. Los resultados de `reproducir/` frente a sus soluciones exactas y frente
a la versión anterior.

**Consistente pero no probado.** La evolución larga con autogravedad en el régimen del artículo
(no hay solución de referencia); la convergencia de h_k con autogravedad; el origen del piso de
energía con autogravedad.

**Sospechoso.** Cualquier resultado que dependa de las órbitas de L pequeño con el paso por
omisión (N1): su energía no se conserva y el código no lo avisa.

**Incorrecto.** N2 (r = 0 con L = 0 da NaN y contamina los diagnósticos globales); N3 (fondos
`nfw`, `burkert`, `isotrun` cerca de r = 0); D2′ (masa perdida para r ≤ −1.5Δr); D3 (Newton de
Kepler); D4 (mapa AA con L = 0 en `utils.f90`). Ninguno alcanza las corridas del artículo.

**Recomendación.** En este orden: (1) N1, un criterio de paso que vea el pericentro, o al menos
un aviso con el error de energía por partícula; (2) N2 y D2′, que son de pocas líneas (término
centrífugo solo si L ≠ 0; lista de celdas por |r|); (3) D3 y D4 con una sola rutina del mapa
ángulo-acción; (4) N3 con series de Taylor. Después, volver a correr `verificacion/correr.sh`:
las cuatro fallas deberían pasar sin tocar las pruebas.
