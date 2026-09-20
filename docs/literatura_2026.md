# Estado del arte y temas de investigación con este código

Búsqueda bibliográfica del 2026-09-20, hecha para decidir qué se puede publicar con
`VlasovPoisson_PIC_sp` (Vlasov–Poisson esférico, perturbaciones ℓ=0, distribución en el
momento angular L, PIC de cuadratura en las variables ángulo–acción).

Todo lo que sigue sobre trabajos ajenos viene de resúmenes y páginas de acceso abierto,
no de una lectura completa. Antes de afirmar novedad en un artículo hay que leer enteros,
como mínimo, Polyachenko y Shukhman (2026), Straub (2024), Fouvry y Prunet (2022) y
Chiba et al. (2025).

## 1. Qué hay publicado

### 1.1 Marco general

- **Hamilton y Fouvry (2024)**, *Kinetic theory of stellar systems: a tutorial*,
  Phys. Plasmas 31, 120901, [arXiv:2402.13322](https://arxiv.org/abs/2402.13322).
  Referencia moderna para el vocabulario común con plasmas: phase mixing como
  cizallamiento en los ángulos, respuesta lineal, amortiguamiento de Landau, relajación
  a dos cuerpos. Es el marco en el que conviene escribir cualquier resultado nuestro.

### 1.2 Phase mixing riguroso

- **Chaturvedi y Luk (2026)**, *Linear and nonlinear phase mixing for the gravitational
  Vlasov–Poisson system under an external Kepler potential*, Arch. Ration. Mech. Anal.,
  [arXiv:2409.14626](https://arxiv.org/abs/2409.14626). Demuestra phase mixing lineal y
  no lineal en simetría esférica con un Kepler externo, y para ello construye variables
  ángulo–acción *definidas dinámicamente* a partir del potencial autoconsistente. Es el
  análogo riguroso de lo que aquí se hace numéricamente con el mapa numérico.
- **Taylor y Velozo Ruiz (2025)**, *Phase mixing and the Vlasov equation in cosmology*,
  [arXiv:2512.04214](https://arxiv.org/abs/2512.04214). Contexto: tasas de mezcla en
  espacios FLRW; sin autogravedad esférica.
- **Nonlinear phase mixing in Einstein–Vlasov near Schwarzschild** (2026),
  [arXiv:2609.13394](https://arxiv.org/abs/2609.13394). Muestra que el tema está activo
  también en el lado relativista.

### 1.3 Modos y amortiguamiento en sistemas esféricos

- **Weinberg (1994)**, ApJ 421, 481: modos débilmente amortiguados en sistemas esféricos
  estables.
- **Fouvry y Prunet (2022)**, MNRAS 509, 2443,
  [arXiv:2105.01371](https://arxiv.org/abs/2105.01371). Método matricial con la
  continuación analítica necesaria para los modos amortiguados; encuentran un modo ℓ=1
  débilmente amortiguado en el isócrono isótropo y lo contrastan con N-cuerpos.
- **Polyachenko, Shukhman y Borodina (2021)**, MNRAS 503, 660,
  [arXiv:2101.08287](https://arxiv.org/abs/2101.08287): distinguen modos genuinos de
  ondas amortiguadas por Landau y advierten que las desviaciones del amortiguamiento de
  Landau son comunes.
- **Modo de dilatación / modo invariante de escala**, [arXiv:2311.05551](https://arxiv.org/abs/2311.05551)
  y A&A 684 (2024), *Exploring the dynamics of collisionless spherical stellar systems
  using the matrix method: insights from the dilation mode*: el caso ℓ=0 tiene una
  solución exacta asociada a la invariancia de escala, útil como control.
- **Polyachenko y Shukhman (2026)**, *Two sets of potential-density basis pairs for the
  study of radial perturbations in collisionless spherical stellar systems*,
  [arXiv:2609.04012](https://arxiv.org/abs/2609.04012) (3 de septiembre de 2026).
  Construyen bases nuevas precisamente porque las estándar convergen mal en ℓ=0. Es el
  estado del arte teórico del problema que este código ataca por simulación, y el
  interlocutor natural de nuestros resultados.

### 1.4 Simulación del sistema esférico

- **Straub (2024)**, *Numerical experiments on stationary, oscillating, and damped
  spherical galaxy models*, Physica D 470, 134351,
  [arXiv:2405.01235](https://arxiv.org/abs/2405.01235); código `radVP` en
  [GitHub](https://github.com/c-straub/radVP). PIC radial del sistema esférico
  (linealizado y completo) alrededor de equilibrios de soporte compacto: encuentra
  oscilaciones parcialmente no amortiguadas en unos modelos y amortiguamiento
  macroscópico en otros, y relaciona el comportamiento con la función de periodo radial.
  Es el trabajo más cercano al nuestro en método y en régimen.
- **Ramming y Rein (2018)**, Physica D 365, 72: oscilaciones numéricas alrededor de
  equilibrios estables.
- **Steady states como puntos fijos de un algoritmo que preserva la masa** (2024),
  [arXiv:2412.01544](https://arxiv.org/abs/2412.01544): otra forma de construir los
  equilibrios que aquí se construyen con `tools/equilibrio_L.py`.

### 1.5 Ecos

- **Chiba, Ding, Hamilton, Kunz y Tremaine (2025)**, *Galactic echoes*,
  [arXiv:2506.16512](https://arxiv.org/abs/2506.16512). Deducen la teoría de ecos de
  plasma en variables ángulo–acción y la aplican a un modelo **unidimensional** del
  movimiento vertical de la Vía Láctea, con partículas de prueba y con difusión por
  dispersión con nubes moleculares. El eco se desenrolla y se vuelve a enrollar; la
  difusión lo amortigua de forma súper-exponencial.
- En plasmas, la supresión de ecos por colisiones está demostrada rigurosamente
  (Bedrossian, Ann. PDE 2017).

### 1.6 Ruido de discreción y anisotropía

- **Formalismo de partícula vestida**: Fouvry, Chavanis y Pichon (2017),
  [arXiv:1706.06009](https://arxiv.org/abs/1706.06009); Weinberg (1998),
  [astro-ph/9707206](https://arxiv.org/abs/astro-ph/9707206). Predicen la amplificación
  gravitatoria del ruido de Poisson (factores ~6–15 según el modelo).
- **Inestabilidad de órbitas radiales (ROI)**: MNRAS 451, 601 (2015),
  [arXiv:1504.03513](https://arxiv.org/abs/1504.03513); con ruido,
  [arXiv:2110.11026](https://arxiv.org/abs/2110.11026). Relevante como límite: la ROI es
  ℓ≥2 y este código, al ser esférico, no puede verla.

## 2. Dónde está el hueco

Lo que no encontré hecho, y que este código ya puede hacer:

1. **Ecos y desenrollado en un sistema esférico con dispersión en L.** Chiba et al.
   trabajan en 1D vertical; aquí la frecuencia depende de dos acciones, ω(J,L), y la
   dispersión en L no está "preparada" por la perturbación. Medimos que borra el
   desenrollado de una espiral: |h₁| pasa de crecer 2.2 veces (L fija) a 2.2e-3 con
   σ_L=0.2. Es un mecanismo de supresión *sin colisiones ni difusión*, distinto del de
   ese trabajo.
2. **Dependencia del amortiguamiento ℓ=0 con la anisotropía**, medida en un montaje
   controlado: con este código σ_L fija la anisotropía del equilibrio F(J,L) sin cambiar
   la masa ni el perfil de densidad de forma esencial. La literatura ℓ=0 es casi toda
   teórica (bases, continuación analítica) o numérica con modelos isótropos.
3. **Criterio cuantitativo de cuándo sirve el mapa ángulo–acción analítico** en presencia
   de autogravedad. Que las acciones calculadas en un potencial aproximado oscilen es
   conocido (Sanders y Binney 2016; Burger et al. 2021), pero la consecuencia medida aquí
   —un residuo *estático* en h_k proporcional a la masa, que imita un modo que no decae—
   no la encontré descrita.
4. **Validación cruzada entre dos códigos independientes** del mismo sistema: `radVP` de
   Straub y este. Straub publica código y modelos; sus equilibrios entran aquí por
   `state=checkpoint`.

## 3. Temas concretos, en orden de madurez

Cada uno indica el montaje, lo que falta implementar y el costo con la máquina actual
(≈39 ms por paso con 128 mil partículas y autogravedad).

### T1. Supresión de ecos por dispersión en momento angular
- **Qué**: amplitud del desenrollado frente a σ_L y a la forma de C(L); después, ecos
  propiamente dichos (dos perturbaciones en t₁ y t₂) con y sin autogravedad.
- **Cómo**: ya existe `dftype=spiral`, el integrador `analytic` y la solución exacta
  (`tools/hk_exacto.py`). Falta un estado con dos perturbaciones y la predicción lineal
  de la amplitud del eco, exp[−½k²(∂ω/∂L)²var_L t*²].
- **Costo**: horas. **Riesgo**: bajo; el resultado principal ya está medido.
- **Interlocutor**: Chiba et al. (2025).

### T2. Amortiguamiento ℓ=0 en equilibrios F(J,L) frente a la anisotropía
- **Qué**: ω y γ de la cola colectiva en función de σ_L (y de a₀), con σ_L→0 reproduciendo
  el valor de L fija (ω=0.05705, γ=5.08e-3 en el otro código).
- **Cómo**: `tools/equilibrio_L.py` + `state=checkpoint`, corridas con ε=0 y ε=0.1 hasta
  t≈3000, ajuste con matrix pencil, controles de N_J, N_L y Δt.
- **Costo**: ~10 corridas de 1–2 h. **Riesgo**: medio; depende de T4 (piso de ruido).
- **Interlocutor**: Polyachenko y Shukhman (2026), Fouvry y Prunet (2022), Straub (2024).

### T3. Teoría lineal con L, sin partículas
- **Qué**: resolver la ecuación linealizada en (Q,J,L) con transporte exacto en Q y la
  función de Green de capas esféricas, como se hizo con L fija, y comparar con T2.
- **Cómo**: portar `lineal.py` del código de L fija agregando la rejilla en L.
- **Costo**: días de programación, minutos de cómputo. **Riesgo**: bajo.
- **Valor**: es el control que vuelve creíble cualquier medición de T2, y permite la
  continuación analítica para decidir si la cola es un polo.

### T4. Piso de ruido de discreción con autogravedad
- **Qué**: cómo crece la componente que no cancela, frente a N_J, N_Q, N_L y Δr, y
  comparación con la amplificación predicha por el formalismo de partícula vestida.
- **Cómo**: repetir una corrida con resoluciones dobles; el diagnóstico ya existe.
- **Costo**: días. **Riesgo**: bajo, pero es condición para T2.
- **Interlocutor**: Fouvry et al. (2017), Weinberg (1998).

### T5. Cuándo sirve el mapa analítico (artículo de métodos)
- **Qué**: la meseta artificial frente a a₀ y σ_L, y un criterio "el mapa del isócrono
  sirve para señales mayores que α a₀ h₀".
- **Cómo**: sale casi gratis de las corridas de T2 con `tools/hk_numerico.py`.
- **Costo**: horas. **Riesgo**: bajo.

### T6. Validación cruzada con `radVP` (Straub)
- **Qué**: reproducir uno o dos de sus modelos (King, politropos, cáscaras) y comparar la
  frecuencia y el amortiguamiento de las oscilaciones ℓ=0.
- **Cómo**: construir el equilibrio con su prescripción, entrar por `checkpoint`.
- **Costo**: días. **Riesgo**: medio (hay que igualar unidades y condiciones de frontera).
- **Valor**: dos códigos independientes que coinciden pesan más que cualquier prueba
  interna, y conecta directamente con literatura reciente.

### T7. Régimen no lineal
- **Qué**: a₀ ≳ 10⁻², donde el efecto colectivo deja de ser corrección, hasta relajación
  violenta.
- **Requisitos**: T2, T3 y T4 primero.

## 4. Lo que el código no puede hoy

- **La perturbación es solo ℓ=0**, que es independiente de tener distribución en L.
  Cada partícula lleva su propio L y es una cáscara esféricamente simétrica de masa
  8 π² f L dL dr dp: las órbitas son rosetas y σ_L fija la anisotropía del equilibrio,
  pero la densidad y el campo que ve Poisson (dΦ/dr = M/r² en una malla radial) solo
  tienen monopolo. En consecuencia la perturbación depende únicamente del ángulo radial
  y la condición de resonancia es n Ω_r = ω con n_ψ = 0 — pero Ω_r = Ω_r(J,L), así que
  la dispersión en L ensancha la banda de resonancia. Ese es el mecanismo que se mide
  aquí. Lo que queda fuera son las perturbaciones no esféricas: la ROI (ℓ=2), los modos
  ℓ=1 de Fouvry y Prunet y las espirales de fase tipo Gaia.
- **Un solo fondo en producción** (isócrono); los demás fondos están implementados y
  verificados, pero los diagnósticos ángulo–acción suponen el isócrono.
- **El mapa verdadero solo en postproceso** (`tools/hk_numerico.py`).
- **El paso de tiempo se elige a mano**, con prueba de convergencia: ningún criterio
  automático del código conoce la frecuencia orbital.
- **Pendiente antes de publicar números de L fija**: el código `vlasov-poisson_PIC`
  comparte las fallas de acoplamiento partícula–malla corregidas aquí; hay que portarlas y
  repetir su corrida de referencia y la de Landau (ver `BUGS_TODO.md`).
