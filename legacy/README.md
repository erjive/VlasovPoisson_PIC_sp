# Código legado

Archivos que estaban en `src/` pero no se compilan: `analysish`, `poisson`,
`poisson_ps` y `reduce_arrays` (sin extensión `.f90`), y `analysish_old.f90`, que ya
no compila (usa variables que no existen). `poisson_ps` requiere un módulo
`chebyshev` que no está en el repositorio. Se conservan como referencia histórica.

`initial_data_estados_2D.f90` guarda las ramas `Plummer`, `compact` y `compact2` de
`initial_data.f90` (y las de `checkpoint`/`other3`, vacías). Venían del código con L fija:
llenaban solo `Nrc*Npc` partículas y nunca asignaban `l_part`, así que con la distribución
en L las partículas no tenían barrera centrífuga y la densidad, la energía y h_k daban
cero. `Plummer` además usaba la energía del isócrono. Para recuperarlos hay que
decidir su dependencia en L.

