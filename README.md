# VlasovPoisson_PIC_sp
This code solves the Vlasov-Poisson system in spherical symmetry

## Compilar y correr

```bash
make                                   # deja exe/VP_PIC
cd exe
./VP_PIC input_parameters              # o sin argumento: lee input_parameters
./VP_PIC input_parameters Nt=2000 directory=prueba/a   # overrides nombre=valor
./VP_PIC --help                        # lista de parámetros
```

El archivo de parámetros tiene líneas `nombre = valor` (sin importar el orden ni las
mayúsculas; `!` o `#` empiezan un comentario). Lo que no aparece conserva el valor de
omisión de `src/parameters.f90`. Cada corrida escribe en su directorio
`params_usados.par`, la configuración completa ya con los overrides, que sirve para
repetirla. Los archivos del formato posicional anterior (`VP_PIC < archivo`) se
convierten con `python3 tools/posicional_a_par.py viejo nuevo`.

## Salidas de h_k

`analysish` proyecta la distribución sobre dos funciones de prueba
Φ_n = exp(-sin²(Q/2)/sq_n²) · J² exp(-(J-j_n)²/sj_n²) · exp(-(L-lt_n)²/slt_n²)
cada `spatial_output` pasos: `hk1.tl`/`hk2.tl` (t, |h_0| … |h_4|) y
`hk1_complex.tl`/`hk2_complex.tl` (t, Re h_0, Im h_0, …). Las instantáneas de
partículas y campos se escriben cada `field_output` pasos (múltiplo de
`spatial_output`; por omisión, igual).
