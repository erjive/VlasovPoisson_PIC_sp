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

## Herramientas (`tools/`)

| herramienta | para qué |
|---|---|
| `posicional_a_par.py` | convierte archivos de parámetros del formato posicional anterior |
| `hk_exacto.py` | h_k exacto de `state=aa_quad` sin autogravedad (continuo y suma discreta) frente a la corrida |
| `aa_numerico_L.py` | variables ángulo-acción numéricas con L por partícula en isócrono + tabla (directo e inverso) |
| `hk_numerico.py` | h_k en las variables verdaderas desde instantáneas HDF5 (necesario con autogravedad) |
| `equilibrio_L.py` | equilibrio autoconsistente F(J,L) y condición inicial para `state=checkpoint` |
| `delta_phi.py` | δΦ(r,t) desde HDF5, con resta opcional de una corrida de referencia |

Requieren Python 3 con `numpy` (y `h5py` las que leen HDF5).

## Video 3D de la DF gaussiana

```bash
cd exe && ./VP_PIC ../reproducir/video/df_gauss.par && cd ..   # 43 s, instantáneas HDF5 (0.7 GB)
python3 reproducir/video/video_df.py --valida                  # comprueba el ángulo azimutal
python3 reproducir/video/video_df.py --procesos 6              # 601 cuadros 1920x1080 y exe/rep/video/df_gauss_3d.mp4
```

Muestra la misma distribución en $(r,p_r,L)$, $(Q,J,L)$ y $(x,y,z)$ hasta t=2400 (unos
23 periodos radiales). En $(x,y,z)$ cada partícula (una capa esférica) se dibuja con
6 estrellas en planos orbitales aleatorios; su ángulo en el plano se obtiene de
dψ/dt = L/r² (validado contra integración directa a 3e-5 rad).

La versión con autogravedad (a0=1e-2, yoshida4 con dt=0.025) es:

```bash
cd exe && ./VP_PIC ../reproducir/video/df_gauss_sg.par && cd ..   # 25 min
python3 reproducir/video/video_df_sg.py --valida                  # error de ψ frente al caso exacto
python3 reproducir/video/video_df_sg.py --aa --procesos 2         # (Q,J) del potencial real, 601 instantáneas
python3 reproducir/video/video_df_sg.py --procesos 2              # exe/rep/video/df_gauss_sg_3d.mp4
```

Aquí $(Q,J)$ sale del mapa ángulo-acción numérico en el potencial de cada instante (con
el del isócrono aparecería el artefacto de coordenadas de las notas) y ψ se integra entre
instantáneas con Hermite cúbico de r(t) y Gauss-Legendre de 8 nodos: el método, probado en
la corrida sin autogravedad donde ψ se conoce exacto, acumula 1.2e-2 rad en 24 órbitas.
Con `--procesos 2` el render usa poca memoria y deja la máquina utilizable.
