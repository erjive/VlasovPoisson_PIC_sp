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
