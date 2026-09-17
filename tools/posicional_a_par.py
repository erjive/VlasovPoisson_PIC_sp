"""Convierte un archivo de parámetros posicional (el formato anterior, un valor por
línea leído con "VP_PIC < archivo") al formato "nombre = valor" de paramfile.f90.

Uso:  python3 tools/posicional_a_par.py viejo [nuevo]      (sin "nuevo": a la salida estándar)

El orden es el que leía read_initial_param. Se toma el primer campo de cada línea,
como hacía la lectura dirigida por lista; si la línea sigue con un comentario que
empieza con "!", se conserva al final de la línea nueva.
"""
import sys

NOMBRES = ['dr', 'Nrc', 'Npc', 'Nlc', 'courant', 'Nt', 'rmin', 'rmax', 'rminc', 'rmaxc',
           'pminc', 'pmaxc', 'lminc', 'lmaxc', 'reduceparticles', 'Nreduce',
           'time_output', 'spatial_output', 'directory', 'a0', 'r0', 'p0', 'l0',
           'sr', 'sp', 'sl', 'state', 'cutoff', 'bsplineorder', 'integrator',
           'spatialorder', 'forcetype', 'BGtype', 'autointeraction', 'output_format']


def convertir(texto):
    lineas = [l for l in texto.splitlines() if l.strip()]
    if len(lineas) < len(NOMBRES):
        sys.exit(f'se esperaban {len(NOMBRES)} valores, hay {len(lineas)}')
    salida = []
    for n, l in zip(NOMBRES, lineas):
        v = l.split()[0]
        c = l[l.index('!'):].strip() if '!' in l[len(v):] else ''
        salida.append(f'{n:<16} = {v:<10} {c}'.rstrip() + '\n')
    return ''.join(salida)


if __name__ == '__main__':
    if len(sys.argv) not in (2, 3):
        sys.exit(__doc__)
    salida = convertir(open(sys.argv[1]).read())
    if len(sys.argv) == 3:
        open(sys.argv[2], 'w').write(salida)
    else:
        sys.stdout.write(salida)
