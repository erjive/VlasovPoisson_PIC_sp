"""Lector de la salida binaria de VP_PIC (output_format = raw).

El archivo <directory>/vlasov_output.raw no se describe a sí mismo: su formato está en
el comentario inicial de src/raw_io.f90 y aquí, y los dos deben cambiar juntos. Todo son
enteros y reales de 8 bytes, en el orden de bytes de la máquina que lo escribió.

    cabecera (una vez):  Nr:int64, autointeraction:int64 (0/1), r:float64[Nr]
    registro (por instantánea):
        l:int64, time, kinetic_energy, potential_energy, total_energy:float64, Npart:int64,
        rho, avg_rho, curr:float64[Nr]            (cada una multiplicada por r^2)
        force, potential:float64[Nr]              (solo con autointeraction)
        r_part, p_part, fl, l_part:float64[Npart] (fl = l_part*f)

Los registros no tienen tamaño fijo (Npart puede bajar con reduceparticles), así que el
archivo se recorre una vez leyendo solo el Npart de cada registro. Los arreglos se leen
al pedirlos.

Uso como módulo:
    from raw_io import CorridaRaw
    c = CorridaRaw('corrida/vlasov_output.raw')
    c.r, c.tiempos, len(c)
    s = c[3]            # instantánea 3: s['r_part'], s['rho'], ..., s['time'], s['l']

Uso desde la línea de comandos:
    python3 tools/raw_io.py corrida/vlasov_output.raw             # resumen
    python3 tools/raw_io.py corrida/vlasov_output.raw salida.h5   # convierte a HDF5

El HDF5 convertido tiene los grupos, conjuntos y atributos que escribe el código con
output_format = hdf5, así que sirve de entrada a las herramientas que leen HDF5.
"""
import sys
import numpy as np

ESCALARES = ('time', 'kinetic_energy', 'potential_energy', 'total_energy')
PARTICULAS = ('r_part', 'p_part', 'fl', 'l_part')


class CorridaRaw:
    def __init__(self, ruta):
        self.ruta = ruta
        self._f = open(ruta, 'rb')
        self.Nr = int(np.fromfile(self._f, np.int64, 1)[0])
        self.autointeraction = bool(np.fromfile(self._f, np.int64, 1)[0])
        self.r = np.fromfile(self._f, np.float64, self.Nr)
        self.malla = ('rho', 'avg_rho', 'curr') + (('force', 'potential') if self.autointeraction else ())
        self._indice = []
        while True:
            inicio = self._f.tell()
            cab = np.fromfile(self._f, np.int64, 1)
            if cab.size == 0:
                break
            esc = np.fromfile(self._f, np.float64, 4)
            npart = np.fromfile(self._f, np.int64, 1)
            resto = 8*(self.Nr*len(self.malla) + (int(npart[0]) if npart.size else 0)*len(PARTICULAS))
            fin = self._f.seek(resto, 1)
            if esc.size < 4 or npart.size == 0 or fin > self._tamano():
                break                               # registro incompleto (corrida interrumpida)
            e = dict(zip(ESCALARES, map(float, esc)))
            e.update(posicion=inicio, l=int(cab[0]), Npart=int(npart[0]))
            self._indice.append(e)

    def _tamano(self):
        pos = self._f.tell()
        fin = self._f.seek(0, 2)
        self._f.seek(pos)
        return fin

    def __len__(self):
        return len(self._indice)

    @property
    def tiempos(self):
        return np.array([e['time'] for e in self._indice])

    @property
    def pasos(self):
        return np.array([e['l'] for e in self._indice])

    def __getitem__(self, i):
        e = self._indice[i]
        self._f.seek(e['posicion'] + 8*6)
        s = {k: e[k] for k in ESCALARES + ('l', 'Npart')}
        for k in self.malla:
            s[k] = np.fromfile(self._f, np.float64, self.Nr)
        for k in PARTICULAS:
            s[k] = np.fromfile(self._f, np.float64, e['Npart'])
        return s

    def cerrar(self):
        self._f.close()

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.cerrar()


def a_hdf5(ruta_raw, ruta_h5):
    """Escribe la corrida en el formato de output_format = hdf5."""
    import h5py
    with CorridaRaw(ruta_raw) as c, h5py.File(ruta_h5, 'w') as h:
        h.create_group('grid').create_dataset('r', data=c.r)
        for i in range(len(c)):
            s = c[i]
            g = h.create_group('step_%010d' % s['l'])
            for k in ESCALARES:
                g.attrs.create(k, np.array([s[k]]))
            for k in c.malla + PARTICULAS:
                g.create_dataset(k, data=s[k])
    return ruta_h5


if __name__ == '__main__':
    if len(sys.argv) not in (2, 3):
        sys.exit(__doc__)
    if len(sys.argv) == 3:
        print('escrito', a_hdf5(sys.argv[1], sys.argv[2]))
    else:
        with CorridaRaw(sys.argv[1]) as c:
            print(f'{c.ruta}: Nr = {c.Nr}, autointeraction = {c.autointeraction}, {len(c)} instantáneas')
            if len(c):
                t = c.tiempos
                print(f'  pasos {c.pasos[0]} a {c.pasos[-1]}, t = {t[0]:g} a {t[-1]:g},'
                      f' partículas: {c._indice[0]["Npart"]} a {c._indice[-1]["Npart"]}')
