# I5. La salida raw contiene lo mismo que la HDF5. Se corre el mismo caso
# (autogravedad, reduceparticles activo para que cambie el numero de
# particulas) con output_format = hdf5 y raw, y se comparan bit a bit la
# malla, los arreglos de malla y de particulas y los escalares de cada
# instantanea, leyendo el raw con tools/raw_io.py. Tambien se convierte el
# raw a HDF5 con la misma herramienta y se compara con el HDF5 del codigo.
# Uso: salida_raw.py <VP_PIC> <dir_trabajo>
import numpy as np, h5py, subprocess, sys, os
exe, wd = sys.argv[1], sys.argv[2]
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools'))
from raw_io import CorridaRaw, a_hdf5, ESCALARES
os.makedirs(wd, exist_ok=True); os.chdir(wd)
plantilla = """dr = 0.1
Nrc = 30
Npc = 30
Nlc = 6
courant = 0.5
Nt = 60
rmin = 0.0
rmax = 6.0
rminc = 0.5
rmaxc = 10.0
pminc = -0.8
pmaxc = 0.8
lminc = 0.0
lmaxc = 2.4
pmax = 2.0
reduceparticles = .true.
Nreduce = 20
time_output = 100000
spatial_output = 5
field_output = 10
a0 = 0.1
l0 = 1.0
sr = 0.1
sp = 0.3
sl = 0.8
state = aa
cutoff = 0.001
bsplineorder = 2
integrator = yoshida4
BGtype = Isochrone
autointeraction = {sg}
output_format = {fmt}
directory = {d}
"""
malos = []; resumen = []
for sg in ('.true.', '.false.'):
    for fmt in ('hdf5', 'raw'):
        d = '%s_%s' % (fmt, 'sg' if sg == '.true.' else 'libre')
        open(d+'.par', 'w').write(plantilla.format(sg=sg, fmt=fmt, d=d))
        subprocess.run([exe, d+'.par'], stdout=open(d+'.log', 'w'), stderr=subprocess.STDOUT, check=True)
    cola = 'sg' if sg == '.true.' else 'libre'
    h = h5py.File('hdf5_%s/vlasov_output.h5' % cola, 'r'); c = CorridaRaw('raw_%s/vlasov_output.raw' % cola)
    st = sorted(k for k in h if k.startswith('step_'))
    n = 0
    if len(st) != len(c): malos.append('%s: %d instantaneas en hdf5, %d en raw' % (cola, len(st), len(c)))
    if h['grid/r'][:].tobytes() != c.r.tobytes(): malos.append('%s: malla' % cola)
    for i, k in enumerate(st[:len(c)]):
        s = c[i]
        if int(k[5:]) != s['l']: malos.append('%s %s: paso' % (cola, k))
        for a in ESCALARES:
            n += 1
            if np.asarray(h[k].attrs[a]).tobytes() != np.float64(s[a]).tobytes(): malos.append('%s %s: %s' % (cola, k, a))
        for dset in h[k]:
            n += 1
            if dset not in s or h[k][dset][:].tobytes() != s[dset].tobytes(): malos.append('%s %s: %s' % (cola, k, dset))
        if set(h[k]) != set(x for x in s if isinstance(s[x], np.ndarray)): malos.append('%s %s: conjuntos distintos' % (cola, k))
    npart = [c[i]['Npart'] for i in range(len(c))]
    # conversion a HDF5
    g = h5py.File(a_hdf5('raw_%s/vlasov_output.raw' % cola, 'conv_%s.h5' % cola), 'r')
    for k in st:
        for dset in h[k]:
            if g[k][dset][:].tobytes() != h[k][dset][:].tobytes(): malos.append('%s %s: %s (convertido)' % (cola, k, dset))
        for a in ESCALARES:
            if np.asarray(g[k].attrs[a]).tobytes() != np.asarray(h[k].attrs[a]).tobytes(): malos.append('%s %s: %s (convertido)' % (cola, k, a))
    tam = (os.path.getsize('hdf5_%s/vlasov_output.h5' % cola), os.path.getsize('raw_%s/vlasov_output.raw' % cola))
    print('  %-5s %d instantaneas, %d arreglos y escalares comparados, particulas %d -> %d, tamano hdf5 %d B, raw %d B' % (cola, len(c), n, npart[0], npart[-1], tam[0], tam[1]))
for m in malos[:10]: print('   distinto:', m)
print(('PASA' if not malos else 'FALLA')+' I5 salida raw igual a la HDF5')
