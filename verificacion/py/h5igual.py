# Compara dos vlasov_output.h5 bit a bit (datasets y atributos).
# Uso: h5igual.py a.h5 b.h5 ; codigo de salida 0 si son identicos.
import h5py, numpy as np, sys
a = h5py.File(sys.argv[1], 'r'); b = h5py.File(sys.argv[2], 'r')
cnt = {'datasets': 0, 'atributos': 0, 'distintos': 0}; peor = [0.0]
def walk(ga, gb):
    for k in ga:
        if k not in gb:
            cnt['distintos'] += 1; continue
        if isinstance(ga[k], h5py.Group):
            for at in ga[k].attrs:
                cnt['atributos'] += 1
                x = np.asarray(ga[k].attrs[at]); y = np.asarray(gb[k].attrs[at])
                if x.tobytes() != y.tobytes():
                    cnt['distintos'] += 1; peor[0] = max(peor[0], float(np.abs(x-y).max()/max(np.abs(x).max(), 1e-300)))
            walk(ga[k], gb[k])
        else:
            cnt['datasets'] += 1; x = ga[k][:]; y = gb[k][:]
            if x.shape != y.shape or x.tobytes() != y.tobytes():
                cnt['distintos'] += 1
                if x.shape == y.shape:
                    peor[0] = max(peor[0], float(np.abs(x-y).max()/max(np.abs(x).max(), 1e-300)))
walk(a, b)
print('  datasets=%d atributos=%d distintos=%d max_rel=%.3e' % (cnt['datasets'], cnt['atributos'], cnt['distintos'], peor[0]))
sys.exit(0 if cnt['distintos'] == 0 and cnt['datasets'] > 0 else 1)
