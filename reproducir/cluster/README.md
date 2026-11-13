# Running the cloud with distributed angular momenta on the cluster

Scripts to build this code on a Slurm cluster and to run and analyze there the series of
`reproducir/scripts/nube_L.py`, the gas cloud far from equilibrium with a spread of the angular
momentum, together with the reduced model of the same clouds. On the laptop the fourteen PIC
runs and their analysis take 7.7 hours and the reduced model 4.7 hours, one after the other.
On a node with 24 cores they run side by side.

The environment (modules, compiler, number of physical cores) is the one of
`vlasov-poisson_PIC`, file `reproducir/cluster/entorno.sh`, whose README describes the cluster
and its queue. That repository is needed in any case: `nube_L.py` takes the particles and the
measures from it, and the reduced model runs there. It must lie next to this one, or its path
must be in the variable `VP_PIC_FIJO`.

| file | purpose |
|---|---|
| `todo.sh` | everything at once: build and test, initial data, and all the jobs; `-t` sends only the test job |
| `prueba.slurm` | a test of ten minutes on a compute node: the node, the executable, the whole chain on a short run, the reduced model, and the cost of each run |
| `comun.sh` | locates the two repositories and loads the environment; sourced by the others |
| `compilar.sh` | builds `exe/VP_PIC` and runs `probar.sh` |
| `probar.sh` | a run of a few seconds compared with the laptop, the Python modules, and whether the initial data are there |
| `nube.slurm` | runs, analyzes and measures the runs given as arguments |
| `enviar.sh` | submits one job per run: `enviar.sh primera`, `enviar.sh resto convergencia`, or run names |
| `estado.sh` | the jobs in the queue, how far each run got, the progress of the reduced model, and the table of measures |
| `traer.sh` | on the laptop: fetches the small files of the finished runs, and with `-r` those of the reduced model |

All commands are given from the root of the repository.

## 0. Check before running anything

In November 2026 the cluster moved from CentOS 7 to AlmaLinux 9, with a new login node. The
modules, the partition and the paths that `entorno.sh` of `vlasov-poisson_PIC` assumes are
those of October 2026, so check in three rounds, from cheap to less cheap:

1. On the login node, nothing is submitted. Clone both repositories into a new folder (section
   1), then:

   ```bash
   bash ../vlasov-poisson_PIC/reproducir/cluster/diagnostico.sh > diagnostico.txt 2>&1
   bash reproducir/cluster/todo.sh -n > todo_n.txt 2>&1
   ```

   `diagnostico.sh` lists the system, the modules, the compilers, Python, every Slurm
   partition and the storage. `todo.sh -n` builds and tests the executable, writes the initial
   data and only prints the `sbatch` commands. If a module or the partition changed, edit
   `entorno.sh` (`module load` lines and `VP_PARTICION`) and repeat.

2. One test job, with the resources of a job of the series:

   ```bash
   bash reproducir/cluster/todo.sh -t > todo_t.txt 2>&1
   ```

   It sends `prueba.slurm` to the partition of `entorno.sh`. In about ten minutes
   `log.prueba.<job>.out` shows the compute node, the executable against the laptop, the whole
   chain (run, analysis, measures) on a run of 2000 steps with six values of L, the tests of the
   reduced model, and the cost of 2000 steps of each of the fourteen runs. It ends with
   `PASA` or `FALLA` for each step.

3. Only if everything passes: `bash reproducir/cluster/todo.sh > todo.txt 2>&1` (section 1).

## 1. Everything at once

On the login node, clone both repositories, on the branches of this work, into a new folder
(110 MB in all), and run `todo.sh`:

```bash
mkdir -p <new folder> && cd <new folder>
git clone -b notas/mezcla-fuera-equilibrio https://github.com/erjive/vlasov-poisson_PIC.git
git clone -b mejoras/portar https://github.com/erjive/VlasovPoisson_PIC_sp.git
cd VlasovPoisson_PIC_sp
bash reproducir/cluster/todo.sh > todo.txt 2>&1
```

`todo.sh` does four things:

1. builds the code and tests it (`compilar.sh`); if the test fails it stops;
2. writes the initial data of the fourteen runs (`nube_L.py preparar`, about a minute);
3. submits the reduced model of the series `Lb` and `control` (`reducido.slurm` of
   `vlasov-poisson_PIC`, 16 CPUs, about 2.5 hours), first, because it is the longest job and
   the queue starts jobs in order of arrival;
4. submits the fourteen runs (`enviar.sh primera resto convergencia`).

With `-n`, or on a machine without `sbatch`, it does the first two and only prints the `sbatch`
commands. What is done, or already in the queue, is skipped, so it can be run again. While jobs
of this series are queued it does not rebuild the executable, because `make clean` would delete
the one they are about to use.

A copy of `vlasov-poisson_PIC` that is already on the cluster can be used instead of a new
clone, but it must be on the branch `notas/mezcla-fuera-equilibrio`, and switching a copy whose
`exe/VP_PIC` was rebuilt there may require `git stash` first. New clones avoid that and leave the
old copy and its runs alone.

The initial data written on the cluster agree with those of the laptop to rounding. To start
from exactly the same particles, copy them before running `todo.sh`, with `DEST_SP` the folder
of this repository on the cluster:

```bash
mkdir -p exe/rep/12_nube/ic                                         # cluster
rsync -av exe/rep/12_nube/ic/ $DEST_SP/exe/rep/12_nube/ic/          # laptop
```

## 2. The jobs

Each run is one job of 8 CPUs, which on `olin` are four cores, and 8 GB: the partition gives
2 GB by default, and the analysis of a run with 63 000 particles holds about 3 GB. The job runs
`nube_L.py correr`, `analizar` and `medir` for its run; what is already done is skipped, so a
job that was interrupted can be submitted again with `enviar.sh`. The output of Slurm goes to
`exe/rep/12_nube/slurm_<run>_<job>.out`. The reduced model writes its log to
`log.reducido.<job>.out` in `vlasov-poisson_PIC` and its results to
`exe/demo_eta/reducido/` there.

To see how far everything got:

```bash
bash reproducir/cluster/estado.sh > estado.txt 2>&1
```

On the laptop the runs take between 1 and 47 minutes and the analysis about 17 minutes per
42 000 particles; four cores of the node of `olin` were 1.4 times faster than the laptop with
the other code. While the reduced model holds 8 cores, four runs fit on the rest of the node.
If the node is free, everything should end in about three hours; the queue, which does not
start a later job before an earlier one, can add hours to that. The snapshots
(`vlasov_output.raw`, 13 GB for the fourteen runs) are only needed by `analizar` and can be
deleted afterwards; the maps (`bulto.npz`, 6.5 GB) are needed to measure anything new.

## 3. Back on the laptop

```bash
export VP_CLUSTER_SP=$DEST_SP                    # with user@ in front if the user names differ
reproducir/cluster/traer.sh -r primera resto convergencia
python3 reproducir/scripts/nube_L.py tabla
```

`traer.sh` copies `medidas.npz`, `capas.npz`, `campo.npz`, `resumen.npz`, `nucleo.npz`,
`params_usados.par` and the `.ok`, `.meta`, `.log` and Slurm output of each run, about 0.2 MB
per run. It refuses to overwrite a run that exists on the laptop unless `-f` is given. With `-r`
it also copies the measures of the reduced model (`<run>_diag.npz` and `<run>_capas.npz`) into
`vlasov-poisson_PIC/exe/demo_eta/reducido/`, never over a file that is already there. `-g` adds
the large files: `bulto.npz` of each PIC run and, with `-r`, the runs of the model. `tabla`
only reads the small files, and puts next to each PIC run the run of the reduced model with the
same particles.
