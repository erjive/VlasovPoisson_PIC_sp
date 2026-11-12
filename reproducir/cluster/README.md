# Running the cloud with distributed angular momenta on the cluster

Scripts to build this code on a Slurm cluster and to run and analyze there the series of
`reproducir/scripts/nube_L.py`: the gas cloud far from equilibrium with a spread of the angular
momentum. On the laptop the fourteen runs and their analysis take 7.7 hours, one after the
other. The runs are independent, so on a node with 24 cores six of them run at the same time.

The environment (modules, compiler, number of physical cores) is the one of
`vlasov-poisson_PIC`, file `reproducir/cluster/entorno.sh`, whose README describes the cluster
and its queue. That repository is needed in any case: `nube_L.py` takes the particles and the
measures from it. It must lie next to this one, or its path must be in the variable
`VP_PIC_FIJO`.

| file | purpose |
|---|---|
| `comun.sh` | locates the two repositories and loads the environment; sourced by the others |
| `compilar.sh` | builds `exe/VP_PIC` and runs `probar.sh` |
| `probar.sh` | a run of a few seconds compared with the laptop, the Python modules, and whether the initial data are there |
| `nube.slurm` | runs, analyzes and measures the runs given as arguments |
| `enviar.sh` | submits one job per run: `enviar.sh primera`, `enviar.sh resto convergencia`, or run names |
| `estado.sh` | the jobs in the queue, how far each run got, and the table of measures |
| `traer.sh` | on the laptop: fetches the small files of the finished runs |

All commands are given from the root of the repository.

## 1. First time

On the login node, clone both repositories, on the branches of this work, into a new folder
(110 MB in all). A copy of `vlasov-poisson_PIC` that is already on the cluster can be used
instead, but it must be on the branch `notas/mezcla-fuera-equilibrio`, and switching a copy
whose `exe/VP_PIC` was rebuilt there may require `git stash` first; new clones avoid that and
leave the old copy and its runs alone.

```bash
mkdir -p <new folder> && cd <new folder>
git clone -b notas/mezcla-fuera-equilibrio https://github.com/erjive/vlasov-poisson_PIC.git
git clone -b mejoras/portar https://github.com/erjive/VlasovPoisson_PIC_sp.git
cd VlasovPoisson_PIC_sp
bash reproducir/cluster/compilar.sh
```

`compilar.sh` ends with `PASA la prueba` when the executable reproduces the result of the
laptop and Python finds what it needs (`numpy`, `scipy`). At that point it still reports that
the initial data are missing. Copy them from the laptop, so that both machines start from the
same particles (30 MB), with `DEST_SP` the folder of this repository on the cluster:

```bash
mkdir -p exe/rep/12_nube/ic                                         # cluster: rsync does not create the folders above ic/
python3 reproducir/scripts/nube_L.py preparar                       # laptop
rsync -av exe/rep/12_nube/ic/ $DEST_SP/exe/rep/12_nube/ic/          # laptop
```

They can also be generated on the cluster with the same command; the particles then agree
with those of the laptop to rounding.

## 2. The runs

```bash
reproducir/cluster/enviar.sh -n primera          # shows the sbatch commands
reproducir/cluster/enviar.sh primera             # five jobs
reproducir/cluster/enviar.sh resto convergencia  # the other nine
bash reproducir/cluster/estado.sh > estado.txt 2>&1
```

Each job asks for 8 CPUs, which on `olin` are four cores, and 8 GB: the partition gives 2 GB by
default, and the analysis of a run with 63 000 particles holds about 3 GB. A job runs
`nube_L.py correr`, `analizar` and `medir` for its run. What is already done is skipped, so a
job that was interrupted can be submitted again. The output of Slurm goes to
`exe/rep/12_nube/slurm_<run>_<job>.out`.

On the laptop the runs take between 1 and 47 minutes and the analysis about 17 minutes per
42 000 particles; four cores of the node of `olin` were 1.4 times faster than the laptop with
the other code. The snapshots (`vlasov_output.raw`, 13 GB for the fourteen runs) are only
needed by `analizar` and can be deleted afterwards; the maps (`bulto.npz`, 6.5 GB) are needed
to measure anything new.

## 3. Back on the laptop

```bash
export VP_CLUSTER_SP=$DEST_SP                    # with user@ in front if the user names differ
reproducir/cluster/traer.sh primera              # 0.2 MB per run
python3 reproducir/scripts/nube_L.py tabla
```

`traer.sh` copies `medidas.npz`, `capas.npz`, `campo.npz`, `resumen.npz`, `nucleo.npz`,
`params_usados.par` and the `.ok`, `.meta`, `.log` and Slurm output of each run; `-g` adds `bulto.npz`.
It refuses to overwrite a run that exists on the laptop unless `-f` is given. `tabla` only reads
those small files, and puts next to each run the run of the reduced model with the same
particles if it is in `vlasov-poisson_PIC/exe/demo_eta/reducido/`.

The reduced model on the same lattices is the series `Lb` of `reducido.py`; on the cluster it is
run with `reproducir/cluster/reducido.slurm` of `vlasov-poisson_PIC`.
