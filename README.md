# Network Recovery and Bayesian Clustering for Populations of Neural Networks

ST980 MSc Dissertation — Akeel Shah, Department of Statistics, University of Warwick.
Supervisors: Dr Anastasia Mantziou and Dr Massimiliano Tamborrino.

A two-stage pipeline for studying directed brain connectivity across epilepsy
patients from scalp EEG:

1. **Recovery.** For each recording, a directed network between four EEG
   channels is recovered by likelihood-free inference (an adapted SMC-ABC
   scheme) for a stochastic multi-population Jansen–Rit neural mass model.
2. **Clustering.** The resulting population of networks is analysed with a
   Bayesian mixture of measurement-error models, grouping recordings by
   connectivity structure and summarising each group by a representative
   network.

Data: the CHB-MIT Scalp EEG Database (23 cases, paediatric epilepsy patients).

---

## Repository structure

```
.
├── README.md
├── R/                          # inference and data-preparation code (R)
│   ├── required_packages.R         # package dependencies
│   ├── matrices_SplittingJRNMM.R   # exponential / covariance matrices for the splitting scheme
│   ├── functions_SMC_ABC_JRNMM.R   # SMC-ABC functions (kernels, summaries, distance, iterations)
│   ├── prepare_EEG_functions.R     # read .edf, parse seizures, window, scale, grid-align, cache
│   ├── make_jobs.R                 # build the SLURM task list (jobs.txt) + selection report
│   └── main_SMC_ABC_JRNMM.R        # run one recovery: Rscript main_...R <record> <period> [idx]
├── hpc/                        # cluster submission scripts (SLURM, Warwick Avon)
│   ├── run_array.sh                # job array: one task per line of jobs.txt
│   ├── run_patient.sh              # per-patient: download .edf + submit before/during jobs
│   └── run_edf.sbatch              # single-recording sbatch (48 cores)
├── clustering/                 # stage-two Bayesian clustering code   ### ADD ###
├── analysis/                   # aggregation + figures
│   └── aggregate_rho.R             # build all_patients_rho.csv from ABC_Results/   ### ADD ###
├── results/
│   ├── all_patients_rho.csv        # recovered networks: per-edge posterior mode + mean, before/during
│   ├── SUBJECTINFO                 # subject age and sex
│   ├── jobs.txt                    # example task list
│   └── jobs_report.csv             # kept/skipped seizures with reasons
└── data/
    └── README.md                   # how to download the CHB-MIT .edf files (not stored here)
```

`### ADD ###` marks folders/files to drop in from your machine (see "To finish", below).

---

## Data

The raw `.edf` recordings are **not** included (tens of gigabytes, and freely
available). Download them from PhysioNet into `data/` — see `data/README.md`.
Each patient's directory also contains a `chbNN-summary.txt` file holding the
seizure annotations parsed by `prepare_EEG_functions.R`.

---

## Environment

- **R 4.4.2.** Install dependencies with `source("R/required_packages.R")`.
- The path simulator is the compiled package **`SplittingJRNMM`**.
- On the Warwick **Avon** cluster: `module load GCC/13.3.0 R/4.4.2`, with
  packages installed into a personal library (`R_LIBS_USER=$HOME/R/library`)
  on a login node, since the compute nodes have no internet access.

---

## Reproducing the results

1. **Download the data** (see `data/README.md`).
2. **Build the task list** for one or more patients:
   ```
   Rscript R/make_jobs.R chb01 chb02 ...
   ```
   This writes `jobs.txt` (one run per line: `record period seizure_index L`)
   and `jobs_report.csv` (every seizure, kept or skipped, with the reason).
3. **Run the recovery** as a SLURM job array sized to the task list:
   ```
   N=$(grep -cve '^[[:space:]]*$' jobs.txt)
   sbatch --array=1-${N} hpc/run_array.sh
   ```
   Each task calls `Rscript R/main_SMC_ABC_JRNMM.R <record> <period> <idx>` and
   writes an `ABC_Results_<record>_<period>/` folder of posterior particles.
   (`hpc/run_patient.sh` is an alternative that also downloads the needed
   `.edf` files and submits a before/during job per seizure recording.)
4. **Aggregate** the per-recording posteriors into the network population:
   ```
   Rscript analysis/aggregate_rho.R      ### ADD this script ###
   ```
   → `results/all_patients_rho.csv` (per edge: posterior mode and mean).
5. **Cluster** the population with the stage-two model (`clustering/`).
6. **Figures/tables** are produced by the scripts in `analysis/`.

---

## Key outputs

- `results/all_patients_rho.csv` — the recovered directed networks: for each
  recording (patient, record, period) the 12 directed edges, each with a
  posterior **mode** (0/1 point estimate) and **mean** (inclusion probability).
- `results/jobs_report.csv` — audit trail of which seizures entered the study.
- `results/SUBJECTINFO` — subject age and sex, for external validation.

---

## To finish

- Add your **stage-two clustering code** under `clustering/`.
- Add the **aggregation script** (`analysis/aggregate_rho.R`) that builds
  `all_patients_rho.csv` from the `ABC_Results_*` folders.
- Copy the authoritative `hpc/` scripts from the cluster (this repo's copies
  should match those in `~/nJRNMM_edf/` on Avon).
- Optionally record exact package versions (`sessionInfo()` or an `renv.lock`).

## Authors and licence

Code by Akeel Shah. The SMC-ABC scheme adapts the method of Ditlevsen,
Tamborrino and Tubikanec (2025); `functions_SMC_ABC_JRNMM.R` and
`matrices_SplittingJRNMM.R` are based on code by Irene Tubikanec. The
clustering follows Mantziou, Lunagómez and Mitra (2024).
