# Network Recovery and Bayesian Model-Based Clustering for Populations of Neural Networks

Technical appendix to the ST980 MSc dissertation of Akeel Shah (2211111),
Department of Statistics, University of Warwick.
Supervisors: Dr Anastasia Mantziou and Dr Massimiliano Tamborrino.

A two-stage pipeline for directed brain connectivity across epilepsy patients
from scalp EEG (CHB-MIT database, 23 cases, 263 recovered networks):

1. **Recovery** (Chapters 3–6). For each recording, a directed network between
   four EEG channels is recovered by an adapted SMC-ABC scheme for a stochastic
   multi-population Jansen–Rit neural mass model, following Ditlevsen,
   Tamborrino and Tubikanec (2025), run at database scale on the Warwick Avon
   cluster.
2. **Clustering** (Chapters 7–9). The population of networks is analysed with
   a Bayesian mixture of measurement-error models (Mantziou, Lunagómez and
   Mitra, 2024), adapted to directed networks.

---

## Layout

```
R/            stage one: EEG preparation and nSMC-ABC recovery (R)
hpc/          SLURM submission scripts for the Avon cluster
analysis/     stage-one post-processing: aggregation, Chapter 6 figures and tests
clustering/   stage two: simulation studies (Ch. 8) and EEG clustering (Ch. 9)
  └─ mantziou2024_supplement/   unmodified helper files from the published supplement
figures/      per-case figures (all cases) and the figures used in the text
results/      recovered networks and derived tables (CSV)
data/         instructions for obtaining the CHB-MIT recordings (not stored here)
```

---

## Provenance

| Component | Source | Modified here? |
|---|---|---|
| `SplittingJRNMM` (path simulator, C++/R) | Ditlevsen et al. (2025), Supplement B | no |
| `R/functions_SMC_ABC_JRNMM.R`, `R/matrices_SplittingJRNMM.R` | Ditlevsen et al. (2025), Supplement C (code by I. Tubikanec) | essentially unchanged |
| `R/prepare_EEG_functions.R`, `R/make_jobs.R`, `R/main_SMC_ABC_JRNMM.R`, `hpc/*` | this work | — |
| `clustering/mantziou2024_supplement/*` | Mantziou et al. (2024), Supplement A | no |
| `clustering/MCMC_mixture_e_sbm_psi.R` | Mantziou et al. (2024) sampler | **yes**: directed edge index set (12 ordered pairs, no diagonal); node neighbourhood extended to `2(n-1)` entries; `B == 1` Erdős–Rényi branch; `tau_conc` argument (Dirichlet concentration on cluster weights) |
| `clustering/MCMC_mixture_e_sbm_diag.R` | as above, instrumented to record acceptance rates | yes |
| `clustering/*.R` drivers, `analysis/*` | this work | — |

---

## Environment

- R 4.4.2 (`source("R/required_packages.R")` installs the R dependencies);
  the compiled package `SplittingJRNMM` from Supplement B of Ditlevsen et al.
- Python 3.11 with `numpy`, `pandas`, `scipy`, `matplotlib` for `analysis/`.
- Avon: `module load GCC/13.3.0 R/4.4.2`; packages installed into
  `R_LIBS_USER=$HOME/R/library` on a login node (compute nodes are offline).

Scripts in `clustering/` are run with `clustering/` as the working directory.

---

## Data

The `.edf` recordings are not stored here. See `data/README.md` for download
from PhysioNet (CHB-MIT v1.0.0). Case `chb17` produced no usable network;
`chb24` has no entry in `SUBJECTINFO`.

---

## Reproducing the results

### Stage one — recovery (Chapters 5–6)

1. `Rscript R/make_jobs.R chb01 chb02 ...` → `jobs.txt`, `jobs_report.csv`
   (every seizure, kept or skipped, with reason).
2. `N=$(grep -cve '^[[:space:]]*$' jobs.txt); sbatch --array=1-${N} hpc/run_array.sh`
   Each task runs `Rscript R/main_SMC_ABC_JRNMM.R <record> <period> <idx>`
   and writes `ABC_Results_<record>_<period>/`. 16 cores, 48 h walltime;
   5–23 h per recording.
3. `analysis/rho_networks_all_patients.ipynb` → `results/all_patients_rho.csv`
   (posterior mode and mean per directed edge) and `figures/chbXX_networks_all.png`.
4. `analysis/abc_posterior_all_patients.ipynb` (uses `abc_posterior_compare.py`)
   → `results/all_patients_continuous_posteriors.csv`,
   `results/continuous_shift_*.csv`, `figures/chbXX_posteriors_before_vs_during.png`,
   `figures/continuous_shift_summary.pdf` (§6.3–6.4, Tables A.3).
5. `analysis/population_connectivity.R` → `figures/conn_heatmap.pdf`,
   `figures/networks_examples.pdf`, Tables 6.1, 6.2 and A.1 (§6.1–6.2).

### Stage two — clustering (Chapters 8–9)

| Dissertation | Script | Seed(s) |
|---|---|---|
| §8.1 reproduction of Mantziou et al. §5.1, random vs true initialisation | `clustering/Simulations_Section_5.1.R` | as in supplement |
| §8.2 directed four-node study | `clustering/four_node_cluster_sim.R` | 1001, 2002 |
| §9.1 chb01, Erdős–Rényi prior | `clustering/eeg_cluster_ER.R` | 2024; 3000 + r |
| §9.1 chb01, two-block prior; §9.3.2 per-patient fits | `clustering/batch_cluster_all_patients.R` | as in file |
| §9.2 distance-based check (Hamming, k-medoids, LOO classifier) | `clustering/eeg_distance_compare.R` | deterministic |
| §9.3.1 MDS under Hamming / Jaccard / spectral | `clustering/eda_mds.R` | deterministic |
| §9.3.3 paired sign tests, BH correction | `clustering/paired_edge_tests.R` | deterministic |
| §9.3.4 clustering of transitions | `clustering/transition_clusters.R` | deterministic |
| §7.6 acceptance rates, R-hat across restarts | `clustering/convergence_diagnostics.R` | as in file |

Settings common to the EEG fits: 5,000 iterations, 1,000 burn-in, no thinning;
`p ~ Beta(1,11)`, `q ~ Beta(1,4)`, `theta ~ Beta(0.5,0.5)`, `tau_conc = 5`.
Simulation studies use `Beta(0.5,0.5)` throughout and `tau_conc = 1`.

### Figures in Chapter 2

`figures/tenwenty.py` (10–20 electrode layout) and `figures/plot_eeg_trace.R`
(four-channel trace of `chb01_03` around seizure onset).

---

## Key outputs

- `results/all_patients_rho.csv` — the 263 recovered directed networks: per
  recording (case, record, period) the 12 edges, each with posterior mode
  (0/1) and mean (inclusion probability). Appendix Table A.2.
- `results/all_patients_continuous_posteriors.csv` — posterior means of the
  ten continuous parameters per recording. Appendix Table A.3.
- `results/jobs_report.csv` — which seizures entered the study and why.
- `results/SUBJECTINFO` — age and sex per case (PhysioNet).

---

## References

Ditlevsen, S., Tamborrino, M. and Tubikanec, I. (2025). Network inference in a
stochastic multi-population neural mass model via approximate Bayesian
computation. *Annals of Applied Statistics*.

Mantziou, A., Lunagómez, S. and Mitra, R. (2024). Bayesian model-based
clustering for populations of network data. *Annals of Applied Statistics*
18(1), 266–302.

Shoeb, A. (2009). *Application of Machine Learning to Epileptic Seizure Onset
Detection and Treatment*. PhD thesis, MIT. Data via PhysioNet (Goldberger et
al., 2000).
