# =====================================================================
# BATCH before/during SBM clustering across all patients
#
# Runs the Section-4.3 measurement-error SBM mixture (C = 2: before vs
# during) SEPARATELY on every  *_rho.csv  in a folder, then writes two
# CSVs:
#
#   (1) <out>/memberships_per_network.csv
#       one row per observed network, all patients stacked:
#       patient, recording, period, true_id, cluster, confidence
#
#   (2) <out>/summary_per_patient.csv
#       one row per patient:
#       patient, N, n_before, n_during, agreement_pct, mean_conf,
#       p_before, p_during, q_before, q_during,
#       during_signature  (edges present>0.8 / absent<0.2 in the
#                           during representative, i.e. the rewiring
#                           you compare ACROSS patients)
#
# WHY per-patient (not pooled): robust to between-patient baseline-wiring
# differences. The during_signature column is the thing to eyeball across
# patients to decide whether a pooled C=2 model would ever split cleanly.
#
# NOTES / CAVEATS
#  * Uses the psi-enabled sampler (MCMC_mixture_e_sbm_psi.R) so the tau
#    balance prior (psi=5) is active -- this is what stopped the 1-vs-12
#    collapse on chb01.
#  * Informative priors: p~Beta(1,11) mean 1/12, q~Beta(1,4) mean 0.2,
#    theta~Beta(0.5,0.5) Jeffreys. IDENTICAL across patients for
#    comparability. For thin-N patients (e.g. N=6) these priors dominate
#    -- treat their p/q as prior-driven and lean on the representative,
#    not the agreement %.
#  * This file format has SEPARATE rho_XX_mode and rho_XX_mean columns.
#    We take the _mode columns only. (The old driver's grep("^rho_")
#    would grab both and break -- fixed in load_patient() below.)
#  * Label-switch aware: random inits mean "cluster 1" may be during;
#    we align to the better of the two labellings before scoring.
# =====================================================================

# ---- 0. Setup -------------------------------------------------------
setwd("C:/Users/123/Documents/aoas1789suppa/Code for Bayesian model-based clustering for populations of netowrk data")
source("simul_fun.R")                 # prob_vec_fun_gen
source("label_fun.R")
source("MCMC_mixture_e_sbm_psi.R")    # psi-enabled sampler (place alongside originals)
source("vec_to_graph.R")

library(LaplacesDemon)   # rdirichlet, p.interval
library(igraph)

# ---- knobs ----------------------------------------------------------
CSV_DIR <- "All_results_cvs"    # subfolder, relative to the code folder above
OUT_DIR   <- "."          # where the two output CSVs go
C         <- 2            # before / during
B         <- 2            # SBM blocks (formal at 4 nodes; do not report block structure)
N_RESTART <- 4            # hclust init + (N_RESTART-1) random restarts, keep best
PSI       <- 5            # tau balance prior (5 = favour ~50/50, mild)
n_iter    <- 5000
burn      <- 1000

# priors (identical for every patient)
mk <- function(v) rep(v, C)
A0 <- mk(1);  B0 <- mk(11)   # p ~ Beta(1,11),  mean 1/12
C0 <- mk(1);  D0 <- mk(4)    # q ~ Beta(1,4),   mean 0.2
E0 <- mk(0.5);F0 <- mk(0.5)  # theta ~ Beta(0.5,0.5) Jeffreys
pert    <- c(0.03, 0.02, 0.04, 0.05, 0.06)
rw_stps <- c(0.18, 0.18, 0.12, 0.12, 0.12, 0.30)

# ---- 1. Directed 4-node edge layout (shared by all patients) --------
n_nodes <- 4
toy_mat <- matrix(nrow = n_nodes, ncol = n_nodes)
lab     <- which(row(toy_mat) != col(toy_mat), arr.ind = TRUE)  # 12 directed
colnames(lab) <- c("row", "col")
slot_names <- paste0("rho_", lab[, "row"], lab[, "col"])         # rho_<src><tgt>
c_init  <- list(c(1,1,2,2), c(1,1,2,2))

# ---- 2. Loader FIXED for the _mode / _mean column format ------------
# Returns list(Sim_data, z_true, raw, N) or NULL if the file is malformed.
load_patient <- function(path) {
  raw <- tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
                  error = function(e) NULL)
  if (is.null(raw)) { warning("could not read ", path); return(NULL) }
  if (!all(c("recording","period") %in% names(raw))) {
    warning("missing recording/period in ", path); return(NULL)
  }
  # take ONLY the _mode columns, strip suffix so names become rho_12 etc.
  mode_cols <- grep("^rho_.*_mode$", names(raw), value = TRUE)
  if (length(mode_cols) == 0) { warning("no _mode columns in ", path); return(NULL) }
  modes <- raw[, mode_cols, drop = FALSE]
  colnames(modes) <- sub("_mode$", "", mode_cols)                # rho_12_mode -> rho_12
  # reorder into lab's slot order; check all 12 present
  if (!all(slot_names %in% colnames(modes))) {
    warning("missing edge cols in ", path, ": ",
            paste(setdiff(slot_names, colnames(modes)), collapse = ", "))
    return(NULL)
  }
  modes <- modes[, slot_names, drop = FALSE]
  # coerce to 0/1 integers defensively
  modes[] <- lapply(modes, function(x) as.integer(round(as.numeric(x))))
  if (any(is.na(modes))) { warning("non-numeric mode cell in ", path); return(NULL) }
  Sim_data <- lapply(seq_len(nrow(modes)), function(r) as.integer(unlist(modes[r, ])))
  if (!all(c("before","during") %in% raw$period)) {
    warning("patient in ", path, " lacks both before and during; skipping")
    return(NULL)
  }
  z_true <- ifelse(raw$period == "before", 1L, 2L)
  list(Sim_data = Sim_data, z_true = z_true, raw = raw, N = length(Sim_data))
}

# ---- 3. One SBM fit (single init) -----------------------------------
fit_once <- function(Sim_data, N_data, z_start, seed) {
  set.seed(seed)
  if (length(unique(z_start)) < C) z_start <- sample(1:C, N_data, replace = TRUE)
  pv <- prob_vec_fun_gen(z_start, Sim_data, C)
  MCMC_mix_err_sbm(C, B, n_iter, 0, n_nodes, N_data,
                   Sim_data, pv, A0, B0, C0, D0, E0, F0,
                   z_start, c_init, 0.05, 0.10, 0.3, pert, rw_stps,
                   psi = PSI)
}

# tightness score (lower within-cluster Hamming = better); -Inf if a cluster empties
score_fit <- function(res, D, keep) {
  zc <- apply(res[[2]][keep, , drop = FALSE], 2,
              function(col) as.integer(names(which.max(table(col)))))
  if (any(tabulate(zc, nbins = C) == 0)) return(-Inf)
  -sum(sapply(1:C, function(k) {
    idx <- which(zc == k); if (length(idx) < 2) return(0)
    mean(D[idx, idx][lower.tri(D[idx, idx])])
  }))
}

# ---- 4. Cluster ONE patient (restarts + label-switch-aware scoring) -
cluster_patient <- function(pat, patient_id) {
  Sim_data <- pat$Sim_data; z_true <- pat$z_true; raw <- pat$raw; N_data <- pat$N
  keep <- (burn + 1):n_iter
  D <- as.matrix(dist(do.call(rbind, Sim_data), method = "manhattan"))  # Hamming

  # build inits: hclust first, then random
  inits <- vector("list", N_RESTART)
  hc <- hclust(as.dist(D), method = "average")
  inits[[1]] <- cutree(hc, k = C)
  if (N_RESTART > 1) {
    set.seed(2024)
    for (r in 2:N_RESTART) inits[[r]] <- sample(1:C, N_data, replace = TRUE)
  }

  runs <- vector("list", N_RESTART); sc <- rep(-Inf, N_RESTART)
  for (r in seq_len(N_RESTART)) {
    runs[[r]] <- tryCatch(fit_once(Sim_data, N_data, inits[[r]], seed = 3000 + r),
                          error = function(e) { warning("fit failed (", patient_id,
                            " restart ", r, "): ", conditionMessage(e)); NULL })
    if (!is.null(runs[[r]])) sc[r] <- score_fit(runs[[r]], D, keep)
  }
  if (all(!is.finite(sc))) { warning("all restarts failed for ", patient_id); return(NULL) }
  res <- runs[[which.max(sc)]]

  # consensus labels + confidence
  z_draws <- res[[2]][keep, , drop = FALSE]
  z_hat  <- apply(z_draws, 2, function(col) as.integer(names(which.max(table(col)))))
  z_conf <- apply(z_draws, 2, function(col) max(table(col)) / length(col))

  # label-switch align to true before/during (better of two labellings)
  if (mean((3L - z_hat) == z_true) > mean(z_hat == z_true)) {
    z_hat <- 3L - z_hat; rep_order <- c(2, 1)
  } else rep_order <- c(1, 2)

  p_draws <- res[[3]][keep, rep_order, drop = FALSE]
  q_draws <- res[[4]][keep, rep_order, drop = FALSE]
  rep_before <- colMeans(res[[5]][[rep_order[1]]][keep, , drop = FALSE])
  rep_during <- colMeans(res[[5]][[rep_order[2]]][keep, , drop = FALSE])
  names(rep_before) <- names(rep_during) <- slot_names

  # during signature: edges clearly present / absent in the during representative
  present <- slot_names[rep_during > 0.8]
  absent  <- slot_names[rep_during < 0.2]
  sig <- paste0("+{", paste(present, collapse = ","), "} -{",
                paste(absent, collapse = ","), "}")

  memberships <- data.frame(
    patient    = patient_id,
    recording  = raw$recording,
    period     = raw$period,
    true_id    = z_true,
    cluster    = z_hat,
    confidence = round(z_conf, 3),
    stringsAsFactors = FALSE)

  summary_row <- data.frame(
    patient       = patient_id,
    N             = N_data,
    n_before      = sum(z_true == 1),
    n_during      = sum(z_true == 2),
    agreement_pct = round(mean(z_hat == z_true) * 100, 1),
    mean_conf     = round(mean(z_conf), 3),
    p_before = round(mean(p_draws[,1]),3), p_during = round(mean(p_draws[,2]),3),
    q_before = round(mean(q_draws[,1]),3), q_during = round(mean(q_draws[,2]),3),
    during_signature = sig,
    stringsAsFactors = FALSE)

  list(memberships = memberships, summary = summary_row)
}

# ---- 5. Loop over all *_rho.csv in the folder -----------------------
files <- list.files(CSV_DIR, pattern = "_rho\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No *_rho.csv files found in ", normalizePath(CSV_DIR))
cat("Found", length(files), "patient files.\n")

all_mem <- list(); all_sum <- list()
for (f in files) {
  # derive patient id: strip any leading upload prefix, keep chbNN
  pid <- sub(".*?(chb[0-9]+)_rho\\.csv$", "\\1", basename(f))
  cat("\n==============================================================\n")
  cat("PATIENT:", pid, " file:", basename(f), "\n")
  pat <- load_patient(f)
  if (is.null(pat)) { cat("  -> skipped (see warning)\n"); next }
  cat("  N =", pat$N, " (before", sum(pat$z_true==1), "/ during", sum(pat$z_true==2), ")\n")
  out <- tryCatch(cluster_patient(pat, pid),
                  error = function(e) { warning("patient ", pid, " failed: ",
                    conditionMessage(e)); NULL })
  if (is.null(out)) { cat("  -> failed\n"); next }
  all_mem[[pid]] <- out$memberships
  all_sum[[pid]] <- out$summary
  cat("  agreement:", out$summary$agreement_pct, "%   during sig:",
      out$summary$during_signature, "\n")
}

# ---- 6. Write the two output CSVs -----------------------------------
if (length(all_mem) == 0) stop("No patients produced results.")
mem_df <- do.call(rbind, all_mem)
sum_df <- do.call(rbind, all_sum)
mem_path <- file.path(OUT_DIR, "memberships_per_network.csv")
sum_path <- file.path(OUT_DIR, "summary_per_patient.csv")
write.csv(mem_df, mem_path, row.names = FALSE)
write.csv(sum_df, sum_path, row.names = FALSE)

cat("\n==============================================================\n")
cat("DONE.", nrow(sum_df), "patients clustered.\n")
cat("  ", mem_path, "  (", nrow(mem_df), "networks)\n")
cat("  ", sum_path, "\n\n")
cat("Cross-patient read: compare the during_signature column. If the\n")
cat("same edges recur across patients, a pooled C=2 model would split\n")
cat("cleanly; if each patient rewires different edges, keep per-patient.\n")
print(sum_df[, c("patient","N","agreement_pct","during_signature")], row.names = FALSE)
