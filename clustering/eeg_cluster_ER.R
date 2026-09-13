# =====================================================================
# EEG before/during clustering — ER representatives + balance prior
#
# Two changes from eeg_cluster_real.R, aimed at the degenerate "1-vs-12"
# collapse we diagnosed (all 20 restarts converged to the same bad split):
#
#   (1) ERDOS-RENYI representatives  -> set B = 1.
#       With one block, every edge in a cluster's representative shares a
#       single probability theta. This is the ER simplification you asked
#       for, achieved structurally (no edits to the block-internal helper
#       files, which we can't see). One theta PER CLUSTER, so the two
#       clusters can still differ in density.
#
#   (2) BALANCE PRIOR on the cluster weights -> tau_conc > 1.
#       The original tau update was rdirichlet(1, 1+h), which happily
#       allows a lone-outlier cluster. Raising the concentration makes
#       lopsided splits costly and pushes toward balanced clusters. This
#       is the change that actually targets the collapse; ER alone would
#       likely reproduce it.
#
# Uses the modified sampler MCMC_mixture_ER.R (same as the original but
# with a tau_conc argument and a B==1 guard in the node-membership step).
#
# Requires the same helper files the original sampler sources
# (Log_post_e_sbm_rev.R, MH_updates.R, etc.) to be on the path.
# =====================================================================

# ---- 0. Setup -------------------------------------------------------
setwd("C:/Users/123/Documents/aoas1789suppa/Code for Bayesian model-based clustering for populations of netowrk data")
source("simul_fun.R")
source("label_fun.R")
source("MCMC_mixture_ER.R")     # <-- the MODIFIED sampler (place it alongside the originals)
source("vec_to_graph.R")

library(LaplacesDemon)
library(igraph)

CSV_PATH <- "rho_table_mode_mean.csv"   # <-- full path if not in getwd()

# ---- knobs you can tune ---------------------------------------------
C_CLUST  <- 2     # 2 = before/during (set 3-4 if you switch back to discovery).
TAU_CONC <- 5     # psi: symmetric Dirichlet concentration on cluster weights.
# 1 = flat/outlier-friendly (original). 5 favours a ~50/50 split and
# discourages the lone-outlier collapse, while staying mild enough that a
# GENUINE outlier can still separate. Raise to 10 if the outlier persists.
B        <- 1     # 1 = Erdos-Renyi. Set 2 to go back to the SBM.
N_RESTART <- 20   # random restarts, as before


# ---- 1. Directed representation -------------------------------------
n_nodes <- 4
toy_mat <- matrix(nrow = n_nodes, ncol = n_nodes)
lab     <- which(row(toy_mat) != col(toy_mat), arr.ind = TRUE)
colnames(lab) <- c("row", "col")
n_edges <- nrow(lab)                                        # 12


# ---- 2. Load real networks (modes) ----------------------------------
raw <- read.csv(CSV_PATH, stringsAsFactors = FALSE, check.names = FALSE)
rho_cols <- grep("^rho_", names(raw), value = TRUE)
modes <- sapply(rho_cols, function(cc) as.integer(sub("^([01]).*", "\\1", raw[[cc]])))
colnames(modes) <- rho_cols
slot_names <- paste0("rho_", lab[, "row"], lab[, "col"])
if (!all(slot_names %in% rho_cols))
  stop("CSV missing edge columns: ", paste(setdiff(slot_names, rho_cols), collapse = ", "))
modes_ord <- modes[, slot_names, drop = FALSE]
Sim_data  <- lapply(seq_len(nrow(modes_ord)), function(r) as.integer(modes_ord[r, ]))
N_data    <- length(Sim_data)
z_true    <- ifelse(raw$period == "before", 1L, 2L)   # held out; scoring only

cat("Loaded", N_data, "networks (before =", sum(z_true==1),
    ", during =", sum(z_true==2), ")   B =", B, "  tau_conc =", TAU_CONC, "\n\n")


# ---- 2b. DIAGNOSTIC: is "before" even a coherent group? -------------
# Runs before any MCMC. This is the check that decides whether a
# prototype-based model can ever recover before/during.
D <- as.matrix(dist(do.call(rbind, Sim_data), method = "manhattan"))  # Hamming
before_idx <- which(z_true == 1); during_idx <- which(z_true == 2)
wb <- mean(D[before_idx, before_idx][lower.tri(D[before_idx, before_idx])])
wd <- mean(D[during_idx, during_idx][lower.tri(D[during_idx, during_idx])])
cat("=== Coherence check (mean within-group Hamming distance) ===\n")
cat("  within-BEFORE:", round(wb, 2), " edges\n")
cat("  within-DURING:", round(wd, 2), " edges\n")
cat("  (if BEFORE >> DURING, 'before' is not a tight cluster -> a\n",
    "   two-prototype model cannot recover the before/during split)\n")
c04b <- which(raw$recording == "chb01_04" & raw$period == "before")
if (length(c04b) == 1) {
  cat("  chb01_04-before mean dist to during-group:",
      round(mean(D[c04b, during_idx]), 2), "\n")
  cat("  chb01_04-before mean dist to other before:",
      round(mean(D[c04b, setdiff(before_idx, c04b)]), 2), "\n")
}
# per-before-network: closer to before-group or during-group?
cat("  each BEFORE network - is it nearer the before crowd or during crowd?\n")
for (i in before_idx) {
  db <- mean(D[i, setdiff(before_idx, i)]); dd <- mean(D[i, during_idx])
  cat("   ", raw$recording[i], ": to-before", round(db,1),
      " to-during", round(dd,1),
      if (dd < db) "  <-- looks DURING-like" else "", "\n")
}
cat("============================================================\n\n")


# ---- 3. Fixed settings ----------------------------------------------
C <- C_CLUST
n_iter <- 5000
burn_in_keep <- 0
burn <- 1000
keep <- (burn + 1):n_iter

# ---- PRIOR HYPERPARAMETERS (matched to the SBM driver) --------------
# p (false positive): Beta(1,11), mean 1/12 = 0.083, weak. "Spurious edges rare."
a_0 <- rep(1, C); b_0 <- rep(11, C)      # p ~ Beta(1,11), mean 1/12
# q (false negative): Beta(1,4), mean 0.2 -- higher than p (nets noisier here).
c_0 <- rep(1, C); d_0 <- rep(4, C)       # q ~ Beta(1,4),  mean 0.2
# theta (ER representative density): kept at Jeffreys Beta(0.5,0.5).
e_0 <- rep(0.5, C); f_0 <- rep(0.5, C)   # theta ~ Beta(0.5,0.5) (unchanged)
pert    <- c(0.01, 0.001, 0.02, 0.005, 0.03)
rw_stps <- c(0.01, 0.005, 0.001, 0.0001, 0.0005, 0.05)

# With B=1 every node is in block 1, so c_init MUST be all 1s. One entry per cluster.
c_init <- if (B == 1) {
  lapply(1:C, function(x) rep(1, n_nodes))
} else {
  lapply(1:C, function(x) c(1,1,2,2))
}


# ---- 4. Restarts with the modified sampler --------------------------
# (D already built in the diagnostic section above)

run_once <- function(z_start, seed) {
  set.seed(seed)
  if (length(unique(z_start)) < C) z_start <- sample(1:C, N_data, replace = TRUE)
  pv <- prob_vec_fun_gen(z_start, Sim_data, C)
  MCMC_mix_err_sbm(C, B, n_iter, burn_in_keep, n_nodes, N_data,
                   Sim_data, pv, a_0, b_0, c_0, d_0, e_0, f_0,
                   z_start, c_init, 0.05, 0.10, 0.3, pert, rw_stps,
                   tau_conc = TAU_CONC)
}

score_run <- function(res_i) {
  zc <- apply(res_i[[2]][keep, , drop=FALSE], 2,
              function(col) as.integer(names(which.max(table(col)))))
  sizes <- tabulate(zc, nbins = C)
  if (any(sizes == 0)) return(-Inf)
  wss <- sum(sapply(1:C, function(k){
    idx <- which(zc==k); if (length(idx)<2) return(0)
    mean(D[idx,idx][lower.tri(D[idx,idx])]) }))
  # for discovery (tau_conc=1) don't reward balance; just tight clusters
  -wss
}

inits <- vector("list", N_RESTART)
hc <- hclust(as.dist(D), method = "average")
inits[[1]] <- cutree(hc, k = C)
set.seed(2024)
for (r in 2:N_RESTART) inits[[r]] <- sample(1:C, N_data, replace = TRUE)

cat("Running", N_RESTART, "restarts (C =", C, ")...\n")
runs <- vector("list", N_RESTART); scores <- numeric(N_RESTART)
for (r in seq_len(N_RESTART)) {
  cat("  restart", r, "/", N_RESTART, "\n")
  runs[[r]] <- run_once(inits[[r]], seed = 3000 + r)
  scores[r] <- score_run(runs[[r]])
}
best <- which.max(scores)
cat("\nBest restart:", best, " score", round(scores[best],3), "\n")
cat("All scores:", round(scores,2), "\n")
res <- runs[[best]]


# ---- 5. Consensus labels (no label-switch fix: exploratory, C clusters) --
z_draws <- res[[2]][keep, , drop=FALSE]
z_hat  <- apply(z_draws, 2, function(col) as.integer(names(which.max(table(col)))))
z_conf <- apply(z_draws, 2, function(col) max(table(col))/length(col))

cat("\n--- Cluster x true-period crosstab -------------------------\n")
print(table(true_period = raw$period, cluster = z_hat))
cat("\ncluster sizes:", tabulate(z_hat, nbins = C), "\n")


# ---- 6. Per-cluster composition + representative sharpness ---------
# A cluster is 'real' if it recurs across restarts AND its representative
# is SHARP (edge probs near 0 or 1, not a mush of ~0.5).
cat("\n--- Per-cluster detail ------------------------------------\n")
for (k in 1:C) {
  members <- which(z_hat == k)
  if (length(members) == 0) { cat("cluster", k, ": (empty)\n\n"); next }
  repk <- colMeans(res[[5]][[k]][keep,,drop=FALSE]); names(repk) <- slot_names
  # sharpness: mean distance of edge-probs from 0.5, scaled to [0,1]
  sharp <- mean(abs(repk - 0.5)) * 2
  cat("cluster", k, " (n =", length(members), ") sharpness =", round(sharp,2),
      if (sharp < 0.4) "  <-- MUSHY: likely noise, not a real state" else "  (sharp)", "\n")
  cat("  members:\n")
  for (m in members)
    cat("     ", raw$recording[m], "-", raw$period[m],
        " (conf", round(z_conf[m],2), ")\n")
  strong <- repk[repk > 0.8 | repk < 0.2]
  cat("  defining edges (prob>0.8 present / <0.2 absent):\n     ",
      paste(names(strong), round(strong,2), collapse="  "), "\n\n")
}


# ---- 7. Per-network --------------------------------------------------
cat("--- Per-network assignments -------------------------------\n")
print(data.frame(recording=raw$recording, period=raw$period,
                 true_period_id=z_true, cluster=z_hat, confidence=round(z_conf,3)),
      row.names=FALSE)

cat("\nINTERPRETATION GUIDE:\n")
cat(" - If DURING recordings all land in one sharp cluster -> during is a\n")
cat("   coherent seizure-onset state (expected, reportable).\n")
cat(" - If BEFORE recordings scatter across clusters or fill a mushy one ->\n")
cat("   pre-seizure is heterogeneous (also reportable, and consistent with\n")
cat("   the coherence check at the top).\n")
cat(" - Trust a cluster only if it is SHARP and RECURS across restarts.\n")