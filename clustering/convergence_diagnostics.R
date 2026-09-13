# =====================================================================
# CONVERGENCE DIAGNOSTICS for the chb01 SBM mixture
#
# Answers three questions your supervisor's perturbation query led to:
#   1. Is the REPRESENTATIVE actually moving? (move-fraction + acceptance)
#   2. Are p and q moving? (acceptance + trace)
#   3. Do independent chains AGREE? (Gelman-Rubin R-hat)
#
# The key claim to test: perturbation size is NOT the bottleneck; the
# representative's MH acceptance is near zero across step sizes because
# its posterior is very peaked. If so, the fix is multi-chain / tempering,
# not tuning `perturb`. R-hat across restarts quantifies non-convergence:
# the fact that restarts give 84.6% vs 61.5% should show up as R-hat >> 1.
#
# HOW TO USE
#   - source the INSTRUMENTED sampler: MCMC_mixture_e_sbm_diag.R
#   - run a few chains from different seeds/inits (helper below)
#   - call report_convergence(chains)
# =====================================================================

source("simul_fun.R"); source("label_fun.R")
source("MCMC_mixture_e_sbm_diag.R")   # instrumented sampler (prints [ACCEPT], returns res[[9]])
source("vec_to_graph.R")

# ---- move-fraction: fraction of kept iters where a chain changed -----
move_fraction <- function(mat) {
  mat <- as.matrix(mat)
  if (nrow(mat) < 2) return(NA_real_)
  mean(rowSums(abs(diff(mat))) != 0)
}

# ---- Gelman-Rubin R-hat for one scalar across m chains ---------------
# Each element of chain_list is a numeric vector (post-burn draws of the scalar).
rhat <- function(chain_list) {
  chain_list <- lapply(chain_list, as.numeric)
  m <- length(chain_list)
  n <- min(sapply(chain_list, length))
  if (m < 2 || n < 2) return(NA_real_)
  X <- sapply(chain_list, function(x) x[1:n])   # n x m
  chain_means <- colMeans(X)
  grand_mean  <- mean(chain_means)
  B <- n / (m - 1) * sum((chain_means - grand_mean)^2)   # between
  W <- mean(apply(X, 2, var))                            # within
  if (W == 0) return(NA_real_)
  var_hat <- (n - 1) / n * W + B / n
  sqrt(var_hat / W)
}

# ---- run one instrumented chain -------------------------------------
# Returns res (with res[[9]] = acceptance list). z_start lets you vary inits.
run_chain <- function(Sim_data, N_data, C, B, n_nodes, prob_vec,
                      A0,B0,C0,D0,E0,F0, z_start, c_init,
                      pert, rw_stps, psi, n_iter, burn, seed) {
  set.seed(seed)
  MCMC_mix_err_sbm(C, B, n_iter, burn, n_nodes, N_data, Sim_data, prob_vec,
                   A0,B0,C0,D0,E0,F0, z_start, c_init,
                   0.05, 0.10, 0.3, pert, rw_stps, psi = psi)
}

# ---- main report -----------------------------------------------------
# chains: list of res objects (each from the instrumented sampler),
#         all with the SAME burn already applied inside (burn_in arg).
# keep_idx: row indices into the stored chains to treat as post-burn.
report_convergence <- function(chains, C = 2) {
  cat("\n==================  CONVERGENCE REPORT  ==================\n")
  m <- length(chains)
  cat("chains:", m, "\n\n")

  # ---- acceptance rates (from instrumentation) ----
  cat("--- MH acceptance rates (per chain) ---------------------\n")
  cat(sprintf("%-6s %10s %10s %8s %8s\n","chain","rep_local","rep_indep","p","q"))
  for (k in seq_len(m)) {
    a <- chains[[k]][[9]]
    cat(sprintf("%-6d %10.3f %10s %8.3f %8.3f\n", k,
                a$rep_local,
                ifelse(is.na(a$rep_indep), " NA", sprintf("%.3f", a$rep_indep)),
                a$p, a$q))
  }
  cat("Target ~0.2-0.4. rep_local near 0 across step sizes = STRUCTURAL\n",
      "(peaked representative posterior), NOT a step-size problem.\n\n")

  # ---- move fractions (representative & params) ----
  cat("--- Move fractions (per chain) --------------------------\n")
  cat(sprintf("%-6s %10s %8s %8s\n","chain","rep[1]","p[1]","q[1]"))
  for (k in seq_len(m)) {
    r <- chains[[k]]
    cat(sprintf("%-6d %10.3f %8.3f %8.3f\n", k,
                move_fraction(r[[5]][[1]]),
                move_fraction(r[[3]][,1,drop=FALSE]),
                move_fraction(r[[4]][,1,drop=FALSE])))
  }
  cat("\n")

  # ---- R-hat across chains for p and q (both clusters) ----
  cat("--- Gelman-Rubin R-hat across chains --------------------\n")
  cat("(R-hat > ~1.1 = chains DISAGREE = not converged)\n")
  for (cl in 1:C) {
    p_lists <- lapply(chains, function(r) r[[3]][, cl])
    q_lists <- lapply(chains, function(r) r[[4]][, cl])
    cat(sprintf("  cluster %d :  R-hat(p) = %.3f    R-hat(q) = %.3f\n",
                cl, rhat(p_lists), rhat(q_lists)))
  }

  # ---- cluster-solution agreement across chains ----
  cat("\n--- Do chains find the same clustering? -----------------\n")
  zhat <- lapply(chains, function(r) {
    apply(r[[2]], 2, function(col) as.integer(names(which.max(table(col)))))
  })
  # align each to chain 1 by the better of identity/swap
  ref <- zhat[[1]]
  for (k in seq_along(zhat)) {
    if (mean((3L - zhat[[k]]) == ref) > mean(zhat[[k]] == ref)) zhat[[k]] <- 3L - zhat[[k]]
  }
  agree_mat <- outer(seq_along(zhat), seq_along(zhat),
                     Vectorize(function(i,j) mean(zhat[[i]] == zhat[[j]])))
  cat("pairwise assignment agreement (after alignment):\n")
  print(round(agree_mat, 2))
  cat("\nIf off-diagonal entries are well below 1.0, chains are landing in\n",
      "DIFFERENT modes -- the bistability (84.6% vs 61.5%) is confirmed as\n",
      "non-convergence, and step-size tuning cannot fix it.\n")
  cat("=========================================================\n")
  invisible(list(agree = agree_mat, zhat = zhat))
}

# ---- EXAMPLE DRIVER (adapt paths/objects to your session) -----------
# Assumes you have Sim_data, z_true, prob_vec_init, c_init, and the prior
# vectors already built exactly as in chb01_clustering_priors.R.
#
# n_iter <- 5000; burn <- 1000; C <- 2; B <- 2; PSI <- 5
# pert    <- c(0.03,0.02,0.04,0.05,0.06)          # mid-range
# rw_stps <- c(0.02,0.02,0.01,0.01,0.01,0.05)
#
# # four chains from different inits (hclust + 3 random)
# D  <- as.matrix(dist(do.call(rbind, Sim_data), method="manhattan"))
# hc <- hclust(as.dist(D), method="average")
# inits <- list(cutree(hc,2),
#               sample(1:2,length(Sim_data),TRUE),
#               sample(1:2,length(Sim_data),TRUE),
#               sample(1:2,length(Sim_data),TRUE))
# chains <- lapply(seq_along(inits), function(k)
#   run_chain(Sim_data, length(Sim_data), C, B, n_nodes, prob_vec_init,
#             a_0,b_0,c_0,d_0,e_0,f_0, inits[[k]], c_init,
#             pert, rw_stps, PSI, n_iter, burn, seed = 1000+k))
#
# report_convergence(chains, C = 2)
