# =====================================================================
# SIMULATION RECOVERY STUDY (simple, single demonstration)
#
# Purpose: feed the SBM mixture data with a KNOWN ground truth and check
# whether it recovers what we put in. On real chb01 data we can't tell
# whether difficulty comes from the METHOD or from the DATA violating the
# model's assumptions. Here the data is generated FROM the model, so if
# recovery fails the method is at fault; if it succeeds, the method works
# and any real-data difficulty is a data/assumption issue.
#
# Design (the choices we agreed):
#   * BINARY true representatives  -> unambiguous ground truth.
#     Taken from the chb01 fit, thresholded at 0.5.
#   * LOW noise everywhere (p = q = 0.05) -> the "does it work at all?"
#     ideal-conditions check. Expect near-100% recovery.
#   * 20 before + 20 during -> balanced, larger than chb01 so recovery
#     isn't at the mercy of tiny-sample noise.
#   * ONE generated dataset, clustered ONCE (demonstration, not a
#     50-replicate study). Structured so replicates can be added later.
#
# The noise process MATCHES the sampler's generative model exactly:
#   present edge (rep=1) is DROPPED with prob q  (false negative)
#   absent  edge (rep=0) is ADDED   with prob p  (false positive)
# =====================================================================

# ---- 0. Setup -------------------------------------------------------
setwd("C:/Users/123/Documents/aoas1789suppa/Code for Bayesian model-based clustering for populations of netowrk data")
source("simul_fun.R")               # prob_vec_fun_gen
source("label_fun.R")
source("MCMC_mixture_e_sbm_psi.R")  # tuned psi sampler (same one chb01 uses)
source("vec_to_graph.R")
library(LaplacesDemon)
library(igraph)

set.seed(1)   # reproducible demonstration

# ---- 1. Edge layout (identical to chb01_clustering.R) ---------------
n_nodes <- 4
toy_mat <- matrix(nrow = n_nodes, ncol = n_nodes)
lab     <- which(row(toy_mat) != col(toy_mat), arr.ind = TRUE)   # 12 directed
colnames(lab) <- c("row", "col")
n_edges    <- nrow(lab)
slot_names <- paste0("rho_", lab[, "row"], lab[, "col"])

# ---- 2. TRUE binary representatives ---------------------------------
# The two chb01 recovered representatives, thresholded to 0/1 at 0.5.
# Named by edge so order is explicit; reordered into slot_names below.
#
# BEST PRACTICE: take these straight from your live chb01 fit so the
# simulation uses exactly what the model recovered. After running
# chb01_clustering.R, do:
#     true_rep_before <- as.integer(rep_hat_1 >= 0.5)
#     true_rep_during <- as.integer(rep_hat_2 >= 0.5)
#     names(true_rep_before) <- names(true_rep_during) <- slot_names
# and skip the hard-coded block below.
#
# Hard-coded fallback (chb01 reps thresholded at 0.5): before is dense,
# during is sparse (edges switched off) -- matching the BEFORE/DURING plots.
.tb <- c(rho_21=0, rho_31=1, rho_41=1, rho_12=1, rho_32=1, rho_42=1,
         rho_13=1, rho_23=1, rho_43=1, rho_14=1, rho_24=0, rho_34=0)
.td <- c(rho_21=1, rho_31=0, rho_41=1, rho_12=1, rho_32=1, rho_42=0,
         rho_13=0, rho_23=1, rho_43=0, rho_14=1, rho_24=0, rho_34=0)
true_rep_before <- as.integer(.tb[slot_names]); names(true_rep_before) <- slot_names
true_rep_during <- as.integer(.td[slot_names]); names(true_rep_during) <- slot_names
true_reps <- list(before = true_rep_before, during = true_rep_during)

cat("True representative edge counts  before:", sum(true_rep_before),
    " during:", sum(true_rep_during), "\n")
cat("They differ on", sum(true_rep_before != true_rep_during), "of 12 edges.\n\n")

# ---- 3. Generate ONE synthetic dataset ------------------------------
N_before <- 20
N_during <- 20
p_true   <- 0.05   # false-positive rate (absent edge added)
q_true   <- 0.05   # false-negative rate (present edge dropped)

# corrupt a binary representative once -> one observed network
simulate_network <- function(rep_vec, p, q) {
  out <- rep_vec
  present <- rep_vec == 1
  absent  <- rep_vec == 0
  # drop present edges with prob q
  out[present] <- rbinom(sum(present), 1, 1 - q)
  # add absent edges with prob p
  out[absent]  <- rbinom(sum(absent),  1, p)
  as.integer(out)
}

sim_before <- lapply(seq_len(N_before), function(i) simulate_network(true_rep_before, p_true, q_true))
sim_during <- lapply(seq_len(N_during), function(i) simulate_network(true_rep_during, p_true, q_true))

Sim_data <- c(sim_before, sim_during)
N_data   <- length(Sim_data)
z_true   <- c(rep(1L, N_before), rep(2L, N_during))   # TRUE labels (known!)

cat("Generated", N_data, "networks (before =", N_before, ", during =", N_during, ")\n")
cat("mean #edges  before:", round(mean(sapply(sim_before, sum)), 2),
    " during:",             round(mean(sapply(sim_during, sum)), 2), "\n\n")

# ---- 4. MCMC settings (identical to the tuned chb01 driver) ---------
C <- 2; B <- 2
n_iter <- 5000; burn_in_keep <- 0
a_0 <- rep(1, C);   b_0 <- rep(11, C)     # p ~ Beta(1,11)
c_0 <- rep(1, C);   d_0 <- rep(4, C)      # q ~ Beta(1,4)
e_0 <- rep(0.5, C); f_0 <- rep(0.5, C)    # theta ~ Beta(0.5,0.5)
psi     <- 5
pert    <- c(0.03, 0.02, 0.04, 0.05, 0.06)
rw_stps <- c(0.18, 0.18, 0.12, 0.12, 0.12, 0.30)   # tuned values

# RANDOM init (do NOT use z_true to initialise -- that would cheat)
z_init <- sample(seq_len(C), N_data, replace = TRUE)
while (length(unique(z_init)) < C)
  z_init <- sample(seq_len(C), N_data, replace = TRUE)
c_init <- list(c(1,1,2,2), c(1,1,2,2))
prob_vec_init <- prob_vec_fun_gen(z_init, Sim_data, C)

# ---- 5. Run the sampler on the synthetic data -----------------------
res <- MCMC_mix_err_sbm(C, B, n_iter, burn_in_keep, n_nodes, N_data,
                        Sim_data, prob_vec_init,
                        a_0, b_0, c_0, d_0, e_0, f_0,
                        z_init, c_init,
                        0.05, 0.10, 0.3, pert, rw_stps, psi = psi)

# ---- 6. Post-process + score against the KNOWN truth ----------------
keep    <- (burn_in_keep + 1):n_iter
burn    <- 1000
keep    <- keep[keep > burn]
z_draws <- res[[2]][keep, , drop = FALSE]
z_hat   <- apply(z_draws, 2, function(col) as.integer(names(which.max(table(col)))))
z_conf  <- apply(z_draws, 2, function(col) max(table(col)) / length(col))

# label-switch align to truth
if (mean((3L - z_hat) == z_true) > mean(z_hat == z_true)) {
  z_hat <- 3L - z_hat; rep_order <- c(2, 1)
} else rep_order <- c(1, 2)

cat("=====================  RECOVERY  =====================\n")
cat("TRUE labels known by construction (this is a simulation).\n\n")
print(table(true = z_true, recovered = z_hat))
cat("\nassignment recovery:", round(mean(z_hat == z_true) * 100, 1), "%\n")
cat("mean confidence:", round(mean(z_conf), 3), "\n\n")

# ---- 7. Parameter recovery: did we get p, q, reps back? -------------
p_draws <- res[[3]][keep, rep_order, drop = FALSE]
q_draws <- res[[4]][keep, rep_order, drop = FALSE]
cat("--- Noise-parameter recovery ------------------------\n")
cat(sprintf("p  true=%.3f   recovered: before=%.3f  during=%.3f\n",
            p_true, mean(p_draws[,1]), mean(p_draws[,2])))
cat(sprintf("q  true=%.3f   recovered: before=%.3f  during=%.3f\n\n",
            q_true, mean(q_draws[,1]), mean(q_draws[,2])))

rep_hat_1 <- colMeans(res[[5]][[rep_order[1]]][keep, , drop = FALSE])
rep_hat_2 <- colMeans(res[[5]][[rep_order[2]]][keep, , drop = FALSE])
names(rep_hat_1) <- names(rep_hat_2) <- slot_names

rec <- data.frame(edge        = slot_names,
                  true_before = true_rep_before,
                  rec_before  = round(rep_hat_1, 2),
                  true_during = true_rep_during,
                  rec_during  = round(rep_hat_2, 2))
cat("--- Representative recovery (true vs recovered prob) ---\n")
print(rec, row.names = FALSE)

# how many edges did we get right (recovered prob on correct side of 0.5)?
rep_err_before <- sum((rep_hat_1 >= 0.5) != (true_rep_before == 1))
rep_err_during <- sum((rep_hat_2 >= 0.5) != (true_rep_during == 1))
cat(sprintf("\nrepresentative edge errors  before: %d/12   during: %d/12\n",
            rep_err_before, rep_err_during))

# ---- 8. Visual: true vs recovered representatives -------------------
sq_layout <- matrix(c(0,1, 1,1, 1,0, 0,0), ncol = 2, byrow = TRUE)
graph_reform_directed <- function(vec, n_nodes, lab, thresh = 0.5) {
  A <- matrix(0, n_nodes, n_nodes)
  for (k in which(vec >= thresh)) A[lab[k,"row"], lab[k,"col"]] <- 1
  graph_from_adjacency_matrix(A, mode = "directed")
}
par(mfrow = c(2, 2))
plot(graph_reform_directed(true_rep_before, n_nodes, lab), layout = sq_layout,
     edge.arrow.size = .5, vertex.color = "lightblue",  main = "TRUE before")
plot(graph_reform_directed(rep_hat_1,      n_nodes, lab), layout = sq_layout,
     edge.arrow.size = .5, vertex.color = "lightblue",  main = "RECOVERED before")
plot(graph_reform_directed(true_rep_during, n_nodes, lab), layout = sq_layout,
     edge.arrow.size = .5, vertex.color = "lightpink",  main = "TRUE during")
plot(graph_reform_directed(rep_hat_2,      n_nodes, lab), layout = sq_layout,
     edge.arrow.size = .5, vertex.color = "lightpink",  main = "RECOVERED during")
par(mfrow = c(1, 1))

cat("\nDONE. Under ideal (low-noise) conditions, recovery should be near\n")
cat("100%% and p/q posteriors should sit near the true 0.05. If so, the\n")
cat("method works when its assumptions hold -- so real-data difficulty is\n")
cat("a data/assumption issue, not a broken sampler. Next: raise p/q to\n")
cat("chb01-realistic levels and re-run to see where recovery degrades.\n")
