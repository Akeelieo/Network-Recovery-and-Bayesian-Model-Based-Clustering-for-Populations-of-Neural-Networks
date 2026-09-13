# =====================================================================
# Directed 4-node networks, 2 clusters
#   Cluster 1 representative: complete digraph (all 12 edges present)
#   Cluster 2 representative: empty graph    (no edges)
#   50 noisy observations per cluster -> N = 100 networks
# Goal: run the measurement-error mixture MCMC and recover z
#       (which network belongs to which cluster)
# =====================================================================

# ---- 0. Setup -------------------------------------------------------
setwd("C:/Users/123/Documents/aoas1789suppa/Code for Bayesian model-based clustering for populations of netowrk data")

source("simul_fun.R")
source("label_fun.R")
source("MCMC_mixture_e_sbm_rev.R")   # NOTE: must be edited for directed case, see README notes
source("vec_to_graph.R")

library(LaplacesDemon)
library(igraph)


# ---- 1. Network representation -------------------------------------
# 4 nodes, DIRECTED, no self-loops  ->  4*3 = 12 possible edges
n_nodes <- 4
toy_mat <- matrix(nrow = n_nodes, ncol = n_nodes)
lab     <- which(row(toy_mat) != col(toy_mat), arr.ind = TRUE)
colnames(lab) <- c("row", "col")
n_edges <- nrow(lab)          # 12
n_edges


# ---- 2. The two "true" representatives ------------------------------
# These are fixed by hand (not drawn from an SBM), so no seed needed.
rep_full  <- rep(1, n_edges)   # cluster 1: every directed edge present
rep_empty <- rep(0, n_edges)   # cluster 2: no edges

repres_true <- list(rep_full, rep_empty)


# ---- 3. Simulate the 100 observed networks --------------------------
# Measurement-error model, applied independently to each of the 12 slots:
#   representative has 0  ->  becomes 1 with prob p_c   (false positive)
#   representative has 1  ->  becomes 0 with prob q_c   (false negative)

simulate_noisy <- function(representative, p_c, q_c, n_obs, seed) {
  set.seed(seed)
  out <- vector("list", n_obs)
  for (m in seq_len(n_obs)) {
    g <- representative
    flip_fp <- (g == 0) & (runif(length(g)) < p_c)   # 0 -> 1
    flip_fn <- (g == 1) & (runif(length(g)) < q_c)   # 1 -> 0
    g[flip_fp] <- 1
    g[flip_fn] <- 0
    out[[m]] <- g
  }
  out
}

p_true <- c(0.10, 0.10)   # false-positive rate, one per cluster
q_true <- c(0.20, 0.20)   # false-negative rate, one per cluster
n_per_cluster <- 50

sim_c1 <- simulate_noisy(repres_true[[1]], p_true[1], q_true[1], n_per_cluster, seed = 1001)
sim_c2 <- simulate_noisy(repres_true[[2]], p_true[2], q_true[2], n_per_cluster, seed = 2002)

Sim_data <- c(sim_c1, sim_c2)              # list of 100 length-12 vectors
z_true   <- c(rep(1, n_per_cluster), rep(2, n_per_cluster))

length(Sim_data)                            # 100
mean(sapply(sim_c1, sum))                   # ~ 12 * (1 - q)  = ~9.6 edges
mean(sapply(sim_c2, sum))                   # ~ 12 * p        = ~1.2 edges


# ---- 4. MCMC settings -----------------------------------------------
C <- 2      # number of clusters
B <- 2      # number of blocks in the SBM prior on the representatives
N_data <- length(Sim_data)
n_iter <- 5000
burn_in_keep <- 0        # keep everything, trim afterwards

# Beta(0.5,0.5) Jeffreys priors, one entry per cluster
a_0 <- rep(0.5, C); b_0 <- rep(0.5, C)   # prior for p
c_0 <- rep(0.5, C); d_0 <- rep(0.5, C)   # prior for q
e_0 <- rep(0.5, C); f_0 <- rep(0.5, C)   # prior for theta

# Metropolis proposal tuning (same values the paper used)
pert    <- c(0.01, 0.001, 0.02, 0.005, 0.03)
rw_stps <- c(0.01, 0.005, 0.001, 0.0001, 0.0005, 0.05)

# Initial block memberships: one length-4 vector per cluster
c_init <- list(c(1, 1, 2, 2),
               c(1, 1, 2, 2))

# Warm-start edge probabilities from the (true) labels
prob_vec_init <- prob_vec_fun(z_true, Sim_data)


# ---- 5. Run the sampler ---------------------------------------------
res <- MCMC_mix_err_sbm(C, B, n_iter, burn_in_keep, n_nodes, N_data,
                        Sim_data, prob_vec_init,
                        a_0, b_0, c_0, d_0, e_0, f_0,
                        z_true, c_init,
                        0.01, 0.01, 0.3,
                        pert, rw_stps)

#res[[1]] tau, [[2]] z, [[3]] p, [[4]] q,
#res[[5]] representatives, [[6]] c, [[7]] w, [[8]] theta


# ---- 6. Posterior summaries -----------------------------------------
burn <- 1000                       # discard first 1000 draws
keep <- (burn + 1):n_iter

# --- noise parameters
p_draws <- res[[3]][keep, , drop = FALSE]
q_draws <- res[[4]][keep, , drop = FALSE]

post_mean_p <- apply(p_draws, 2, mean)
post_mean_q <- apply(q_draws, 2, mean)

cat("posterior mean p:", round(post_mean_p, 3), " (true:", p_true, ")\n")
cat("posterior mean q:", round(post_mean_q, 3), " (true:", q_true, ")\n")

apply(p_draws, 2, p.interval)      # 95% credible intervals
apply(q_draws, 2, p.interval)

# --- mixing weights
apply(res[[1]][keep, , drop = FALSE], 2, mean)   # should be ~0.5, 0.5


# ---- 7. THE MAIN THING: cluster labels z ----------------------------
z_draws <- res[[2]][keep, , drop = FALSE]        # (kept draws) x 100

# consensus label = most frequent cluster across draws, per network
z_hat <- apply(z_draws, 2, function(col) as.integer(names(which.max(table(col)))))

# how confident is each assignment? (proportion of draws in the modal cluster)
z_conf <- apply(z_draws, 2, function(col) max(table(col)) / length(col))

# cross-tabulate against the truth
confusion <- table(true = z_true, estimated = z_hat)
print(confusion)

cat("\nagreement with truth:",
    round(mean(z_hat == z_true) * 100, 1), "%\n")
cat("mean assignment confidence:", round(mean(z_conf), 3), "\n")

# per-network detail
z_summary <- data.frame(network    = seq_len(N_data),
                        true       = z_true,
                        estimated  = z_hat,
                        confidence = round(z_conf, 3))
head(z_summary, 10)
tail(z_summary, 10)


# ---- 8. Recovered representatives ------------------------------------
# each res[[5]][[c]] is (draws x 12); average post-burn-in to get the
# posterior probability that each directed edge is present
rep_hat_1 <- colMeans(res[[5]][[1]][keep, , drop = FALSE])
rep_hat_2 <- colMeans(res[[5]][[2]][keep, , drop = FALSE])

round(rep_hat_1, 2)     # expect all near 1  (complete graph)
round(rep_hat_2, 2)     # expect all near 0  (empty graph)


# ---- 9. Optional: visualise one recovered representative -------------
# requires graph_reform in vec_to_graph.R to use mode = "directed"
 g1 <- graph_reform(round(rep_hat_1), n_nodes)
 plot(g1, edge.arrow.size = 0.5, vertex.color = "lightblue")
 