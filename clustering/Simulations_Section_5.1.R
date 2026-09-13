# =====================================================================
# REPRODUCTION STUDY (Mantziou et al. 2024, Section 5.1) -- clean driver
#
# Runs the 12 regimes at REDUCED chain length (1000 iters, 200 burn-in),
# saves each regime's draws, then builds two result tables:
#   (A) posterior means of p, q, theta   -> tab:repro-p and companions
#   (B) clustering recovery: ARI, purity, entropy  -> tab:repro-clust
#
# All paths go through ONE variable, OUT. Set it once, below.
# Burn-in is 200 to match the 1000-iteration run (NOT 150000).
# =====================================================================

# ---- 0. setup -------------------------------------------------------
setwd("C:/Users/123/Downloads/aoas1789suppa (3)/Code for Bayesian model-based clustering for populations of netowrk data")
source("simul_fun.R")
source("label_fun.R")
source("MCMC_mixture_e_sbm_rev.R")   # ORIGINAL undirected sampler
source("vec_to_graph.R")

library(LaplacesDemon)   # p.interval
library(fst)
library(mclust)          # adjustedRandIndex  (install.packages("mclust") if needed)

OUT   <- "C:/Users/123/Documents/results"   # <-- one output dir for everything
N_ITER <- 1000
BURN   <- 200
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- 1. representatives (exactly as Mantziou) -----------------------
toy_mat <- matrix(nrow = 21, ncol = 21)
lab <- which(upper.tri(toy_mat), arr.ind = TRUE)
lab_theta <- c("11", "12", "22")

create_repres <- function(memb, thet, labl, sd) {
  a <- lapply(1:3, function(x) vector(length = 210))
  for (i in 1:3) {
    a[[i]] <- lab_adj(a[[i]], memb[[i]], labl)
    a[[i]] <- gener_edges_sbm(a[[i]], thet[[i]], sd)
    sd <- sd + 5000
  }
  a
}

# SBM structure 1
c_1 <- c(rep(1,5),rep(2,5),rep(1,5),rep(2,6))
c_2 <- c(rep(2,10),rep(1,11))
c_3 <- c(rep(2,5),rep(1,10),rep(2,6))
c_study_1 <- list(c_1,c_2,c_3)
theta_study_1 <- rep(list(setNames(c(0.8,0.2,0.8), lab_theta)), 3)
repres_study_1 <- create_repres(c_study_1, theta_study_1, lab, 5000)

# SBM structure 2
c_1 <- c(rep(1,7),rep(2,3),rep(1,8),rep(2,3))
c_2 <- c(rep(2,10),rep(1,11))
c_3 <- c(rep(2,7),rep(1,6),rep(2,8))
c_study_2 <- list(c_1,c_2,c_3)
theta_study_2 <- rep(list(setNames(c(0.7,0.05,0.8), lab_theta)), 3)
repres_study_2 <- create_repres(c_study_2, theta_study_2, lab, 15000)

repres <- list(repres_study_1, repres_study_2)
c_all  <- list(c_study_1,       c_study_2)

# ---- 2. priors + MCMC tuning ----------------------------------------
a_0 <- b_0 <- c_0 <- d_0 <- e_0 <- f_0 <- rep(0.5, 3)
pert    <- c(0.01,0.001,0.02,0.005,0.03)
rw_stps <- c(0.01,0.005,0.001,0.0001,0.0005,0.05)

# ---- 3. generate + fit the 12 regimes -------------------------------
# regime layout matches the triple loop: p in {.1,.2,.3}, q != p, 2 structures
z_tr_mix <- c(rep(1,60), rep(2,60), rep(3,60))
t_vals   <- seq(0.1, 0.3, 0.1)
data_seed  <- 1
regime_ind <- 1
regime_key <- data.frame()   # record (regime, p, q, structure) for the tables

for (i in t_vals) {
  for (k in t_vals[t_vals != i]) {
    for (j in 1:2) {
      p <- rep(i, 3); q <- rep(k, 3)
      Sim_mix <- simulate_data_e(repres[[j]], p, q, 180, 3, data_seed)
      
      # random starting partition (NOT the truth) -- forces the chain to find the clusters
      set.seed(data_seed)
      z_init <- sample(1:3, 180, replace = TRUE)
      while (length(unique(z_init)) < 3)
        z_init <- sample(1:3, 180, replace = TRUE)
      
      prob_vec_mix <- prob_vec_fun(z_init, Sim_mix)   # build from z_init, not the truth
      res <- MCMC_mix_err_sbm(3, 2, N_ITER, 0, 21, 180, Sim_mix, prob_vec_mix,
                              a_0,b_0,c_0,d_0,e_0,f_0, z_init, c_all[[j]],
                              0.01, 0.01, 0.3, pert, rw_stps)
      
      rdir <- file.path(OUT, paste0("regime_", regime_ind))
      dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(res[[2]], file.path(rdir, "z.rds"),     version = 2)
      saveRDS(res[[3]], file.path(rdir, "p.rds"),     version = 2)
      saveRDS(res[[4]], file.path(rdir, "q.rds"),     version = 2)
      saveRDS(res[[8]], file.path(rdir, "theta.rds"), version = 2)
      
      regime_key <- rbind(regime_key,
                          data.frame(regime = regime_ind, p = i, q = k, structure = j))
      data_seed  <- data_seed + 50000
      regime_ind <- regime_ind + 1
      cat("done regime", regime_ind - 1, ": p =", i, " q =", k, " struct =", j, "\n")
    }
  }
}
write.csv(regime_key, file.path(OUT, "regime_key.csv"), row.names = FALSE)

# ---- 4. TABLE A: posterior means of p, q, theta ---------------------
triple <- function(v) sprintf("(%.2f,%.2f,%.2f)", v[1], v[2], v[3])
drop_burn <- function(m) m[-(1:BURN), , drop = FALSE]

tabA <- data.frame(regime = 1:12, true_p = NA, true_q = NA,
                   p_hat = NA, q_hat = NA,
                   theta1 = NA, theta2 = NA, theta3 = NA)
for (i in 1:12) {
  rdir <- file.path(OUT, paste0("regime_", i))
  p  <- drop_burn(readRDS(file.path(rdir, "p.rds")))
  q  <- drop_burn(readRDS(file.path(rdir, "q.rds")))
  th <- readRDS(file.path(rdir, "theta.rds"))
  tabA$true_p[i] <- regime_key$p[i]; tabA$true_q[i] <- regime_key$q[i]
  tabA$p_hat[i]  <- triple(colMeans(p))
  tabA$q_hat[i]  <- triple(colMeans(q))
  tabA$theta1[i] <- triple(colMeans(drop_burn(th[[1]])))
  tabA$theta2[i] <- triple(colMeans(drop_burn(th[[2]])))
  tabA$theta3[i] <- triple(colMeans(drop_burn(th[[3]])))
}
cat("\n===== TABLE A: posterior means (p, q per cluster) =====\n")
print(tabA[, c("regime","true_p","true_q","p_hat","q_hat")], row.names = FALSE)
write.csv(tabA, file.path(OUT, "table_pqtheta.csv"), row.names = FALSE)

# ---- 5. TABLE B: clustering recovery (ARI, purity, entropy) ---------
# modal label per network from the post-burn z draws
modal_labels <- function(zmat) apply(zmat, 2, function(col)
  as.integer(names(which.max(table(col)))))

purity_fn <- function(est, tru) {
  ct <- table(est, tru)
  sum(apply(ct, 1, max)) / length(tru)
}
entropy_fn <- function(est, tru) {   # normalised conditional entropy, lower = purer
  ct <- table(est, tru); n <- sum(ct)
  H <- 0
  for (r in 1:nrow(ct)) {
    rs <- sum(ct[r, ]); if (rs == 0) next
    pr <- ct[r, ] / rs; pr <- pr[pr > 0]
    H <- H + (rs / n) * (-sum(pr * log(pr)))
  }
  H
}

tabB <- data.frame(regime = 1:12, true_p = NA, true_q = NA,
                   ARI = NA, purity = NA, entropy = NA, mean_conf = NA)
for (i in 1:12) {
  rdir <- file.path(OUT, paste0("regime_", i))
  zmat <- drop_burn(readRDS(file.path(rdir, "z.rds")))   # (iters-burn) x 180
  est  <- modal_labels(zmat)
  conf <- mean(apply(zmat, 2, function(col) max(table(col)) / length(col)))
  tabB$true_p[i]    <- regime_key$p[i]; tabB$true_q[i] <- regime_key$q[i]
  tabB$ARI[i]       <- round(adjustedRandIndex(est, z_tr_mix), 3)
  tabB$purity[i]    <- round(purity_fn(est, z_tr_mix), 3)
  tabB$entropy[i]   <- round(entropy_fn(est, z_tr_mix), 3)
  tabB$mean_conf[i] <- round(conf, 3)
}
cat("\n===== TABLE B: clustering recovery =====\n")
print(tabB, row.names = FALSE)
write.csv(tabB, file.path(OUT, "table_clustering.csv"), row.names = FALSE)

cat("\nDONE. Wrote regime_key.csv, table_pqtheta.csv, table_clustering.csv to",
    OUT, "\nTable B (ARI/purity/entropy/confidence) is your tab:repro-clust.\n")


