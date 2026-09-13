# =====================================================================
# Self-contained distance-based comparison for the 13 chb01 networks.
# Loads the CSV itself -- depends on NO other script. Only needs the
# 'cluster' package (for PAM / k-medoids).
#
#   install.packages("cluster")   # once, if not already installed
# =====================================================================

CSV_PATH <- "rho_table_mode_mean__1_.csv"   # <-- full path if not in getwd()

# ---- load data (mirrors the loader in the other drivers) ------------
raw <- read.csv(CSV_PATH, stringsAsFactors = FALSE, check.names = FALSE)
rho_cols <- grep("^rho_", names(raw), value = TRUE)
modes <- sapply(rho_cols, function(cc) as.integer(sub("^([01]).*", "\\1", raw[[cc]])))
colnames(modes) <- rho_cols

# column-major directed edge order (rho_21, rho_31, ... ) to match the model
n_nodes <- 4
toy <- matrix(nrow = n_nodes, ncol = n_nodes)
lab <- which(row(toy) != col(toy), arr.ind = TRUE); colnames(lab) <- c("row","col")
slot_names <- paste0("rho_", lab[,"row"], lab[,"col"])
X <- modes[, slot_names, drop = FALSE]          # 13 x 12 binary matrix
Sim_data <- lapply(seq_len(nrow(X)), function(r) as.integer(X[r, ]))
z_true   <- ifelse(raw$period == "before", 1L, 2L)
lab_txt  <- paste0(raw$recording, "-", raw$period)
truth    <- z_true

D <- as.matrix(dist(X, method = "manhattan"))   # Hamming distance
agree_best <- function(pred, truth) {
  a <- mean(pred == truth); b <- mean((3L - pred) == truth); max(a, b)
}

cat("N =", nrow(X), " (before =", sum(truth==1), ", during =", sum(truth==2), ")\n\n")

# ---- (1) Hierarchical clustering ------------------------------------
cat("=== (1) Hierarchical clustering (average linkage) ===\n")
hc <- hclust(as.dist(D), method = "average")
for (k in c(2,3)) {
  cl <- cutree(hc, k = k)
  cat("\n-- k =", k, "--\n"); print(table(period = raw$period, cluster = cl))
  if (k == 2) cat("best-match accuracy:", round(agree_best(cl, truth)*100,1), "%\n")
}

# ---- (2) Leave-one-out nearest-centroid classifier ------------------
cat("\n=== (2) Leave-one-out nearest-centroid classifier ===\n")
pred_loo <- integer(nrow(X))
for (i in seq_len(nrow(X))) {
  tr <- setdiff(seq_len(nrow(X)), i)
  cen1 <- colMeans(X[tr[truth[tr]==1], , drop=FALSE])
  cen2 <- colMeans(X[tr[truth[tr]==2], , drop=FALSE])
  pred_loo[i] <- if (sum(abs(X[i,]-cen1)) <= sum(abs(X[i,]-cen2))) 1L else 2L
}
cat("LOO accuracy:", round(mean(pred_loo==truth)*100,1), "%\n")
print(table(true = raw$period, predicted = ifelse(pred_loo==1,"before","during")))
cat("misclassified:", paste(lab_txt[pred_loo != truth], collapse=", "), "\n")

# ---- (3) k-medoids / PAM (k = 2) ------------------------------------
cat("\n=== (3) k-medoids / PAM (k = 2) ===\n")
if (requireNamespace("cluster", quietly = TRUE)) {
  pam2 <- cluster::pam(as.dist(D), k = 2)
  print(table(period = raw$period, cluster = pam2$clustering))
  cat("best-match accuracy:", round(agree_best(pam2$clustering, truth)*100,1), "%\n")
  cat("medoids:", paste(lab_txt[pam2$id.med], collapse=", "), "\n")
} else cat("  install.packages('cluster') to enable PAM\n")

# ---- (4) Outlier ranking --------------------------------------------
cat("\n=== (4) Outlier check (mean Hamming distance to all others) ===\n")
meand <- rowSums(D) / (nrow(D)-1)
for (i in order(meand, decreasing = TRUE))
  cat("   ", sprintf("%-18s", lab_txt[i]), "mean dist", round(meand[i],2),
      if (meand[i] > mean(meand)+sd(meand)) "  <-- OUTLIER" else "", "\n")