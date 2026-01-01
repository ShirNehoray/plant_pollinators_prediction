set.seed(42)

# =========================
# Packages
# =========================
if (!requireNamespace("softImpute", quietly = TRUE)) install.packages("softImpute")
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")

library(softImpute)
library(pROC)

# =========================
# 1) Create a 30x30 toy edge list of 1s
# =========================
n_rows <- 30
n_cols <- 30

rows <- paste0("R", 1:n_rows)
cols <- paste0("C", 1:n_cols)

n_ones <- 177  # density control


all_pairs <- expand.grid(node_from = rows, node_to = cols, stringsAsFactors = FALSE)
one_idx <- sample(seq_len(nrow(all_pairs)), n_ones)

aggregated_df <- data.frame(
  node_from   = all_pairs$node_from[one_idx],
  node_to     = all_pairs$node_to[one_idx],
  interaction = 1
)

# =========================
# 2) Build ONE interaction matrix P (30x30)
# =========================
P <- matrix(0, nrow = length(rows), ncol = length(cols),
            dimnames = list(rows, cols))

idx <- cbind(match(aggregated_df$node_from, rows),
             match(aggregated_df$node_to,   cols))
P[idx] <- aggregated_df$interaction

cat("P dimensions:", dim(P), "\n")
cat("Total 1s:", sum(P == 1), "Total 0s:", sum(P == 0), "\n")

# =========================
# 3) Create 5 disjoint folds
# =========================
K <- 5
make_folds <- function(n, K) sample(rep(seq_len(K), length.out = n))

ones_in_P  <- which(P == 1, arr.ind = TRUE)
zeros_in_P <- which(P == 0, arr.ind = TRUE)

pos_folds  <- make_folds(nrow(ones_in_P),  K)
zero_folds <- make_folds(nrow(zeros_in_P), K)

cat("\n1s per fold:\n"); print(table(pos_folds))
cat("\n0s per fold:\n"); print(table(zero_folds))

# store coverage check
edge_names <- function(idx_mat, Pmat) {
  apply(idx_mat, 1, function(rc) paste0(rownames(Pmat)[rc[1]], "->", colnames(Pmat)[rc[2]]))
}
held_out_pos_all  <- character(0)
held_out_zero_all <- character(0)

# results storage
all_fold_results <- data.frame()

# =========================
# 4) CV loop: mask + softImpute + evaluate AUC
# =========================
k_values <- c(2, 5, 10)

for (fold in 1:K) {
  
  # Start from the full matrix each fold
  P_fold <- P
  
  remove_pos  <- ones_in_P[pos_folds == fold,  , drop = FALSE]
  remove_zero <- zeros_in_P[zero_folds == fold, , drop = FALSE]
  
  # Mask held-out edges
  P_fold[remove_pos]  <- NA
  P_fold[remove_zero] <- NA
  
  cat("\n====================\n")
  cat("Fold", fold, "\n")
  cat("====================\n")
  cat("Held-out 1s:", nrow(remove_pos),  "\n")
  cat("Held-out 0s:", nrow(remove_zero), "\n")
  
  # coverage check collection
  held_out_pos_all  <- c(held_out_pos_all,  edge_names(remove_pos, P))
  held_out_zero_all <- c(held_out_zero_all, edge_names(remove_zero, P))
  
  # ---- Center matrix (softImpute expects NA for missing) ----
  C <- biScale(P_fold, row.center = TRUE, col.center = TRUE,
               row.scale = FALSE, col.scale = FALSE)
  
  # lambda grid (include lambda0(C))
  lam0 <- lambda0(C)
  # lambda_values <- c(1, 5, 50, 100, lam0)
  lambda_values <- lam0
  
  # Build evaluation labels once
  y_true <- c(rep(1, nrow(remove_pos)), rep(0, nrow(remove_zero)))
  rc_all <- rbind(remove_pos, remove_zero)
  
  for (k in k_values) {
    for (lambda in lambda_values) {
      
      fit <- softImpute(C, rank.max = k, lambda = lambda, type = "svd")
      C_hat <- complete(C, fit)  # filled matrix
      
      # pull predictions for held-out edges
      preds <- C_hat[rc_all]
      
      # AUC (skip if degenerate)
      auc <- as.numeric(pROC::auc(pROC::roc(y_true, preds, quiet = TRUE)))
      
      all_fold_results <- rbind(
        all_fold_results,
        data.frame(
          fold = fold,
          k = k,
          lambda = lambda,
          auc = auc,
          n_holdout_1 = nrow(remove_pos),
          n_holdout_0 = nrow(remove_zero)
        )
      )
      
      cat(sprintf("  k=%2d | lambda=%10.4f | AUC=%.6f\n", k, lambda, auc))
    }
  }
}

# =========================
# 5) Coverage checks (no repeats across folds)
# =========================
stopifnot(all(table(held_out_pos_all) == 1))
stopifnot(all(table(held_out_zero_all) == 1))
cat("\n✅ Coverage OK: every 1-edge and every 0-edge was masked exactly once across folds.\n")

# =========================
# 6) Save results
# =========================
out_file <- "toy_cv_softimpute_results.csv"
write.csv(all_fold_results, out_file, row.names = FALSE)
cat("✅ Saved:", out_file, "\n")

# Show best setting
best <- all_fold_results[which.max(all_fold_results$auc), ]
print(best)
