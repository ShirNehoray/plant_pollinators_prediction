# =========================
# Packages (only for convenience printing; SVD is base R)
# =========================
# install.packages("dplyr") # if needed
library(dplyr)

set.seed(42)

# =========================
# 1) Build a 50x50 toy edge list of 1s
# =========================
n_rows <- 30
n_cols <- 30
rows <- paste0("R", 1:n_rows)
cols <- paste0("C", 1:n_cols)

n_ones <- 130  # density ~500/2500 = 20%

all_pairs <- expand.grid(node_from = rows, node_to = cols, stringsAsFactors = FALSE)
one_idx <- sample(seq_len(nrow(all_pairs)), n_ones)

aggregated_df <- data.frame(
  node_from   = all_pairs$node_from[one_idx],
  node_to     = all_pairs$node_to[one_idx],
  interaction = 1
)

# =========================
# 2) Build ONE matrix P (train=test in your case)
#    NOTE: rows = node_to, cols = node_from (like your original function)
# =========================
build_interaction_matrix <- function(data, all_rows, all_cols) {
  M <- matrix(0, nrow = length(all_rows), ncol = length(all_cols),
              dimnames = list(all_rows, all_cols))
  
  idx <- cbind(match(data$node_from, all_rows),
               match(data$node_to,   all_cols))
  
  M[idx] <- data$interaction
  M
}

P <- build_interaction_matrix(aggregated_df, rows, cols)


cat("P dim:", dim(P), "\n")
cat("Total 1s:", sum(P == 1), "Total 0s:", sum(P == 0), "\n")


# =========================
# 3) 5-fold edge withholding: ~20% of 1s + ~20% of 0s per fold
#    Disjoint folds => no repeats
# =========================
K <- 5
set.seed(123)

make_folds <- function(n, K) sample(rep(seq_len(K), length.out = n))

ones_in_P  <- which(P == 1, arr.ind = TRUE)
zeros_in_P <- which(P == 0, arr.ind = TRUE)

pos_folds  <- make_folds(nrow(ones_in_P),  K)
zero_folds <- make_folds(nrow(zeros_in_P), K)

cat("\n1s per fold:\n"); print(table(pos_folds))
cat("\n0s per fold:\n"); print(table(zero_folds))


# =========================
# 4) SVD + soft-threshold prediction loop (NO implement_impute)
# =========================
k_values <- 2
lambda_values <- c(1, 5, 20)  # keep small in toy; you can add more

all_results <- data.frame()

for (fold in 1:K) {
  
  # --- mask fold edges ---
  P_fold <- P
  
  remove_pos  <- ones_in_P[pos_folds == fold,  , drop = FALSE]
  remove_zero <- zeros_in_P[zero_folds == fold, , drop = FALSE]
  
  P_fold[remove_pos]  <- NA
  P_fold[remove_zero] <- NA
  
  cat("\n====================\n")
  cat("Fold", fold, "\n")
  cat("====================\n")
  cat("Held-out 1s:", nrow(remove_pos), "Held-out 0s:", nrow(remove_zero), "\n")
  
  # --- build C ---
  # Because you said train=test (same network), we use C = P_fold.
  # (If you had A + P, you'd combine them here.)
  C <- P_fold
  
  # --- center (explicit version of "biScale row.center=TRUE col.center=TRUE") ---
  # We must handle NAs: use means computed from observed values only.
  row_means <- apply(C, 1, function(x) mean(x, na.rm = TRUE))
  col_means <- apply(C, 2, function(x) mean(x, na.rm = TRUE))
  grand_mean <- mean(C, na.rm = TRUE)
  
  # Create centered matrix (still with NAs)
  C_centered <- C
  for (i in 1:nrow(C_centered)) {
    for (j in 1:ncol(C_centered)) {
      if (!is.na(C_centered[i, j])) {
        C_centered[i, j] <- C_centered[i, j] - row_means[i] - col_means[j] + grand_mean
      }
    }
  }
  
  # For SVD we need no NAs. A simple choice: replace NAs with 0 in centered space.
  # (There are other choices; this keeps the toy simple and transparent.)
  C_svd <- C_centered
  C_svd[is.na(C_svd)] <- 0
  
  # --- SVD ---
  sv <- svd(C_svd)
  d <- sv$d
  
  for (k in k_values) {
    k_use <- min(k, length(d))  # safety
    
    U_k <- sv$u[, 1:k_use, drop = FALSE]
    V_k <- sv$v[, 1:k_use, drop = FALSE]
    d_k <- d[1:k_use]
    
    for (lambda in lambda_values) {
      
      # --- soft-threshold singular values (nuclear norm shrinkage) ---
      d_shrunk <- pmax(d_k - lambda, 0)
      
      # --- reconstruct centered prediction ---
      C_hat_centered <- U_k %*% diag(d_shrunk, nrow = k_use, ncol = k_use) %*% t(V_k)
      
      # --- uncenter back to original scale ---
      # C_hat = centered + row_mean + col_mean - grand_mean
      C_hat <- C_hat_centered
      for (i in 1:nrow(C_hat)) {
        for (j in 1:ncol(C_hat)) {
          C_hat[i, j] <- C_hat[i, j] + row_means[i] + col_means[j] - grand_mean
        }
      }
      
      # --- extract predictions for held-out edges only ---
      # build a results df for held-out positives
      df_pos <- data.frame(
        fold = fold,
        k = k_use,
        lambda = lambda,
        type = "heldout_pos",
        node_to = rownames(P)[remove_pos[, "row"]],
        node_from = colnames(P)[remove_pos[, "col"]],
        true = 1,
        pred = C_hat[remove_pos]
      )
      
      # held-out zeros
      df_zero <- data.frame(
        fold = fold,
        k = k_use,
        lambda = lambda,
        type = "heldout_zero",
        node_to = rownames(P)[remove_zero[, "row"]],
        node_from = colnames(P)[remove_zero[, "col"]],
        true = 0,
        pred = C_hat[remove_zero]
      )
      
      all_results <- rbind(all_results, df_pos, df_zero)
    }
  }
}

# =========================
# 5) Quick look at results
# =========================
cat("\nAll results rows:", nrow(all_results), "\n")
print(head(all_results, 10))

# Example: summary by (k, lambda) of mean predicted score for pos vs zero
summary_tbl <- all_results %>%
  group_by(k, lambda, type) %>%
  summarise(mean_pred = mean(pred), .groups = "drop")

cat("\nMean predictions by k/lambda/type:\n")
print(summary_tbl)

