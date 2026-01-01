# =========================
# 0) Packages
# =========================
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("softImpute", quietly = TRUE)) install.packages("softImpute")
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")

library(dplyr)
library(softImpute)
library(pROC)

# Study_id_chose <- "1_Bartomeus"


Interaction_data <- read.csv("./csv/Interaction_edges.csv", encoding = "latin1", stringsAsFactors = FALSE)
one_study_interactions <- Interaction_data %>% 
  filter(Study_id == Study_id_chose)
A_l <-  one_study_interactions



# Create lookup table for network IDs
layer_names <- unique(A_l$layer_from)
layer_map <- data.frame(
  original = layer_names,
  new_id = paste0("layer_", seq_along(layer_names))
)

# Apply the mapping
aggregated_df <- A_l %>%
  left_join(layer_map, by = c("layer_from" = "original")) %>%
  mutate(layer_from = new_id,
         layer_to   = new_id) %>%
  select(layer_from, node_from, layer_to, node_to, weight, type)

# =========================
# 1) Your matrix builder (unchanged)
# =========================
build_interaction_matrix <- function(data, layers_to_filter) {
  # Step 1: Filter rows based on specified layers
  layers <- paste0("layer_", layers_to_filter)
  filtered_data <- subset(data, layer_from %in% layers)
  
  # Step 2: Aggregate weights for identical species pairs
  aggregated_data <- filtered_data %>%
    group_by(node_from, node_to) %>%
    summarise(weight = sum(weight), .groups = "drop")
  
  # Step 3: Create the matrix with specific row and column species
  species_from <- unique(aggregated_data$node_from)  # Columns
  species_to   <- unique(aggregated_data$node_to)    # Rows
  
  # Initialize an empty matrix
  interaction_matrix <- matrix(
    0,
    nrow = length(species_to),
    ncol = length(species_from),
    dimnames = list(species_to, species_from)
  )
  
  # Populate the matrix with aggregated weights
  for (i in 1:nrow(aggregated_data)) {
    row <- aggregated_data$node_to[i]    # Rows represent 'node_to'
    col <- aggregated_data$node_from[i]  # Cols represent 'node_from'
    interaction_matrix[row, col] <- aggregated_data$weight[i]
  }
  
  return(interaction_matrix)
}

# =========================
# 2) Your implement_impute (unchanged)
#    NOTE: uses global variables: P, P_original, remove_indices, zeros_to_remove_indices
# =========================
implement_impute <- function(C, k, lambda) { 
  fit <- softImpute(C, rank.max = k, lambda = lambda, type = "svd", maxit = 600)
  
  C_reconstructed <- softImpute::complete(C, fit)
  
  # Extract reconstructed P
  P_reconstructed <- C_reconstructed[rownames(P), colnames(P)]
  
  # Combine indices of removed ones and zeros
  if (is.null(dim(remove_indices))) {
    test_indices <- rbind(
      data.frame(row = remove_indices["row"], col = remove_indices["col"], label = 1),
      data.frame(row = zeros_to_remove_indices[, "row"], col = zeros_to_remove_indices[, "col"],
                 label = rep(0, nrow(zeros_to_remove_indices)))
    )
  } else {
    test_indices <- rbind(
      data.frame(row = remove_indices[, "row"], col = remove_indices[, "col"], label = rep(1, nrow(remove_indices))),
      data.frame(row = zeros_to_remove_indices[, "row"], col = zeros_to_remove_indices[, "col"],
                 label = rep(0, nrow(zeros_to_remove_indices)))
    )
  }
  
  test_rows <- rownames(P)[test_indices$row]
  test_cols <- colnames(P)[test_indices$col]
  
  original_links   <- P_original[cbind(test_rows, test_cols)]
  predicted_values <- P_reconstructed[cbind(test_rows, test_cols)]
  
  results <- data.frame(
    k = k,
    lambda = lambda,
    original_links = original_links,
    predicted_values = predicted_values,
    node_to = test_rows,
    node_from = test_cols,
    removed = 1
  )
  
  # Enumerate all edges
  all_edges <- expand.grid(
    node_to = rownames(P),
    node_from = colnames(P),
    k = k,
    lambda = lambda,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  all_edges$original_links <- mapply(
    function(r, c) P_original[r, c],
    all_edges$node_to,
    all_edges$node_from
  )
  
  removed_edges_idx <- data.frame(
    node_to = test_rows,
    node_from = test_cols,
    stringsAsFactors = FALSE
  )
  
  not_removed <- all_edges[
    !paste(all_edges$node_to, all_edges$node_from) %in%
      paste(removed_edges_idx$node_to, removed_edges_idx$node_from),
  ]
  not_removed$removed <- 0
  not_removed$predicted_values <- NA
  not_removed$k <- k
  not_removed$lambda <- lambda
  
  list(results = results, not_removed = not_removed)
  
}

sigmoid <- function(x) 1 / (1 + exp(-x))



# =========================
# 3) Choose ONE train layer and ONE predict layer
# =========================
layers_to_train  <- "1"   # <-- change to your layer id
layer_to_predict <- "1"   # <-- change to your layer id

A <- build_interaction_matrix(aggregated_df, layers_to_train)
P <- build_interaction_matrix(aggregated_df, layer_to_predict)

cat("Raw dims: A =", paste(dim(A), collapse="x"),
    " | P =", paste(dim(P), collapse="x"), "\n")

# =========================
# # 4) Align A and P to the same row/col universe (IMPORTANT for A+P)
# # =========================
# all_row_ids <- sort(unique(c(rownames(A_raw), rownames(P_raw))))
# all_col_ids <- sort(unique(c(colnames(A_raw), colnames(P_raw))))
# 
# A <- matrix(0, nrow = length(all_row_ids), ncol = length(all_col_ids),
#             dimnames = list(all_row_ids, all_col_ids))
# P <- matrix(0, nrow = length(all_row_ids), ncol = length(all_col_ids),
#             dimnames = list(all_row_ids, all_col_ids))
# 
# A[rownames(A_raw), colnames(A_raw)] <- A_raw
# P[rownames(P_raw), colnames(P_raw)] <- P_raw
# 
# # Save original P before masking (your implement_impute uses this)
# P_original <- P

# =========================
# 5) Create fold assignment (5-fold) and run only fold = 1
#    Hold out ~20% of ALL 1s and ~20% of ALL 0s (disjoint folds)
# =========================
set.seed(42)
K <- 5
fold <- 1

# ---- 1) Make A and P binary ----
A <- (A > 0) * 1
P <- (P > 0) * 1
P_original <- P  # keep original before masking

# ---- 2) Get ALL positive and ALL negative indices in P ----
ones_in_P  <- which(P_original == 1, arr.ind = TRUE)
zeros_in_P <- which(P_original == 0, arr.ind = TRUE)

n_pos  <- nrow(ones_in_P)
n_zero <- nrow(zeros_in_P)

cat("\n===== TOTAL COUNTS IN P =====\n")
cat("Total 1s:", n_pos, "\n")
cat("Total 0s:", n_zero, "\n")
cat("Total edges:", n_pos + n_zero, "\n")

# ---- 3) Balanced fold IDs (sizes differ by at most 1) ----
make_folds <- function(n, K) sample(rep(seq_len(K), length.out = n))

pos_folds  <- make_folds(n_pos,  K)   # assigns EACH 1-edge to exactly one fold
zero_folds <- make_folds(n_zero, K)   # assigns EACH 0-edge to exactly one fold

cat("\n===== PER-FOLD COUNTS =====\n")
cat("1s per fold:\n"); print(table(pos_folds))
cat("0s per fold:\n"); print(table(zero_folds))

# ---- 4) Select held-out edges for the chosen fold (here fold=1) ----
remove_pos  <- ones_in_P[pos_folds == fold,  , drop = FALSE]
remove_zero <- zeros_in_P[zero_folds == fold, , drop = FALSE]

cat("\n===== HOLDOUT (fold 1) =====\n")
cat("Holdout 1s:", nrow(remove_pos), " (~", round(100*nrow(remove_pos)/n_pos, 2), "% )\n")
cat("Holdout 0s:", nrow(remove_zero), " (~", round(100*nrow(remove_zero)/n_zero, 2), "% )\n")

# ---- 5) Mask them in P ----
P_masked <- P_original
P_masked[remove_pos]  <- NA
P_masked[remove_zero] <- NA

# IMPORTANT: variables your implement_impute() expects
P <- P_masked
remove_indices <- remove_pos
zeros_to_remove_indices <- remove_zero

# =========================
# 6) Build C and run ONE (k, lambda) (no loops)
# =========================
# If you truly want A+P (your pipeline):
C <- A
C[rownames(P), colnames(P)] <- ifelse(is.na(P), NA, C[rownames(P), colnames(P)] + P)

# Centering like your pipeline
C <- biScale(C, row.center = TRUE, col.center = TRUE,
             row.scale = FALSE, col.scale = FALSE)

# Pick ONE k and ONE lambda
k <- 2
lambda <- lambda0(C)   

r <- implement_impute(C, k, lambda)



final_removed <- data.frame(
  Study_id         = Study_id_chose,
  K_folder         = fold,
  train_layer      = layers_to_train,
  test_layer       = layer_to_predict,
  node_from        = r$results$node_from,
  node_to          = r$results$node_to,
  predicted_values = r$results$predicted_values,
  original_links         = r$results$original_links,
  removed          = r$results$removed,     # 1
  lambda           = r$results$lambda,
  k                = r$results$k,
  stringsAsFactors = FALSE
)

final_not_removed <- data.frame(
  Study_id         = Study_id_chose,
  K_folder         = fold,
  train_layer      = layers_to_train,
  test_layer       = layer_to_predict,
  node_from        = r$not_removed$node_from,
  node_to          = r$not_removed$node_to,
  predicted_values = r$not_removed$predicted_values, # NA in your function
  original_links   = r$not_removed$original_links,
  removed          = r$not_removed$removed,          # 0
  lambda           = r$not_removed$lambda,
  k                = r$not_removed$k,
  stringsAsFactors = FALSE
)

final_df <- rbind(final_removed, final_not_removed)

# sanity check
table(final_df$removed)


df <- final_df %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values))

plant_species <- unique(df$node_from)           # get unique names
pollinator_species <- unique(df$node_to)
length(plant_species)
length(pollinator_species)


thresholds <- seq(0, 1, by = 0.1)

# 1) filter & prep
df_prepped <- df %>%
  filter(removed == 1) %>%
  mutate(
    predicted_prob   = sigmoid(predicted_values),
    original_binary  = if_else(original_links > 0, 1, 0)
  )

thresholds <- seq(0, 1, by = 0.1)

# 1) filter & prep
df_prepped <- df %>%
  filter(removed == 1) %>%
  mutate(
    predicted_prob   = sigmoid(predicted_values),
    original_binary  = if_else(original_links > 0, 1, 0)
  )




#--------
# 2) expand to one row per threshold
df_thresh <- df_prepped %>%
  tidyr::expand_grid(threshold = thresholds) %>%  # <-- switch here
  mutate(
    predicted_bin = if_else(predicted_prob > threshold, 1, 0)
  ) %>%
  group_by( train_layer, test_layer, K_folder, threshold) %>%
  summarise(
    TP = sum(original_binary == 1 & predicted_bin == 1),
    FN = sum(original_binary == 1 & predicted_bin == 0),
    TN = sum(original_binary == 0 & predicted_bin == 0),
    FP = sum(original_binary == 0 & predicted_bin == 1),
    specificity      = TN / (TN + FP),
    precision        = TP / (TP + FP),
    recall           = TP / (TP + FN),
    f1_score         = 2 * (precision * recall) / (precision + recall),
    balanced_accuracy= (recall + specificity) / 2,
    mcc = (TP * TN - FP * FN) /
      sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)),
    mse  = mean((predicted_values - original_links)^2, na.rm = TRUE),
    rmse = sqrt(mse)#,
    #.groups = "drop"
  ) %>%
  ungroup() %>%
  group_by(train_layer, test_layer, threshold) %>%
  summarise(
    TP = mean(TP, na.rm = TRUE),
    FN = mean(FN, na.rm = TRUE),
    TN = mean(TN, na.rm = TRUE),
    FP = mean(FP, na.rm = TRUE),
    specificity = mean(specificity, na.rm = TRUE),
    precision = mean(precision, na.rm = TRUE),
    recall = mean(recall, na.rm = TRUE),
    f1_score = mean(f1_score, na.rm = TRUE),
    balanced_accuracy = mean(balanced_accuracy, na.rm = TRUE),
    mcc = mean(mcc, na.rm = TRUE),
    mse = mean(mse, na.rm = TRUE),
    rmse = mean(rmse, na.rm = TRUE)
  ) %>%
  ungroup()
