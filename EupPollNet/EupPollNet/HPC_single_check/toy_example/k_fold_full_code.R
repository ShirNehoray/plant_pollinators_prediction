# =========================
# 0) Libraries
# =========================
packages <- c("dplyr", "softImpute", "pROC")
missing <- setdiff(packages, rownames(installed.packages()))
if (length(missing)) install.packages(missing, dependencies = TRUE)
invisible(lapply(packages, library, character.only = TRUE))

## ---- parameters ----
Study_id_chose <- "1_Bartomeus"


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
# 1) Functions (same spirit as HPC)
# =========================

build_interaction_matrix <- function(data, layers_to_filter) {
  layers <- paste0("layer_", layers_to_filter)
  filtered_data <- subset(data, layer_from %in% layers)
  
  aggregated_data <- filtered_data %>%
    dplyr::group_by(node_from, node_to) %>%
    dplyr::summarise(weight = sum(weight), .groups = "drop")
  
  species_from <- unique(aggregated_data$node_from)  # cols
  species_to   <- unique(aggregated_data$node_to)    # rows
  
  interaction_matrix <- matrix(
    0,
    nrow = length(species_to),
    ncol = length(species_from),
    dimnames = list(species_to, species_from)
  )
  
  for (i in 1:nrow(aggregated_data)) {
    r <- aggregated_data$node_to[i]
    c <- aggregated_data$node_from[i]
    interaction_matrix[r, c] <- aggregated_data$weight[i]
  }
  
  interaction_matrix
}

# IMPORTANT: implement_impute uses globals: P, P_original, remove_indices, zeros_to_remove_indices
implement_impute <- function(C, k, lambda) {
  
  fit <- softImpute::softImpute(C, rank.max = k, lambda = lambda, type = "svd", maxit = 600)
  C_reconstructed <- softImpute::complete(C, fit)
  
  P_reconstructed <- C_reconstructed[rownames(P), colnames(P)]
  
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
  

  original_links <- P_original_weighted[cbind(test_rows, test_cols)]
  
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
  
  removed_edges_idx <- data.frame(node_to = test_rows, node_from = test_cols, stringsAsFactors = FALSE)
  
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

make_folds <- function(n, K) sample(rep(seq_len(K), length.out = n))

# =========================
# 2) USER PARAMETERS
# =========================
K <- 5
set.seed(42)

# same as your HPC ranks / lambdas
k_values <- c(2, 5, 10)

# pick lambda grid. (lambda0 depends on C; we append it inside the loop)
lambda_grid_fixed <- c(1, 5, 50, 100)

# if TRUE: make matrices binary before CV
binary_network <- TRUE

# =========================
# 3) INPUT: you already created aggregated_df for ONE study (as in your code)
#    aggregated_df must have: layer_from, node_from, node_to, weight
#
#    If instead you want multiple studies: wrap the entire block below in a loop over Study_id,
#    and set aggregated_df each time like your HPC script does.
# =========================

# Example placeholders (make sure these exist in your session)
# Study_id_chose <- "1_Bartomeus"
# aggregated_df <- ...

# total layers in this study dataset:
layer_names <- sort(unique(aggregated_df$layer_from))
num_layers <- length(layer_names)

cat("\n=============================\n")
cat("Study:", Study_id_chose, "\n")
cat("Layers found:", num_layers, "\n")
cat("=============================\n")

combined_results <- data.frame()

# =========================
# 4) MAIN LOOP: each layer predicts itself
# =========================
for (layer_id in 1:num_layers) {
  
  cat("\n------------------------------------------\n")
  cat("Self-prediction on layer:", layer_id, "\n")
  cat("------------------------------------------\n")
  
  # Build P (and A if you still want A+P; for self-to-self you can set A=P)
  P0 <- build_interaction_matrix(aggregated_df, layers_to_filter = layer_id)
  # Build weighted matrix
  P0_weighted <- build_interaction_matrix(aggregated_df, layers_to_filter = layer_id)

  # Save weighted version forever
  P_original_weighted <- P0_weighted

  # Binarize if requested
  if (binary_network) {
    P0 <- (P0_weighted > 0) * 1
  } else {
    P0 <- P0_weighted
  }

  # Binary matrix used for CV & evaluation
  P_original <- P0

  # Optional: for self-to-self, A can simply be P0 (or you can keep A separate if you prefer)
  A0 <- P0
  
  # Indices of all 1s and all 0s
  ones_in_P  <- which(P_original == 1, arr.ind = TRUE)
  zeros_in_P <- which(P_original == 0, arr.ind = TRUE)
  
  n_pos  <- nrow(ones_in_P)
  n_zero <- nrow(zeros_in_P)
  
  cat("Matrix dim:", paste(dim(P_original), collapse = "x"),
      "| #1s:", n_pos, "| #0s:", n_zero, "\n")
  
  if (n_pos < K) {
    cat("Skipping layer", layer_id, "- too few positives for", K, "folds.\n")
    next
  }
  
  # Disjoint fold assignment: covers ALL 1s and ALL 0s exactly once over K folds
  pos_folds  <- make_folds(n_pos,  K)
  zero_folds <- make_folds(n_zero, K)
  
  cat("1s per fold: "); print(table(pos_folds))
  cat("0s per fold: "); print(table(zero_folds))
  
  # ---- fold loop ----
  for (fold in 1:K) {
    
    cat(sprintf(">>> START fold %d/%d | Study=%s | Layer=%d\n",
                fold, K, as.character(Study_id_chose), layer_id))
    
    # Choose held-out sets (20%-ish each fold; disjoint)
    remove_pos  <- ones_in_P[pos_folds  == fold, , drop = FALSE]
    remove_zero <- zeros_in_P[zero_folds == fold, , drop = FALSE]
    
    # Mask P
    P <- P_original
    P[remove_pos]  <- NA
    P[remove_zero] <- NA
    
    # globals expected by implement_impute()
    remove_indices <- remove_pos
    zeros_to_remove_indices <- remove_zero
    
    # ---- Build combined C (HPC-style) ----
    # For self-prediction, A and P share the same universe, but we keep the same structure.
    all_row_ids <- unique(c(rownames(A0), rownames(P)))
    all_col_ids <- unique(c(colnames(A0), colnames(P)))
    
    C <- matrix(0,
                nrow = length(all_row_ids),
                ncol = length(all_col_ids),
                dimnames = list(all_row_ids, all_col_ids))
    
    # place A
    C[rownames(A0), colnames(A0)] <- A0
    
    # place P while preserving NA positions from P (THIS is the correct masking logic)
    subC <- C[rownames(P), colnames(P)]
    C[rownames(P), colnames(P)] <- ifelse(is.na(P), NA, subC + P)
    
    # SAFETY: ensure C is really a matrix before calling lambda0(C)
    stopifnot(is.matrix(C))
    storage.mode(C) <- "double"
    
    # center
    C <- softImpute::biScale(C, row.center = TRUE, col.center = TRUE,
                             row.scale = FALSE, col.scale = FALSE)
    
    cat("NAs in C:", sum(is.na(C)), "\n")
    
    # lambda values (append lambda0)
    lam0 <- softImpute::lambda0(C)
    lambda_values <- c(lambda_grid_fixed, lam0)
    
    # ---- k/lambda loop (same as HPC) ----
    results <- data.frame()
    not_removed_all <- NULL
    
    for (k in k_values) {
      for (lambda in lambda_values) {
        
        r <- implement_impute(C, k, lambda)
        
        r$results$input_lambda <- lambda
        r$not_removed$input_lambda <- lambda
        
        results <- rbind(results, r$results)
        not_removed_all <- rbind(not_removed_all, r$not_removed)
      }
    }
    
    complete_edges_all <- rbind(results, not_removed_all)
    complete_edges_all$itr <- fold
    complete_edges_all$cv_type <- "Kfold_edges"
    
    combined_results <- rbind(
      combined_results,
      cbind(
        data.frame(
          Study_id = Study_id_chose,
          train_layer = layer_id,
          test_layer  = layer_id,
          K = K,
          fold = fold
        ),
        complete_edges_all
      )
    )
  }
}

cat("\nDONE. combined_results rows:", nrow(combined_results), "\n")

# =========================
# 5) OPTIONAL: Make your final "CSS-style" table (removed edges only)
# =========================


df <- combined_results %>%
  mutate(predicted_values = if_else(predicted_values < 0, 0, predicted_values))

# 1) filter & prep
df_prepped <- df %>%
  filter(removed == 1) %>%
  mutate(
    predicted_prob   = sigmoid(predicted_values),
    original_binary  = if_else(original_links > 0, 1, 0)
  )

final_df <- df_prepped %>%
  dplyr::filter(removed == 1) %>%
  dplyr::transmute(
    Study_id         = Study_id,
    K_folder         = fold,
    train_layer      = train_layer,
    test_layer       = test_layer,
    node_from        = node_from,
    node_to          = node_to,
    predicted_values = predicted_prob,
    original_links   = original_links,
    original_binary  = original_binary,
    removed          = removed,   # always 1 here
    lambda           = lambda,
    k                = k
  )

#few checks
all(final_df$train_layer == final_df$test_layer) #see that train layer and test are the same 

write.csv(final_df, "fina_df_one_study.csv")
