# --- Load required packages ---
library(dplyr)
library(tidyr)
library(ggplot2)

# --- Step 1: Read the dataset ---
# Make sure "Interaction_data.csv" is in your R working directory.
data <- read.csv("Interaction_data.csv")

# --- Step 2: Define bioregions and initialize a list ---
bioregions <- c("Alpine", "Atlantic", "Boreal", "Continental")
metaweb_list <- list()

# --- Step 3: Loop through each bioregion to create metawebs ---
for (bio in bioregions) {
  # Filter rows for this bioregion
  df_bio <- data %>%
    filter(Bioregion == bio) %>%
    select(Plant_accepted_name, Pollinator_accepted_name, Interaction)
  
  # Sum interactions for each unique plant-pollinator pair
  df_sum <- df_bio %>%
    group_by(Plant_accepted_name, Pollinator_accepted_name) %>%
    summarise(total_interactions = sum(Interaction, na.rm = TRUE), .groups = 'drop')
  
  # Create an interaction matrix
  meta_matrix <- df_sum %>%
    pivot_wider(names_from = Pollinator_accepted_name,
                values_from = total_interactions,
                values_fill = 0)
  
  # Turn the matrix into a data frame
  meta_df <- as.data.frame(meta_matrix)
  
  # Store both versions in the list
  metaweb_list[[bio]] <- list(
    matrix = meta_matrix,
    dataframe = meta_df
  )
}

# --- Step 4: Prepare a list of unique interaction pairs for each bioregion ---
interaction_sets <- list()

for (bio in names(metaweb_list)) {
  df <- metaweb_list[[bio]]$dataframe
  
  # Convert matrix back to a long format to get pairs
  df_long <- df %>%
    pivot_longer(-Plant_accepted_name, names_to = "Pollinator_accepted_name", values_to = "value") %>%
    filter(value > 0) %>%
    mutate(pair = paste(Plant_accepted_name, Pollinator_accepted_name, sep = "_"))
  
  interaction_sets[[bio]] <- unique(df_long$pair)
}

# --- Step 5: Calculate pairwise Jaccard similarities ---
bioregion_names <- names(interaction_sets)
n <- length(bioregion_names)
jaccard_matrix <- matrix(NA, nrow = n, ncol = n,
                         dimnames = list(bioregion_names, bioregion_names))

for (i in 1:n) {
  for (j in 1:n) {
    set_i <- interaction_sets[[i]]
    set_j <- interaction_sets[[j]]
    
    inter <- length(intersect(set_i, set_j))
    union_val <- length(union(set_i, set_j))
    
    # Avoid division by zero if a bioregion has no interactions
    if (union_val == 0) {
      jaccard_matrix[i, j] <- 0
    } else {
      jaccard_matrix[i, j] <- inter / union_val
    }
  }
}

# --- Step 6: Convert to a tidy data frame for plotting ---
jaccard_df <- as.data.frame(jaccard_matrix) %>%
  mutate(Bioregion1 = rownames(.)) %>%
  pivot_longer(-Bioregion1, names_to = "Bioregion2", values_to = "Jaccard")

# --- Step 7: Plot the heatmap with Jaccard values displayed ---
ggplot(jaccard_df, aes(x = Bioregion1, y = Bioregion2, fill = Jaccard)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Jaccard, 2)), color = "black", size = 4) +
scale_fill_viridis_c(option = "viridis", name = "Jaccard\nSimilarity", limits = c(0, 1)) +
  theme_minimal() +
  labs(title = "Jaccard Similarity of Plant-Pollinator Metawebs",
       x = "Bioregion", y = "Bioregion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))


