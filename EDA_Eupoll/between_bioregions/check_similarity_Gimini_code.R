# Load packages
library(dplyr)
library(ggplot2)
library(tidyr)

# --- Step 1: Load Data ---

# Read the interaction data
data <- read.csv("Interaction_data.csv", encoding = "latin1")

# --- Step 2: Check Required Columns ---

# Check for essential columns used in the script
required_cols <- c(
  "Bioregion", "Plant_accepted_name", "Pollinator_family",
  "Interaction", "Pollinator_accepted_name", "EuPPollNet_habitat"
)

# --- Step 3: Aggregate Data by Bioregion ---

# for Apidae

# Get all plants that interact with Apidae
cat("Identifying plants that interact with Apidae...\n")
plants_with_apidae <- data %>%
  filter(Pollinator_family == "Apidae", Interaction > 0) %>%
  select(Plant_accepted_name) %>%
  distinct() %>%
  pull(Plant_accepted_name) # Get a vector of plant names

cat(sprintf("Found %d unique plant species that interact with Apidae.\n", length(plants_with_apidae)))

# Aggregate only plants that interact with Apidae, grouped by Bioregion
plant_species_per_bioregion <- data %>%
  filter(Plant_accepted_name %in% plants_with_apidae) %>% # Filter for relevant plants 
  group_by(Bioregion) %>%
  summarise(
    Plants = list(unique(Plant_accepted_name)),
    .groups = "drop"
  )

# Filter for Apidae and aggregate pollinator species by Bioregion
filtered_pollinator_species_per_bioregion <- data %>%
  filter(Pollinator_family == "Apidae", Interaction > 0) %>% #to insure that there is an interaction
  group_by(Bioregion) %>%
  summarise(
    Pollinators = list(unique(Pollinator_accepted_name)),
    .groups = "drop"
  )

# Create a data frame with all unique bioregions
all_bioregions_df <- data %>%
  select(Bioregion) %>%
  distinct()

# Join the plant and pollinator data
species_per_bioregion <- all_bioregions_df %>%
  left_join(plant_species_per_bioregion, by = "Bioregion") %>%
  left_join(filtered_pollinator_species_per_bioregion, by = "Bioregion")

# Replace NA/NULL list entries with empty lists for the Jaccard calculation
for (i in 1:nrow(species_per_bioregion)) {
  if (is.null(species_per_bioregion$Pollinators[[i]])) {
    # Use list() to represent an empty set
    species_per_bioregion$Pollinators[i] <- list(list())
  }
  # Add check for plants
  if (is.null(species_per_bioregion$Plants[[i]])) {
    species_per_bioregion$Plants[i] <- list(list())
  }
}


# --- Step 4: Aggregate Data by Habitat  ---

# Aggregate only plants that interact with Apidae, grouped by Habitat
plant_species_per_habitat <- data %>%
  filter(Plant_accepted_name %in% plants_with_apidae) %>% # Filter for relevant plants
  group_by(EuPPollNet_habitat) %>%
  summarise(
    Plants = list(unique(Plant_accepted_name)),
    .groups = "drop"
  )

# Filter for Apidae with interaction > 0 and aggregate pollinator species by Habitat
filtered_pollinator_species_per_habitat <- data %>%
  filter(Pollinator_family == "Apidae", Interaction > 0) %>%
  group_by(EuPPollNet_habitat) %>%
  summarise(
    Pollinators = list(unique(Pollinator_accepted_name)),
    .groups = "drop"
  )

# Create a data frame with all unique habitats
all_habitats_df <- data %>%
  select(EuPPollNet_habitat) %>%
  distinct()

# Join the plant and pollinator data (per habitat)
species_per_habitat <- all_habitats_df %>%
  left_join(plant_species_per_habitat, by = "EuPPollNet_habitat") %>%
  left_join(filtered_pollinator_species_per_habitat, by = "EuPPollNet_habitat")

# Replace NA/NULL list entries with empty lists for the Jaccard calculation
for (i in 1:nrow(species_per_habitat)) {
  if (is.null(species_per_habitat$Pollinators[[i]])) {
    species_per_habitat$Pollinators[i] <- list(list())
  }
  if (is.null(species_per_habitat$Plants[[i]])) { # Also add check for plants, just in case
    species_per_habitat$Plants[i] <- list(list())
  }
}


# --- Step 5: Similarity Functions ---

#  Jaccard similarity function
jaccard_similarity <- function(set1, set2) {
  # Find the number of species common to both sets
  intersection_size <- length(intersect(set1, set2))
  # Find the total number of unique species across both sets
  union_size <- length(union(set1, set2))
  
  # Avoid division by zero if both sets are empty
  if (union_size == 0) {
    # If both sets are empty, they are 100% similar, but Jaccard is 0/0.
    # Let's return 1.0 if intersection is 0 and union is 0.
    if(intersection_size == 0) return(1.0) 
    return(NA) # Should not be reached if union_size is 0
  } else {
    return(intersection_size / union_size)
  }
}

# create a similarity matrix
create_similarity_matrix <- function(species_data, species_type, grouping_col) {
  num_groups <- nrow(species_data)
  group_names <- species_data[[grouping_col]]
  
  # Create an empty matrix to store the similarity scores
  sim_matrix <- matrix(NA,
                       nrow = num_groups,
                       ncol = num_groups,
                       dimnames = list(group_names, group_names))
  
  # Loop through each pair of groups to calculate similarity
  for (i in 1:num_groups) {
    for (j in 1:num_groups) {
      # Get the species lists for the two groups
      species_set1 <- unlist(species_data[[species_type]][i])
      species_set2 <- unlist(species_data[[species_type]][j])
      
      # Calculate and store the Jaccard similarity score
      sim_matrix[i, j] <- jaccard_similarity(species_set1, species_set2)
    }
  }
  return(sim_matrix)
}

# --- Function to create a cross-similarity matrix (e.g., Bioregion 1 vs Bioregion 2) ---
create_cross_similarity_matrix <- function(data_group1, data_group2, species_type, grouping_col) {
  
  # Get the names of the habitats in each group
  group1_names <- data_group1[[grouping_col]]
  group2_names <- data_group2[[grouping_col]]
  
  # Create an empty matrix
  sim_matrix <- matrix(NA,
                       nrow = length(group1_names),
                       ncol = length(group2_names),
                       dimnames = list(group1_names, group2_names))
  
  # Loop through each habitat in group 1
  for (i in 1:length(group1_names)) {
    # Loop through each habitat in group 2
    for (j in 1:length(group2_names)) {
      # Get the species lists for the two habitats
      species_set1 <- unlist(data_group1[[species_type]][i])
      species_set2 <- unlist(data_group2[[species_type]][j])
      
      # Calculate and store the Jaccard similarity score
      sim_matrix[i, j] <- jaccard_similarity(species_set1, species_set2)
    }
  }
  return(sim_matrix)
}

# --- Plotting Functions  ---
# Make a function to plot a heatmap
plot_heatmap <- function(similarity_matrix, title, x_label, y_label) {
  
  # Convert the matrix to a long-format data frame, which is required for ggplot2
  long_format_data <- as.data.frame(as.table(similarity_matrix))
  # Use generic names for the columns
  names(long_format_data) <- c("Group1", "Group2", "Similarity")
  
  # Create the heatmap plot
  heatmap_plot <- ggplot(long_format_data, aes(x = Group1, y = Group2, fill = Similarity)) +
    geom_tile(color = "white") + # Add white lines between tiles for better separation
    # Add text labels for the Jaccard index values, rounded to 2 decimal places
    geom_text(aes(label = round(Similarity, 2)), color = "black", size = 3) +
    # Use the 'viridis' color palette as requested
    scale_fill_viridis_c(option = "viridis", name = "Jaccard\nSimilarity", na.value = "grey80") + # Added na.value
    labs(
      title = title,
      x = x_label,
      y = y_label
    ) +
    theme_minimal() + # Use a clean, minimal theme
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10), # Rotate x-axis labels
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), # Center the title
      panel.grid.major = element_blank(), # Remove grid lines
      panel.grid.minor = element_blank()
    ) +
    coord_fixed() # Ensure the tiles are square
  
  return(heatmap_plot)
}


# --- Step 6: Create Matrices ---

# Create for BIOREGION comparison 
plant_similarity_matrix_bioregion <- create_similarity_matrix(species_per_bioregion, "Plants", "Bioregion")
pollinator_similarity_matrix_bioregion <- create_similarity_matrix(species_per_bioregion, "Pollinators", "Bioregion")

# Create for HABITAT comparison
plant_similarity_matrix_habitat <- create_similarity_matrix(species_per_habitat, "Plants", "EuPPollNet_habitat")
pollinator_similarity_matrix_habitat <- create_similarity_matrix(species_per_habitat, "Pollinators", "EuPPollNet_habitat")


# --- Step 7: Find and Print Most Similar Bioregions ---

cat("\n--- Analysis of Most Similar Bioregions ---\n")
# Initialize variables for the new plots
pollinator_cross_heatmap <- NULL
plant_cross_heatmap <- NULL

# For Apidae Pollinators
pollinator_matrix_copy <- pollinator_similarity_matrix_bioregion
# Set the diagonal to NA (to exclude self-similarity, which is always 1)
diag(pollinator_matrix_copy) <- NA

# Find the maximum similarity value
max_pollinator_sim <- max(pollinator_matrix_copy, na.rm = TRUE)

# Find the row indices of this max value
max_pollinator_indices <- which(pollinator_matrix_copy == max_pollinator_sim, arr.ind = TRUE)

# Take the first pair, since the matrix is symmetrical
if(nrow(max_pollinator_indices) > 0) {
  bio1_poll_name <- rownames(pollinator_matrix_copy)[max_pollinator_indices[1, 1]]
  bio2_poll_name <- colnames(pollinator_matrix_copy)[max_pollinator_indices[1, 2]]
  
  cat(sprintf("Highest Apidae Pollinator Similarity (Bioregions):\n"))
  cat(sprintf("  Bioregion 1: %s\n", bio1_poll_name))
  cat(sprintf("  Bioregion 2: %s\n", bio2_poll_name))
  cat(sprintf("  Jaccard Similarity: %.4f\n", max_pollinator_sim))
  
  cat(sprintf("\n  Comparison of these two bioregions (Apidae Pollinators):\n"))
  
  # Get the species lists for these two bioregions
  species_list_1 <- unlist(species_per_bioregion$Pollinators[species_per_bioregion$Bioregion == bio1_poll_name])
  species_list_2 <- unlist(species_per_bioregion$Pollinators[species_per_bioregion$Bioregion == bio2_poll_name])
  
  # Calculate intersection (shared species)
  shared_species <- intersect(species_list_1, species_list_2)
  shared_count <- length(shared_species)
  
  # Calculate union (total unique species)
  total_species <- union(species_list_1, species_list_2)
  total_count <- length(total_species)
  
  # # Print the comparison details
  # cat(sprintf("    Total unique Apidae species in %s: %d\n", bio1_poll_name, length(species_list_1)))
  # cat(sprintf("    Total unique Apidae species in %s: %d\n", bio2_poll_name, length(species_list_2)))
  # cat(sprintf("    Shared Apidae species: %d\n", shared_count))
  # cat(sprintf("    Total unique Apidae species (Union): %d\n", total_count))
  # 
  # Verify the calculation
  if(total_count > 0) {
    cat(sprintf("    Calculated Jaccard (Shared/Union): %.4f\n", shared_count / total_count))
  } else {
    cat("    Calculated Jaccard (Shared/Union): NA (no species)\n")
  }
  # 
  # # Print the actual list of shared species, if any
  # if(shared_count > 0) {
  #   cat("    Shared Species List:\n")
  #   # Use paste() to format the list nicely
  #   cat(paste("      -", shared_species, collapse = "\n"))
  #   cat("\n\n")
  # } else {
  #   cat("    No shared Apidae species found between these two regions.\n\n")
  # }
  
  # --- CROSS-HABITAT HEATMAP (POLLINATORS) ---
  cat(sprintf("  Generating Habitat Cross-Similarity Heatmap for %s vs %s (Apidae Pollinators)...\n", bio1_poll_name, bio2_poll_name))
  
  # Get species lists grouped by habitat, only for these two bioregions
  data_for_poll_cross_sim <- data %>%
    filter(Bioregion %in% c(bio1_poll_name, bio2_poll_name), Pollinator_family == "Apidae", Interaction > 0) %>%
    group_by(Bioregion, EuPPollNet_habitat) %>%
    summarise(Pollinators = list(unique(Pollinator_accepted_name)), .groups = "drop")
  
  # Separate the data into two groups
  data_group1_poll <- data_for_poll_cross_sim %>% filter(Bioregion == bio1_poll_name)
  data_group2_poll <- data_for_poll_cross_sim %>% filter(Bioregion == bio2_poll_name)
  
  # Only create plot if both bioregions have data
  if (nrow(data_group1_poll) > 0 && nrow(data_group2_poll) > 0) {
    # Create the cross-similarity matrix
    pollinator_cross_sim_matrix <- create_cross_similarity_matrix(data_group1_poll, data_group2_poll, "Pollinators", "EuPPollNet_habitat")
    
    # Generate the heatmap
    pollinator_cross_heatmap <- plot_heatmap(
      pollinator_cross_sim_matrix,
      paste("Apidae Habitat Similarity:", bio1_poll_name, "vs.", bio2_poll_name),
      paste("Habitats in", bio2_poll_name),
      paste("Habitats in", bio1_poll_name)
    )
  } else {
    cat(sprintf("    Skipping cross-habitat plot: not enough data in one or both bioregions.\n\n"))
  }
  
  
} else {
  cat("Could not find max pollinator similarity pair for Bioregions.\n\n")
}

# --- For Plants  ---
# Make a copy of the plant matrix
plant_matrix_copy <- plant_similarity_matrix_bioregion
diag(plant_matrix_copy) <- NA

max_plant_sim <- max(plant_matrix_copy, na.rm = TRUE)
max_plant_indices <- which(plant_matrix_copy == max_plant_sim, arr.ind = TRUE)

if(nrow(max_plant_indices) > 0) {
  bio1_plant_name <- rownames(plant_matrix_copy)[max_plant_indices[1, 1]]
  bio2_plant_name <- colnames(plant_matrix_copy)[max_plant_indices[1, 2]]
  
  cat(sprintf("Highest Plant (Apidae-Interacting) Species Similarity (Bioregions):\n"))
  cat(sprintf("  Bioregion 1: %s\n", bio1_plant_name))
  cat(sprintf("  Bioregion 2: %s\n", bio2_plant_name))
  cat(sprintf("  Jaccard Similarity: %.4f\n", max_plant_sim))
  
  # --- NEW DETAILED COMPARISON (PLANTS) ---
  cat(sprintf("\n  Comparison of these two bioregions (Plants - Apidae-Interacting):\n"))
  
  # Get the species lists for these two bioregions
  species_list_1 <- unlist(species_per_bioregion$Plants[species_per_bioregion$Bioregion == bio1_plant_name])
  species_list_2 <- unlist(species_per_bioregion$Plants[species_per_bioregion$Bioregion == bio2_plant_name])
  
  # Calculate intersection (shared species)
  shared_species <- intersect(species_list_1, species_list_2)
  shared_count <- length(shared_species)
  
  # Calculate union (total unique species)
  total_species <- union(species_list_1, species_list_2)
  total_count <- length(total_species)
  
  # Print the comparison details
  cat(sprintf("    Total unique Plant (Apidae-Interacting) species in %s: %d\n", bio1_plant_name, length(species_list_1)))
  cat(sprintf("    Total unique Plant (Apidae-Interacting) species in %s: %d\n", bio2_plant_name, length(species_list_2)))
  cat(sprintf("    Shared Plant (Apidae-Interacting) species: %d\n", shared_count))
  cat(sprintf("    Total unique Plant (Apidae-Interacting) species (Union): %d\n", total_count))
  
  # Verify the calculation
  if(total_count > 0) {
    cat(sprintf("    Calculated Jaccard (Shared/Union): %.4f\n", shared_count / total_count))
  } else {
    cat("    Calculated Jaccard (Shared/Union): NA (no species)\n")
  }
  
  # Print the actual list of shared species, if any
  if(shared_count > 0) {
    cat("    Shared Species List:\n")
    # Use paste() to format the list nicely
    cat(paste("      -", shared_species, collapse = "\n"))
    cat("\n\n")
  } else {
    cat("    No shared Plant (Apidae-Interacting) species found between these two regions.\n\n")
  }
  
  # --- NEW: CROSS-HABITAT HEATMAP (PLANTS) ---
  cat(sprintf("  Generating Habitat Cross-Similarity Heatmap for %s vs %s (Plants - Apidae-Interacting)...\n", bio1_plant_name, bio2_plant_name))
  
  # Get species lists grouped by habitat, only for these two bioregions
  data_for_plant_cross_sim <- data %>%
    filter(Bioregion %in% c(bio1_plant_name, bio2_plant_name), 
           Plant_accepted_name %in% plants_with_apidae) %>% # Use the filter
    group_by(Bioregion, EuPPollNet_habitat) %>%
    summarise(Plants = list(unique(Plant_accepted_name)), .groups = "drop")
  
  # Separate the data into two groups
  data_group1_plant <- data_for_plant_cross_sim %>% filter(Bioregion == bio1_plant_name)
  data_group2_plant <- data_for_plant_cross_sim %>% filter(Bioregion == bio2_plant_name)
  
  # Only create plot if both bioregions have data
  if (nrow(data_group1_plant) > 0 && nrow(data_group2_plant) > 0) {
    # Create the cross-similarity matrix
    plant_cross_sim_matrix <- create_cross_similarity_matrix(data_group1_plant, data_group2_plant, "Plants", "EuPPollNet_habitat")
    
    # Generate the heatmap
    plant_cross_heatmap <- plot_heatmap(
      plant_cross_sim_matrix,
      paste("Plant (Apidae-Interacting) Habitat Similarity:", bio1_plant_name, "vs.", bio2_plant_name),
      paste("Habitats in", bio2_plant_name),
      paste("Habitats in", bio1_plant_name)
    )
  } else {
    cat(sprintf("    Skipping cross-habitat plot: not enough data in one or both bioregions.\n\n"))
  }
  
} else {
  cat("Could not find max plant similarity pair for Bioregions.\n\n")
}


# --- Step 7.5: Create Heatmaps ---
# cat("\n--- Generating Presence/Absence Heatmaps ---\n")
pa_heatmap_pollinators <- NULL
pa_heatmap_plants <- NULL

# --- For the Apidae Pollinator Pair ---
if(nrow(max_pollinator_indices) > 0) {
  cat(sprintf("  Generating Presence/Absence plot for: %s vs %s (Apidae)\n", bio1_poll_name, bio2_poll_name))
  
  # Get presence/absence data for just these two bioregions
  pa_data_poll <- data %>%
    filter(Bioregion %in% c(bio1_poll_name, bio2_poll_name), Pollinator_family == "Apidae", Interaction > 0) %>%
    select(Bioregion, Pollinator_accepted_name) %>%
    distinct() %>%
    mutate(Presence = 1)
  
  if (nrow(pa_data_poll) > 0) {
    pa_heatmap_pollinators <- ggplot(pa_data_poll, aes(x = Bioregion, y = Pollinator_accepted_name, fill = as.factor(Presence))) +
      geom_tile(color = "white") +
      scale_fill_manual(values = c("1" = "steelblue"), guide = "none") + # Single color for presence
      labs(
        title = paste("Apidae Species Presence:", bio1_poll_name, "vs.", bio2_poll_name),
        x = "Bioregion",
        y = "Pollinator Species"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
        axis.text.y = element_text(size = 5), # May need to be smaller if many species
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  } else {
    cat("    No Apidae data found for one or both bioregions to create Presence/Absence plot.\n")
  }
}

# --- For the Plant Pair ---
if(nrow(max_plant_indices) > 0) {
  cat(sprintf("  Generating Presence/Absence plot for: %s vs %s (Plants - Apidae-Interacting)\n", bio1_plant_name, bio2_plant_name))
  
  # Get presence/absence data for just these two bioregions
  pa_data_plant <- data %>%
    filter(Bioregion %in% c(bio1_plant_name, bio2_plant_name), 
           Plant_accepted_name %in% plants_with_apidae) %>% # Keep the Apidae-interacting filter
    select(Bioregion, Plant_accepted_name) %>%
    distinct() %>%
    mutate(Presence = 1)
  
  if(nrow(pa_data_plant) > 0) {
    pa_heatmap_plants <- ggplot(pa_data_plant, aes(x = Bioregion, y = Plant_accepted_name, fill = as.factor(Presence))) +
      geom_tile(color = "darkgreen") + # Use color for the tile border
      scale_fill_manual(values = c("1" = "lightgreen"), guide = "none") + # Single color for presence
      labs(
        title = paste("Plant (Apidae-Interacting) Species Presence:", bio1_plant_name, "vs.", bio2_plant_name),
        x = "Bioregion",
        y = "Plant Species"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 16),
        axis.text.y = element_text(size = 1), # May need to be smaller if many species
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
  } else {
    cat("    No Plant (Apidae-Interacting) data found for one or both bioregions to create Presence/Absence plot.\n")
  }
}


# --- Step 8: Generate and Print Plots (from your code) ---

# --- Generate and print the BIOREGION heatmaps ---
plant_heatmap_bioregion <- plot_heatmap(
  plant_similarity_matrix_bioregion, 
  "Plant (Apidae-Interacting) Species Similarity Between Bioregions",
  "Bioregion", "Bioregion"
)
pollinator_heatmap_bioregion <- plot_heatmap(
  pollinator_similarity_matrix_bioregion, 
  "Apidae Pollinator Species Similarity Between Bioregions",
  "Bioregion", "Bioregion"
)

# --- Generate and print the new HABITAT heatmaps ---
plant_heatmap_habitat <- plot_heatmap(
  plant_similarity_matrix_habitat,
  "Heatmap of Plant (Apidae-Interacting) Species Similarity Between Habitats",
  "Habitat", "Habitat"
)
pollinator_heatmap_habitat <- plot_heatmap(
  pollinator_similarity_matrix_habitat,
  "Apidae Pollinator Species Similarity Between Habitats",
  "Habitat", "Habitat"
)

# Print all plots 
cat("--- Generating Plots ---\n")
print(plant_heatmap_bioregion)
print(pollinator_heatmap_bioregion)
print(plant_heatmap_habitat)
print(pollinator_heatmap_habitat)

# --- NEW: Print the cross-similarity plots ---
if (!is.null(pollinator_cross_heatmap)) {
  print(pollinator_cross_heatmap)
}
if (!is.null(plant_cross_heatmap)) {
  print(plant_cross_heatmap)
}

# --- NEW: Print the Presence/Absence plots ---
if (!is.null(pa_heatmap_pollinators)) {
  print(pa_heatmap_pollinators)
}
if (!is.null(pa_heatmap_plants)) {
  print(pa_heatmap_plants)
}

cat("--- Script Finished ---\n")

