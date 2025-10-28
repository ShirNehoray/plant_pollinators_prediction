# Load necessary packages
# Ensure you have these installed: install.packages(c("dplyr", "ggplot2", "tidyr"))
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)

setwd("/Users/shirnehoray/EDA/plant_pollinators_prediction/EDA_Eupoll")


# Load the data from the CSV file
data <- read_csv("Interaction_data.csv")



# --- Step 2: Aggregate Species by Bioregion (All Plants, Apidae Pollinators with Interaction > 0) ---

# Aggregate all plant species from the original dataset
plant_species_per_bioregion <- data %>%
  group_by(Bioregion) %>%
  summarise(
    Plants = list(unique(Plant_accepted_name)),
    .groups = "drop"
  )

# Filter FOR Apidae with interaction > 0 and aggregate pollinator species
filtered_pollinator_species_per_bioregion <- data %>%
  filter(Pollinator_family == "Apidae") %>%
  group_by(Bioregion) %>%
  summarise(
    Pollinators = list(unique(Pollinator_accepted_name)),
    .groups = "drop"
  )

# Create a data frame with all unique bioregions to ensure matrices are comparable
all_bioregions_df <- data %>%
  select(Bioregion) %>%
  distinct()

# Join the plant and pollinator data. This ensures that bioregions without qualifying
# pollinators are included, but their 'Pollinators' list will be NA after the join.
species_per_bioregion <- all_bioregions_df %>%
  left_join(plant_species_per_bioregion, by = "Bioregion") %>%
  left_join(filtered_pollinator_species_per_bioregion, by = "Bioregion")

# Replace NA/NULL list entries with empty lists for the Jaccard calculation.
# A NULL entry means a bioregion had no qualifying pollinators after filtering.
for (i in 1:nrow(species_per_bioregion)) {
  if (is.null(species_per_bioregion$Pollinators[[i]])) {
    species_per_bioregion$Pollinators[[i]] <- list()
  }
}


# Get the list of unique bioregions
bioregions <- species_per_bioregion$Bioregion


# --- Step 3: Calculate Jaccard Similarity Matrix ---

# Define the Jaccard similarity function
jaccard_similarity <- function(set1, set2) {
  # Find the number of species common to both sets
  intersection_size <- length(intersect(set1, set2))
  # Find the total number of unique species across both sets
  union_size <- length(union(set1, set2))
  
  # Avoid division by zero if both sets are empty
  if (union_size == 0) {
    return(NA)
  } else {
    return(intersection_size / union_size)
  }
}

# Function to create a similarity matrix for a given species type ("Plants" or "Pollinators")
create_similarity_matrix <- function(species_data, species_type) {
  num_bioregions <- nrow(species_data)
  
  # Create an empty matrix to store the similarity scores
  sim_matrix <- matrix(NA,
                       nrow = num_bioregions,
                       ncol = num_bioregions,
                       dimnames = list(species_data$Bioregion, species_data$Bioregion))
  
  # Loop through each pair of bioregions to calculate similarity
  for (i in 1:num_bioregions) {
    for (j in 1:num_bioregions) {
      # Get the species lists for the two bioregions
      species_set1 <- unlist(species_data[[species_type]][i])
      species_set2 <- unlist(species_data[[species_type]][j])
      
      # Calculate and store the Jaccard similarity score
      sim_matrix[i, j] <- jaccard_similarity(species_set1, species_set2)
    }
  }
  return(sim_matrix)
}

# Create the similarity matrices for both plants and filtered pollinators
plant_similarity_matrix <- create_similarity_matrix(species_per_bioregion, "Plants")
pollinator_similarity_matrix <- create_similarity_matrix(species_per_bioregion, "Pollinators")


# --- Step 4: Visualize the Heatmaps using ggplot2 ---

# Function to plot a heatmap from a similarity matrix
plot_heatmap <- function(similarity_matrix, title) {
  
  # Convert the matrix to a long-format data frame, which is required for ggplot2
  # This creates three columns: Var1 (Bioregion 1), Var2 (Bioregion 2), and value (Similarity)
  long_format_data <- as.data.frame(as.table(similarity_matrix))
  names(long_format_data) <- c("Bioregion1", "Bioregion2", "Similarity")
  
  # Create the heatmap plot
  heatmap_plot <- ggplot(long_format_data, aes(x = Bioregion1, y = Bioregion2, fill = Similarity)) +
    geom_tile(color = "white") + # Add white lines between tiles for better separation
    # Add text labels for the Jaccard index values, rounded to 2 decimal places
    geom_text(aes(label = round(Similarity, 2)), color = "black", size = 3) +
    # Use the 'viridis' color palette as requested
    scale_fill_viridis_c(option = "viridis", name = "Jaccard\nSimilarity") +
    labs(
      title = title,
      x = "Bioregion",
      y = "Bioregion"
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

# Generate and print the two heatmaps
plant_heatmap <- plot_heatmap(plant_similarity_matrix, " Plant Species Similarity Between Bioregions")
pollinator_heatmap <- plot_heatmap(pollinator_similarity_matrix, "Apidae Pollinator Species Similarity")

# Print the plots to the console/viewer
print(plant_heatmap)
print(pollinator_heatmap)





# Aggregate all unique habitats from the original dataset
habitat_per_bioregion <- data %>%
  group_by(Bioregion) %>%
  summarise(
    Habitats = list(unique(EuPPollNet_habitat)),
    .groups = "drop"
  )

# Create a data frame with all unique bioregions to ensure matrices are comparable
all_bioregions_df <- data %>%
  select(Bioregion) %>%
  distinct()

# Join the plant, pollinator, and habitat data.
species_per_bioregion <- all_bioregions_df %>%
  left_join(plant_species_per_bioregion, by = "Bioregion") %>%
  left_join(filtered_pollinator_species_per_bioregion, by = "Bioregion") %>%
  left_join(habitat_per_bioregion, by = "Bioregion")

# Replace NA/NULL list entries with empty lists for the Jaccard calculation.
for (i in 1:nrow(species_per_bioregion)) {
  if (is.null(species_per_bioregion$Pollinators[[i]])) {
    species_per_bioregion$Pollinators[[i]] <- list()
  }
  if (is.null(species_per_bioregion$Habitats[[i]])) {
    species_per_bioregion$Habitats[[i]] <- list()
  }
}

# Get the list of unique bioregions
bioregions <- species_per_bioregion$Bioregion


# --- Step 3: Calculate Jaccard Similarity Matrix ---

# Define the Jaccard similarity function
jaccard_similarity <- function(set1, set2) {
  # Find the number of items common to both sets
  intersection_size <- length(intersect(set1, set2))
  # Find the total number of unique items across both sets
  union_size <- length(union(set1, set2))
  
  # Avoid division by zero if both sets are empty
  if (union_size == 0) {
    return(NA)
  } else {
    return(intersection_size / union_size)
  }
}

# Function to create a similarity matrix for a given type ("Plants", "Pollinators", or "Habitats")
create_similarity_matrix <- function(species_data, data_type) {
  num_bioregions <- nrow(species_data)
  
  # Create an empty matrix to store the similarity scores
  sim_matrix <- matrix(NA,
                       nrow = num_bioregions,
                       ncol = num_bioregions,
                       dimnames = list(species_data$Bioregion, species_data$Bioregion))
  
  # Loop through each pair of bioregions to calculate similarity
  for (i in 1:num_bioregions) {
    for (j in 1:num_bioregions) {
      # Get the lists for the two bioregions
      set1 <- unlist(species_data[[data_type]][i])
      set2 <- unlist(species_data[[data_type]][j])
      
      # Calculate and store the Jaccard similarity score
      sim_matrix[i, j] <- jaccard_similarity(set1, set2)
    }
  }
  return(sim_matrix)
}

# Create the similarity matrices
plant_similarity_matrix <- create_similarity_matrix(species_per_bioregion, "Plants")
pollinator_similarity_matrix <- create_similarity_matrix(species_per_bioregion, "Pollinators")
habitat_similarity_matrix <- create_similarity_matrix(species_per_bioregion, "Habitats")


# --- Step 4: Visualize the Heatmaps using ggplot2 ---

# Function to plot a heatmap from a similarity matrix
plot_heatmap <- function(similarity_matrix, title) {
  
  # Convert the matrix to a long-format data frame, which is required for ggplot2
  long_format_data <- as.data.frame(as.table(similarity_matrix))
  names(long_format_data) <- c("Bioregion1", "Bioregion2", "Similarity")
  
  # Create the heatmap plot
  heatmap_plot <- ggplot(long_format_data, aes(x = Bioregion1, y = Bioregion2, fill = Similarity)) +
    geom_tile(color = "white") + # Add white lines between tiles for better separation
    # Add text labels for the Jaccard index values, rounded to 2 decimal places
    geom_text(aes(label = round(Similarity, 2)), color = "black", size = 3) +
    # Use the 'viridis' color palette as requested
    scale_fill_viridis_c(option = "viridis", name = "Jaccard\nSimilarity") +
    labs(
      title = title,
      x = "Bioregion",
      y = "Bioregion"
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

# Generate and print the three heatmaps
plant_heatmap <- plot_heatmap(plant_similarity_matrix, "Heatmap of Plant Species Similarity Between Bioregions")
pollinator_heatmap <- plot_heatmap(pollinator_similarity_matrix, "Heatmap of Apidae (Interaction > 0) Pollinator Species Similarity")
habitat_heatmap <- plot_heatmap(habitat_similarity_matrix, "Heatmap of Habitat Similarity Between Bioregions")

# Print the plots to the console/viewer
print(plant_heatmap)
print(pollinator_heatmap)
print(habitat_heatmap)


