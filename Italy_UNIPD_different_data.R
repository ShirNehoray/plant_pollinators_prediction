setwd("/Users/shirnehoray/Downloads")
setwd("/Users/shirn/Documents/Master/data")


# Load required packages
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(bipartite)  # For network analysis
library(vegan)
library(tibble)
library(lubridate)
library(igraph)
library(parallel)
library(foreach)
library(doParallel)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(geosphere)
library(readr)

# Load the data from the CSV file
data <- read_csv("Interaction_data_italy.csv")

# --- Data Preparation for Matrix Creation ---

data_clean <- data %>%
  mutate(
    Interaction = as.numeric(Interaction),
    Interaction = replace_na(Interaction, 0)
  ) %>%
  group_by(Site_id, Plant_species, Pollinator_species) %>%
  summarise(Interaction = sum(Interaction, na.rm = TRUE), .groups = "drop") %>%
  # Filter out pairs with no interactions after aggregation
  filter(Interaction > 0)


# --- Create a list of matrices, one for each network ---

interaction_matrices <- data_clean %>%
  # Split the dataframe into a list based on Site_id
  group_split(Site_id) %>%
  # Set the name of each list element to be its Site_id
  setNames(sapply(., function(df) first(df$Site_id))) %>%
  # Apply the matrix creation logic to each dataframe in the list
  lapply(function(network_df) {
    network_df %>%
      # Remove the Site_id column 
      select(-Site_id) %>%
      # Pivot the data: plants as rows, pollinators as columns
      pivot_wider(
        names_from = Pollinator_species,
        values_from = Interaction,
        values_fill = 0 # If a pair has no interaction, fill with 0
      ) %>%
      # Convert the 'Plant_species' column into row names, creating a true matrix
      column_to_rownames(var = "Plant_species") %>%
      # Officially convert the resulting data frame to a matrix object
      as.matrix()
  })

# --- Verification ---

# Print the names of all created matrices (the Site_ids)
cat("Created", length(interaction_matrices), "matrices. The network IDs are:\n")
print(names(interaction_matrices))

# # 2. Display the first few rows and columns of a single matrix as an example
# cat("\n--- Example Matrix: UNIPD01_C ---\n")
# # Using head() to keep the output manageable
# print(head(interaction_matrices$UNIPD01_C))

# --- Step 2: Identify and Remove Networks with the Invasive Plant ---

# First, ensure the matrix list from the previous step exists.
# If not, run the code from the previous response to create 'interaction_matrices'.

# 1. Identify which network names need to be removed.
#    The conditions are:
#    a) The network name (the key in the list) starts with "UNIPD01".
#    b) The plant "Buddleja_davidii" exists as a row name in that network's matrix.

networks_to_remove <- names(interaction_matrices)[
  # Condition (a): Check which names start with "UNIPD01"
  grepl("^UNIPD01", names(interaction_matrices)) &
    # Condition (b): Check which matrices contain the plant - Buddleja_davidii
    sapply(interaction_matrices, function(mat) "Buddleja_davidii" %in% rownames(mat))
]

# Print the names of the networks that will be removed
cat("The following networks are from the UNIPD01 survey and contain 'Buddleja_davidii':\n")
if (length(networks_to_remove) > 0) {
  print(networks_to_remove)
} else {
  cat("None.\n")
}

# Create a new, filtered list of matrices that excludes the identified networks
interaction_matrices_filtered <- interaction_matrices[
  !names(interaction_matrices) %in% networks_to_remove
]


# --- Verification ---
# Print the counts before and after to confirm the operation was successful.
cat("\n--- Summary of Filtering ---\n")
cat("Original number of matrices:", length(interaction_matrices), "\n")
cat("Number of matrices removed:", length(networks_to_remove), "\n")
cat("Final number of matrices after filtering:", length(interaction_matrices_filtered), "\n")


final_network_names <- names(interaction_matrices_filtered)

data_for_summary <- data %>%
  filter(Site_id %in% final_network_names)

# Filter out all networks that start with "UNIPD03"
data_without_unipd03 <- data_for_summary %>%
  filter(!grepl("^UNIPD03", Site_id))


# create the summary_table using your original code 

# Data Cleaning
data_without_unipd03$Interaction <- as.numeric(as.character(data_without_unipd03$Interaction))
data_without_unipd03$Interaction[is.na(data_without_unipd03$Interaction)] <- 0

# Rename Site_id to network_id
data_without_unipd03 <- data_without_unipd03 %>% rename(network_id = Site_id)

# Aggregate duplicates
data_aggregated <- data_without_unipd03 %>%
  group_by(network_id, Plant_species, Pollinator_species, Sampling_method,
           Sampling_area_square_meters, Habitat, Country, Latitude, Longitude, Year) %>%
  summarise(Interaction = sum(Interaction, na.rm = TRUE), .groups = "drop") # calculationg number of interactions 

# Create the final summary table
summary_table <- data_aggregated %>%
  group_by(network_id) %>%
  summarise(
    total_interaction = sum(Interaction, na.rm = TRUE),
    total_pollinators = n_distinct(Pollinator_species, na.rm = TRUE),
    total_plants = n_distinct(Plant_species, na.rm = TRUE),
    Sampling_method = first(na.omit(Sampling_method)),
    Sampling_area_square = first(na.omit(Sampling_area_square_meters)),
    Habitat = first(na.omit(Habitat)),
    Country = first(na.omit(Country)),
    Latitude = first(na.omit(Latitude)),
    Longitude = first(na.omit(Longitude)),
    year = first(na.omit(Year))
  ) %>%
  ungroup()


# --- Display the final, correct summary table ---
cat("--- Final Summary Table (after excluding networks with Buddleja_davidii) ---\n")
print(summary_table)


# --- Plant Similarity Analysis ---

# Load required packages if not already loaded
library(vegan)
library(pheatmap)

# --- Step 1.1: Create a Site-by-PLANT Matrix ---
# The logic keeps every element WHERE the name does NOT start with "UNIPD03".
interaction_matrices_filtered <- interaction_matrices_filtered[
  !grepl("^UNIPD03", names(interaction_matrices_filtered))
]

# Get a list of all unique plant species from our filtered matrices
all_plants <- unique(unlist(sapply(interaction_matrices_filtered, rownames)))

# Get the names of the networks we are analyzing
network_names <- names(interaction_matrices_filtered)

# Create an empty matrix for plants
site_plant_matrix <- matrix(0,
                            nrow = length(network_names),
                            ncol = length(all_plants),
                            dimnames = list(network_names, all_plants))

# Fill the matrix with presence/absence data (1 if plant is in the network)
for (net_name in network_names) {
  current_plants <- rownames(interaction_matrices_filtered[[net_name]])
  site_plant_matrix[net_name, current_plants] <- 1
}

# --- Step 1.2: Calculate Dissimilarity and Create Heatmap ---

# Calculate Bray-Curtis dissimilarity for plants
plant_bray_dist <- vegdist(site_plant_matrix, method = "bray")

# Convert to a similarity matrix for visualization
plant_similarity_matrix <- 1 - as.matrix(plant_bray_dist)

# Generate the heatmap for PLANTS
pheatmap(plant_similarity_matrix,
         main = "Heatmap of PLANT Similarity Between Networks",
         clustering_distance_rows = plant_bray_dist,
         clustering_distance_cols = plant_bray_dist,
         color = viridis::viridis(100),
         border_color = "grey60",
         fontsize = 8)
