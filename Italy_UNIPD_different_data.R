setwd("/Users/shirnehoray/Downloads")

# Load required packages
library(readxl)
library(dplyr)

# Load the data from the Excel file
# Adjust the path if the file is not in your working directory
data <- read_excel("Raw (3).xlsx", sheet = "Interaction_data")

# Data Cleaning
# Ensure Interaction is numeric
data$Interaction <- as.numeric(as.character(data$Interaction))

# Handle NA values in Interaction (replace with 0 for summation)
data$Interaction[is.na(data$Interaction)] <- 0

# Rename Site_id to network_id
data <- data %>% rename(network_id = Site_id)

# Aggregate duplicates to avoid double-counting interactions for the same plant-pollinator pair
data <- data %>%
  group_by(network_id, Plant_species, Pollinator_species, Sampling_method, 
           Sampling_area_square_meters, Habitat, Country, Latitude, Longitude, Year) %>%
  summarise(Interaction = sum(Interaction, na.rm = TRUE), .groups = "drop")

# Create the summary table, ensuring unique species counts
summary_table <- data %>%
  group_by(network_id) %>%
  summarise(
    total_interaction = sum(Interaction, na.rm = TRUE),
    total_pollinators = n_distinct(Pollinator_species, na.rm = TRUE),  # Count unique pollinator species
    total_plants = n_distinct(Plant_species, na.rm = TRUE),  # Count unique plant species
    Sampling_method = first(na.omit(Sampling_method)),  # First non-NA value
    Sampling_area_square = first(na.omit(Sampling_area_square_meters)),  # First non-NA value
    Habitat = first(na.omit(Habitat)),  # First non-NA value
    Country = first(na.omit(Country)),  # First non-NA value
    Latitude = first(na.omit(Latitude)),  # First non-NA value
    Longitude = first(na.omit(Longitude)),  # First non-NA value
    year = first(na.omit(Year))  # First non-NA value
  ) %>%
  ungroup()

# Display the table
print(summary_table)

