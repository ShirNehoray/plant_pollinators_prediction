# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(bipartite)  # For network analysis
library(vegan)

# Read the data
interactions <- read.csv("Interaction_data.csv")


# Number of unique networks
n_networks <- interactions %>%
  distinct(Network_id) %>%
  count()

print(n_networks)

inter_per_net <- interactions %>%
  group_by(Network_id) %>%
  summarise(Total_Interactions = sum(Interaction))

interactions$Network_id <- as.character(interactions$Network_id)

inter_per_net <- interactions %>%
  group_by(Network_id) %>%
  summarise(Total_Interactions = sum(Interaction)) %>%
  arrange(desc(Total_Interactions))


# Get unique Network_id - Bioregem pairs
region_per_net <- interactions %>%
  select(Network_id, Bioregion) %>%
  distinct()

# Join with interaction summary
inter_per_net_with_region <- inter_per_net %>%
  left_join(region_per_net, by = "Network_id")

## Get the data for North Italy 
# # Filter all rows where the network country is Italy
# italy_data <- interactions %>%
#   filter(Country %in% c("Italy"))
# 
# # Save to CSV
# write.csv(italy_data, "italy_networks.csv", row.names = FALSE)

# Count networks per bioregion
networks_per_bioregion <- interactions %>%
  select(Network_id, Bioregion) %>%
  distinct() %>%
  count(Bioregion)

# Plot
ggplot(networks_per_bioregion, aes(x = reorder(Bioregion, -n), y = n)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Number of Networks per Bioregion",
       x = "Bioregion", y = "Number of Networks") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Count networks per country
networks_per_country <- interactions %>%
  select(Network_id, Country) %>%
  distinct() %>%
  count(Country)

# Plot
ggplot(networks_per_country, aes(x = reorder(Country, -n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Number of Networks per Country",
       x = "Country", y = "Number of Networks") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Get unique network-country-bioregion combinations
country_bioregion <- interactions %>%
  select(Network_id, Country, Bioregion) %>%
  distinct()

# Count networks per country and attach bioregion
networks_per_country <- country_bioregion %>%
  distinct(Network_id, Country, Bioregion) %>%
  count(Country, Bioregion)

# Plot
ggplot(networks_per_country, aes(x = reorder(Country, -n), y = n, fill = Bioregion)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of Networks per Country (Colored by Bioregion)",
       x = "Country", y = "Number of Networks") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")  # or use "Dark2", "Paired", etc.


# Species distribution 
# Number of networks each plant appears in
plant_distribution <- interactions %>%
  select(Network_id, Plant_accepted_name) %>%
  distinct() %>%
  count(Plant_accepted_name) %>%
  arrange(desc(n))

# Plot top 40
ggplot(plant_distribution %>% slice(1:40),
       aes(x = reorder(Plant_accepted_name, n), y = n)) +
  geom_bar(stat = "identity", fill = "seagreen") +
  labs(title = "Top 20 Plant Species by Number of Networks",
       x = "Plant Species", y = "Number of Networks") +
  coord_flip() +
  theme_minimal()


# Histogram
ggplot(plant_distribution, aes(x = n)) +
  geom_histogram(binwidth = 1, fill = "seagreen") +
  labs(title = "Distribution of Plant Species Across Networks",
       x = "Number of Networks", y = "Number of Plant Species") +
  theme_minimal()

# Number of networks each pollinator appears in
pollinator_distribution <- interactions %>%
  select(Network_id, Pollinator_accepted_name) %>%
  distinct() %>%
  count(Pollinator_accepted_name) %>%
  arrange(desc(n))

# Plot top 40
ggplot(pollinator_distribution %>% slice(1:40),
       aes(x = reorder(Pollinator_accepted_name, n), y = n)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  labs(title = "Top 20 Pollinator Species by Number of Networks",
       x = "Pollinator Species", y = "Number of Networks") +
  coord_flip() +
  theme_minimal()


# Histogram
ggplot(pollinator_distribution, aes(x = n)) +
  geom_histogram(binwidth = 1, fill = "darkorange") +
  labs(title = "Distribution of Plant Species Across Networks",
       x = "Number of Networks", y = "Number of Plant Species") +
  theme_minimal()


#------
install.packages("vegan")
install.packages("iNEXT")  # Optional, for extrapolation like panel c
library(dplyr)
library(vegan)
library(ggplot2)
library(tibble)
library(tidyr)


# Create interaction matrix: rows = networks, columns = species
mat_pollinators <- interactions %>%
  select(Network_id, Pollinator_accepted_name) %>%
  distinct() %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Pollinator_accepted_name, values_from = present, values_fill = 0) %>%
  column_to_rownames("Network_id") %>%
  as.matrix()


mat_plants <- interactions %>%
  select(Network_id, Plant_accepted_name) %>%
  distinct() %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Plant_accepted_name, values_from = present, values_fill = 0) %>%
  column_to_rownames("Network_id") %>%
  as.matrix()

spec_curve <- specaccum(mat_pollinators, method = "random")

plot(spec_curve, 
     ci.type = "poly", col = "black", lwd = 2, 
     ci.col = "grey", xlab = "Networks", ylab = "Pollinator species",
     main = "Pollinator Accumulation Curve")


# Calculate occurrence proportion for each pollinator
poll_occurrence <- interactions %>%
  select(Network_id, Pollinator_accepted_name) %>%
  distinct() %>%
  count(Pollinator_accepted_name) %>%
  mutate(freq = n / n_distinct(interactions$Network_id)) %>%
  arrange(desc(freq)) %>%
  mutate(rank = row_number())

# Plot
ggplot(poll_occurrence, aes(x = rank, y = freq)) +
  geom_area(fill = "orange", alpha = 0.8) +
  labs(title = "Pollinator Species Occurrence Across Networks",
       x = "Ranked Pollinator Species",
       y = "Occurrence across networks") +
  theme_minimal()
