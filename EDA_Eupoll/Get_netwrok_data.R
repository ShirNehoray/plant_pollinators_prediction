setwd("/Users/shirn/OneDrive/Documents/master/data")

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(bipartite)  # For network analysis
library(vegan)
library(tibble)



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



#------

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

# Create interaction matrix: rows = plants species, columns = pollinators species
mat_plant_x_poll <- interactions %>%
  select(Plant_accepted_name, Pollinator_accepted_name) %>%
  distinct() %>%                              
  mutate(present = 1L) %>%
  pivot_wider(names_from = Pollinator_accepted_name,
              values_from = present, values_fill = 0) %>%
  column_to_rownames("Plant_accepted_name") %>%
  as.matrix()


# accumulation curve for pollinatiors 
spec_curve <- specaccum(mat_pollinators, method = "random")

df_curve <- data.frame(
  Sites = spec_curve$sites,
  Richness = spec_curve$richness,
  SD = spec_curve$sd
)

# graphing 
ggplot(df_curve, aes(x = Sites, y = Richness)) +
  geom_line(color = "black", size = 1.2) +                       
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), 
              fill = "grey", alpha = 0.5) +                       
  geom_line(data = df_curve, aes(x = Sites, y = Richness))+
  labs(x = "Networks", y = "Pollinator species",
       title = "Pollinator Accumulation Curve") +
  theme_minimal(base_size = 14)



# accumulation curve for plants
spec_curve_plants <- specaccum(mat_plants, method = "random")

df_curve_plants <- data.frame(
  Sites = spec_curve_plants$sites,
  Richness = spec_curve_plants$richness,
  SD = spec_curve_plants$sd
)

# graphing 
ggplot(df_curve_plants, aes(x = Sites, y = Richness)) +
  geom_line(color = "black", size = 1.2) +                       
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), 
              fill = "grey", alpha = 0.5) +                       
  geom_line(data = df_curve_plants, aes(x = Sites, y = Richness))+
  labs(x = "Networks", y = "Plants species",
       title = "Plants Accumulation Curve") +
  theme_minimal(base_size = 14)



# accumulation curve for plants and pollinators 
spec_curve_d <- specaccum(mat_plant_x_poll, method = "random", permutations = 200)

df_d <- data.frame(
  Plants      = spec_curve_d$sites,     
  Pollinators = spec_curve_d$richness,   
  SD          = spec_curve_d$sd
)


#graphing
ggplot(df_d, aes(x = Plants, y = Pollinators)) +
  geom_ribbon(aes(ymin = Pollinators - SD, ymax = Pollinators + SD),
              fill = "grey", alpha = 0.5) +
  geom_line(color = "black", size = 1.2) +
  geom_point(data = df_d[nrow(df_d), , drop = FALSE], size = 2, color = "black") +
  labs(x = "Plant species", y = "Pollinator species",
       title = "Aaccumulation curve of Pollinator across Plant") +
  theme_minimal(base_size = 14)


#------ Species distribution ----

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
  geom_area(fill = "purple", alpha = 0.8) +
  labs(title = "Distribution of Pollinators Species Across Networks",
       x = "Ranked Pollinator Species",
       y = "Pollinator Species") +
  theme_minimal()


# Calculate occurrence proportion for each plant
plant_occurrence <- interactions %>%
  select(Network_id, Plant_accepted_name) %>%
  distinct() %>%
  count(Plant_accepted_name) %>%
  mutate(freq = n / n_distinct(interactions$Network_id)) %>%
  arrange(desc(freq)) %>%
  mutate(rank = row_number())

# Plot
ggplot(plant_occurrence, aes(x = rank, y = freq)) +
  geom_area(fill = "seagreen", alpha = 0.8) +
  labs(title = "Distribution of Plant Species Across Networks",
       x = "Ranked Pollinator Species",
       y = "Pollinator Species") +
  theme_minimal()



# Plot top 40
ggplot(plant_distribution %>% slice(1:40),
       aes(x = reorder(Plant_accepted_name, n), y = n)) +
  geom_bar(stat = "identity", fill = "seagreen") +
  labs(title = "Top 20 Plant Species by Number of Networks",
       x = "Plant Species", y = "Number of Networks") +
  coord_flip() +
  theme_minimal()


# # Histogram
# ggplot(plant_distribution, aes(x = n)) +
#   geom_histogram(binwidth = 1, fill = "seagreen") +
#   labs(title = "Distribution of Plant Species Across Networks",
#        x = "Number of Networks", y = "Number of Plant Species") +
#   theme_minimal()

# # Number of networks each pollinator appears in
# pollinator_distribution <- interactions %>%
#   select(Network_id, Pollinator_accepted_name) %>%
#   distinct() %>%
#   count(Pollinator_accepted_name) %>%
#   arrange(desc(n))

# # Histogram
# ggplot(pollinator_distribution, aes(x = n)) +
#   geom_histogram(binwidth = 1, fill = "darkorange") +
#   labs(title = "Distribution of Plant Species Across Networks",
#        x = "Number of Networks", y = "Number of Plant Species") +
#   theme_minimal()


# Plot top 40
ggplot(pollinator_distribution %>% slice(1:40),
       aes(x = reorder(Pollinator_accepted_name, n), y = n)) +
  geom_bar(stat = "identity", fill = "darkorange") +
  labs(title = "Top 20 Pollinator Species by Number of Networks",
       x = "Pollinator Species", y = "Number of Networks") +
  coord_flip() +
  theme_minimal()


# Habitats 
# Count unique networks per habitat
habitat_counts <- interactions %>%
  select(Network_id, EuPPollNet_habitat) %>%
  distinct() %>%
  count(EuPPollNet_habitat, name = "Num_networks")

# Plot
ggplot(habitat_counts, aes(x = reorder(EuPPollNet_habitat, Num_networks), y = Num_networks)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  labs(title = "Number of Networks per Habitat",
       x = "Habitat", y = "Number of Networks") +
  theme_minimal()

habitat_country_counts <- interactions %>%
  select(Network_id, EuPPollNet_habitat, Country) %>%
  distinct() %>%
  count(EuPPollNet_habitat, Country)

ggplot(habitat_country_counts, aes(x = reorder(EuPPollNet_habitat, -n), y = n, fill = Country)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Networks per Habitat and Country",
       x = "Habitat", y = "Number of Networks") +
  theme_minimal()


ggplot(habitat_country_counts, aes(x = reorder(EuPPollNet_habitat, -n), y = n)) +
  geom_bar(stat = "identity", fill = "darkseagreen") +
  coord_flip() +
  facet_wrap(~ Country, scales = "free_x") +
  labs(title = "Networks per Habitat by Country",
       x = "Habitat", y = "Number of Networks") +
  theme_minimal(base_size = 6)



# ---summarizing species per habitat and pollinator order----
pollinator_habitat <- interactions %>%
  select(Network_id, Habitat = EuPPollNet_habitat, Pollinator_order, Pollinator_accepted_name) %>%
  distinct() %>%
  group_by(Habitat, Pollinator_order) %>%
  summarise(Species_count = n_distinct(Pollinator_accepted_name), .groups = "drop")

# Total species per habitat
total_species_per_habitat <- pollinator_habitat %>%
  group_by(Habitat) %>%
  summarise(Total_species = sum(Species_count))

# Merge and calculate proportion
pollinator_props <- pollinator_habitat %>%
  left_join(total_species_per_habitat, by = "Habitat") %>%
  mutate(Proportion = Species_count / Total_species)

# Add number of studies per habitat (for side bars)
study_counts <- interactions %>%
  select(Network_id, Habitat = EuPPollNet_habitat) %>%
  distinct() %>%
  count(Habitat, name = "Num_studies")

pollinator_props <- pollinator_props %>%
  left_join(study_counts, by = "Habitat")


# Order habitats by number of studies
pollinator_props$Habitat <- factor(pollinator_props$Habitat, 
                                   levels = pollinator_props %>%
                                     group_by(Habitat) %>%
                                     summarise(n = mean(Num_studies)) %>%
                                     arrange(n) %>%
                                     pull(Habitat))

# Stacked bar: proportion of species
p1 <- ggplot(pollinator_props, aes(x = Proportion, y = Habitat, fill = Pollinator_order)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(
    "Hymenoptera" = "#4B0082",
    "Diptera" = "#00BFFF",
    "Lepidoptera" = "#90EE90",
    "Coleoptera" = "#FFD700"
  )) +
  theme_minimal() +
  labs(x = "Proportion of species", y = NULL, fill = "Pollinator Order") +
  theme(legend.position = "bottom")

# Side bar: number of studies
p2 <- ggplot(pollinator_props %>% distinct(Habitat, Num_studies), 
             aes(x = Num_studies, y = Habitat)) +
  geom_bar(stat = "identity", fill = "gray30") +
  theme_minimal() +
  labs(x = "Number of studies", y = NULL) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

install.packages("ggplot2")

install.packages("cowplot")
library(cowplot)

plot_grid(p1, p2, nrow = 1, rel_widths = c(3, 1))



#Network Size


net_species_counts <- interactions %>%
  select(Network_id, Plant_accepted_name, Pollinator_accepted_name) %>%
  distinct() %>%
  group_by(Network_id) %>%
  summarise(
    N_plants = n_distinct(Plant_accepted_name),
    N_pollinators = n_distinct(Pollinator_accepted_name),
    .groups = "drop"
  )

# Plot
ggplot(net_species_counts, aes(x = N_plants, y = N_pollinators)) +
  geom_point(alpha = 0.6, color = "darkblue") +
  labs(title = "Network Size: Number of Species per Network",
       x = "Number of Plant Species", y = "Number of Pollinator Species") +
  theme_minimal()


# Summarise total interaction count per network and country
net_interactions_country <- interactions %>%
  group_by(Network_id, Country) %>%
  summarise(Total_interactions = sum(Interaction), .groups = "drop")

# Plot boxplot of interaction counts by country
ggplot(net_interactions_country, aes(x = Country, y = Total_interactions)) +
  geom_boxplot(fill = "lightblue") +
  labs(title = "Number of Interactions per Network by Country",
       x = "Country", y = "Total Number of Interactions") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(bipartite)
library(tidyr)

# Create presence column
interactions$present <- ifelse(interactions$Interaction > 0, 1, 0)

# Create a named list of matrices (one per network)
network_mats <- interactions %>%
  select(Network_id, Plant_accepted_name, Pollinator_accepted_name, present) %>%
  distinct() %>%
  pivot_wider(
    names_from = Pollinator_accepted_name,
    values_from = present,
    values_fill = 0
  ) %>%
  split(.$Network_id)

# Convert each to matrix
network_matrices <- lapply(network_mats, function(df) {
  mat <- as.matrix(df[,-1])
  rownames(mat) <- df$Plant_accepted_name
  mat
})

# Example: calculate connectance for each matrix
network_stats <- lapply(network_matrices, function(m) {
  m_bin <- (m > 0) * 1
  data.frame(
    Connectance = networklevel(m_bin, index = "connectance"),
    Nestedness  = networklevel(m_bin, index = "nestedness"),
    H2          = networklevel(m_bin, index = "H2")
  )
}) %>% bind_rows(.id = "Network_id")

#box plot for the connectance
network_stats_full <- network_stats %>%
  left_join(interactions %>% distinct(Network_id, Country, Bioregion, EuPPollNet_habitat), by = "Network_id")
ggplot(network_stats_full, aes(x = Country, y = Connectance)) +
  geom_boxplot(fill = "lightgreen") +
  theme_minimal() +
  labs(title = "Connectance Across Countries", y = "Connectance", x = "Country") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(network_stats_full, aes(x = Country, y = Nestedness)) +
  geom_boxplot(fill = "salmon") +
  theme_minimal() +
  labs(title = "Nestedness Across Countries", y = "Nestedness", x = "Country") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
                    