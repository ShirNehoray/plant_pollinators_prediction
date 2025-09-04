setwd("/Users/shirnehoray/EDA/EDA_Eupoll")

# Pacages 
library(dplyr)
library(ggplot2)
library(bipartite)   # for network analysis
library(tidyr)
library(readr)
library(ggplot2)
library(vegan)


# Read in the Italy-only dataset you exported
italy_data <- read_csv("italy_networks.csv")

n_networks <- italy_data %>%
  distinct(Network_id) %>%
  nrow()

n_networks # number of networks in Italy 

# Get in how many studies explored this country
unique_studies <- italy_data %>%
  distinct(Study_id)

print(unique_studies)

networks_per_study <- italy_data %>%
  distinct(Study_id, Network_id) %>%
  count(Study_id, name = "N_networks")

print(networks_per_study)

study_bioregion <- italy_data %>%
  distinct(Study_id, Network_id, Bioregion) %>%
  group_by(Study_id, Bioregion) %>%
  summarise(N_networks = n_distinct(Network_id), .groups = "drop") %>%
  arrange(Study_id, Bioregion)

print(study_bioregion)

ggplot(study_bioregion, aes(x = Study_id, y = N_networks, fill = Bioregion)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of Networks per Study and Bioregion (Italy)",
       x = "Study ID",
       y = "Number of Networks") +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 1))


# Get what bioregions is includes 
unique_bioregion <- italy_data %>%
  distinct(Bioregion)

print(unique_bioregion)

networks_per_bioregion <- italy_data %>%
  distinct(Network_id, Bioregion) %>%   # each network counted once per bioregion
  count(Bioregion)

print(networks_per_bioregion)

# Number of network per region 
ggplot(networks_per_bioregion, aes(x = Bioregion, y = n, fill = Bioregion)) +
  geom_bar(stat = "identity") +
  labs(title = "Number of Networks per Bioregion (Italy)",
       x = "Bioregion",
       y = "Number of Networks") +
  theme(axis.text.x = element_text( hjust = 1),
        legend.position = "none")

# Network analysis 
analyze_network <- function(df, net_id) {
  sub_df <- df %>% filter(Network_id == net_id)

  # Create incidence matrix
  mat <- sub_df %>%
    group_by(Plant_accepted_name, Pollinator_accepted_name) %>%
    summarise(freq = sum(Interaction), .groups = "drop") %>%
    pivot_wider(names_from = Pollinator_accepted_name, 
                values_from = freq, 
                values_fill = 0) %>%
    column_to_rownames("Plant_accepted_name") %>%
    as.matrix()

  bin_mat <- (mat > 0) * 1   # binary version

  data.frame(
    Network_id   = net_id,
    Plants       = nrow(bin_mat),
    Pollinators  = ncol(bin_mat),
    Links        = sum(bin_mat),
    Connectance  = networklevel(bin_mat, index = "connectance"),
    Nestedness   = networklevel(bin_mat, index = "nestedness")
  )
}

# Apply to all Italian networks
all_ids <- unique(italy_data$Network_id)
network_summary <- do.call(rbind, lapply(all_ids, function(id) analyze_network(italy_data, id)))

print(network_summary)


ggplot(network_summary, aes(x = Connectance)) +
  geom_histogram(binwidth = 0.01, fill = "tomato", color = "white") +
  labs(title = "Distribution of Connectance",
       x = "Connectance",
       y = "Count of Networks") +
  theme_minimal()



# ---- Species counts per network ----
species_counts <- df %>%
  group_by(Network_id) %>%
  summarise(
    N_plants       = n_distinct(Plant_accepted_name, na.rm = TRUE),
    N_pollinators  = n_distinct(Pollinator_accepted_name, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(N_total = N_plants + N_pollinators)

# make an output folder
dir.create("species_distribution", showWarnings = FALSE)

# save the table too
write.csv(species_counts, "species_distribution/species_counts_per_network.csv", row.names = FALSE)

# ---- Histograms: species distribution across networks ----

# 1) Plant species per network
p1 <- ggplot(species_counts, aes(x = N_plants)) +
  geom_histogram(fill = "lightgreen", bins = 30) +
  labs(title = "Distribution of PLANT species per network",
       x = "Number of plant species (per network)",
       y = "Number of networks") +
  theme_minimal() 
ggsave("species_distribution/hist_plants_per_network.png", p1, width = 8, height = 5, dpi = 300)

# 2) Pollinator species per network
p2 <- ggplot(species_counts, aes(x = N_pollinators)) +
  geom_histogram(fill = "salmon", bins = 30) +
  labs(title = "Distribution of POLLINATOR species per network",
       x = "Number of pollinator species (per network)",
       y = "Number of networks") +
  theme_minimal()
ggsave("species_distribution/hist_pollinators_per_network.png", p2, width = 8, height = 5, dpi = 300)

# 3) Total species per network (plants + pollinators)
p3 <- ggplot(species_counts, aes(x = N_total)) +
  geom_histogram(fill = "#1f78b4", bins = 30) +
  labs(title = "Distribution of TOTAL species per network",
       x = "Total species per network (plants + pollinators)",
       y = "Number of networks") +
  theme_minimal()
ggsave("species_distribution/hist_total_species_per_network.png", p3, width = 8, height = 5, dpi = 300)




