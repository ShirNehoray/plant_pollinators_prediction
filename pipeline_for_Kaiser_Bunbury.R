# Load necessary libraries
library(readxl)     # For reading Excel
library(dplyr)
library(tidyr)
library(ggplot2)
library(bipartite)  # For network analysis
library(vegan)
library(tibble)
library(lubridate)  # Though limited use here since no full dates
library(igraph)
library(parallel)
library(foreach)
library(doParallel)
# Note: Skipping spatial libraries (sf, rnaturalearth, geosphere) for simplicity as no per-network coords available

# Read the data
net_pa <- read_excel("Kaiser_Bunbury_et_al_2017.xlsx", sheet = "64 networks_no.visits")      # number of visits
net_freq <- read_excel("Kaiser_Bunbury_et_al_2017.xlsx", sheet = "64 networks_visitfreq")    # frequency (not used here, but available)
plants <- read_excel("Kaiser_Bunbury_et_al_2017.xlsx",
                     sheet = "Plant species",
                     skip = 1,      # skips the note row
                     col_names = TRUE)

pollinators <- read_excel("Kaiser_Bunbury_et_al_2017.xlsx",
                          sheet = "Pollinator species",
                          skip = 1,  # skips the note row
                          col_names = TRUE)
info <- read_excel("Kaiser_Bunbury_et_al_2017.xlsx", sheet = "info")


# Prepare interactions data in long format (using number of visits)
interactions <- net_pa %>%
  pivot_longer(cols = C1:YW, names_to = "Pollinator_id", values_to = "Interaction")
  # filter(Interaction > 0)  # Keep only positive interactions to mimic presence of links

interactions <- interactions %>%
  rename(Network_ID = `Network ID`, Plant_species_ID = 'Plant species ID')



# # Join with plant and pollinator names
# interactions <- interactions %>%
#   left_join(plants %>% select(`Plant species ID`, `Plant species name`), by = c("Plant species ID" = "Plant species ID")) %>%
#   left_join(pollinators %>% select(`Pollinator species ID`, `Pollinator species name`), by = c("Pollinator_id" = "Pollinator species ID")) %>%
#   rename(Plant_accepted_name = `Plant species name`,
#          Pollinator_accepted_name = `Pollinator species name`,
#          Network_id = `Network ID`,
#          Plant_original_name = `Plant species ID`,  # Use ID as placeholder if needed
#          Pollinator_original_name = Pollinator_id,
#          Habitat = Treatment)  # Use Treatment as proxy for Habitat/Bioregion


# Summarize network-level information, counting species only for Interaction > 1
network_summary <- interactions %>%
  group_by(Network_ID, Site, Treatment) %>%
  summarise(
    n_plants = n_distinct(`Plant species ID`[Interaction > 1]),
    n_pollinators = n_distinct(Pollinator_id[Interaction > 1]),
    .groups = "drop"
  )


# Number of unique networks (should be 64)
n_networks <- network_summary %>%
  distinct(Network_ID) %>%
  count()

print(n_networks)

# Interactions per network 
inter_per_net <- interactions %>%
  group_by(Network_ID) %>%
  summarise(Total_Interactions = sum(Interaction)) %>%
  arrange(desc(Total_Interactions))

# Get unique Network_id - Bioregion pairs (using Treatment as proxy)
treatmant_per_net <- interactions %>%
  select(Network_ID, Treatment) %>%
  distinct()

# Join with interaction summary
inter_per_net_with_treatmant <- inter_per_net %>%
  left_join(treatmant_per_net, by = "Network_ID")

# Count networks per bioregion (proxy)
networks_per_treatmant <- network_summary %>%
  select(Network_ID, Treatment) %>%
  distinct() %>%
  count(Treatment)


# Plot
ggplot(networks_per_treatmant, aes(x = reorder(Treatment, -n), y = n)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Number of Networks per Bioregion (Proxy)",
       x = "Treatment", y = "Number of Networks") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#---- Interactions analysis ---- 

# Create interaction matrix: rows = networks, columns = species
# Pollinators
mat_pollinators <- interactions %>%
  select(Network_ID, Pollinator_id) %>%
  distinct() %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Pollinator_id, values_from = present, values_fill = 0) %>%
  column_to_rownames("Network_ID") %>%
  as.matrix()



# Plants
mat_plants <- interactions %>%
  select(Network_ID, Plant_species_ID) %>%
  distinct() %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Plant_species_ID, values_from = present, values_fill = 0) %>%
  column_to_rownames("Network_ID") %>%
  as.matrix()

# Plants x Pollinators (global)
mat_plant_x_poll <- interactions %>%
  filter(Interaction > 1) %>%  # Filter for interactions with more than 1 visit
  select(Plant_species_ID, Pollinator_id) %>%
  distinct() %>%  # Keep unique plant-pollinator pairs
  mutate(present = 1L) %>%
  pivot_wider(names_from = Pollinator_id,
              values_from = present,
              values_fill = 0) %>%
  column_to_rownames("Plant_species_ID") %>%
  as.matrix()

# Accumulation curves 
# Accumulation curve for pollinators 
spec_curve <- specaccum(mat_pollinators, method = "random")

df_curve <- data.frame(
  Sites = spec_curve$sites,
  Richness = spec_curve$richness,
  SD = spec_curve$sd
)

# Graphing 
ggplot(df_curve, aes(x = Sites, y = Richness)) +
  geom_line(color = "black", size = 1.2) +                       
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), 
              fill = "grey", alpha = 0.5) +                       
  geom_line(data = df_curve, aes(x = Sites, y = Richness))+
  labs(x = "Networks", y = "Pollinator species",
       title = "Pollinator Accumulation Curve") +
  theme_minimal(base_size = 14)

# Accumulation curve for plants
spec_curve_plants <- specaccum(mat_plants, method = "random")

df_curve_plants <- data.frame(
  Sites = spec_curve_plants$sites,
  Richness = spec_curve_plants$richness,
  SD = spec_curve_plants$sd
)

# Graphing 
ggplot(df_curve_plants, aes(x = Sites, y = Richness)) +
  geom_line(color = "black", size = 1.2) +                       
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), 
              fill = "grey", alpha = 0.5) +                       
  geom_line(data = df_curve_plants, aes(x = Sites, y = Richness))+
  labs(x = "Networks", y = "Plant species",
       title = "Plant Accumulation Curve") +
  theme_minimal(base_size = 14)

#Accumulation curve for pollinators across plants 
spec_curve_d <- specaccum(mat_plant_x_poll, method = "random", permutations = 200)

df_d <- data.frame(
  Plants = spec_curve_d$sites,     
  Pollinators = spec_curve_d$richness,   
  SD = spec_curve_d$sd
)

# Graphing
ggplot(df_d, aes(x = Plants, y = Pollinators)) +
  geom_ribbon(aes(ymin = Pollinators - SD, ymax = Pollinators + SD),
              fill = "grey", alpha = 0.5) +
  geom_line(color = "black", size = 1.2) +
  labs(x = "Plant species", y = "Pollinator species",
       title = "Accumulation curve of Pollinator across Plant") +
  theme_minimal(base_size = 14)

#---- Species distribution ----

# Calculate occurrence proportion for each pollinator
poll_occurrence <- interactions %>%
  filter(Interaction > 1) %>%  # Filter for interactions with more than 1 visit
  select(Network_ID, Pollinator_id) %>%
  distinct() %>%
  count(Pollinator_id) %>%
  mutate(freq = n / n_distinct(interactions$Network_ID)) %>%
  arrange(desc(freq)) %>%
  mutate(rank = row_number())

# Plots
ggplot(poll_occurrence, aes(x = n)) +
  geom_histogram( fill = "purple", color = "black") +
  labs(title = "Distribution of Number of Networks per Number of Pollinators",
       x = "Number of networks",
       y = "Number of pollinator species") +
  theme_minimal()

#Frequency 
ggplot(poll_occurrence, aes(x = freq)) +
  geom_histogram(binwidth = 0.05, fill = "purple", color = "black") +
  labs(title = "Distribution of Pollinator Frequencies",
       x = "Frequency across networks",
       y = "Number of species") +
  theme_minimal()

# Calculate occurrence proportion for each plant
plant_occurrence <- interactions %>%
  select(Network_ID, Plant_species_ID) %>%
  distinct() %>%
  count(Plant_species_ID) %>%
  mutate(freq = n / n_distinct(interactions$Network_ID)) %>%
  arrange(desc(freq)) %>%
  mutate(rank = row_number())

# Plot
ggplot(plant_occurrence, aes(x = n)) +
  geom_histogram(binwidth = 1, fill = "seagreen", color = "black") +
  labs(title = "Distribution of Number of Networks per Number of Plants",
       x = "Number of networks",
       y = "Number of Plant species") +
  theme_minimal()

#fraquency 
ggplot(plant_occurrence, aes(x = freq)) +
  geom_histogram(binwidth = 0.05, fill = "seagreen", color = "black") +
  labs(title = "Distribution of Plant Frequencies",
       x = "Frequency across networks",
       y = "Number of species") +
  theme_minimal()



#---- Network Analysis ----

# Plot geom point graph for the number of plants and pollinators 
ggplot(network_summary, aes(x = n_plants, y = n_pollinators)) +
  geom_point(alpha = 0.6, color = "darkblue") +
  scale_x_continuous(breaks = seq(min(network_summary$n_plants),
                                    max(network_summary$n_plants), 
                                    by = 1)) +
  labs(title = "Network Size: Number of Species per Network",
       x = "Number of Plant Species", y = "Number of Pollinator Species") +
  theme_minimal()

# # Filter and plot large networks (e.g., N_plants and N_pollinators >= 10, adjust threshold as needed)
# large_net_threshold <- 10  # Lowered threshold for this dataset
# large_net_counts <- network_summary %>%
#   filter(n_plants >= large_net_threshold & n_pollinators >= large_net_threshold)
# 
# ggplot(large_net_counts, aes(x = n_plants, y = n_pollinators, color = Country)) +  # Country all same, or use Site
#   geom_point(alpha = 0.6, size = 3) +
#   labs(title = paste0("Large Networks (≥ ", large_net_threshold, " Species)"),
#        x = "Number of Plant Species", y = "Number of Pollinator Species") +
#   theme_minimal() +
#   theme(legend.position = "right")


# Summaries total interaction count per network and site (proxy for country)
net_interactions_site <- interactions %>%
  group_by(Network_ID, Site) %>%
  summarise(Total_interactions = sum(Interaction), .groups = "drop")

# Plot boxplot of interaction counts by site
ggplot(net_interactions_site, aes(x = Site, y = Total_interactions)) +
  geom_boxplot(fill = "lightblue") +
  labs(title = "Number of Interactions per Network by Site",
       x = "Site", y = "Total Number of Interactions") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# For indexes calculations

# Create presence column (1 = present)
interactions$present <- ifelse(interactions$Interaction > 0)

# Create a named list of matrices (one per network)
network_mats <- interactions %>%
  select(Network_ID, Plant_species_ID, Pollinator_id, present) %>%
  distinct() %>%
  pivot_wider(
    names_from = Pollinator_id,
    values_from = present,
    values_fill = 0
  ) %>%
  split(.$Network_ID)

# Convert each to matrix
network_matrices <- lapply(network_mats, function(df) {
  mat <- as.matrix(df[,-1])  # Exclude Network_id
  rownames(mat) <- df$Plant_accepted_name
  mat
})

# Set up parallel backend
cores <- detectCores() - 1
registerDoParallel(cores)

# Step 1: Get network sizes to filter small networks
network_sizes <- lapply(names(network_matrices), function(net_id) {
  m <- network_matrices[[net_id]]
  data.frame(Network_ID = net_id, N_plants = nrow(m), N_pollinators = ncol(m))
}) %>% bind_rows()

# Filter large networks (e.g., >=3 plants and >=3 pollinators)
large_networks <- network_sizes %>%
  filter(N_plants >= 3, N_pollinators >= 3) %>%
  pull(Network_ID)

# Step 2: Calculate network metrics for large networks
network_stats_large <- foreach(net_id = large_networks, .combine = rbind, 
                               .packages = c("igraph", "bipartite")) %dopar% {
                                 m_bin <- (network_matrices[[net_id]] > 0) * 1  # Binary matrix
                                 g <- graph_from_incidence_matrix(m_bin, directed = FALSE)
                                 g_poll <- bipartite.projection(g)$proj2  # Pollinator projection
                                 
                                 # Network metrics
                                 net_levels <- networklevel(m_bin, index = c("connectance", "NODF"))
                                 connectance <- net_levels["connectance"]
                                 nodf <- net_levels["NODF"]
                                 evenness <- networklevel(m_bin, index = "interaction evenness")
                                 
                                 # Modularity (skip if projected network too small)
                                 mod <- if (vcount(g_poll) >= 5 && ecount(g_poll) >= 5) {
                                   comm <- cluster_fast_greedy(g_poll)
                                   modularity(g_poll, membership(comm))
                                 } else {
                                   NA
                                 }
                                 
                                 data.frame(
                                   Network_ID = net_id,
                                   Connectance = connectance,
                                   Nestedness_NODF = nodf,
                                   Modularity = mod,
                                   Interaction_Evenness = evenness
                                 )
                               }

# Stop parallel backend
stopImplicitCluster()

# # Step 3: Assign NA to small networks
#  small_networks <- setdiff(names(network_matrices), large_networks)
# network_stats_small <- data.frame(
#   Network_id = small_networks,
#   Connectance = NA,
#   Nestedness_NODF = NA,
#   Modularity = NA,
#   Interaction_Evenness = NA
# )

# # Combine results
# network_stats <- bind_rows(network_stats_large, network_stats_small)

# Join with metadata
network_stats_full <- network_stats_large %>%
  left_join(interactions %>% distinct(Network_ID, Site), by = "Network_ID") %>%
  left_join(network_sizes, by = "Network_ID")

# Visualizations
# Connectance
ggplot(network_stats_full, aes(x = Site, y = Connectance)) +  # Country all same; adapt to Site if needed
  geom_boxplot(fill = "lightgreen") +
  theme_minimal() +
  labs(title = "Connectance Across Sites", y = "Connectance", x = "Site") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Nestedness NODF
ggplot(network_stats_full, aes(x = Site, y = Nestedness_NODF)) +
  geom_boxplot(fill = "salmon") +
  theme_minimal() +
  labs(title = "Nestedness (NODF) Across Sites", y = "NODF", x = "Site") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Modularity
ggplot(network_stats_full %>% filter(!is.na(Modularity)), aes(x = Site, y = Modularity)) +
  geom_boxplot(fill = "lightpink") +
  theme_minimal() +
  labs(title = "Modularity Across Sites", y = "Modularity", x = "Site") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Interaction Evenness
ggplot(network_stats_full, aes(x = Site, y = Interaction_Evenness)) +
  geom_boxplot(fill = "lightblue") +
  theme_minimal() +
  labs(title = "Interaction Evenness Across Sites", y = "Evenness", x = "Site") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# Species-level analysis - degree distribution 
species_roles <- lapply(network_matrices, function(m) {
  m_bin <- (m > 0) * 1
  roles <- specieslevel(m_bin, index = c("degree", "normalized degree"))
  # Combine plant + pollinator into one dataframe
  df <- bind_rows(
    as.data.frame(roles$`higher level`) %>%
      rownames_to_column("Species") %>%
      mutate(Level = "Plant"),
    as.data.frame(roles$`lower level`) %>%
      rownames_to_column("Species") %>%
      mutate(Level = "Pollinator")
  )
  df
}) %>%
  bind_rows(.id = "Network_ID")


# Join with metadata (using majority or first occurrence)
meta_by_net <- interactions %>%
  distinct(Network_ID, Site, Treatment)

# Safe join
species_roles_full <- species_roles %>%
  left_join(meta_by_net, by = "Network_ID")

# Explore top high degree
top_generalists <- species_roles_full %>%
  arrange(desc(degree)) %>%
  group_by(Level) %>%
  slice_head(n = 10)

print(top_generalists)

# Plots
ggplot(species_roles_full, aes(x = degree, fill = Level)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  theme_minimal() +
  labs(title = "Degree distribution of plants and pollinators",
       x = "Degree", y = "Number of species")


