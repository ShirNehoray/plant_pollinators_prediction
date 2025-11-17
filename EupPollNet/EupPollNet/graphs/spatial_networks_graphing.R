library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(vegan)        
library(purrr)

tme <-  theme(axis.text = element_text(size = 10, color = "black"),
              axis.title = element_text(size = 16, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
theme_set(theme_bw())

# Read the interaction data
data <- read.csv("Interaction_data.csv", encoding = "latin1", stringsAsFactors = FALSE)
metadata <- read_csv("metadata.csv", show_col_types = FALSE)


# Number of networks per each study ----- 
# classify the studies 
classification_col <- grep("clasification|classification", names(metadata), ignore.case = TRUE, value = TRUE)[1]

metadata <- metadata %>%
  mutate(!!sym(classification_col) := str_trim(!!sym(classification_col))) %>%
  mutate(!!sym(classification_col) := case_when(
    str_detect(!!sym(classification_col), "spatial.*temporal|temporal.*spatial") ~ "spatial and temporal",
    !!sym(classification_col) == "spatial" ~ "spatial",
    !!sym(classification_col) == "temporal" ~ "temporal",
    TRUE ~ !!sym(classification_col)
  ))


# summary the classify the studies 
study_summary <- metadata %>%
  group_by(Study_id) %>%
  summarise(
    total_networks = sum(number_of_networks, na.rm = TRUE), 
    classification = case_when(
      any(!!sym(classification_col) == "spatial and temporal") ~ "spatial and temporal",
      any(!!sym(classification_col) == "temporal") ~ "temporal",
      any(!!sym(classification_col) == "spatial") ~ "spatial",
      TRUE ~ "unknown"
    ),
    .groups = "drop"
  ) %>%
  filter(total_networks > 0) %>%
  arrange(desc(total_networks)) %>%
  mutate(Study_id = factor(Study_id, levels = Study_id))  # keep sorted order

# graphing the number of network per study and classify the study
ggplot(study_summary, aes(x = Study_id, y = total_networks, fill = classification)) +
  geom_col(width = 0.9) +
  geom_text(aes(label = total_networks), 
            vjust = -0.5, size = 2, color = "black") +
  scale_fill_manual(
    values = c(
      "spatial" = "plum4",
      "temporal" = "lightblue3",
      "spatial and temporal" = "lightcoral"
    ),
    name = "Classification"
  ) +
  labs(
    title = " Number of Networks per Study",
    x = "Study_id",
    y = "Number of Networks"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8.5),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  tme


# Networks I chose for the research -----


spatial_df <- study_summary %>%
  filter(classification %in% c("spatial", "spatial and temporal")) %>%
  arrange(desc(total_networks)) %>%
  mutate(Study_id = factor(Study_id, levels = Study_id))   # keep sorted order

# take studies with more then 4 networks 
spatial_df_filtered <- spatial_df %>%
  filter(total_networks >= 4) %>%
  mutate(Study_id = factor(Study_id, levels = Study_id))   # keep the new order

interactions<- data %>%
  semi_join(spatial_df_filtered, by = "Study_id")

#  Calculate pollinator, plant, and total species for each network
network_sizes <- interactions %>%
  dplyr::group_by(Study_id, Network_id) %>%
  dplyr::summarise(
    n_pollinators = dplyr::n_distinct(Pollinator_accepted_name), 
    n_plants = dplyr::n_distinct(Plant_accepted_name),           
    .groups = 'drop' 
  ) %>%
  # Calculate total species for the network
  dplyr::mutate(
    total_species = n_pollinators + n_plants
  )

# Apply a threshold: Keep only networks with 20 or more species
chosen_networks <- network_sizes %>%
  dplyr::filter(total_species >= 15) %>%
  dplyr::arrange(desc(total_species))

metadata_chosen <- metadata %>%
  # Make sure the key columns exist in metadata
  # (most metadata files have Study_id and Network_id – adjust if yours uses a different name)
  semi_join(chosen_networks, by = c("Study_id"))


study_summary <- metadata_chosen %>%
  group_by(Study_id) %>%
  summarise(
    total_networks = sum(number_of_networks, na.rm = TRUE),
    classification = case_when(
      any(!!sym(classification_col) == "spatial and temporal") ~ "spatial and temporal",
      any(!!sym(classification_col) == "temporal")            ~ "temporal",
      any(!!sym(classification_col) == "spatial")             ~ "spatial",
      TRUE                                                    ~ "unknown"
    ),
    .groups = "drop"
  ) %>%
  filter(total_networks > 10) %>%
  arrange(desc(total_networks)) %>%
  mutate(Study_id = factor(Study_id, levels = Study_id))


ggplot(study_summary,
       aes(x = Study_id, y = total_networks, fill = classification)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.4) +   # ← linewidth
  geom_text(aes(label = total_networks),
            vjust = -0.5, size = 3, colour = "black", fontface = "bold") +
  
scale_fill_manual(
  values = c(
    "spatial"              = "plum4",
    "spatial and temporal" = "lightcoral"
  ),
  name = "Classification"
) +
  
  labs(
    title = "Number of Networks per Study (spatial & spatial-and-temporal only)",
    x = "Study_id",
    y = "Number of Networks"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.title    = element_text(face = "bold"),
    axis.text.x   = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
    legend.position = "right",
    legend.title  = element_text(face = "bold"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  tme



#take the data i need
selected_study_ids <- unique(study_summary$Study_id)

data_filtered <- data %>%
  filter(Study_id %in% selected_study_ids)



interactions_binary <- data_filtered %>%
  select(Study_id, Network_id, Plant_accepted_name, Pollinator_accepted_name) %>%
  distinct() %>%
  unite("interaction", Plant_accepted_name, Pollinator_accepted_name, sep = "_") %>%
  mutate(presence = 1) %>%
  pivot_wider(
    names_from = Network_id,
    values_from = presence,
    values_fill = 0
  )


jaccard_results <- interactions_binary %>%
  group_by(Study_id) %>%
  nest() %>%
  mutate(
    jaccard_pairs = map(data, ~ {
      mat <- .x %>%
        select(-interaction) %>%
        as.matrix() %>%
        t()  
      
      net_names <- rownames(mat)
      
      if (nrow(mat) < 2) return(tibble())
      
      dist_mat <- 1-vegdist(mat, method = "jaccard")
      
      as.data.frame(as.matrix(dist_mat)) %>%
        rownames_to_column("Network1") %>%
        pivot_longer(-Network1, names_to = "Network2", values_to = "Jaccard") %>%
        filter(Network1 < Network2) %>%
        mutate(
          Network1 = net_names[as.numeric(Network1)],
          Network2 = net_names[as.numeric(Network2)]
        )
    })
  ) %>%
  select(Study_id, jaccard_pairs) %>%
  unnest(jaccard_pairs)

print(jaccard_results)



#Plot the Jaccard
jaccard_clean <- jaccard_results %>% 
  filter(!is.na(Jaccard))

ggplot(jaccard_clean, aes(x = Study_id, y = Jaccard)) +
  geom_boxplot(
    fill = "lightblue3", 
    color = "black", 
    outlier.colour = "darkseagreen", 
    outlier.shape = 16,
    outlier.alpha = 0.7
  ) +
  labs(
    title = "Jaccard index Between Networks Within Each Study",
    x = "Study ID",
    y = "Jaccard Index"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "gray50"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )+
  tme








#--------


# Spices per network

interactions <- data_filtered %>%
  semi_join(spatial_df_filtered, by = "Study_id")
interactions$Plant_accepted_name


#  Calculate pollinator, plant, and total species for each network
network_sizes <- interactions %>%
  dplyr::group_by(Study_id, Network_id) %>%
  dplyr::summarise(
    n_pollinators = dplyr::n_distinct(Pollinator_accepted_name), 
    n_plants = dplyr::n_distinct(Plant_accepted_name),           
    .groups = 'drop' 
  ) %>%
  # Calculate total species for the network
  dplyr::mutate(
    total_species = n_pollinators + n_plants
  )

# Apply a threshold: Keep only networks with 20 or more species
large_networks <- network_sizes %>%
  dplyr::filter(total_species >= 20) %>%
  dplyr::arrange(desc(total_species))


# Plotting a box plot to show the min, max, median and quartiles for each study
ggplot(large_networks, aes(x = reorder(Study_id, -total_species, FUN = median), y = total_species)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5, color = "orange") +
  theme_minimal() +
  labs(
    title = " Network Sizes (>= 20 species) for each Study",
    x = "Study ID",
    y = "Total Species (Pollinators + Plants)"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(), # Remove grid lines
    panel.grid.minor = element_blank()
  ) +
  tme


# calculate connectence -----

connectance_results <- data_filtered %>%
  select(Study_id, Network_id, Plant_accepted_name, Pollinator_accepted_name) %>%
  distinct() %>%
  group_by(Study_id, Network_id) %>%
  summarise(
    n_interactions = n(),
    n_plants = n_distinct(Plant_accepted_name),
    n_pollinators = n_distinct(Pollinator_accepted_name),
    .groups = "drop"
  ) %>%
  mutate(
    Connectance = n_interactions / (n_plants * n_pollinators)
  ) %>%
  select(Study_id, Network_id, Connectance)

ggplot(connectance_results, 
       aes(x = reorder(Study_id, -Connectance, FUN = median),  
           y = Connectance)) +
  geom_boxplot(fill = "white", color = "black") +
  geom_jitter(width = 0.1, alpha = 0.5, color = "lightcyan3") + 
  theme_minimal() +
  labs(
    title = "Network Connectance for each Study",
    x = "Study ID",
    y = "Connectance"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  tme +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1)
  )














