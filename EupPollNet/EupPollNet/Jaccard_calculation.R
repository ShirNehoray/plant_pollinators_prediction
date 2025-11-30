library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(dplyr)
library(lubridate)
library(dplyr)
library(readr)
library(ggplot2)


Interacion_data_filtered <- read.csv("./csv/Interactions_data_filtered_studies.csv", encoding = "latin1", stringsAsFactors = FALSE)
tme <-  theme(axis.text = element_text(size = 10, color = "black"),
              axis.title = element_text(size = 16, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
theme_set(theme_bw())


studies <- unique(Interacion_data_filtered$Study_id)
res_list <- vector("list", length(studies))


for (i in seq_along(studies)) {
  s <- studies[i]
  
  dat_s <- Interacion_data_filtered %>%
    filter(Study_id == s)
  
  # set interaction_binary for this single study
  interactions_binary_s <- dat_s %>%
    select(Network_id, Plant_accepted_name, Pollinator_accepted_name) %>%
    distinct() %>%
    unite("interaction", Plant_accepted_name, Pollinator_accepted_name, sep = "_") %>%
    mutate(presence = 1L) %>%
    pivot_wider(
      names_from  = Network_id,
      values_from = presence,
      values_fill = 0
    )
  
  #  sets of interactions per network
  net_sets <- interactions_binary_s %>%
    pivot_longer(
      cols      = -interaction,
      names_to  = "Network",
      values_to = "value"
    ) %>%
    filter(value == 1L) %>%
    group_by(Network) %>%
    summarise(interactions = list(interaction), .groups = "drop")
  
  if (nrow(net_sets) < 2) {
    res_list[[i]] <- tibble(
      Study_id = s,
      Network1 = character(0),
      Network2 = character(0),
      Jaccard  = numeric(0)
    )
    next
  }
  
  pairs <- combn(net_sets$Network, 2, simplify = FALSE)
  
  jaccard_s <- map_dfr(pairs, function(p) {
    set1 <- net_sets$interactions[[which(net_sets$Network == p[1])]]
    set2 <- net_sets$interactions[[which(net_sets$Network == p[2])]]
    
    inter <- length(intersect(set1, set2))
    uni   <- length(union(set1, set2))
    
    tibble(
      Study_id = s,
      Network1 = p[1],
      Network2 = p[2],
      Jaccard  = ifelse(uni == 0, NA_real_, inter / uni)
    )
  })
  
  res_list[[i]] <- jaccard_s
}

jaccard_results <- bind_rows(res_list)

ggplot(jaccard_results, aes(x = Study_id, y = Jaccard)) +
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


# Jaccard for species ------

compute_jaccard_species <- function(species_df) {
  # species_df has: Network_id, species
  net_sets <- species_df %>%
    group_by(Network_id) %>%
    summarise(species = list(unique(species)), .groups = "drop")
  
  if(nrow(net_sets) < 2) return(tibble())
  
  pairs <- combn(net_sets$Network_id, 2, simplify = FALSE)
  
  map_dfr(pairs, function(p) {
    set1 <- net_sets$species[[which(net_sets$Network_id == p[1])]]
    set2 <- net_sets$species[[which(net_sets$Network_id == p[2])]]
    
    inter <- length(intersect(set1, set2))
    uni   <- length(union(set1, set2))
    
    tibble(
      Network1 = p[1],
      Network2 = p[2],
      Jaccard  = ifelse(uni == 0, NA, inter / uni)
    )
  })
}
species_jaccard_results <- Interacion_data_filtered %>%
  group_by(Study_id) %>%
  group_modify(~{
    
    # --- Plants ---
    plant_df <- .x %>%
      select(Network_id, species = Plant_accepted_name) %>%
      distinct()
    
    plant_jaccard <- compute_jaccard_species(plant_df) %>%
      mutate(Type = "Plant")
    
    # --- Pollinators ---
    poll_df <- .x %>%
      select(Network_id, species = Pollinator_accepted_name) %>%
      distinct()
    
    poll_jaccard <- compute_jaccard_species(poll_df) %>%
      mutate(Type = "Pollinator")
    
    bind_rows(plant_jaccard, poll_jaccard)
  }) %>%
  ungroup()


ggplot(species_jaccard_results, aes(x = Study_id, y = Jaccard, fill = Type)) +
  geom_boxplot(
    color = "black",
    outlier.colour = "darkseagreen",
    outlier.shape = 16,
    outlier.alpha = 0.7
  ) +
  labs(
    title = "Species Jaccard Similarity Between Networks",
    subtitle = "Comparison of Plant and Pollinator overlap",
    x = "Study ID",
    y = "Jaccard Similarity"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  scale_fill_manual(values = c("Plant" = "lightblue3", "Pollinator" = "lightcoral")) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "gray50"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "top"
  )+
  tme



