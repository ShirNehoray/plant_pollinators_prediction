# Network analysis after the aggregation 

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# --- summary table for the amount of species ---
sizes_df <- lapply(names(site_mats_freq), function(s){
  m <- site_mats_freq[[s]]
  data.frame(Site = s,
             Plants = nrow(m),
             Pollinators = ncol(m),
             Links = sum((m > 0)))
}) %>% bind_rows()

# --- long table ---
sizes_long <- sizes_df %>%
  select(Site, Plants, Pollinators) %>%
  pivot_longer(c(Plants, Pollinators), names_to = "Group", values_to = "Count")

# --- histogram bar ---
ggplot(sizes_long, aes(x = reorder(Site, Count), y = Count, fill = Group)) +
  geom_col(position = "dodge") +          # # שתי עמודות זו-לצד-זו לכל אתר
  coord_flip() +
  labs(title = "Network size per site", x = "Site", y = "Number of species", fill = "") +
  theme_minimal()



# ============= Species site-frequency =============
# Count in how many sites each species occurs

library(dplyr)
library(ggplot2)
library(tidyr)

# --- Plants ---
plant_sets <- lapply(site_mats_freq, function(m){
  if (is.null(m) || length(m) == 0) return(character(0))
  rn <- rownames(m)
  rn[rowSums(m > 0, na.rm = TRUE) > 0]
})
plant_counts <- table(unlist(plant_sets))                  # count sites per plant
plant_df <- data.frame(Species = names(plant_counts),
                       n_sites = as.integer(plant_counts),
                       row.names = NULL) %>%
  arrange(desc(n_sites))

# --- Pollinators ---
poll_sets <- lapply(site_mats_freq, function(m){
  if (is.null(m) || length(m) == 0) return(character(0))
  cn <- colnames(m)
  cn[colSums(m > 0, na.rm = TRUE) > 0]
})
poll_counts <- table(unlist(poll_sets))                    # count sites per pollinator
poll_df <- data.frame(Species = names(poll_counts),
                      n_sites = as.integer(poll_counts),
                      row.names = NULL) %>%
  arrange(desc(n_sites))

# --- Plot plants ---
ggplot(plant_df, aes(x = reorder(Species, n_sites), y = n_sites)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Plants: number of sites per species",
       x = "Plant species", y = "Number of Sites") +
  theme_minimal()

# --- Top-30 plants ---
topN <- 30
ggplot(head(plant_df, topN),
       aes(x = reorder(Species, n_sites), y = n_sites)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = paste("Plants: number of sites per species (Top", topN, ")"),
       x = "Plant species", y = "Number of Sites") +
  theme_minimal()

# --- Plot pollinators ---
ggplot(poll_df, aes(x = reorder(Species, n_sites), y = n_sites)) +
  geom_col(fill = "darkorange") +
  coord_flip() +
  labs(title = "Pollinators: number of sites per species",
       x = "Pollinator species", y = "Number of Sites") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 3))   # אפשר לשחק עם המספר

# --- Top-30 pollinators ---
ggplot(head(poll_df, topN),
       aes(x = reorder(Species, n_sites), y = n_sites)) +
  geom_col(fill = "darkorange") +
  coord_flip() +
  labs(title = paste("Pollinators: number of sites per species (Top", topN, ")"),
       x = "Pollinator species", y = "Number of Sites") +
  theme_minimal()

