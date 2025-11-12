setwd("/Users/shirnehoray/EDA")


# 1. Load packages and data -----

library(readxl)      # to read Excel files
library(dplyr)       # for data wrangling
library(igraph)      # for network analysis
library(bipartite)   # specialized for ecological networks
library(ggplot2)     # for plots
library(tibble)
library(tidyr)


# Read the data sheets ----
net_pa <- read_excel("Kaiser_Bunbury_et_al_2017.xlsx", sheet = "64 networks_no.visits")      # number of visits
net_freq <- read_excel("Kaiser_Bunbury_et_al_2017.xlsx", sheet = "64 networks_visitfreq")    # frequency
plants <- read_excel("Kaiser_Bunbury_et_al_2017.xlsx",
                     sheet = "Plant species",
                     skip = 1,      # skips the note row
                     col_names = TRUE)

pollinators <- read_excel("Kaiser_Bunbury_et_al_2017.xlsx",
                          sheet = "Pollinator species",
                          skip = 1,  # skips the note row
                          col_names = TRUE)
info <- read_excel("Kaiser_Bunbury_et_al_2017.xlsx", sheet = "info")



# Basic data set inspection----

head(net_pa)
unique(net_pa$Site)   # sites names (should be 8 sites)
unique(net_pa$Month)  # months

# Check is there are missing values
sum(is.na(net_pa))
sum(is.na(net_freq))



# Species amount ----

# Function that calaulate the number of pollinators and palantsin every web

get_richness <- function(df){
  # מוצאים את תחילת עמודות המאביקים בצורה דינמית (ללא "קסם 7")
  poll_start <- which(names(df) == "Floral abundance") + 1
  stopifnot(length(poll_start) == 1)  # ביטחון: חייב להיות עמודה אחת כזו
  
  # בונים מטריצת אינטראקציות בלבד (צמחים × מאביקים)
  mat <- as.matrix(df[, poll_start:ncol(df)])
  storage.mode(mat) <- "numeric"        # לוודא שהכול מספרי
  rownames(mat) <- df$`Plant species ID`
  
  # עושר "ברשת" = כמה צמחים קיבלו ביקור לפחות פעם אחת
  n_plants_active <- sum(rowSums(mat, na.rm = TRUE) > 0)
  # עושר "פוטנציאלי" = כמה צמחים מופיעים בטבלה (גם אם לא ביקרו בהם)
  n_plants_present <- dplyr::n_distinct(df$`Plant species ID`)
  
  # עושר מאביקים (כמה מיני מאביקים הופיעו בפועל)
  n_polls <- sum(colSums(mat, na.rm = TRUE) > 0)
  
  data.frame(
    Site        = unique(df$Site),
    Month       = unique(df$Month),
    NetworkID   = unique(df$`Network ID`),
    Plants_active  = n_plants_active,   # צמחים שקיבלו ביקור
    Plants_present = n_plants_present,  # כל הצמחים המופיעים ברשת (גם 0 ביקורים)
    Pollinators    = n_polls
  )
}


# calculate the richness for all the networks 
richness_results <- net_pa %>%
  group_by(Site, Month, `Network ID`) %>%
  group_split() %>%
  lapply(get_richness) %>%
  bind_rows()


richness_long <- richness_results %>%
  tidyr::pivot_longer(cols = c("Plants_active","Pollinators"),
                      names_to = "Group", values_to = "Richness")

# Genom point
ggplot(richness_long, aes(x = Month, y = Richness, color = Group)) +
  geom_line(aes(group = interaction(Site, Group))) +
  geom_point() +
  facet_wrap(~Site) +
  labs(title = "Number of Plants and Pollinators",
       y = "Number of Species") +
  theme_minimal()

# Histograms
ggplot(richness_long, aes(x = factor(Month), y = Richness, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Site) +
  labs(x = "Month", y = "Number of Species", fill = "Group") +
  theme_minimal()

# Heatmap
ggplot(richness_long, aes(x = factor(Month), y = Site, fill = Richness)) +
  geom_tile(color = "white") +
  facet_wrap(~Group) +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(x = "Month", y = "Site", title = "Species Richness Per Site") +
  theme_minimal()


# 4. species simillarities between sites 

# pairwise matrix function 
pairwise_shared <- function(incidence){
  mat <- as.matrix(incidence[ , -1, drop = FALSE])   # remove site columns
  rownames(mat) <- incidence$Site
  inter <- mat %*% t(mat)                            # pairwise number of shared species between sites
  tot   <- rowSums(mat)                              # number of total species in קש site 
  uni   <- outer(tot, tot, "+") - inter              # union
  jacc  <- inter / pmax(uni, 1)                      # Jaccard
  list(shared = inter, jaccard = jacc)
}

# create an heatmap
plot_heat <- function(M, title){
  as.data.frame(M) |>
    rownames_to_column("Site1") |>
    pivot_longer(-Site1, names_to = "Site2", values_to = "value") |>
    ggplot(aes(Site1, Site2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    coord_equal() +
    labs(title = title, x = "Site", y = "Site", fill = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Get pollinators ID
poll_start <- which(names(net_pa) == "Floral abundance") + 1
poll_cols  <- names(net_pa)[poll_start:ncol(net_pa)]

#  Active plants that got visit in at least one of the months
plants_active_by_site <- net_pa %>%
  group_by(Site, `Plant species ID`) %>%
  summarise(active = as.numeric(
    any(rowSums(across(all_of(poll_cols)), na.rm = TRUE) > 0)),
    .groups = "drop") %>%
  pivot_wider(names_from = `Plant species ID`,
              values_from = active, values_fill = 0) %>%
  arrange(Site)

# Jaccard
P <- pairwise_shared(plants_active_by_site)

# Heatmap plot
plot_heat(P$jaccard, "Plants — Jaccard similarity (aggregated across months)")

# Pollinator presents 
polls_active_by_site <- net_pa %>%
  group_by(Site) %>%
  summarise(across(all_of(poll_cols), ~ as.numeric(sum(.x, na.rm = TRUE) > 0)),
            .groups = "drop") %>%
  arrange(Site)

# Jaccard
Q <- pairwise_shared(polls_active_by_site)

# Heatmap plot
plot_heat(Q$jaccard, "Pollinators — Jaccard similarity (aggregated across months)")







