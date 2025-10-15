## ============================
## Parameters
## ============================

library(dplyr)
library(readr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf) # for coord_sf CRS handling and sf objects
library(tidyr)
library(purrr)

bio <- "Continental"     # bioregion to analyze
min_total <- 15     # minimum total species per network
top_n_pairs <- 100   # how many top pairs to keep (optional)

## ============================
## Read data 
## ============================
interactions <- read.csv("Interaction_data.csv", stringsAsFactors = FALSE)

## ============================
## Utilities (base R)
## ============================
normalize_name <- function(x) {
  # Optional: standardize species names (trim, collapse spaces, lower-case)
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)
  tolower(x)
}

jaccard_sets <- function(a, b) {
  a <- unique(a); b <- unique(b)
  ia <- length(intersect(a, b))
  ua <- length(union(a, b))
  if (ua == 0) NA_real_ else ia / ua
}

## ============================
## 0) Basic column checks
## ============================
req_cols <- c("Bioregion", "Network_id", "Plant_accepted_name", "Pollinator_accepted_name")
missing_cols <- setdiff(req_cols, names(interactions))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

## ============================
## 1) Filter to chosen bioregion
## ============================
dat_bio <- interactions[interactions$Bioregion == bio, , drop = FALSE]

if (nrow(dat_bio) == 0) stop(paste("No rows for Bioregion:", bio))

## ============================
## 2) Compute species counts per network (to filter small networks)
##    We need distinct counts of plants & pollinators per Network_id
## ============================
# Distinct plants per network
plants_by_net <- aggregate(Plant_accepted_name ~ Network_id, data = dat_bio,
                           FUN = function(x) length(unique(x)))
names(plants_by_net)[2] <- "n_plants"

# Distinct pollinators per network
polls_by_net <- aggregate(Pollinator_accepted_name ~ Network_id, data = dat_bio,
                          FUN = function(x) length(unique(x)))
names(polls_by_net)[2] <- "n_pollinators"

# Merge counts and compute total
counts <- merge(plants_by_net, polls_by_net, by = "Network_id", all = TRUE)
counts$n_plants[is.na(counts$n_plants)] <- 0
counts$n_pollinators[is.na(counts$n_pollinators)] <- 0
counts$total_species <- counts$n_plants + counts$n_pollinators

# Keep networks meeting the threshold
nets_keep <- counts$Network_id[counts$total_species >= min_total]
if (length(nets_keep) == 0) {
  stop(paste("No networks with total_species >=", min_total, "in", bio))
}

## ============================
## 3) Build species sets per network (plants + pollinators together)
## ============================
# Subset to kept networks
dat_bio_sub <- dat_bio[dat_bio$Network_id %in% nets_keep, , drop = FALSE]

# Normalize names (optional but recommended)
plant_names <- normalize_name(dat_bio_sub$Plant_accepted_name)
poll_names  <- normalize_name(dat_bio_sub$Pollinator_accepted_name)

# Tag guilds to avoid accidental collisions
plant_tags <- paste0("PLANT::", plant_names)
poll_tags  <- paste0("POLL::",  poll_names)

# Split row indices by Network_id so we can collect species per network
idx_by_net <- split(seq_len(nrow(dat_bio_sub)), dat_bio_sub$Network_id)

# For each network, take the union of unique plant & pollinator tags
species_items <- lapply(idx_by_net, function(idx) {
  unique(c(plant_tags[idx], poll_tags[idx]))
})
# 'species_items' is a named list: names(species_items) are Network_id values

net_ids <- names(species_items)
n <- length(net_ids)

## ============================
## 4) Pairwise Jaccard matrix (upper triangle fill)
## ============================
sim_mat <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(net_ids, net_ids))

for (i in seq_len(n)) {
  # diagonal
  sim_mat[i, i] <- NA_real_
  # upper triangle only
  if (i < n) {
    for (j in (i + 1L):n) {
      val <- jaccard_sets(species_items[[i]], species_items[[j]])
      sim_mat[i, j] <- val
      sim_mat[j, i] <- val
    }
  }
}

## ============================
## 5) Tidy table of unique pairs, sorted high-to-low, with Rank
## ============================
# Convert matrix to data.frame of pairs
# We'll keep only i<j to avoid duplicates
pairs_list <- list()
k <- 0L
for (i in seq_len(n)) {
  if (i < n) {
    for (j in (i + 1L):n) {
      k <- k + 1L
      pairs_list[[k]] <- data.frame(
        Net1 = net_ids[i],
        Net2 = net_ids[j],
        Jaccard = sim_mat[i, j],
        stringsAsFactors = FALSE
      )
    }
  }
}
sim_df <- do.call(rbind, pairs_list)

# Remove NA (shouldn't be any except degenerate cases) and sort descending
sim_df <- sim_df[!is.na(sim_df$Jaccard), , drop = FALSE]
ord <- order(sim_df$Jaccard, decreasing = TRUE)
sim_df <- sim_df[ord, , drop = FALSE]

# Add Rank
sim_df$Rank <- seq_len(nrow(sim_df))

# Reorder columns
sim_df <- sim_df[, c("Rank", "Net1", "Net2", "Jaccard")]

# Keep top N pairs (optional)
if (!is.null(top_n_pairs) && is.finite(top_n_pairs)) {
  top_n <- min(top_n_pairs, nrow(sim_df))
  sim_top <- sim_df[seq_len(top_n), , drop = FALSE]
} else {
  sim_top <- sim_df
}

## ============================
## 6) Inspect results
## ============================
print(head(sim_top, 20), row.names = FALSE)

## If you want the unique networks appearing in the top pairs:
unique_nets <- unique(c(sim_top$Net1, sim_top$Net2))

## ============================
## 7) Collect per-network coordinates for the networks in top pairs
## ============================
# unique networks from the top pairs
unique_nets <- unique(c(sim_top$Net1, sim_top$Net2))

# subset interactions to the chosen bioregion and those networks
dat_map <- dat_bio_sub[dat_bio_sub$Network_id %in% unique_nets, , drop = FALSE]

# basic sanity
if (!all(c("Latitude","Longitude") %in% names(dat_map))) {
  stop("Latitude/Longitude columns are missing in the dataset.")
}

# Ensure numeric (safely)
dat_map$Latitude  <- suppressWarnings(as.numeric(dat_map$Latitude))
dat_map$Longitude <- suppressWarnings(as.numeric(dat_map$Longitude))

# remove rows with missing coordinates
dat_map <- dat_map[!is.na(dat_map$Latitude) & !is.na(dat_map$Longitude), , drop = FALSE]
if (nrow(dat_map) == 0) stop("No valid coordinates to plot.")

## For each Network_id, take a single representative lat/lon (first unique non-NA)
first_unique <- function(x) {
  ux <- unique(x)
  ux[!is.na(ux)][1]
}

# aggregate coordinates per network (base R)
lat_by_net <- aggregate(Latitude  ~ Network_id, data = dat_map, FUN = first_unique)
lon_by_net <- aggregate(Longitude ~ Network_id, data = dat_map, FUN = first_unique)

# merge back
net_locs <- merge(lat_by_net, lon_by_net, by = "Network_id", all = TRUE)

# (optional) add total species to size points (we already computed 'counts' earlier)
net_locs <- merge(net_locs, counts[, c("Network_id","total_species")], by = "Network_id", all.x = TRUE)

# keep only networks we actually plot
net_locs <- net_locs[net_locs$Network_id %in% unique_nets, , drop = FALSE]

## ============================
## 8) Determine map extents with a small buffer
## ============================
lat_range <- range(net_locs$Latitude,  na.rm = TRUE)
lon_range <- range(net_locs$Longitude, na.rm = TRUE)
lat_buffer <- 1
lon_buffer <- 2
xlim <- c(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer)
ylim <- c(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer)

## ============================
## 9) Base map (Europe) and plot
## ============================
## ============================
## 9) Base map (Europe) and plot
## ============================
europe <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf", continent = "Europe")

# simple palette for distinct networks
n_cols <- length(unique(net_locs$Network_id))
cols <- grDevices::hcl.colors(n_cols, palette = "Spectral", rev = TRUE)

p_map <- ggplot(europe) +
  geom_sf(fill = "grey90", color = "white") +
  geom_point(
    data = net_locs,
    aes(x = Longitude, y = Latitude, color = Network_id, size = total_species),
    alpha = 0.8
  ) +
  scale_color_manual(values = cols, guide = "none") +  # Remove Network_id legend
  scale_size_continuous(range = c(2, 8), name = "Total Species") +
  # IMPORTANT for newer ggplot2: allow numeric xlim/ylim on lon/lat
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, default_crs = NULL) +
  labs(
    title = sprintf("Top similar networks in %s (species-based Jaccard)", bio),
    subtitle = sprintf("Points sized by total species; Lat range: %.2f–%.2f",
                       lat_range[1], lat_range[2]),
    x = "Longitude", y = "Latitude", color = "Network_id"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_map)
