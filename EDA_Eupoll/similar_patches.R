# ============================
# Libraries

# ============================
library(dplyr)
library(readr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf) # for coord_sf CRS handling and sf objects
library(tidyr)
library(purrr)

# ============================
# Parameters (easy to tweak)
# ============================
bio <- "Mediterranean"      # Bioregion to focus on
min_total <- 15     # Minimum total species per network (plants + pollinators)
top_n_pairs <- 50    # How many top-similar network pairs to keep

cat(sprintf("Processing Bioregion for Pollinators: %s\n", bio))
cat(sprintf("Filtering networks with total species >= %d\n", min_total))

# ============================
# Read data
# Notes:
# 1) Don't over-specify col_types for columns that might not exist in all files.
# 2) Enforce key columns' types after reading.
# ============================
interactions <- read_csv("Interaction_data.csv",
                         show_col_types = FALSE)

# Enforce expected types where applicable (only if columns exist)
force_character <- intersect(c("Network_id", "Bioregion",
                               "Plant_accepted_name", "Pollinator_accepted_name",
                               "Country"), names(interactions))
force_double    <- intersect(c("Latitude", "Longitude"), names(interactions))

interactions <- interactions %>%
  mutate(across(all_of(force_character), as.character)) %>%
  mutate(across(all_of(force_double), as.double))

# Create a robust "Site" column by coalescing any available site-like fields.
# This safely handles missing columns (uses only those that exist).
site_candidates <- intersect(
  c("EuPPollNet_habitat", "Site_name", "Habitat", "Location"),
  names(interactions))



interactions <- read_csv("Interaction_data.csv", show_col_types = FALSE) %>%
    mutate(
      Network_id            = as.character(Network_id),
      Bioregion             = as.character(Bioregion),
      Plant_accepted_name   = as.character(Plant_accepted_name),
      Pollinator_accepted_name = as.character(Pollinator_accepted_name),
      Latitude              = suppressWarnings(as.double(Latitude)),
      Longitude             = suppressWarnings(as.double(Longitude)),
      Country               = as.character(Country)
    )
  
  # ---------- Jaccard ----------
  jaccard_sets <- function(a, b) {
    ia <- length(base::intersect(a, b))
    ua <- length(base::union(a, b))
    if (ua == 0) NA_real_ else ia / ua
  }
  

  
  # ---------- number of species per network ----------
  species_per_network <- interactions %>%
    group_by(Bioregion, Network_id) %>%
    summarise(
      n_plants       = n_distinct(Plant_accepted_name),
      n_pollinators  = n_distinct(Pollinator_accepted_name),
      .groups = "drop"
    ) %>%
    mutate(total_species = n_plants + n_pollinators) %>%
    filter(total_species >= min_total)
  
  sub_data <- species_per_network %>% filter(Bioregion == bio)
  nets <- sub_data$Network_id
  
  if (length(nets) == 0) {
    stop(sprintf("No networks in %s with at least %d total species.", bio, min_total))
  }
  
  # # ---------- Build EDGE sets (and optional weights) per network ----------
  # # Edge key format keeps guild identity explicit to avoid name collisions.
  # edge_key <- function(p, a) paste0("EDGE::PLANT::", p, " | POLL::", a)
  # 
  # dat_bio <- interactions %>%
  #   filter(Bioregion == bio, Network_id %in% nets) %>%
  #   # Drop rows with missing species names (cannot form an edge)
  #   filter(!is.na(Plant_accepted_name), !is.na(Pollinator_accepted_name)) %>%
  #   mutate(edge = edge_key(Plant_accepted_name, Pollinator_accepted_name))
  # 
  # # unique edges per network
  # edges_tbl <- dat_bio %>%
  #   group_by(Network_id) %>%
  #   summarise(items = list(unique(edge)), .groups = "drop")
  
  # ---------- Compute pairwise similarity matrix ----------
  net_ids <- edges_tbl$Network_id
  n <- length(net_ids)
  sim_mat <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(net_ids, net_ids))
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) {
        sim_mat[i, j] <- NA_real_
      } else {
        it_i <- edges_tbl$items[[i]]
        it_j <- edges_tbl$items[[j]]
        wi   <- edges_tbl$weights[[i]]
        wj   <- edges_tbl$weights[[j]]
        
        sim_mat[i, j] <- if (!is.null(wi) && !is.null(wj)) {
          # Weighted Jaccard on edges
          weighted_jaccard(wi, wj)
        } else {
          # Plain Jaccard on edge sets
          jaccard_sets(it_i, it_j)
        }
      }
    }
  }
  
  # ---------- Flatten to unique pairs (upper triangle) and pick top-N ----------
  sim_df <- as.data.frame(as.table(sim_mat), stringsAsFactors = FALSE) %>%
    filter(!is.na(Freq), as.character(Var1) < as.character(Var2)) %>%
    arrange(desc(Freq)) %>%
    slice_head(n = top_n_pairs)
  
  if (nrow(sim_df) == 0) {
    stop("No off-diagonal similarities found (all NA).")
  }
  
  cat("Top similar network pairs by EDGE Jaccard:\n")
  print(head(sim_df, 20), row.names = FALSE)
  
  # ---------- Collect unique networks from top pairs (to reuse in your map, etc.) ----------
  unique_nets <- unique(c(as.character(sim_df$Var1), as.character(sim_df$Var2)))
  cat("Unique networks appearing in top pairs:\n",
      paste(unique_nets, collapse = ", "), "\n")
  
  # ---------- Summarize locations & sizes for those networks (for mapping) ----------
  first_unique <- function(x) {
    ux <- unique(x)
    ux[!is.na(ux)][1] %||% NA
  }
  
  top_network_locs <- interactions %>%
    filter(Bioregion == bio, Network_id %in% unique_nets) %>%
    group_by(Network_id) %>%
    summarise(
      Latitude       = first_unique(Latitude),
      Longitude      = first_unique(Longitude),
      Country        = first_unique(Country),
      n_plants       = n_distinct(Plant_accepted_name),
      n_pollinators  = n_distinct(Pollinator_accepted_name),
      .groups = "drop"
    ) %>%
    mutate(total_species = n_plants + n_pollinators)
  
  if (any(is.na(top_network_locs$Latitude)) || any(is.na(top_network_locs$Longitude))) {
    warning("Some networks lack coordinates; the map may omit them.")
  }
  
  # ---------- Map (same pattern as before) ----------
  lat_range <- range(top_network_locs$Latitude, na.rm = TRUE)
  lon_range <- range(top_network_locs$Longitude, na.rm = TRUE)
  lat_buffer <- 1; lon_buffer <- 2
  xlim <- c(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer)
  ylim <- c(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer)
  
  europe <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf", continent = "Europe")
  
  cols <- grDevices::hcl.colors(length(unique_nets), "Spectral", rev = TRUE)
  
  p_map <- ggplot(europe) +
  geom_sf(fill = "grey90", color = "white") +
  geom_point(
    data = top_network_locs,
    aes(x = Longitude, y = Latitude, color = Network_id, size = total_species),
    alpha = 0.8
  ) +
  scale_color_manual(values = cols) +
  scale_size_continuous(range = c(2, 8), name = "Total Species") +
  coord_sf(
    xlim = xlim,
    ylim = ylim,
    expand = FALSE,
    default_crs = NULL   # <<< הוספתי את זה
  ) +
  labs(
    title = sprintf("Most similar networks in %s (Jaccard on edges)", bio),
    subtitle = sprintf("Similarity computed on interaction sets; Lat range: %.2f–%.2f",
                       lat_range[1], lat_range[2]),
    x = "Longitude", y = "Latitude", color = "Network_id"
  ) +
  theme_minimal() + theme(legend.position = "bottom")

  print(p_map)
  # ggsave("boreal_top_similar_networks_edges.png", p_map, width = 10, height = 7, dpi = 300)
  
