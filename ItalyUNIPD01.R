
##  Packages
library(tidyr)
library(purrr)
library(ggplot2)
library(dplyr)
library(readr)
library(rnaturalearth)
library(tidyverse)

##  Settings
bioregions   <- c("Continental")


##  Read data
interactions <- read_csv("Interaction_data.csv",
                         col_types = cols(Network_id = col_character()))

interactions <- interactions %>%
  filter(str_detect(Network_id, "UNIPD01|UNIPD02|UNIPD03"))


##  Utilities (base R)
normalize_name <- function(x) {
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


##  Loop over bioregions
for (bio in bioregions) {
  cat(sprintf("\n=== Processing Bioregion: %s ===\n", bio))
  
  req_cols <- c("Bioregion", "Network_id",
                "Plant_accepted_name", "Pollinator_accepted_name")
  missing_cols <- setdiff(req_cols, names(interactions))
  if (length(missing_cols) > 0) {
    cat(paste("Missing required columns:", paste(missing_cols, collapse = ", "), "\n"))
    next
  }
  
  ## 1) Filter to chosen bioregion
  dat_bio <- interactions[interactions$Bioregion == bio, , drop = FALSE]
  if (nrow(dat_bio) == 0) {
    cat(paste("No rows for Bioregion:", bio, "\n"))
    next
  }
  

  ## 2) Species counts per network
  plants_by_net <- aggregate(Plant_accepted_name ~ Network_id,
                             data = dat_bio,
                             FUN = function(x) length(unique(x)))
  names(plants_by_net)[2] <- "n_plants"
  
  polls_by_net <- aggregate(Pollinator_accepted_name ~ Network_id,
                            data = dat_bio,
                            FUN = function(x) length(unique(x)))
  names(polls_by_net)[2] <- "n_pollinators"
  
  counts <- merge(plants_by_net, polls_by_net,
                  by = "Network_id", all = TRUE)
  counts$n_plants[is.na(counts$n_plants)]       <- 0
  counts$n_pollinators[is.na(counts$n_pollinators)] <- 0
  counts$total_species <- counts$n_plants + counts$n_pollinators
  
  # NO FILTER: use all networks
  nets_keep <- counts$Network_id
  
  if (length(nets_keep) == 0) {
    cat(paste("No networks in", bio, "\n"))
    next
  }
  
  ## 3) Build species sets per network (plants + pollinators)
  dat_bio_sub <- dat_bio[dat_bio$Network_id %in% nets_keep, , drop = FALSE]
  
  plant_names <- normalize_name(dat_bio_sub$Plant_accepted_name)
  poll_names  <- normalize_name(dat_bio_sub$Pollinator_accepted_name)
  
  plant_tags <- paste0("PLANT::", plant_names)
  poll_tags  <- paste0("POLL::",  poll_names)
  
  idx_by_net <- split(seq_len(nrow(dat_bio_sub)), dat_bio_sub$Network_id)
  
  species_items <- lapply(idx_by_net, function(idx) {
    unique(c(plant_tags[idx], poll_tags[idx]))
  })
  net_ids <- names(species_items)
  n <- length(net_ids)
  
  ## 4) Pairwise Jaccard matrix (upper triangle)
  sim_mat <- matrix(NA_real_, nrow = n, ncol = n,
                    dimnames = list(net_ids, net_ids))
  
  for (i in seq_len(n)) {
    sim_mat[i, i] <- NA_real_
    if (i < n) {
      for (j in (i + 1L):n) {
        val <- jaccard_sets(species_items[[i]], species_items[[j]])
        sim_mat[i, j] <- val
        sim_mat[j, i] <- val
      }
    }
  }
  
  ## 5) Tidy table of unique pairs, sorted high-to-low
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
  
  if (is.null(sim_df) || nrow(sim_df) == 0 ||
      !("Jaccard" %in% names(sim_df)) || all(is.na(sim_df$Jaccard))) {
    cat(sprintf("No valid Jaccard similarities for %s (empty or all NA).\n", bio))
    next
  }
  
  sim_df <- sim_df[!is.na(sim_df$Jaccard), , drop = FALSE]
  ord <- order(sim_df$Jaccard, decreasing = TRUE)
  sim_df <- sim_df[ord, , drop = FALSE]
  
  # NO top_n_pairs: use ALL pairs
  sim_top <- sim_df
  
  ## 6) Inspect results (only top 10 to avoid clutter)
  cat(sprintf("\n=== Top 10 similar network pairs for %s ===\n", bio))
  print(head(sim_top, 10), row.names = FALSE)
  if (nrow(sim_top) > 10) {
    cat(sprintf("... and %d more pairs (total: %d)\n", 
                nrow(sim_top) - 10, nrow(sim_top)))
  }
  
  unique_nets <- unique(c(sim_top$Net1, sim_top$Net2))
  

  ## 7) Coordinates for networks in top pairs
  dat_map <- dat_bio_sub[dat_bio_sub$Network_id %in% unique_nets, , drop = FALSE]
  
  if (!all(c("Latitude","Longitude") %in% names(dat_map))) {
    cat("Latitude/Longitude columns are missing in the dataset.\n")
    next
  }
  
  dat_map$Latitude  <- suppressWarnings(as.numeric(dat_map$Latitude))
  dat_map$Longitude <- suppressWarnings(as.numeric(dat_map$Longitude))
  
  dat_map <- dat_map[!is.na(dat_map$Latitude) & !is.na(dat_map$Longitude), , drop = FALSE]
  if (nrow(dat_map) == 0) {
    cat("No valid coordinates to plot.\n")
    next
  }
  
  first_unique <- function(x) {
    ux <- unique(x)
    ux[!is.na(ux)][1]
  }
  
  lat_by_net <- aggregate(Latitude  ~ Network_id, data = dat_map,
                          FUN = first_unique)
  lon_by_net <- aggregate(Longitude ~ Network_id, data = dat_map,
                          FUN = first_unique)
  
  net_locs <- merge(lat_by_net, lon_by_net, by = "Network_id", all = TRUE)
  net_locs <- merge(net_locs,
                    counts[, c("Network_id","total_species")],
                    by = "Network_id", all.x = TRUE)
  net_locs <- net_locs[net_locs$Network_id %in% unique_nets, , drop = FALSE]
  

  ## 8) Map extents 
  lat_range <- range(net_locs$Latitude,  na.rm = TRUE)
  lon_range <- range(net_locs$Longitude, na.rm = TRUE)
  lat_buffer <- 1
  lon_buffer <- -3
  xlim <- c(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer)
  ylim <- c(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer)
  
  ## 9) Base map (Europe) + plot + NETWORK LABELS 
  europe <- rnaturalearth::ne_countries(scale = "medium",
                                        returnclass = "sf",
                                        continent = "Europe")
  
  n_cols <- length(unique(net_locs$Network_id))
  cols   <- grDevices::hcl.colors(n_cols, palette = "Spectral", rev = TRUE)
  names(cols) <- unique(net_locs$Network_id)
  
  p_map <- ggplot(europe) +
    geom_sf(fill = "beige", color = "black", linewidth = 0.2) +
    geom_point(data = net_locs,
               aes(x = Longitude, y = Latitude,
                   color = Network_id, size = total_species),
               alpha = 0.8) +
    
    # NETWORK LABELS: using base ggplot2 only
    geom_text(
      data = net_locs,
      aes(x = Longitude, y = Latitude, label = Network_id),
      size = 3.2,
      vjust = -1.3,           # label above point
      hjust = 0.5,
      color = "black",
      fontface = "bold",
      check_overlap = TRUE,   # hides overlapping labels
      show.legend = FALSE
    ) +
    
    scale_color_manual(values = cols, guide = "none") +
    scale_size_continuous(range = c(2, 8), name = "Total Species") +
    coord_sf(xlim = xlim, ylim = ylim,
             expand = FALSE, default_crs = NULL) +
    labs(title = sprintf("Top similar networks in %s (species-based Jaccard)", bio),
         subtitle = sprintf("Points sized by total species; Lat range: %.2f–%.2f",
                            lat_range[1], lat_range[2]),
         x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10))
  
  print(p_map)
  
  ## Save the map
  # ggsave("similar_networks_map.png", p_map, width = 12, height = 8, dpi = 300)
  

  # 10) Heat-map of species similarity – better colours + diagonal
  # Networks that appear on the map
  net_ids_map <- net_locs$Network_id
  n_map       <- length(net_ids_map)
  
  if (n_map < 2) {
    cat("Only one network with coordinates – skipping heatmap.\n")
  } else {
    ## 10-b) Sub-matrix for those networks
    sim_sub <- sim_mat[net_ids_map, net_ids_map]
    
    ## 10-c) Fill the diagonal with 1 (self-similarity)
    diag(sim_sub) <- 1
    
    #Order rows/columns alphabetically (by Network_id)
    ord <- order(net_ids_map) # <--- CORRECT: Sort by the existing object
    sim_sub    <- sim_sub[ord, ord]
    net_labels <- net_ids_map[ord]
    
    # Long format for ggplot2
    heat_df <- expand.grid(Net1 = net_labels,
                           Net2 = net_labels,
                           KEEP.OUT.ATTRS = FALSE,
                           stringsAsFactors = FALSE)
    heat_df$Jaccard <- as.vector(sim_sub)
    
    # Heat-map with a high-contrast palette
    p_heat <- ggplot(heat_df, aes(x = Net1, y = Net2, fill = Jaccard)) +
      geom_tile(color = "white", linewidth = 0.2) +
      
    scale_fill_viridis_c(
      option = "magma",               # dark → bright gradient
      limits = c(0, 1),
      na.value = "grey90",
      name = "Jaccard\nsimilarity"
    ) +
      
      labs(
        title = sprintf("Species similarity (Jaccard) – %s", bio),
        subtitle = sprintf("%d networks – diagonal = 1 (self-similarity)", n_map),
        x = "Network (Net1)",
        y = "Network (Net2)"
      ) +
      
      theme_minimal(base_size = 7) +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(),
        panel.grid = element_blank(),
        legend.position = "right"
      )
    
    ## 10-g) Show & save
    print(p_heat)
    # ggsave(
    #   filename = sprintf("heatmap_species_similarity_%s.png", bio),
    #   plot     = p_heat,
    #   width    = max(8, 0.6 * n_map),
    #   height   = max(6, 0.5 * n_map),
    #   dpi      = 300
    # )
    # cat(sprintf("Heat-map saved as 'heatmap_species_similarity_%s.png'\n", bio))
  }
}
