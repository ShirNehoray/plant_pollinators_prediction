# ============================
# EDA Pipeline for Top Jaccard Pairs in Continental Bioregion - Network Comparisons
# ============================
# This pipeline:
# 1. Filters to Continental bioregion.
# 2. Computes Jaccard similarities as in your original code.
# 3. Selects the top N pairs with highest Jaccard (default top_n_pairs = 10).
# 4. Performs Exploratory Data Analysis (EDA) with focus on comparisons between networks:
#    - Summary statistics (species richness, interaction counts, shared species).
#    - Bar graphs for plant and pollinator species frequency across networks (Y: number of networks, X: species, sorted by frequency).
#    - Visualizations (histograms of species richness, shared species bar plots, interaction heatmaps for top pair).
#    - Network metrics (connectance, modularity if using bipartite package).
#    - Spatial analysis (map of top networks).
# 5. Outputs results to console and plots.

# Note: Assumes Interaction_data.csv has columns: Bioregion, Network_id, Plant_accepted_name, Pollinator_accepted_name, Interaction (if available), Latitude, Longitude.
# Install additional packages if needed: install.packages(c("bipartite", "igraph")) for network metrics.

library(dplyr)
library(readr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(tidyr)
library(purrr)
library(bipartite)  # For network metrics like connectance
library(igraph)    # For additional graph analysis

setwd("/Users/shirn/OneDrive/Documents/master/data")

# Parameters (easy to tweak)
bio <- "Continental"  # Focus on this bioregion
min_total <- 15       # Minimum total species per network
top_n_pairs <- 50     # Focus on top 10 pairs with highest Jaccard (adjust as needed)
top_species_n <- 100   # Number of top species to show in bar plots for plants/pollinators

# ============================
# Read data
# ============================
interactions <- read.csv("Interaction_data.csv", stringsAsFactors = FALSE)

# ============================
# Utilities
# ============================
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


# ============================
# Basic column checks
# ============================
req_cols <- c("Bioregion", "Network_id", "Plant_accepted_name", "Pollinator_accepted_name")
missing_cols <- setdiff(req_cols, names(interactions))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
}

# ============================
# 1) Filter to Continental bioregion
# ============================
dat_bio <- interactions[interactions$Bioregion == bio, , drop = FALSE]

if (nrow(dat_bio) == 0) {
  stop(paste("No rows for Bioregion:", bio))
}

# ============================
# 2) Compute species counts per network
# ============================
plants_by_net <- aggregate(Plant_accepted_name ~ Network_id, data = dat_bio,
                           FUN = function(x) length(unique(x)))
names(plants_by_net)[2] <- "n_plants"

polls_by_net <- aggregate(Pollinator_accepted_name ~ Network_id, data = dat_bio,
                          FUN = function(x) length(unique(x)))
names(polls_by_net)[2] <- "n_pollinators"

counts <- merge(plants_by_net, polls_by_net, by = "Network_id", all = TRUE)
counts$n_plants[is.na(counts$n_plants)] <- 0
counts$n_pollinators[is.na(counts$n_pollinators)] <- 0
counts$total_species <- counts$n_plants + counts$n_pollinators

nets_keep <- counts$Network_id[counts$total_species >= min_total]
if (length(nets_keep) == 0) {
  stop(paste("No networks with total_species >=", min_total, "in", bio))
}

# ============================
# 3) Build species sets per network
# ============================
dat_bio_sub <- dat_bio[dat_bio$Network_id %in% nets_keep, , drop = FALSE]

plant_names <- normalize_name(dat_bio_sub$Plant_accepted_name)
poll_names  <- normalize_name(dat_bio_sub$Pollinator_accepted_name)

plant_tags <- paste0("PLANT::", plant_names)
poll_tags  <- paste0("POLL::", poll_names)

idx_by_net <- split(seq_len(nrow(dat_bio_sub)), dat_bio_sub$Network_id)

species_items <- lapply(idx_by_net, function(idx) {
  unique(c(plant_tags[idx], poll_tags[idx]))
})

net_ids <- names(species_items)
n <- length(net_ids)

# ============================
# 4) Pairwise Jaccard matrix
# ============================
sim_mat <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(net_ids, net_ids))

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

# ============================
# 5) Tidy table of unique pairs, sorted high-to-low
# ============================
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

if (is.null(sim_df) || nrow(sim_df) == 0 || !("Jaccard" %in% names(sim_df)) || all(is.na(sim_df$Jaccard))) {
  stop(sprintf("No valid Jaccard similarities for %s (empty or all NA).", bio))
}

sim_df <- sim_df[!is.na(sim_df$Jaccard), , drop = FALSE]
ord <- order(sim_df$Jaccard, decreasing = TRUE)
sim_df <- sim_df[ord, , drop = FALSE]

top_n <- min(top_n_pairs, nrow(sim_df))
sim_top <- sim_df[seq_len(top_n), , drop = FALSE]

# ============================
# 6) Inspect top pairs
# ============================
cat(sprintf("\n=== Top %d similar network pairs in %s ===\n", top_n, bio))
print(sim_top, row.names = FALSE)

# Collect unique networks from top pairs
unique_nets <- unique(c(sim_top$Net1, sim_top$Net2))

# ============================
# EDA Pipeline on Top Networks - Focus on Comparisons
# ============================
# Subset data to top networks
top_data <- dat_bio_sub[dat_bio_sub$Network_id %in% unique_nets, , drop = FALSE]

# 6a) Summary Statistics
# - Species richness per network
species_summary <- counts[counts$Network_id %in% unique_nets, ]
cat("\n=== Species Richness Summary for Top Networks ===\n")
print(species_summary)

# - Shared species for top pairs
shared_species <- data.frame()
for (i in seq_len(nrow(sim_top))) {
  net1 <- sim_top$Net1[i]
  net2 <- sim_top$Net2[i]
  sp1 <- species_items[[net1]]
  sp2 <- species_items[[net2]]
  shared <- length(intersect(sp1, sp2))
  shared_species <- rbind(shared_species, data.frame(Net1 = net1, Net2 = net2, Shared = shared, Jaccard = sim_top$Jaccard[i]))
}
cat("\n=== Shared Species for Top Pairs ===\n")
print(shared_species)

# - Interaction counts (assuming 'Interaction' column exists; adjust if named differently)
if ("Interaction" %in% names(top_data)) {
  interaction_summary <- aggregate(Interaction ~ Network_id, data = top_data, FUN = sum)
  names(interaction_summary)[2] <- "Total_Interactions"
  cat("\n=== Interaction Counts for Top Networks ===\n")
  print(interaction_summary)
} else {
  cat("\nNote: No 'Interaction' column found; skipping interaction counts.\n")
}

# 6b) Species Frequency Across Networks
# - Plants: Count number of networks each plant appears in
plant_freq <- top_data %>%
  group_by(Plant_accepted_name) %>%
  summarise(Num_Networks = n_distinct(Network_id), .groups = "drop") %>%
  arrange(desc(Num_Networks)) %>%
  slice_head(n = top_species_n)  # Top N plants

cat("\n=== Top Plants by Frequency Across Networks ===\n")
print(plant_freq)

# Bar plot for plants
p_plant_freq <- ggplot(plant_freq, aes(x = reorder(Plant_accepted_name, Num_Networks), y = Num_Networks)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  coord_flip() +  # Flip for better readability
  labs(title = "Top Shared Plants ", 
       x = "Plant Species", 
       y = "Number of Networks") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 5))  

print(p_plant_freq)


# - Pollinators: Count number of networks each pollinator appears in
poll_freq <- top_data %>%
  group_by(Pollinator_accepted_name) %>%
  summarise(Num_Networks = n_distinct(Network_id), .groups = "drop") %>%
  arrange(desc(Num_Networks)) %>%
  slice_head(n = top_species_n)  # Top N pollinators

cat("\n=== Top Pollinators by Frequency Across Networks ===\n")
print(poll_freq)

# Bar plot for pollinators
p_poll_freq <- ggplot(poll_freq, aes(x = reorder(Pollinator_accepted_name, Num_Networks), y = Num_Networks)) +
  geom_bar(stat = "identity", fill = "purple") +
  coord_flip() +  # Flip for better readability
  labs(title = "Top Shared Pollinators ", x = "Pollinator Species", y = "Number of Networks") +
  theme(axis.text.y = element_text(size = 5))
print(p_poll_freq)





# 6d) Network Metrics (using bipartite package) (? check if it is wrtitten good!)

# ================================================================
# NETWORK METRICS PER NETWORK: Connectance, Nestedness (NODF),
# Modularity (Q), and H2' Specialization
# ================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(bipartite)   # networklevel(), computeModules(), modularity()
# install.packages("GGally") # if needed for correlation matrix plot
# library(GGally)

# ----------------------------
# 0) Helpers (robust wrappers)
# ----------------------------

# Build a plant × pollinator matrix (counts). Returns a matrix with dimnames.
build_web <- function(df, plant_col, poll_col) {
  # xtabs on character columns; works as counts and supports binary presence
  xtabs(reformulate(c(plant_col, poll_col)), data = df)
  # ================================================================
  # NETWORK METRICS PER NETWORK — SEPARATE PLOTS FOR EACH METRIC
  # ================================================================
  
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(bipartite)
  
  # ----------------------------
  # Helper functions (safe metric calculations)
  # ----------------------------
  calc_connectance_safe <- function(web) {
    if (is.null(dim(web)) || any(dim(web) == 0)) return(NA_real_)
    as.numeric(bipartite::networklevel(web, index = "connectance"))
  }
  
  calc_nestedness_safe <- function(web) {
    if (is.null(dim(web)) || any(dim(web) == 0)) return(NA_real_)
    tryCatch(as.numeric(bipartite::networklevel(web, index = "nestedness")),
             error = function(e) NA_real_)
  }
  
  calc_H2_safe <- function(web) {
    if (is.null(dim(web)) || any(dim(web) == 0)) return(NA_real_)
    tryCatch(as.numeric(bipartite::networklevel(web, index = "H2")),
             error = function(e) NA_real_)
  }
  
  calc_modularity_safe <- function(web) {
    if (is.null(dim(web)) || any(dim(web) == 0)) return(NA_real_)
    tryCatch({
      M <- bipartite::computeModules(web)
      as.numeric(bipartite::modularity(M))
    }, error = function(e) NA_real_)
  }
  
  # ----------------------------
  # Select columns (normalized if available)
  # ----------------------------
  use_norm <- all(c("Plant_norm", "Poll_norm") %in% names(top_data))
  plant_col <- if (use_norm) "Plant_norm" else "Plant_accepted_name"
  poll_col  <- if (use_norm) "Poll_norm"  else "Pollinator_accepted_name"
  
  # ----------------------------
  # Compute metrics for each network
  # ----------------------------
  network_metrics <- top_data %>%
    group_by(Network_id) %>%
    group_modify(~{
      web <- xtabs(reformulate(c(plant_col, poll_col)), data = .x)
      P <- nrow(web)
      A <- ncol(web)
      L <- sum(web > 0)
      tibble(
        Plants      = P,
        Pollinators = A,
        Links       = L,
        Possible    = P * A,
        Connectance = calc_connectance_safe(web),
        Nestedness  = calc_nestedness_safe(web),
        ModularityQ = calc_modularity_safe(web),
        H2_prime    = calc_H2_safe(web)
      )
    }) %>%
    ungroup()
  
  cat("\n=== Network Metrics Summary ===\n")
  print(network_metrics)
  
  # ================================================================
  # INDIVIDUAL PLOTS FOR EACH METRIC
  # ================================================================
  
  # ================================================================
  # Separate plots per metric: BOX + HISTOGRAM (with robust binwidth)
  # Assumes you already have `network_metrics` with columns:
  #   Connectance, Nestedness, ModularityQ, H2_prime
  # ================================================================
  
  library(dplyr)
  library(ggplot2)
  
  plot_metric_hist <- function(df, metric_name, color = "orange") {
    # Extract and validate metric column
    df_metric <- df %>%
      dplyr::select(Network_id, dplyr::all_of(metric_name)) %>%
      dplyr::rename(Value = dplyr::all_of(metric_name))
    
    if (!"Value" %in% names(df_metric) || all(is.na(df_metric$Value))) {
      message(sprintf("No valid (non-NA) values for metric '%s' — skipping plots.", metric_name))
      return(invisible(NULL))
    }
    
    df_metric <- dplyr::filter(df_metric, !is.na(Value))
    if (nrow(df_metric) == 0) {
      message(sprintf("No rows to plot for metric '%s' after filtering NA — skipping.", metric_name))
      return(invisible(NULL))
    }
    
    # Compute robust binwidth (Freedman–Diaconis); fallback if degenerate
    vals <- df_metric$Value
    bw_fd <- 2 * stats::IQR(vals) / (length(vals)^(1/3))
    rng   <- diff(range(vals))
    if (!is.finite(bw_fd) || bw_fd <= 0) {
      # fallback: 10% of range or a small positive default
      bw_fd <- max(rng * 0.1, 0.01)
    }
    # If all values are identical, force a single bin around that value
    single_value <- length(unique(vals)) == 1
    med <- stats::median(vals)
    
    # ---------------------------
    # 1) Boxplot (distribution)
    # ---------------------------
    p_box <- ggplot2::ggplot(df_metric, ggplot2::aes(x = "", y = Value)) +
      ggplot2::geom_boxplot(fill = color, alpha = 0.6, outlier.color = "red") +
      ggplot2::geom_jitter(width = 0.1, alpha = 0.7, size = 2) +
      ggplot2::labs(
        title = paste(metric_name, "Distribution Across Networks"),
        y = metric_name, x = NULL
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold")
      )
    print(p_box)
    
    # ---------------------------
    # 2) Histogram (values distribution)
    # ---------------------------
    # If only one unique value, show a single-bar histogram centered on that value.
    if (single_value) {
      # Create a tiny jitter to make ggplot draw a bar
      p_hist <- ggplot2::ggplot(df_metric, ggplot2::aes(x = Value)) +
        ggplot2::geom_histogram(bins = 1, fill = color, color = "white", alpha = 0.8) +
        ggplot2::geom_vline(xintercept = med, linetype = "dashed") +
        ggplot2::labs(
          title = paste(metric_name, "Histogram (Single Value)"),
          x = metric_name, y = "Count"
        ) +
        ggplot2::theme_minimal()
    } else {
      p_hist <- ggplot2::ggplot(df_metric, ggplot2::aes(x = Value)) +
        ggplot2::geom_histogram(binwidth = bw_fd, boundary = min(vals),
                                closed = "left", fill = color,
                                color = "white", alpha = 0.8) +
        ggplot2::geom_vline(xintercept = med, linetype = "dashed") +
        ggplot2::labs(
          title = paste(metric_name, "Histogram"),
          x = metric_name, y = "Count"
        ) +
        ggplot2::theme_minimal()
    }
    print(p_hist)
  }
  
  # ================================================================
  # Create separate plots for each key metric
  # (Adjust colors if you like)
  # ================================================================
  plot_metric_hist(network_metrics, "Connectance", color = "#E69F00")
  plot_metric_hist(network_metrics, "Nestedness",  color = "#56B4E9")
  # plot_metric_hist(network_metrics, "ModularityQ", color = "#009E73")
  plot_metric_hist(network_metrics, "H2_prime",    color = "#CC79A7")
  

#-----
# 6e) Spatial Analysis (Map of top networks)
net_locs <- aggregate(cbind(Latitude, Longitude) ~ Network_id, data = top_data, FUN = mean, na.rm = TRUE)
net_locs <- merge(net_locs, species_summary[, c("Network_id", "total_species")], by = "Network_id")

europe <- ne_countries(scale = "medium", returnclass = "sf", continent = "Europe")

lat_range <- range(net_locs$Latitude, na.rm = TRUE)
lon_range <- range(net_locs$Longitude, na.rm = TRUE)
lat_buffer <- 1; lon_buffer <- 2
xlim <- c(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer)
ylim <- c(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer)

cols <- grDevices::hcl.colors(nrow(net_locs), "Spectral", rev = TRUE)

p_map <- ggplot(europe) +
  geom_sf(fill = "beige", color = "black", linewidth = 0.2) +
  geom_point(
    data = net_locs,
    aes(x = Longitude, y = Latitude, color = Network_id, size = total_species),
    alpha = 0.8
  ) +
  scale_color_manual(values = cols, guide = "none") +
  scale_size_continuous(range = c(2, 8), name = "Total Species") +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(
    title = sprintf("Top Similar Networks in %s", bio),
    subtitle = sprintf("Points sized by total species; Lat range: %.2f–%.2f", lat_range[1], lat_range[2]),
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p_map)


# End of EDA Pipeline
