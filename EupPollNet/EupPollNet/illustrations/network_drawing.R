library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(scales)


my_data <- read.csv("Interaction_data.csv")
  
study_summary <- my_data %>%
      group_by(Study_id) %>%
      summarize(
        Latitude = first(Latitude), 
        Longitude = first(Longitude), 
        Country = first(Country),
        network_count = n_distinct(Network_id) 
      )
    
    # 3. GET THE BACKGROUND MAP DATA
    world_map <- ne_countries(scale = "medium", returnclass = "sf")
    
    # ---------------------------------------------------------------------
    # 4. CUSTOM LABEL FUNCTIONS
    # ---------------------------------------------------------------------
    lat_labels <- function(breaks) {
      paste0(abs(breaks), "°N")
    }
    lon_labels <- function(breaks) {
      paste0(abs(breaks), "°", ifelse(breaks < 0, "W", ifelse(breaks > 0, "E", "")))
    }
    
    
    # ---------------------------------------------------------------------
    # 5. CREATE THE PLOT (One Size, One Color)
    # ---------------------------------------------------------------------
    
    ggplot() +
      
      # --- Layer 1: The Map ---
      geom_sf(data = world_map, fill = "ivory", color = "black", linewidth = 0.3) +
      
      # --- Layer 2: Your Study Points ---
      # *** KEY CHANGE IS HERE ***
      # 'size' and 'color' are MOVED OUTSIDE of aes()
      # This makes them fixed values, not mapped to data.
      geom_point(
        data = study_summary, 
        aes(
          x = Longitude,
          y = Latitude
          # No size or color mapping inside aes()
        ),
        shape = 21,
        color = "black",
        size = 8,         # Fixed size
        fill = "indianred3",   # Fixed color
        alpha = 0.7 ,
        stroke= 1
      ) +
      
      # --- Layer 3: Coordinate System & Zoom ---
      coord_sf(xlim = c(-12, 32), ylim = c(34, 71), expand = FALSE) +
      
      # --- Layer 4: Axis Labels and Gridlines ---
      # We removed scale_size_continuous() as it is no longer needed
      scale_x_continuous(name = "Longitude", breaks = seq(-10, 30, by = 10), labels = lon_labels) +
      scale_y_continuous(name = "Latitude", breaks = seq(35, 70, by = 5), labels = lat_labels) +
      
      # --- Layer 5: Titles and Theme ---
      labs(
        title = "Geographic Distribution of Studies",
        subtitle = "All study locations"
        # No legend titles ('color', 'size') are needed
      ) +
      
      # --- Layer 6: Theme ---
      theme_bw() +
      theme(
        panel.background = element_rect(fill = "lightcyan3"), 
        # panel.grid.major = element_line(color = "gray80", linetype = "dashed", linewidth = 0.3),
        panel.grid.major = element_blank(), # <-- Hides major gridlines
        panel.grid.minor = element_blank()  # <-- Hides minor gridlines
      )
    
    
    # ----------- For one study -------------------------
    world_map <- ne_countries(scale = "medium", returnclass = "sf")

    lat_labels <- function(breaks) {
      paste0(abs(breaks), "°N")
    }
    lon_labels <- function(breaks) {
      paste0(abs(breaks), "°", ifelse(breaks < 0, "W", ifelse(breaks > 0, "E", "")))
    }
    
    # ---------------------------------------------------------------------
    # 3. FILTER & AGGREGATE DATA (using logic from your example)
    # ---------------------------------------------------------------------
    target_study_id <- "36_Larkin" # <-- THIS IS THE ONLY CHANGE
    
    # First, filter for the study, then aggregate by Network_id
    # to get one location per network. This creates our 'net_locs' table.
    networks_to_plot <- my_data %>%
      filter(Study_id == target_study_id) %>%
      group_by(Network_id) %>%
      summarize(
        Latitude = first(Latitude), # Get the first Lat for this network
        Longitude = first(Longitude), # Get the first Lon for this network
        .groups = 'drop' # Drop the grouping
      )
    
    # Check that we found data
    if(nrow(networks_to_plot) == 0) {
      stop(paste("Error: Study ID '", target_study_id, "' not found, or has no networks.", sep=""))
    }
    
    # ---------------------------------------------------------------------
    # 4. DETERMINE MAP EXTENTS (using logic from your example)
    # ---------------------------------------------------------------------
    # This calculates the zoom boundaries automatically based on the data
    lat_range <- range(networks_to_plot$Latitude,  na.rm = TRUE)
    lon_range <- range(networks_to_plot$Longitude, na.rm = TRUE)
    
    # Set the buffer (padding) around the points
    lat_buffer <- 1 # 1 degree N/S buffer (you can change this)
    lon_buffer <- 2 # 2 degrees E/W buffer (you can change this)
    
    # Define the final zoom limits
    xlim <- c(lon_range[1] - lon_buffer, lon_range[2] + lon_buffer)
    ylim <- c(lat_range[1] - lat_buffer, lat_range[2] + lat_buffer)
    
    
    # ---------------------------------------------------------------------
    # 5. GET BASE MAP (Europe, like your example)
    # ---------------------------------------------------------------------
    europe <- rnaturalearth::ne_countries(
      scale = "medium", 
      returnclass = "sf", 
      continent = "Europe" # Filter map data to Europe
    )
    
    # ---------------------------------------------------------------------
    # 6. CREATE THE PLOT (Zoomed, Blue dots, No grid)
    # ---------------------------------------------------------------------
    
    ggplot() + # We start with ggplot(), then add layers
      
      # --- Layer 1: The Map ---
      # Use 'europe' data for the map layer
      geom_sf(data = europe, fill = "ivory", color = "black", linewidth = 0.3) +
      
      # --- Layer 2: Your Study Points ---
      geom_point(
        data = networks_to_plot, # Use our new aggregated table
        aes(
          x = Longitude,
          y = Latitude
        ),
        shape = 21,
        color = "black",
        size = 8,         # Fixed size
        fill = "aquamarine4",   # Fixed color
        alpha = 0.7 ,
        stroke= 1
      ) +
      
      # --- Layer 3: Coordinate System & Zoom ---
      # Use the xlim/ylim we calculated automatically
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      
      # --- Layer 4: Axis Labels ---
      scale_x_continuous(name = "Longitude", labels = lon_labels) +
      scale_y_continuous(name = "Latitude", labels = lat_labels) +
      
      # --- Layer 5: Titles ---
      labs(
        title = "Geographic Distribution of Networks",
        subtitle = paste("Networks for study:", target_study_id) # Subtitle updates automatically
      ) +
      
      # --- Layer 6: Theme (No Gridlines) ---
      theme_bw() + # Use theme_bw for a clean background
      theme(
        panel.background = element_rect(fill = "lightsteelblue1"), # Ocean color
        panel.grid.major = element_blank(), # Hide gridlines
        panel.grid.minor = element_blank()  # Hide gridlines
      )
    
    