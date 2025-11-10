library(dplyr)
library(ggplot2)

# Read the interaction data
# data <- read.csv("Interaction_data.csv", encoding = "latin1", stringsAsFactors = FALSE)
data <- Interaction_data

# Number of networks per each study ----
# Count unique networks per study 
network_counts <- data %>%
  dplyr::distinct(Study_id, Network_id) %>%
  dplyr::group_by(Study_id) %>%
  dplyr::summarise(unique_networks = dplyr::n()) %>%
  dplyr::arrange(desc(unique_networks))


# Create a bar plot
ggplot(network_counts, aes(x = reorder(Study_id, -unique_networks), y = unique_networks)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() +
  labs(title = "Number of Unique Networks per Study",
       x = "Study ID",
       y = "Number of Unique Networks") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5), 
panel.grid.major = element_blank(), # remove grid lines
panel.grid.minor = element_blank()


)


# Remove the large network -  13_Karise
no_large_network <- network_counts %>% filter(Study_id != "13_Karise")

ggplot(no_large_network, aes(x = reorder(Study_id, -unique_networks), y = unique_networks)) +
  geom_bar(stat = "identity", fill = "#335588") +
  theme_minimal() +
  labs(title = "Number of Unique Networks per Study (excluding 13_Karise)",
       x = "Study ID",
       y = "Number of Unique Networks") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(hjust = 0.5),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank()
)


# Create the distribution data
network_distribution <- network_counts %>%
  dplyr::group_by(unique_networks) %>%
  dplyr::summarise(study_count = dplyr::n()) %>%
  dplyr::arrange(unique_networks)


#  Plot the Bar Chart
ggplot(network_distribution, aes(x = study_count, y = factor(unique_networks))) +
  geom_col(fill = "steelblue") + # geom_col is the same as geom_bar(stat="identity")
  theme_minimal() +
  labs(title = "Distribution of Studies by Network Count",
       x = "Number of Studies",           # Swapped label
       y = "Number of Unique Networks") +  # Swapped label
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()  
  )



# Distribution for no_large_network
network_distribution_no_large_network<- no_large_network %>%
  dplyr::group_by(unique_networks) %>%
  dplyr::summarise(study_count = dplyr::n()) %>%
  dplyr::arrange(unique_networks)

# Plot distribution 
ggplot(network_distribution_no_large_network, aes(x = study_count, y = unique_networks)) +
  geom_line(color = "darkgreen", linewidth = 1) + # 'linewidth' is the modern argument for size
  geom_point(color = "darkgreen", size = 3) +
  theme_minimal() +
  labs(title = "Distribution of Number of Networks per Study",
       x = "Number of Studies",
       y = "Number of Unique Networks") +
  theme(
    plot.title = element_text(hjust = 0.5),
    # panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()  
  ) +
  scale_x_continuous(breaks = seq(0, max(network_distribution$unique_networks, na.rm = TRUE), by = 1))


# Size of network -----

#  Calculate pollinator, plant, and total species for each network
network_sizes <- data %>%
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
  )


# 1. Boxplot for PLANT species count per study
ggplot(large_networks, aes(x = reorder(Study_id, -n_plants, FUN = median), y = n_plants)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5, color = "darkgreen") +
  theme_minimal() +
  labs(
    title = "Plant Species per Network (>= 20 total species) by Study",
    x = "Study ID",
    y = "Number of Plant Species"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# 2. Boxplot for POLLINATOR species count per study
ggplot(large_networks, aes(x = reorder(Study_id, -n_pollinators, FUN = median), y = n_pollinators)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, alpha = 0.5, color = "purple") +
  theme_minimal() +
  labs(
    title = "Pollinator Species per Network (>= 20 total species) by Study",
    x = "Study ID",
    y = "Number of Pollinator Species"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
