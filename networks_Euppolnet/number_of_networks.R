library(dplyr)
library(ggplot2)

# Read the interaction data
data <- read.csv("Interaction_data.csv", encoding = "latin1", stringsAsFactors = FALSE)

# Count unique networks per study 
network_counts <- data %>%
  # Keep only unique study_id + network_id combinations
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
  theme(plot.title = element_text(hjust = 0.5),  # Center the title
# --- REMOVE ALL GRID LINES ---
panel.grid.major = element_blank(),
panel.grid.minor = element_blank()
)


# remove 13_Karise only for plotting
plot_df <- network_counts %>% filter(Study_id != "13_Karise")


ggplot(plot_df, aes(x = reorder(Study_id, -unique_networks), y = unique_networks)) +
  geom_bar(stat = "identity", fill = "#358") +
  theme_minimal() +
  labs(title = "Number of Unique Networks per Study (excluding 13_Karise)",
       x = "Study ID",
       y = "Number of Unique Networks") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(hjust = 0.5),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank()
)
