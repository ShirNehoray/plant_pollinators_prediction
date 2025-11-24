# Remove the large network -  13_Karise
no_large_network <- study_summary %>% filter(Study_id != "13_Karise")

ggplot(no_large_network, aes(x = reorder(Study_id, -unique_networks), y = unique_networks)) +
  geom_bar(stat = "identity", fill = "#335588") +
  theme_minimal() +
  labs(title = "Number of Networks per Study (excluding 13_Karise)",
       x = "Study ID",
       y = "Number of Networks") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  tme


# Create the distribution data
network_distribution <- network_counts %>%
  dplyr::group_by(unique_networks) %>%
  dplyr::summarise(study_count = dplyr::n()) %>%
  dplyr::arrange(unique_networks)


# #  Plot the Bar Chart
# ggplot(network_distribution, aes(x = study_count, y = factor(unique_networks))) +
#   geom_col(fill = "steelblue") + # geom_col is the same as geom_bar(stat="identity")
#   theme_minimal() +
#   labs(title = "Distribution of Studies by Network Count",
#        x = "Number of Studies",           # Swapped label
#        y = "Number Networks") +  # Swapped label
#   theme(plot.title = element_text(hjust = 0.5),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank()  
#   )+
#   tme



# Distribution for no_large_network
network_distribution_no_large_network<- no_large_network %>%
  dplyr::group_by(unique_networks) %>%
  dplyr::summarise(study_count = dplyr::n()) %>%
  dplyr::arrange(unique_networks)

# Plot distribution 
ggplot(network_distribution_no_large_network, aes(x = study_count, y = unique_networks)) +
  geom_line(color = "steelblue", linewidth = 1) + # 'linewidth' is the modern argument for size
  geom_point(color = "steelblue", size = 3) +
  theme_minimal() +
  labs(title = "Distribution of Number of Networks per Study",
       x = "Number of Studies",
       y = "Number of Networks") +
  theme(
    plot.title = element_text(hjust = 0.5),
    # panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()  
  ) +
  scale_x_continuous(breaks = seq(0, max(network_distribution$unique_networks, na.rm = TRUE), by = 1)) +
  tme


networks_per_year_Mullen <- data %>%
  filter(Study_id == "35_Mullen") %>%        # keep only this study
  mutate(Year = year(ymd(Date))) %>%         # extract year from Date
  distinct(Year, Network_id) %>%             # unique networks per year
  group_by(Year) %>%
  summarise(n_networks = n()) %>%
  arrange(Year)

networks_per_year_Mullen