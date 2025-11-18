library(dplyr)
library(readr)
library(ggplot2)
data <- read.csv("Interaction_data.csv", encoding = "latin1", stringsAsFactors = FALSE)
metadata <- read_csv("metadata.csv", show_col_types = FALSE)
tme <-  theme(axis.text = element_text(size = 10, color = "black"),
              axis.title = element_text(size = 16, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
theme_set(theme_bw())


data <- data %>%
  mutate(
    Plant_accepted_name      = gsub(" ", "_", Plant_accepted_name),
    Pollinator_accepted_name = gsub(" ", "_", Pollinator_accepted_name)
  )

network_lookup <- data %>%
  distinct(Study_id, Network_id) %>%
  group_by(Study_id) %>%
  arrange(Network_id, .by_group = TRUE) %>%
  mutate(
    Study_id_clean = gsub("[^A-Za-z0-9]", "", Study_id),
    Network_id_new = paste0(Study_id_clean, "_Network", row_number())
  ) %>%
  ungroup() %>%
  select(Study_id, Network_id, Network_id_new)


data_new <- data%>%
  left_join(network_lookup, by = c("Study_id", "Network_id"))

unique(data_new$Network_id)



data_new <- data %>%
  # For each Study_id, assign a running number to each unique Network_id
  group_by(Study_id) %>%
  mutate(
    Network_number = dense_rank(Network_id),   # assigns 1,2,3,... per Study_id
    Network_id_new = paste0(Study_id, "_Network", Network_number)
  ) %>%
  ungroup()

data_new$Network_id <- data_new$Network_id_new
data_new$Network_id_new <- NULL


study_summary <- data_new %>%   # or data_new if that's your main object
  group_by(Study_id) %>%
  summarise(
    n_networks      = n_distinct(Network_id),
    n_plant         = n_distinct(Plant_accepted_name),
    n_pollinators   = n_distinct(Pollinator_accepted_name),
    n_total_species = n_plant + n_pollinators,
    .groups = "drop"
  )




study_summary_full <- study_summary %>%
  left_join(
    metadata %>%
      select( Study_id, Sampling_method,`Study_classification `),
    by = "Study_id"
  )

metadata$`Study_classification `

filtered_studies <- study_summary_full %>%
  filter(
    n_networks > 10,
    n_total_species > 30,
    `Study_classification ` %in% c("spatial", "spatial and temporal")
  )

filtered_studies <- study_summary_full %>%
  filter(
    n_networks > 10,
    n_total_species > 20,
    `Study_classification ` %in% c("spatial", "spatial and temporal")
  ) %>%
  distinct(Study_id, .keep_all = TRUE)   # ensure no duplicates

study_summary_unique <- filtered_studies %>%
  distinct(Study_id, .keep_all = TRUE)


# write.csv(
#   study_summary_unique,
#   "filtered_studies.csv",
#   row.names = FALSE
# )
# 
ggplot(study_summary_unique,
       aes(x = Study_id, y = n_networks, fill = `Study_classification `)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.4) +   # ← linewidth
  geom_text(aes(label = n_networks),
            vjust = -0.5, size = 3, colour = "black", fontface = "bold") +
  
  scale_fill_manual(
    values = c(
      "spatial"              = "plum4",
      "spatial and temporal" = "lightcoral"
    ),
    name = "Classification"
  ) +
  
  labs(
    title = "Number of Networks per Study (spatial & spatial-and-temporal only)",
    x = "Study_id",
    y = "Number of Networks"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.title    = element_text(face = "bold"),
    axis.text.x   = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
    legend.position = "right",
    legend.title  = element_text(face = "bold"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  tme



interactions_binary <- study_summary_unique %>%
  select(Study_id, Network_id, Plant_accepted_name, Pollinator_accepted_name) %>%
  distinct() %>%
  unite("interaction", Plant_accepted_name, Pollinator_accepted_name, sep = "_") %>%
  mutate(presence = 1) %>%
  pivot_wider(
    names_from = Network_id,
    values_from = presence,
    values_fill = 0
  )
jaccard_results <- interactions_binary %>%
  group_by(Study_id) %>%
  nest() %>%
  mutate(
    jaccard_pairs = map(data, ~ {
      mat <- .x %>%
        select(-interaction) %>%
        as.matrix() %>%
        t()  
      
      net_names <- rownames(mat)
      
      if (nrow(mat) < 2) return(tibble())
      
      dist_mat <- 1-vegdist(mat, method = "jaccard")
      
      as.data.frame(as.matrix(dist_mat)) %>%
        rownames_to_column("Network1") %>%
        pivot_longer(-Network1, names_to = "Network2", values_to = "Jaccard") %>%
        filter(Network1 < Network2) %>%
        mutate(
          Network1 = net_names[as.numeric(Network1)],
          Network2 = net_names[as.numeric(Network2)]
        )
    })
  ) %>%
  select(Study_id, jaccard_pairs) %>%
  unnest(jaccard_pairs)

print(jaccard_results)





library(dplyr)
library(lubridate)

networks_per_year_Mullen <- data %>%
  filter(Study_id == "35_Mullen") %>%        # keep only this study
  mutate(Year = year(ymd(Date))) %>%         # extract year from Date
  distinct(Year, Network_id) %>%             # unique networks per year
  group_by(Year) %>%
  summarise(n_networks = n()) %>%
  arrange(Year)

networks_per_year_Mullen
