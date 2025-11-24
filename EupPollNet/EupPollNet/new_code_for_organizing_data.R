library(dplyr)
library(readr)
library(ggplot2)

data <- read.csv("Interaction_data.csv", encoding = "latin1", stringsAsFactors = FALSE)
metadata <- read_csv("metadata.csv", show_col_types = FALSE)


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

write.csv(data_new, "Interaction_data_fix.csv", row.names = FALSE)


#summary the data
study_summary <- data_new %>%   
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


write.csv(
  study_summary_unique,
  "summary_filtered_studies.csv",
  row.names = FALSE
)


selected_study_ids <- unique(study_summary_unique$Study_id)

interactions_filter <- data_new %>%
  filter(Study_id %in% selected_study_ids)

write.csv(
  interactions_filter,
  "Interactions_data_filtered_studies.csv",
  row.names = FALSE
)






