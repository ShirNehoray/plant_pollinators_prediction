library(dplyr)
data <- read.csv("Interaction_data.csv", encoding = "latin1", stringsAsFactors = FALSE)
metadata <- read_csv("metadata.csv", show_col_types = FALSE)


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


write.csv(
  study_summary_unique,
  "filtered_studies.csv",
  row.names = FALSE
)



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
