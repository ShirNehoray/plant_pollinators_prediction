library(dplyr)

Interacion_data_filtered <- read.csv("./csv/Interactions_data_filtered_studies.csv", encoding = "latin1", stringsAsFactors = FALSE)



A_l <- Interacion_data_filtered %>% 
  transmute(
    Study_id    = Study_id,
    layer_from = Network_id,
    node_from    = Plant_accepted_name,
    layer_to   = Network_id,
    node_to      = Pollinator_accepted_name,
    type         = "pollination",
    weight       = Interaction
  )


A_l_new <- A_l %>%
  group_by(Study_id,layer_from, node_from,layer_to, node_to, type) %>%
  summarise(weight = sum(weight), .groups = "drop")

write.csv(A_l_new, file="Interaction_edges.csv", row.names = FALSE)
