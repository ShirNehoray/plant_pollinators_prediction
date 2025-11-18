study_summary_unique

selected_study_ids <- unique(study_summary_unique$Study_id)

interactions_filter <- data_new %>%
  filter(Study_id %in% selected_study_ids)


interactions_binary <- interactions_filter %>%
  select(Study_id, Network_id, Plant_accepted_name, Pollinator_accepted_name) %>%
  distinct() %>%
  unite("interaction", Plant_accepted_name, Pollinator_accepted_name, sep = "_") %>%
  mutate(presence = 1) %>%
  pivot_wider(
    names_from = Network_id,
    values_from = presence,
    values_fill = 0
  )





library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

jaccard_results <- interactions_binary %>%
  group_by(Study_id) %>%
  nest() %>%
  mutate(
    jaccard_pairs = map(data, ~ {
      
      net_sets <- .x %>%
        pivot_longer(
          cols = -interaction,
          names_to = "Network",
          values_to = "value"
        ) %>%
        filter(value == 1) %>%                  # keep only true interactions
        group_by(Network) %>%
        summarise(interactions = list(interaction), .groups = "drop")
      
      if (nrow(net_sets) < 2) return(tibble())
      
      # Compute pairwise Jaccard manually
      pairs <- combn(net_sets$Network, 2, simplify = FALSE)
      
      map_dfr(pairs, function(p) {
        set1 <- net_sets$interactions[[which(net_sets$Network == p[1])]]
        set2 <- net_sets$interactions[[which(net_sets$Network == p[2])]]
        
        inter <- length(intersect(set1, set2))
        uni   <- length(union(set1, set2))
        
        tibble(
          Network1 = p[1],
          Network2 = p[2],
          Jaccard  = ifelse(uni == 0, NA, inter / uni)
        )
      })
    })
  ) %>%
  unnest(jaccard_pairs)

print(jaccard_results)

jaccard_clean <- jaccard_results %>% 
  filter(!is.na(Jaccard))

ggplot(jaccard_clean, aes(x = Study_id, y = Jaccard)) +
  geom_boxplot(
    fill = "lightblue3", 
    color = "black", 
    outlier.colour = "darkseagreen", 
    outlier.shape = 16,
    outlier.alpha = 0.7
  ) +
  labs(
    title = "Jaccard index Between Networks Within Each Study",
    x = "Study ID",
    y = "Jaccard Index"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "gray50"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )+
  tme





