#Network-level indices

library(bipartite)

net_metrics <- function(m){
  mb <- (m > 0) * 1
  data.frame(
    Plants      = nrow(m),
    Pollinators = ncol(m),
    Links       = sum(mb),
    Connectance = networklevel(mb, index = "connectance"),
    Nestedness  = networklevel(mb, index = "nestedness"),
    Modularity  = networklevel(mb, index = "modularity"),
    Evenness    = networklevel(mb, index = "interaction evenness")
  )
}

metrics_by_site <- lapply(site_mats_freq, net_metrics) %>%
  bind_rows(.id = "Site")

library(tidyr)

metrics_long <- metrics_by_site %>%
  pivot_longer(cols = c(Connectance, Nestedness, Modularity, Evenness),
               names_to = "Metric", values_to = "Value")

ggplot(metrics_long, aes(x = Site, y = Value, fill = Metric)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = round(Value, 2)),
            position = position_dodge(width = 0.9), 
            vjust = -0.3, size = 3) +
  labs(title = "Network metrics per site", x = "Site", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
