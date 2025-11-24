#graphing chosen networks 

summary_data <- read.csv("summary_filtered_studies.csv", encoding = "latin1", stringsAsFactors = FALSE)
tme <-  theme(axis.text = element_text(size = 10, color = "black"),
              axis.title = element_text(size = 16, face = "bold"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
theme_set(theme_bw())



ggplot(summary_data,
       aes(x = Study_id, y = n_networks, fill = Study_classification)) +
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

