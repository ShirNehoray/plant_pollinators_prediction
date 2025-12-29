library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

df <- read_csv("df_eval_links1_Bartomeus.csv")

# Choose which "layer" you want to summarize by:
# layer_col <- "train_layer"
layer_col <- "test_layer"

# 1) Count evaluations per (network, layer, itr)
counts <- df %>%
  mutate(layer = .data[[layer_col]]) %>%
  filter(!is.na(evaluation), !is.na(layer), !is.na(Study_id), !is.na(itr)) %>%
  count(Study_id, layer, itr, evaluation, name = "n_eval")

# 2) Convert to portions within each (network, layer, itr)
props <- counts %>%
  group_by(Study_id, layer, itr) %>%
  mutate(total = sum(n_eval), portion = n_eval / total) %>%
  ungroup() %>%
  mutate(
    evaluation = factor(evaluation, levels = c("TP","TN","FP","FN")),
    layer = factor(layer)
  )

# 3) Boxplot: distribution across iterations, for each layer, within each network
ggplot(props, aes(x = layer, y = portion, fill = evaluation)) +
  geom_boxplot(outlier_alpha = 0.25, position = position_dodge(width = 0.75)) +
  facet_wrap(~ Study_id, scales = "free_x") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = paste0("Portion of link evaluations by ", layer_col, " (per network)"),
    x = layer_col,
    y = "Portion of interactions",
    fill = "Evaluation"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))






