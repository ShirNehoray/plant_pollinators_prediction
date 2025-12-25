# This  script is for generating loops for running different combinations of layer prediction on HPC
library(tidyverse)

filtered_studies <- read.csv("./csv/summary_filtered_studies.csv", encoding = "latin1", stringsAsFactors = FALSE)




# Extract unique study IDs
Study_ids <- unique(filtered_studies$Study_id)

# Output bash file
output_path <- "HPC_bash_jobs_all_selected_studies.sh"

# Create bash commands
bash_lines <- c(
  "#!/bin/bash",
  paste("qsub qsub bash_studies_networks_prediction.sh", shQuote(Study_ids))
)

# Write to file
writeLines(bash_lines, output_path)
