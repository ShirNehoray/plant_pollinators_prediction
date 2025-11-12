# Load required packages
library(readxl)      
library(dplyr)       # for data wrangling (filter, select, group_by, summarise, etc.)
library(igraph)      # for network analysis
library(bipartite)   # specialized tools for bipartite ecological networks
library(ggplot2)     # for visualization
library(tibble)      # for tidy data structures
library(tidyr)       # for data reshaping (pivot_longer, pivot_wider)

# --- Import data from Excel file ---
net_pa <- read_excel("Kaiser_Bunbury_et_al_2017.xlsx", sheet = "64 networks_no.visits")       
net_freq <- read_excel("Kaiser_Bunbury_et_al_2017.xlsx", sheet = "64 networks_visitfreq")    
plants <- read_excel("Kaiser_Bunbury_et_al_2017.xlsx",
                     sheet = "Plant species",
                     skip = 1,        # skip first row with notes
                     col_names = TRUE) # use header names
pollinators <- read_excel("Kaiser_Bunbury_et_al_2017.xlsx",
                          sheet = "Pollinator species",
                          skip = 1,      # skip first row with notes
                          col_names = TRUE)

info <- read_excel("Kaiser_Bunbury_et_al_2017.xlsx", sheet = "info")  

# --- Identify pollinator columns ---
poll_start <- which(names(net_pa) == "Floral abundance") + 1  
# position of the first pollinator column
poll_cols  <- names(net_pa)[poll_start:ncol(net_pa)]          
# names of all pollinator columns

# --- Build quantitative (weighted) network edges ---
edges_freq <- net_freq %>%
  select(Site, `Plant species ID`, all_of(poll_cols)) %>%     # keep only site, plant ID, and pollinators
  pivot_longer(cols = all_of(poll_cols), 
               names_to = "PollinatorID",
               values_to = "weight", 
               values_drop_na = TRUE) %>%                     # convert from wide to long format
  mutate(weight = as.numeric(weight)) %>%                     # ensure numeric weights
  group_by(Site, `Plant species ID`, PollinatorID) %>%        # group by site, plant, and pollinator
  summarise(weight = sum(weight, na.rm = TRUE), .groups = "drop") %>% 
  filter(weight > 0)                                          # keep only interactions with visits

# --- Build binary (presence/absence) network edges ---
edges_bin <- net_pa %>%
  select(Site, `Plant species ID`, all_of(poll_cols)) %>%
  pivot_longer(cols = all_of(poll_cols), 
               names_to = "PollinatorID",
               values_to = "x", 
               values_drop_na = TRUE) %>%
  mutate(x = as.numeric(x)) %>%                               # convert to numeric
  group_by(Site, `Plant species ID`, PollinatorID) %>%
  summarise(weight = as.numeric(any(x > 0, na.rm = TRUE)),     # set weight = 1 if any presence
            .groups = "drop")  

# --- Convert edge lists into matrices per site ---
names(edges_freq) <- trimws(names(edges_freq))   # remove extra spaces in column names

# Function to convert edge list of a site into adjacency matrix
edge_to_matrix <- function(df_site){
  m <- xtabs(weight ~ `Plant species ID` + PollinatorID, data = df_site)  
  as.matrix(m)   # return as standard matrix
}

# Apply function to all sites
site_mats_freq <- split(edges_freq, edges_freq$Site) |>   # split data by site
  lapply(edge_to_matrix)                                  # convert each to matrix

# --- Check results ---
length(site_mats_freq)       # number of networks (sites) -> should be 8
names(site_mats_freq)        # names of the sites
dim(site_mats_freq[[1]])     # dimensions of the first site's matrix (plants × pollinators)

# --- Alternative method: use pivot_wider to build matrices ---
edge_to_matrix_wider <- function(df_site){
  df_site |>
    group_by(`Plant species ID`, PollinatorID) |>
    summarise(weight = sum(weight, na.rm = TRUE), .groups = "drop") |>
    pivot_wider(names_from = PollinatorID, 
                values_from = weight,
                values_fill = 0,                           # fill missing interactions with 0
                values_fn = list(weight = sum)) |>
    as.data.frame() |>
    (\(x) { rownames(x) <- make.unique(x$`Plant species ID`);  
    as.matrix(x[, setdiff(names(x), "Plant species ID"), drop = FALSE]) })()
  # assign plant IDs as rownames, remove column "Plant species ID"
}

site_mats_freq2 <- split(edges_freq, edges_freq$Site) |> 
  lapply(edge_to_matrix_wider)

# --- Another way to split data per site ---
g <- edges_freq |> dplyr::group_by(Site)   
site_list  <- g |> group_split()                           # list of edge lists per site
site_keys  <- g |> group_keys() |> dplyr::pull(Site)       # get site names
names(site_list) <- site_keys                              # name each list element
site_mats_freq <- lapply(site_list, edge_to_matrix)        # convert each to matrix

# --- Summary information ---
n_sites <- dplyr::n_distinct(edges_freq$Site)   # number of unique sites
n_sites                           

length(site_mats_freq)             # number of networks built

# Get number of plants and pollinators in each network
sizes <- sapply(site_mats_freq, dim)            # dimensions of each matrix
rownames(sizes) <- c("Plants", "Pollinators")   # label rows
sizes



# Creating the 8 networks (network per site) ---- 
#quantitive 
site_mats_freq <- split(edges_freq, edges_freq$Site) |>
  lapply(edge_to_matrix)

# Binary
site_mats_bin  <- split(edges_bin,  edges_bin$Site) |>
  lapply(edge_to_matrix)


# Save as csv each network matrix
dir.create("aggregated_by_site", showWarnings = FALSE)

invisible(lapply(names(site_mats_freq), function(s){
  write.csv(site_mats_freq[[s]],
            file = file.path("aggregated_by_site", paste0("network_freq_", s, ".csv")))
}))

invisible(lapply(names(site_mats_bin), function(s){
  write.csv(site_mats_bin[[s]],
            file = file.path("aggregated_by_site", paste0("network_bin_", s, ".csv")))
}))



# --- Graphing bipartite networks ---

# --- Safe CSV reader: convert each CSV to a numeric matrix ---
read_net_csv <- function(path){
  # Assumes the first column contains Plant species IDs
  x <- read.csv(path, check.names = FALSE, row.names = 1)  # keep original colnames, use first col as rownames
  m <- as.matrix(x)                                        # convert to matrix
  storage.mode(m) <- "numeric"                             # ensure numeric entries
  m[is.na(m)] <- 0                                         # replace missing values with 0
  m[m < 0] <- 0                                            # prevent negative values
  return(m)
}

# --- Collect all frequency (weighted) network files ---
freq_files <- list.files("aggregated_by_site", 
                         pattern = "^network_freq_.*\\.csv$", 
                         full.names = TRUE)                # find all matching CSVs
site_names <- sub("^network_freq_(.*)\\.csv$", "\\1", basename(freq_files))  
# extract site names from file names
site_mats_freq <- setNames(lapply(freq_files, read_net_csv), site_names)  
# read each file → matrix, store in list with site names

# --- Collect all binary (presence/absence) network files, if they exist ---
bin_files  <- list.files("aggregated_by_site", 
                         pattern = "^network_bin_.*\\.csv$", 
                         full.names = TRUE)
if(length(bin_files) > 0){
  site_names_bin <- sub("^network_bin_(.*)\\.csv$", "\\1", basename(bin_files))
  site_mats_bin  <- setNames(lapply(bin_files, read_net_csv), site_names_bin)
}

# --- Plot weighted networks into a PDF ---
pdf("plotweb_8_sites_freq_from_csv.pdf", width = 10, height = 12)  
par(mfrow = c(4,2), mar = c(2,2,2,2))   # arrange 8 plots (4 rows × 2 cols), small margins

for(s in names(site_mats_freq)){
  plotweb(site_mats_freq[[s]], method = "normal", text.rot = 90)  
  title(main = s)    # add site name as plot title
}

dev.off() 

# --- Plot binary networks into a PDF (if available) ---
if(exists("site_mats_bin")){
  pdf("plotweb_8_sites_bin_from_csv.pdf", width = 10, height = 12)
  par(mfrow = c(4,2), mar = c(2,2,2,2))   # 8 plots layout
  
  for(s in names(site_mats_bin)){
    plotweb(site_mats_bin[[s]], method = "normal", text.rot = 90)  
    title(main = s)   # site name as title
  }
  dev.off()
}


