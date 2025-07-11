# ---------------------------------------------
# Overview:
# This script compares morphometric, genetic and geographical
# geographical distances of mosquitoes using landmark data,
# genetic distance matrices and environmental variables.
# It performs Mantel tests to analyse  association between
# morphological, genetic and environmental factors.
# ---------------------------------------------
# Input:
# - ‘wings_data.rds’: Pre-processed morphometric wing data
# - ‘genetic_distance.csv’: Genetic distance matrix
# - ‘data_collection_bachelor.csv’: Metadata with GPS coordinates
# - ‘df_wide_export.csv’: Prepared data with environmental variables
# ---------------------------------------------
# Output:
# - No direct file output; results of the Mantel tests
# are output in the console
# ---------------------------------------------



# Load required packages
library(tidyverse)    # Data manipulation and visualization
library(geomorph)     # Morphometric analysis
library(vegan)        # Ecological analysis, including Mantel test
library(Morpho)       # Additional morphometric tools

# Load morphometric data (assumed to include 3D landmark coordinates and metadata)
wings <- readRDS("wings_data.rds")

# Perform Generalized Procrustes Analysis (GPA) on landmark data
Y.gpa <- gpagen(wings$landmarks)  # Aligns shapes to remove size, orientation, and position

# Extract aligned coordinates
coords <- Y.gpa$coords  # Array: [Landmark, Dimension, Specimen_ID]

# Convert 3D coordinates to 2D matrix (specimens as rows)
coords_2d <- two.d.array(coords)

# Compute pairwise morphological distances (Euclidean) between specimens
dist_matrix <- dist(coords_2d)


### GEOGRAPHICAL DISTANCE AND MORPHOLOGICAL DISTANCE ###

# Extract longitude and a geographic grouping variable
geo_coords <- as.matrix(data.frame(
  Longitude = wings$Längengrad,
  Landscape_Group = wings$landscape_no_so
))

# Compute Euclidean distance matrix for geographic coordinates
dist_matrix_geo <- dist(geo_coords)

# Compute morphological (phenotypic) distance matrix again (optional if already done above)
phen_dist <- dist(coords_2d)

# Perform Mantel test between phenotype and geographic distance
mantel_phen_geo <- mantel(phen_dist, dist_matrix_geo)


### GENETIC DISTANCE MATRIX PREPARATION ###

# Load genetic distance matrix from CSV
data_gen <- read.csv("genetic_distance.csv", row.names = 1, check.names = FALSE)
data_gen[is.na(data_gen)] <- 0  # Replace NAs with zeros
gen_mat <- as.matrix(data_gen)

# Clean row/column names: keep only the numeric part before the first underscore
old_labels <- rownames(gen_mat)
new_labels <- sub("_.*", "", old_labels)

rownames(gen_mat) <- new_labels
colnames(gen_mat) <- new_labels

# Convert to dist object for analysis
genetic_dist_clean <- as.dist(gen_mat)


### MATCH GENETIC IDS WITH GEOGRAPHIC COORDINATES ###

# Load full metadata (with lat/lon and filename)
data <- read.table("data_collection_bachelor.csv", sep=";", dec=",", header=TRUE)

# Extract individual ID from filename (keep number before 'r' or 'l')
data$id <- sub("^(\\d+)[rl].*", "\\1", data$File.Name)

# Remove duplicates by ID
data_unique <- data[!duplicated(data$id), ]

# Filter rows where ID is in genetic distance labels
filtered_df <- data_unique[data_unique$id %in% new_labels, 
                           c("id", "Capture.location.LAT", "Capture.location.LONG")]
names(filtered_df) <- c("id", "Latitude", "Longitude")

# Calculate geographic distance matrix
coords2 <- as.matrix(filtered_df[, c("Longitude", "Latitude")])
dist_matrix2_geo <- dist(coords2)

# Match genetic distance to geographic coordinates
common_ids <- intersect(new_labels, filtered_df$id)
gen_mat_sub <- gen_mat[common_ids, common_ids]
gen_dis_geo_sub <- as.dist(gen_mat_sub)

# Perform Mantel test: genetic vs. geographic distance
mantel_gen_geo <- mantel(gen_dis_geo_sub, dist_matrix2_geo)


### COMPARISON: GENETIC VS PHENOTYPIC DISTANCE ###

# Ensure consistent labeling for both genetic and morphometric distances
df_wide <- read.table("df_wide_export.csv", sep=",", dec=".", header=TRUE) 
# Assign correct labels to morphometric distances
phen_dist <- dist(coords_2d)  # Recompute if needed
new_labels <- as.character(df_wide$id_base)  # Assumes df_wide contains cleaned IDs

phen_mat <- as.matrix(phen_dist)
rownames(phen_mat) <- new_labels
colnames(phen_mat) <- new_labels
phen_dist_named <- as.dist(phen_mat)

# Clean genetic matrix again to match phenotypic IDs
genetic_dist_clean <- as.dist(gen_mat)

# Check label mismatches
phen_labels <- labels(phen_dist_named)
gen_labels <- labels(genetic_dist_clean)

setdiff(phen_labels, gen_labels)  # In phen, not in gen
setdiff(gen_labels, phen_labels)  # In gen, not in phen
common_ids <- intersect(phen_labels, gen_labels)

# Subset both matrices to include only common individuals
phen_mat <- as.matrix(phen_dist_named)[common_ids, common_ids]
phen_dist_final <- as.dist(phen_mat)

gen_mat <- as.matrix(genetic_dist_clean)[common_ids, common_ids]
gen_dist_final <- as.dist(gen_mat)

# Mantel test: genetic vs. morphometric distances
mantel_gen_phen <- mantel(gen_dist_final, phen_dist_final)


### GENETIC DISTANCE VS. ENVIRONMENTAL VARIABLES ###

### Environmental Variable: Latitude ###
env_lat <- filtered_df["Latitude"]
rownames(env_lat) <- filtered_df$id
common_ids <- intersect(rownames(gen_mat), rownames(env_lat))

# Subset matrices
gen_mat_sub <- gen_mat[common_ids, , drop = FALSE]
env_mat_lat <- env_lat[common_ids, , drop = FALSE]

# Mantel test: genetic distance vs. latitude
mantel_gen_lat <- mantel(dist(gen_mat_sub), dist(env_mat_lat))


### Environmental Variable: Tree Cover ###
env_tree <- df_wide["tree_gradient"]
rownames(env_tree) <- df_wide$id_base
common_ids <- intersect(rownames(gen_mat), rownames(env_tree))

gen_mat_sub <- gen_mat[common_ids, , drop = FALSE]
env_mat_tree <- env_tree[common_ids, , drop = FALSE]

# Mantel test: genetic distance vs. tree cover
mantel_gen_tree <- mantel(dist(gen_mat_sub), dist(env_mat_tree))


### Environmental Variable: Crop Cover ###
env_crop <- df_wide["crop_gradient"]
rownames(env_crop) <- df_wide$id_base
common_ids <- intersect(rownames(gen_mat), rownames(env_crop))

gen_mat_sub <- gen_mat[common_ids, , drop = FALSE]
env_mat_crop <- env_crop[common_ids, , drop = FALSE]

# Mantel test: genetic distance vs. crop cover
mantel_gen_crop <- mantel(dist(gen_mat_sub), dist(env_mat_crop))
