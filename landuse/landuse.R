# ---------------------------------------------
# Overview:
# This script processes landscape data, performs a PCA (principal component analysis)
# for dimension reduction and identifies relevant gradients
# tree crop land and urban.
# It also checks whether there are significant correlations between these gradients.
# ---------------------------------------------
# INPUT:
# - ‘landuse_raw.csv’: Raw data on land use
# ---------------------------------------------
# OUTPUT:
# - ‘PCA_biplot_landuse.pdf’: Biplot of the PCA
# - ‘landuse_needed_gradients.csv’: Table with relevant gradients for subsequent analyses
# ---------------------------------------------

# load packages
library(ggplot2)
library(ggfortify)
library(vegan)
library(dplyr)

# read landscape data from CSV 
land_use <- read.table("landuse_raw.csv", header = TRUE, sep = ";", dec = ",")

# check for rows with any NAs
missing_rows <- which(rowSums(is.na(land_use)) > 0)

# view first few rows of the dataset
head(land_use)

# remove first three columns 
landuse_all <- land_use[, -c(1:3)]

# perform PCA 
PCA <- princomp(landuse_all)

# visualize PCA results
pdf("PCA_biplot_landuse.pdf", width = 8, height = 6)
biplot(PCA)
dev.off()

# summary of PCA components
summary(PCA)

## interpretation from PCA: crop_500m, tree_500m, and builtup_500m are the most important variables

# test correlations between gradients to check for redundancy
cor.test(land_use$tree_500m, land_use$builtup_500m)   # correlation
cor.test(land_use$tree_500m, land_use$crop_500m)      # no correlation
cor.test(land_use$builtup_500m, land_use$crop_500m)   # strong significant negative correlation 
##remove builtup

# new dataframe with only needed variables for further analysis
landuse <- data.frame(
  site = land_use$coord.Location,
  tree_gradient = land_use$tree_500m,
  crop_gradient = land_use$crop_500m,
  builtup_gradient = land_use$builtup_500m,
  noso_gradient = land_use$WGS84_Y,   # North-South gradient (latitude)
  longitude = land_use$WGS84_X        
)

# dataframe summary
write.table(
  landuse, 
  file = "landuse_needed_gradients.csv", 
  sep = ";", 
  dec = ",", 
  row.names = FALSE
)