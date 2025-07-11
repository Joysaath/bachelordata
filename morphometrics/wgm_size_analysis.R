# ---------------------------------------------
# Overview:
# This script loads preprocessed morphometric wing data,
# performs Generalised Procrustes Analysis (GPA), calculates centroid sizes
# centroid sizes, tests size differences between sites,
# investigates the influence of environmental gradients on size and
# analyses allometric effects between size and shape.
# ---------------------------------------------
# Input:
# - ‘wings_data.rds’: preprocessed landmark data of the wings
# ---------------------------------------------
# Output:
# - ‘procrustes_variances_size.csv’: Procrustes variances based
# on centroid sizes per site to analyse morphological variability
# ---------------------------------------------

### Load packages
library(geomorph)
library(tidyverse)
library(Morpho)

### Load preprocessed morphometric data
wings <- readRDS("wings_data.rds")

# Generalized Procrustes Analysis (GPA)
Y.gpa <- gpagen(wings$landmarks)

# Calculate Centroid Sizes
centroid_sizes <- Y.gpa$Csize

### ANOVA: Does centroid size differ between sites?
anova_model <- aov(centroid_sizes ~ wings$site)
summary(anova_model)

### Linear regression: Centroid size ~  gradients
regression_result3 <- lm(
  centroid_sizes ~ wings$gradient_tree + wings$gradient_crop + wings$landscape_no_so
)
summary(regression_result3)

### Morphological disparity based on centroid size
groups <- as.factor(wings$site)
res <- morphol.disparity(
  centroid_sizes ~ 1,
  groups = groups,
  iter = 10000,
  distance = "mahal"  # Mahalanobis distance used for size disparity
)

# Extract Procrustes variances
vars <- res$Procrustes.var

# Convert to data frame
var_df <- data.frame(
  group = names(vars),
  procrustes_variance = as.numeric(vars)
)

# Save to CSV
write.csv(var_df, "procrustes_variances_size.csv", row.names = FALSE)

### Allometry test
allom_test <- procD.lm(Y.gpa$coords ~ log(centroid_sizes), iter = 999)
summary(allom_test)