# ---------------------------------------------
# Overview:
# This script loads preprocessed morphometric wing data,
# performs Generalised Procrustes Analysis (GPA), Principal Component
# Analysis (PCA) and Procrustes-ANOVA and analyses
# morphological variability (disparity) between sites.
# ---------------------------------------------
# Input:
# - ‘wings_data.rds’: pre-processed landmark data of the wings
# ---------------------------------------------
# Output:
# - ‘procrustes_variances_shape.csv’: Procrustes variants per site
# for the representation of morphological variability
# ---------------------------------------------

### load packages
library(geomorph)
library(tidyverse)
library(Morpho)

### load preprocessed morphometric data
wings <- readRDS("wings_data.rds")

### 1. visualize raw landmark configurations
plotAllSpecimens(wings$landmarks, label = FALSE)
title("Raw Specimens")

### 2. perform GPA
Y.gpa <- gpagen(wings$landmarks)
str(Y.gpa)  # Inspect structure of GPA result

### 3. visualize aligned landmark data
plotAllSpecimens(Y.gpa$coords, label = FALSE)
title("Aligned Specimens (GPA)")

### 4. PCA
analyse <- gm.prcomp(Y.gpa$coords)

### 5. summarize PCA results
summary(analyse)

### 6. prepare shape data for linear modeling
shape_data <- two.d.array(Y.gpa$coords)

### 7. Procrustes ANOVA 
model <- procD.lm(
  shape_data ~ wings$landscape_no_so + wings$gradient_tree + wings$gradient_crop,
  iter = 99  
)

summary(model)

result_site <- procD.lm(shape_data ~ wings$site, iter = 99)
summary(result_site)

### 8. variance of locations
group <- as.factor(wings$site)
gdf <- geomorph.data.frame(coords = Y.gpa$coords, group = group)

# run morphological disparity analysis
res <- morphol.disparity(coords ~ 1, groups = ~ group, data = gdf, iter = 499)

# extract procrustes variances 
vars <- res$Procrustes.var

# convert to data frame
var_df <- data.frame(
  group = names(vars),
  procrustes_variance = as.numeric(vars)
)

# save as CSV 
write.csv(var_df, "procrustes_variances_shape.csv", row.names = FALSE)