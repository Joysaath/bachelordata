# -----------------------------------------------
# Overview:
# This script performs a morphometric analysis of wing data.
# It includes a Generalised Procrustes Analysis (GPA),
# principal component analysis (PCA), visualisation of the PCA results
# with environmental gradients and the investigation of shape differences
# between groups based on tree cover.
# In addition, variances of individual landmarks between groups are analysed
# analysed and statistically compared.
# -----------------------------------------------
# Input:
# - ‘wings_data.rds’: Pre-processed 3D landmark coordinates
# and associated environmental data (tree cover, landscape variables)
# -----------------------------------------------
# Output:
# - PDF files with PCA plots (‘PCA_blanco2.pdf’, ‘PCA_tree.pdf’,
# ‘PCA_crop.pdf’, ‘PCA_sn.pdf’)
# - PDF with shape differences between groups (‘predicted_wingshape2.pdf’)
# - PDF with landmark variance plots (‘landmark_variance_plot.pdf’)
# - CSV files with variances and Levene test results
# -----------------------------------------------



# Load required libraries
# -----------------------------------------------
library(geomorph)
library(Morpho)
library(vegan)
library(dplyr)
library(ggplot2)
library(car)
# -----------------------------------------------
# Load and prepare data
# -----------------------------------------------
wings <- readRDS("wings_data.rds")

# -----------------------------------------------
# Generalized Procrustes Analysis (GPA) and PCA
# -----------------------------------------------
Y.gpa <- gpagen(wings$landmarks)
analyse <- gm.prcomp(Y.gpa$coords)

# Save PCA plot
pdf("PCA_blanco2.pdf", width = 7, height = 5)
plot(analyse, 
     axis1 = 1, 
     axis2 = 2, 
     main = "PCA of Procrustes-Aligned Shapes",
     asp = 1,
     pch = 1,          # open circle
     cex = 1.5,
     cex.lab = 1.5,
     cex.axis = 1.5)
grid(nx = NA, ny = NA)  # suppress grid lines
dev.off()

# -----------------------------------------------
# PCA with environmental gradients (ordisurf)
# -----------------------------------------------
y <- two.d.array(Y.gpa$coords)
pca_data <- as.data.frame(y)
vegan_pca <- rda(pca_data)

env_data <- data.frame(
  gradient_tree = wings$gradient_tree,
  gradient_sn = wings$landscape_no_so,
  gradient_crop = wings$gradient_crop
)

# Fit environmental surfaces
surf1 <- ordisurf(vegan_pca ~ gradient_tree, data = env_data, knots = 10, isotropic = TRUE)
surf2 <- ordisurf(vegan_pca ~ gradient_crop, data = env_data, knots = 10, isotropic = TRUE)
surf3 <- ordisurf(vegan_pca ~ gradient_sn, data = env_data, knots = 10, isotropic = TRUE)

# Helper function to plot PCA with surface
plot_pca_surf <- function(surf, title, file) {
  pdf(file, width = 7, height = 5)
  pca_scores <- scores(vegan_pca, display = "sites")
  plot(pca_scores[, 1], pca_scores[, 2],
       main = title,
       asp = 1,
       pch = 1,
       cex = 1.5,
       xlab = "PC1",
       ylab = "PC2",
       cex.lab = 1.5,
       cex.axis = 1.5)
  plot(surf, add = TRUE)
  dev.off()
}

# Save the plots
plot_pca_surf(surf1, "Tree Cover Across PCA Space", "PCA_tree.pdf")
plot_pca_surf(surf2, "Crop Cover Across PCA Space", "PCA_crop.pdf")
plot_pca_surf(surf3, "Latitude Across PCA Space", "PCA_sn.pdf")

# -----------------------------------------------
# Mean wing shape by tree cover group
# -----------------------------------------------
tree_groups <- factor(paste(wings$categorytree))
coords_split <- coords.subset(A = Y.gpa$coords, group = tree_groups)
mean_shapes <- lapply(coords_split, mshape)

vein_links <- matrix(c(
  1, 2, 2, 3, 3, 4, 4, 5, 5, 6,
  6, 7, 7, 8, 8, 9, 9, 10, 10, 11,
  11, 12, 12, 13, 13, 14, 14, 15, 15, 16,
  16, 17, 17, 18
), ncol = 2, byrow = TRUE)

# Plot wing shape difference between groups
pdf("predicted_wingshape2.pdf", width = 6, height = 4.5) 
GP1 <- gridPar(
  pt.bg = "black", pt.size = 1.3,
  tar.pt.bg = "red", tar.pt.size = 1.3,
  link.col = "red"
)

plotRefToTarget(
  mean_shapes$`viel baum`,   # high tree cover
  mean_shapes$`wenig baum`,  # low tree cover
  method = "points",
  mag = 1,
  links = vein_links,
  gridPars = GP1
)

title(main = "Predicted Wing Shape", cex.main = 1.5, font.main = 2)

legend(
  "topright",
  legend = c("High tree cover", "Low tree cover"),
  pch = 21,
  pt.bg = c("black", "red"),
  col = "black",
  pt.cex = 1.5,
  cex = 1.5,
  bty = "n"
)

coords <- mean_shapes$`viel baum`
text(coords, labels = 1:nrow(coords), pos = 3, cex = 1.3, col = "black")
dev.off()

# -----------------------------------------------
# Visualize and export variance in key landmarks
# -----------------------------------------------
# Function to extract landmark coordinates by group
extract_lm_data <- function(coords, group_name, lm_ids) {
  do.call(rbind, lapply(lm_ids, function(lm) {
    lm_data <- t(coords[lm, , ])
    data.frame(
      X = lm_data[, 1],
      Y = lm_data[, 2],
      Landmark = as.factor(paste0("LM ", lm)),
      Group = group_name
    )
  }))
}

lm_ids <- c(1, 2, 3, 10, 17, 18)
viel_coords <- coords_split[["viel baum"]]
wenig_coords <- coords_split[["wenig baum"]]

df_high <- extract_lm_data(viel_coords, "High tree cover", lm_ids)
df_low <- extract_lm_data(wenig_coords, "Low tree cover", lm_ids)
df_all <- rbind(df_high, df_low)

# Plot variance
ggplot(df_all, aes(x = X, y = Y, color = Group)) +
  geom_point(alpha = 0.7, linewidth = 2.5) +
  stat_ellipse(type = "norm", level = 0.68, linetype = "dashed", size = 0.8) +
  facet_wrap(~ Landmark, scales = "free") +
  theme_minimal(base_size = 15) +
  labs(
    title = "Variance in Landmarks 1, 2, 3, 10, 17, 18",
    x = "X-Coordinate",
    y = "Y-Coordinate",
    color = NULL
  ) +
  scale_color_manual(
    values = c("black", "red"),
    labels = c("High tree cover", "Low tree cover")
  ) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 17),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 20),
    legend.position = "top"
  )

ggsave("landmark_variance_plot.pdf", plot = last_plot(), width = 8, height = 6, dpi = 300)

# -----------------------------------------------
# Calculate and export variance by landmark & group
# -----------------------------------------------
variances <- df_all %>%
  group_by(Landmark, Group) %>%
  summarise(
    Variance_X = var(X),
    Variance_Y = var(Y),
    .groups = "drop"
  )

write.csv(variances, file = "landmark_variances.csv", row.names = FALSE)

levene_results <- list()

for (lm in unique(df_all$Landmark)) {
  subset_data <- df_all %>% filter(Landmark == lm)
  
  # Levene-Test for X
  levene_x <- leveneTest(X ~ Group, data = subset_data)
  
  # Levene-Test for Y
  levene_y <- leveneTest(Y ~ Group, data = subset_data)
  
  # p values
  levene_results[[as.character(lm)]] <- data.frame(
    Landmark = lm,
    p_X = levene_x$`Pr(>F)`[1],
    p_Y = levene_y$`Pr(>F)`[1]
  )
}

# summarize data
levene_df <- do.call(rbind, levene_results)

# save csv
write.csv(levene_df, file = "levene_test_results.csv", row.names = FALSE)
