# ---------------------------------------------
# Overview:
# This script generates a clustered heatmap of pairwise genetic distances.
# Missing values are replaced with zeros. The color scale is customized 
# to highlight different distance ranges.
# ---------------------------------------------
#Input:
  # - "genetic_distance.csv": CSV file with pairwise genetic distances,
  #   rows and columns must match in order and naming.
# ---------------------------------------------
# Output:
# - "heatmap_plot.pdf": PDF file containing the heatmap.
# ---------------------------------------------

# loading libraries
library(pheatmap)

# loading matrix: genetic pairwise distances 
data <- read.csv("genetic_distance.csv", row.names = 1, check.names = FALSE)

# replace all NAs with 0
data[is.na(data)] <- 0

# subdivide value areas
breaks <- c(
  seq(0, 0.01, length.out = 21),     
  seq(0.011, 0.02, length.out = 21),  
  seq(0.021, 0.03, length.out = 11),  
  seq(0.031, 0.0315, length.out = 6)  
)
# colour gradient according to the breaks
colors <- c(
  colorRampPalette(c("lightblue", "blue"))(20),
  colorRampPalette(c("blue", "purple4"))(20),
  colorRampPalette(c("purple4", "black"))(10),
  rep("black", 5)  
)


# heatmap
pdf("heatmap_plot.pdf", width = 8, height = 6)

pheatmap(data,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colors,
         breaks = breaks,
         display_numbers = FALSE,
         legend = TRUE,
         show_rownames = FALSE,
         show_colnames = FALSE)

dev.off()
