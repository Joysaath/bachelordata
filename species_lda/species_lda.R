# ---------------------------------------------
# Overview:
# This script analyses landmark data of the wings of mosquitoes
# using Generalised Procrustes Analysis (GPA) and Linear Discriminant Analysis (LDA).
# The aim is to check whether the Culex species identified using the COI gene are also morphometrically distinguishable # on the basis of their wing shape (WGM).
# ---------------------------------------------
# INPUT:
# - ‘Landmark_Koo_zusammen.csv’: Landmark data in long format (incl. X/Y coordinates, label etc.)
# ---------------------------------------------
# OUTPUT:
# - ‘accuracy_per_class.csv’: LDA recognition accuracy per species
# - ‘species_lda_plot.pdf’: Plot of the LDA projection (LD1 vs. LD2)
# ---------------------------------------------

# load packages
library(geomorph)
library(tidyr)
library(dplyr)
library(knitr)
library(ggplot2)
library(MASS)
library(caret)

### 1. Load and Prepare Data

# load CSV data 
data <- read.table("Landmark_Koo_zusammen.csv", sep=";", dec=".", header=TRUE)
# remove first column
data <- data[, -1]  

# convert coordinate columns and scale 
data$X <- as.numeric(gsub(",", ".", data$X)) / 10
data$Y <- as.numeric(gsub(",", ".", data$Y)) / 10

# add landmark numbers (1 to 18, repeated)
data$number <- rep_len(1:18, nrow(data))

# initialize species as Cx. pipiens
data$species <- rep_len("Cx.pipiens", nrow(data))

# assign Cx. pipiens molestus based on COI label
molestus <- c("26", "107", "111", "472", "477", "766", "526", "474")
data$species[rowSums(sapply(molestus, function(x) grepl(x, data$Label))) > 0] <- "Cx. pipiens molestus"

# assign Cx. torrentium based on COI label
torrentium <- c("411", "544", "747", "654", "539", "642", "678", "740", "656",
                "792", "664", "509", "714", "698", "737", "767", "545", "733", 
                "553", "745", "673", "653", "547", "511", "651", "681", "741", 
                "781", "464", "546", "542")

data$species[rowSums(sapply(torrentium, function(x) grepl(x, data$Label))) > 0] <- "Cx. torrentium"
data$species[rowSums(sapply(molestus, function(x) grepl(x, data$Label))) > 0] <- "Cx. pipiens molestus"

# remove duplicates 
data_unique <- data %>%
  group_by(Label, number) %>%
  slice(1) %>%
  ungroup()

# convert to wide format
df_wide <- data_unique %>%
  pivot_wider(names_from = number, values_from = c(X, Y))

# clean up labels
df_wide <- df_wide %>%
  rename(id = Label)
df_wide$id <- sub("_.*", "", df_wide$id)

# remove right wings if left wing exists for same individual
df_wide <- df_wide %>%
  group_by(id_base = gsub("[lr]", "", id)) %>%
  filter(!(grepl("r$", id) & any(grepl("l$", id)))) %>%
  ungroup()

# remove rows with any missing coordinates
df_wide <- df_wide[!apply(df_wide[, grep("^X_|^Y_", colnames(df_wide))], 1, function(row) any(is.na(row))), ]

### 2. Build Landmark Matrix and Convert to Array

# extract only landmark columns
landmarks <- as.matrix(df_wide[, grep("^X_|^Y_", colnames(df_wide))])

# split X and Y columns
X_data <- landmarks[, 1:18]
Y_data <- landmarks[, 19:36]

# interleave X and Y to match landmark structure
landmarks_reordered <- as.matrix(
  do.call(cbind, lapply(1:18, function(i) cbind(X_data[, i], Y_data[, i])))
)
colnames(landmarks_reordered) <- as.vector(
  sapply(1:18, function(i) c(paste0("X", i), paste0("Y", i)))
)

# convert to 3D array: [landmarks, dimensions, specimens]
array_3D <- arrayspecs(landmarks_reordered, 18, 2)

### 3. Build Final Object with Needed Data Only

wings <- list()
wings$landmarks <- array_3D
wings$species <- df_wide$species

### 4. Generalized Procrustes Analysis (GPA)
# perform GPA 
Y.gpa <- gpagen(wings$landmarks)

### 5. LDA Analysis with Cross-Validation

# convert 3D GPA coords to 2D data frame
coords_flat <- two.d.array(Y.gpa$coords)
df <- as.data.frame(coords_flat)
df$species <- as.factor(wings$species)

# LDA with cross-validation
lda_model <- lda(species ~ ., data = df, CV = TRUE)

# confusion matrix and class-wise recall
confusion_matrix <- table(lda_model$class, df$species)

recall_per_class <- diag(confusion_matrix) / rowSums(confusion_matrix)
recall_percent <- round(recall_per_class * 100, 2)  # output recall in percent
recall_df <- data.frame(
  Species = names(recall_percent),
  RecallPercent = recall_percent
)
write.csv(recall_df, "accuracy_per_class.csv", row.names = FALSE)

### 6. LDA Projection for Visualization

# split into training (80%) and test (20%)
set.seed(123)
train_index <- createDataPartition(df$species, p = 0.8, list = FALSE)
train_data <- df[train_index, ]
test_data <- df[-train_index, ]

# fit LDA on training data
lda_model_vis <- lda(species ~ ., data = train_data)
predictions <- predict(lda_model_vis, newdata = test_data)

# prepare for plotting
lda_data <- as_tibble(predictions$x)
lda_data$class <- predictions$class

# plot LD1 and LD2
species_lda <- ggplot(lda_data, aes(x = LD1, y = LD2, color = class)) +
  geom_jitter(height = 0.1, size = 2, alpha = 0.7) +
  labs(title = "LDA Projection", x = "LD1", y = "LD2") +
  theme_minimal()

ggsave("species_lda_plot.pdf", plot = species_lda, device = "pdf", width = 8, height = 6)
