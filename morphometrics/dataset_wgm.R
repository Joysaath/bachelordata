# ---------------------------------------------
# Overview:
# This script loads and prepares landmark data of the wings of mosquitoes.
# It cleans the data, assigns species, formats the landmark data,
# creates connection structures between landmarks and saves a
# data object for further morphometric analyses.
# ---------------------------------------------
# Input:
# - ‘Landmark_Koo_zusammen.csv’: Landmark data in long format
# with X/Y coordinates, labels and environmental variables
# ---------------------------------------------
# Output:
# - ‘wings_data.rds’: list with processed landmark data, types,
# connection structure and environmental variables for morphometric analyses
# ---------------------------------------------


# load packages
library(dplyr)
library(tidyr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tidyverse)

##### 1. Load and Prepare the Data #####

# load dataset 
data <- read.table("Landmark_Koo_zusammen.csv", sep = ";", dec = ".", header = TRUE)

# remove the first column
data <- data[, -1]

# convert decimal commas to points and scale X, Y coordinates
data$X <- as.numeric(gsub(",", ".", data$X)) / 10
data$Y <- as.numeric(gsub(",", ".", data$Y)) / 10
data$S_N_gradient <- as.numeric(gsub(",", ".", data$S_N_gradient))
data$tree_gradient <- as.numeric(gsub(",", ".", data$tree_gradient))
data$crop_gradient <- as.numeric(gsub(",", ".", data$crop_gradient))
data$builtup_gradient <- as.numeric(gsub(",", ".", data$builtup_gradient))
data$gradient <- as.numeric(gsub(",", ".", data$gradient))
data$Längengrad <- as.numeric(gsub(",", ".", data$Längengrad))

# add landmark number
data$number <- rep_len(1:18, nrow(data))

# assign default species
data$species <- rep_len("Cx.pipiens", nrow(data))

# reassign to Cx. pipiens molestus 
molestus <- c("26", "107", "111", "472", "477", "766", "526", "474")
data$species[rowSums(sapply(molestus, function(x) grepl(x, data$Label))) > 0] <- "Cx. pipiens molestus"

# reassign to Cx. torrentium 
torrentium <- c("411", "544", "747", "654", "539", "642", "678", "740", "656",
                "792", "664", "509", "714", "698", "737", "767", "545", "733", 
                "553", "745", "673", "653", "547", "511", "651", "681", "741", 
                "781", "464", "546", "542")

data$species[rowSums(sapply(torrentium, function(x) grepl(x, data$Label))) > 0] <- "Cx. torrentium"

# remove Cx. torrentium from the dataset
data <- data[!grepl("Cx\\. torrentium", data$species), ]

##### 2. Pivot Table (Reshape Data) #####

# remove duplicate rows: keep only the first occurrence per Label and landmark number
data_unique <- data %>%
  group_by(Label, number) %>%
  slice(1) %>%
  ungroup()

# reshape data to wide format: one row per individual, X and Y columns for each landmark
df_wide <- data_unique %>%
  pivot_wider(
    names_from = number,
    values_from = c(X, Y)
  )

# all reshaped landmark columns are numeric
stopifnot(all(sapply(df_wide[, grep("^X_|^Y_", colnames(df_wide))], is.numeric)))

# rename Label column to id and remove suffix after underscore 
df_wide <- df_wide %>%
  rename(id = Label) %>%
  mutate(id = sub("_.*", "", id))

# keep only left wing if both left ("l") and right ("r") exist; otherwise, keep whichever exists
df_wide <- df_wide %>%
  group_by(id_base = gsub("[lr]", "", id)) %>%
  filter(!(grepl("r$", id) & any(grepl("l$", id)))) %>%
  ungroup()

# remove any individuals with missing landmark coordinates
df_wide <- df_wide[!apply(df_wide[, grep("^X_|^Y_", colnames(df_wide))], 1, function(row) any(is.na(row))), ]

# Create connection structure (landmark links)
links_data <- data.frame(
  Spalte1 = c(1,1,2,3,4,4,5,5,6,6,7,7,8,8,9,9,10,11,12,13,14,14,15,16),
  Spalte2 = c(3,16,3,4,5,17,6,17,7,15,8,18,9,18,10,12,11,12,13,14,15,18,16,17)
)

### 2. Extract and reorder landmark coordinates
landmarks <- as.matrix(df_wide[, 11:46])

# Split X and Y
X_data <- landmarks[, 1:18]
Y_data <- landmarks[, 19:36]

# Interleave X and Y: X1,Y1,X2,Y2,...
landmarks_reordered <- as.matrix(
  do.call(cbind, lapply(1:18, function(i) cbind(X_data[, i], Y_data[, i])))
)

# Rename columns for clarity 
colnames(landmarks_reordered) <- as.vector(
  sapply(1:18, function(i) c(paste0("X", i), paste0("Y", i)))
)

# Create 3D landmark array (p, k, n): 18 landmarks, 2 dimensions, n specimens
array_3D <- arrayspecs(landmarks_reordered, 18, 2)

# Check array dimensions
dim(array_3D)  # Should return 18, 2, n

### 3. Add derived variables
df_wide$baumkategorie <- ifelse(df_wide$tree_gradient <= 19, "wenig baum", "viel baum")

### 4. Create connection structure (landmark links)
links_data <- data.frame(
  Spalte1 = c(1,1,2,3,4,4,5,5,6,6,7,7,8,8,9,9,10,11,12,13,14,14,15,16),
  Spalte2 = c(3,16,3,4,5,17,6,17,7,15,8,18,9,18,10,12,11,12,13,14,15,18,16,17)
)

### 5. Create list object for morphometric analysis
wings <- list(
  landmarks = array_3D,
  species = df_wide$species,
  links = links_data,
  site = df_wide$site,
  landscape_category = df_wide$category,
  landscape_gradient = df_wide$gradient,
  landscape_no_so = df_wide$S_N_gradient,
  gradient_tree = df_wide$tree_gradient,
  gradient_crop = df_wide$crop_gradient,
  gradient_builtup = df_wide$builtup_gradient,
  Längengrad = df_wide$Längengrad,
  categorytree = df_wide$baumkategorie
)

#save data
saveRDS(wings, file = "wings_data.rds")
