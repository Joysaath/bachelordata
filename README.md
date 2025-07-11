# Phylogenetic and Morphometric Study of _Culex pipiens_ s.s. in Germany

## Project overview
This project investigates the associations between environmental factors (land use, temperature, latitude), genetic variability, and wing characteristics (shape and size) of _Culex pipiens_ s.s. mosquito populations in Germany.

Mosquitoes of the _Culex pipiens_ complex are relevant vectors of the West Nile virus, which has been detected locally in Germany since 2018. Environmental changes, including land-use changes and climate change, may affect genetic and morphometric traits of mosquito populations, which could alter their vector competence.

In this study, up to 20 female mosquitoes from 23 sites in Germany were analyzed. Wing shape and size were analyzed using wing geometric morphometrics, and genetic variability was determined using mitochondrial COI gene sequencing.

## Repository Structure
- `/genetics`: Scripts and sequence data for the genetic analysis.
- `/landuse`: Scripts and data for the analysis of environmental factors.
- `/morphometrics`: Analyses and data on wing shape and size using morphometric methods.
- `/species_lda`: Script for linear discriminant analysis for species identification.
- `/germanymap`: Script for visualizing coordinates on a map of Germany.
- `/output`: Subfolder for plots, tables and results from the analyses.

## Dependecies
The following R packages are required for the analyses:

- `ape`
- `car`
- `caret`
- `devtools`
- `dplyr`
- `geomorph`
- `ggfortify`
- `ggplot2`
- `ggspatial`
- `knitr`
- `MASS`
- `MORPHO`
- `neutralitytestr`
- `pegas`
- `pheatmap`
- `PopGenome`
- `rnaturalearth`
- `rnaturalearthdata`
- `rnaturalearthhires`
- `sf`
- `stringr`
- `tidyverse`
- `tidyr`
- `vegan`
  
Install the packages in R with :
```r
install.packages(c(
  "ape",
  "car",
  "caret",
  "devtools",
  "dplyr",
  "geomorph",
  "ggfortify",
  "ggplot2",
  "ggspatial",
  "knitr",
  "MASS",
  "Morpho",
  "neutralitytestr",
  "pegas",
  "pheatmap",
  "PopGenome",
  "rnaturalearth",
  "rnaturalearthdata",
  "rnaturalearthhires",
  "sf"
  "stringr",
  "tidyverse",
  "tidyr",
  "vegan"
))
