# ---------------------------------------------
# Overview:
# This script analyses genetic diversity and neutrality statistics
# Tajima's D, nucleotide diversity and haplotype diversity
# for an entire data set as well as population-specific.
# The input consists of alignments in FASTA format.
# The results are exported as CSV files.
# ---------------------------------------------
# INPUT:
# - FASTA file in the folder ‘alignment_fasta/’
# - Contain genetic sequence data for Culex mosquitoes
# ---------------------------------------------
# OUTPUT:
# - ‘diversity_summary_all.csv’: Results for entire dataset
# - ‘diversity_by_population.csv’: Results for each population
# ---------------------------------------------

# loading packages
library(devtools)
library(PopGenome)
library(car)
library(ggplot2)
library(stringr)
library(pegas)
library(neutralitytestr)
library(ape)

# -----------------------------------------
# ANALYSIS OF THE WHOLE DATASET
# -----------------------------------------

# loading alignments (FASTA-data from "coiculex/" file)
mydata <- readData("alignment_fasta/")

# summary
get.sum.data(mydata)

# neutrality stats (Tajima's D)
mydata <- neutrality.stats(mydata)

# diversity stats (pi, Hd)
mydata <- diversity.stats(mydata)

# extract values 
pi_raw_all    <- mydata@nuc.diversity.within
pi_norm_all   <- pi_raw_all / mydata@n.sites
haplo_div_all <- mydata@hap.diversity.within
tajima_d_all  <- mydata@Tajima.D

# show results
cat("\n--- Gesamtanalyse ---\n")
cat("Pi (roh):", round(pi_raw_all, 6), "\n")
cat("Pi (normalisiert):", round(pi_norm_all, 6), "\n")
cat("Haplotype Diversity (Hd):", round(haplo_div_all, 6), "\n")
cat("Tajima's D:", round(tajima_d_all, 6), "\n")

# export as CSV 
summary_all_df <- data.frame(
  Pi_raw             = round(as.numeric(pi_raw_all[1]), 6),
  Pi_normalized      = round(as.numeric(pi_norm_all[1]), 6),
  HaplotypeDiversity = round(as.numeric(haplo_div_all[1]), 6),
  TajimaD            = round(as.numeric(tajima_d_all[1]), 6)
)

write.csv(summary_all_df, "diversity_summary_all.csv", row.names = FALSE)

# -----------------------------------------
# ANALYSIS FOR POPULATIONS
# -----------------------------------------

# extract data
culex <- get.individuals(mydata)[[1]]
culex <- culex[1:349]

# clean labels
new_Labels <- sub("^(\\d+_[^A-Za-z0-9]*[A-Za-z0-9]+).*Alignmentof.*", "\\1", culex)

# assign individuals to populations
pop.wm  <- grep("Wm", new_Labels)
pop.la  <- grep("La", new_Labels)
pop.ha1 <- grep("Ha1", new_Labels)
pop.ha2 <- grep("Ha2", new_Labels)
pop.ha3 <- grep("Ha3", new_Labels)
pop.ha4 <- grep("Ha4", new_Labels)
pop.ru  <- grep("Ru", new_Labels)
pop.kr  <- grep("Kr", new_Labels)
pop.ko  <- grep("Ko", new_Labels)
pop.va  <- grep("Va", new_Labels)
pop.va2 <- grep("Va2", new_Labels)
pop.inn <- grep("In", new_Labels)
pop.bu  <- grep("Bu", new_Labels)
pop.kt  <- grep("Kt", new_Labels)
pop.ra  <- grep("Ra", new_Labels)
pop.ma  <- grep("Ma", new_Labels)
pop.li  <- grep("Li", new_Labels)
pop.fu  <- grep("Fu", new_Labels)
pop.el  <- grep("El", new_Labels)
pop.re  <- grep("Re", new_Labels)
pop.ba  <- grep("Ba", new_Labels)

# assign populations
mydata.pop <- set.populations(mydata, list(
  wm = pop.wm, la = pop.la,
  ha1 = pop.ha1, ha2 = pop.ha2, ha3 = pop.ha3, ha4 = pop.ha4,
  ru = pop.ru, kr = pop.kr, ko = pop.ko,
  va = pop.va, va2 = pop.va2,
  inn = pop.inn, bu = pop.bu,
  kt = pop.kt, ra = pop.ra,
  ma = pop.ma, li = pop.li,
  fu = pop.fu, el = pop.el,
  re = pop.re, ba = pop.ba
))

# neutrality stats per population
mydata.pop <- neutrality.stats(mydata.pop)

# diversity per population
mydata.pop <- diversity.stats(mydata.pop)

# extract values
pi_raw        <- as.vector(mydata.pop@nuc.diversity.within)
pi_norm       <- as.vector(pi_raw / mydata.pop@n.sites)
haplo_div     <- as.vector(mydata.pop@hap.diversity.within)
tajima_d      <- as.vector(mydata.pop@Tajima.D)

# get populationnames 
pop_names <- names(mydata.pop@populations)

# summary
populations_stats <- data.frame(
  Population         = pop_names,
  Pi_raw             = round(pi_raw, 6),
  Pi_normalized      = round(pi_norm, 6),
  HaplotypeDiversity = round(haplo_div, 6),
  TajimaD            = round(tajima_d, 6)
)

# show results
cat("\n--- Populationsweise Analyse ---\n")
print(populations_stats)

# export as CSV
write.csv(results_df, "diversity_by_population.csv", row.names = FALSE)