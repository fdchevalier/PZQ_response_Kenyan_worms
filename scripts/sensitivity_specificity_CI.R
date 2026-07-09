#!/usr/bin/env Rscript
# Title: Sensitivity_Specificity_CI.R
# Version: 0.0
# Author: Tim JC Anderson <tanderso@txbiomed.org>
# Created in: 2026-06-03
# Modified in:
# Licence: GPL v3



#==========#
# Comments #
#==========#

"Estimate the number of genotype per infected snail."



#==========#
# Versions #
#==========#

# v0.0 - 2026-06-03: creation



#==========#
# Packages #
#==========#

suppressMessages({
    library("magrittr")
	library("binom")
})



#===========#
# Variables #
#===========#

# System
options(width = 150)

# Working directory
setwd(file.path(getwd(), "scripts"))

# Folders
data_fd  <- "../data/"
res_fd   <- "../results/"

# Files
geno_fl <- paste0(data_fd, "lab_worm_genotypes.tsv")
mvt_fl  <- paste0(res_fd, "2-Analysis/lab/Plate_mean.tsv")

myseed <- 7



#=================#
# Data processing #
#=================#

# Load data
geno <- read.delim(geno_fl, stringsAsFactors = FALSE)
mvt  <- read.delim(mvt_fl, stringsAsFactors = FALSE)

mvt_wm    <- mvt[mvt[,6] >= 1, ] %>% {apply(.[, 1:4], 1, paste0, collapse = "")}
no_mvt_wm <- mvt[mvt[,6] < 1, ]  %>% {apply(.[, 1:4], 1, paste0, collapse = "")}
res_wm    <- geno[ geno[, 5] == "R/R", ] %>% {apply(.[, 1:4], 1, paste0, collapse = "")}


# Data
TP <- sum(mvt_wm %in% res_wm)
FN <- sum(no_mvt_wm %in% res_wm)
TN <- nrow(mvt) - length(mvt_wm %in% res_wm)
FP <- sum(! mvt_wm %in% res_wm)

# Sensitivity
sens <- TP / (TP + FN)
sens.ci <- binom.confint(
  x = TP,
  n = TP + FN,
  methods = "exact"
)

# Specificity
spec <- TN / (TN + FP)
spec.ci <- binom.confint(
  x = TN,
  n = TN + FP,
  methods = "exact"
)

cat("Sensitivity =", round(sens, 4), "\n")
cat("95% CI =", round(sens.ci$lower, 4), "-", round(sens.ci$upper, 4), "\n\n")

cat("Specificity =", round(spec, 4), "\n")
cat("95% CI =", round(spec.ci$lower, 4), "-", round(spec.ci$upper, 4), "\n")
