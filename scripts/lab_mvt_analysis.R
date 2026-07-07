#!/usr/bin/env Rscript
# Title: lab_mvt_analysis.R
# Version: 1.1
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2022-08-20
# Modified in: 2026-07-06
# Licence: GPL v3



#==========#
# Comments #
#==========#

"Analyze movement of individual worms from pools with different proportion of PZQ resistant and sensitive worms."



#==========#
# Versions #
#==========#

# v1.1 - 2026-07-06: change movement normalization / remove legend
# v1.0 - 2026-05-26: rewrite to focus on PZQ only
# v0.0 - 2022-08-20: creation



#==========#
# Packages #
#==========#

suppressMessages({
    library("magrittr")
})



#===========#
# Variables #
#===========#

# Working directory
setwd(file.path(getwd(), "scripts"))

# Folders
data_fd  <- "../data/"
res_fd   <- "../results/"
graph_fd <- "../graphs/"

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

# # Normalize data for each proportion
# my_col <- ncol(mvt) + 1
# for (i in unique(mvt[, 2])) {
# 	mvt[ mvt[, 2] == i, my_col] <- mvt[ mvt[,2] == i, ] %>% {.[, 6] / max(.[, 6])}
# }

# Normalize data overall
mvt[, ncol(mvt) + 1] <- mvt[, 6] / max(mvt[, 6])

# Identify genotype and associate color
mvt[, ncol(mvt) + 1] <- NA
for (i in 1:nrow(mvt)) {
	y <- geno[ apply(geno[, 2:4], 1, paste0, collapse="") == paste0(mvt[i, 2:4], collapse=""), 5] %>% {ifelse(. == "S/S", "blue", "red")}
	# if (length(y) == 0) y <- "black"
	if (length(y) == 0) y <- "blue"
	mvt[i, ncol(mvt)] <- y
}

# Reordering
for (i in unique(mvt[, 2])) {
	mvt[ mvt[, 2] == i, ] %<>% .[ order(.[,8], .[,6]), ]
}



#=========#
# Figures #
#=========#

if( ! dir.exists(graph_fd)) { dir.create(graph_fd, recursive = TRUE) }

gp_pc <- unique(mvt[, 2]) %>% paste0(., "%")

pdf(paste0(graph_fd, "Fig. 3 - movement_genotype.pdf"), width = 3, height = 4)
# par(mar = c(4, 4, 2.1, 0) + 0.1)
par(mar = c(4, 4, 0, 0) + 0.1)
pch <- 16

set.seed(myseed)
stripchart(mvt[,7] ~ mvt[,2], method = "jitter", pch = NA, vertical = TRUE, group.names = gp_pc, frame = FALSE, ylab = "Normalized index movement", xlab = "Percentage of resistant\nworms in pools", xaxt = "n", xlim = c(0.5, 2.5))

# Recreate the jittered x positions manually
groups <- as.numeric(as.factor(mvt[,2]))
set.seed(myseed)
x_jitter <- jitter(groups)

points(x_jitter, mvt[,7], pch = pch, col = mvt[, ncol(mvt)])

# X axis
axis(1, at = 1:length(gp_pc), labels = gp_pc, lwd = 0, lwd.ticks = 1)
x_diff <- (min(x_jitter) %>% round(., 1) - 1) %>% abs()
for (i in 1:length(gp_pc)) segments(i - x_diff, par("usr")[3], i + x_diff, , lwd = 2)

# mid_x <- grconvertX(0.5, from = "ndc", to = "user")
# legend(mid_x, par("usr")[4], legend = c("PZQ-ER", "PQZ-ES", "ND") , col = c("red", "blue", "black"), pch = pch, bty = "n", horiz = TRUE, xjust = 0.4, yjust = -0.2, xpd = TRUE)

dev.off()

