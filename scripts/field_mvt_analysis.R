#!/usr/bin/env Rscript
# Title: field_mvt_analysis.R
# Version: 1.1
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2024-11-08
# Modified in: 2026-07-06
# Licence: GPL v3



#==========#
# Comments #
#==========#

"Analyze movement of individual worms from pools generated from field schistosome parasites."



#==========#
# Versions #
#==========#

# v1.1 - 2026-07-06: remove legend
# v1.0 - 2026-05-27: polish code
# v0.0 - 2024-11-08: creation



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
res_fd   <- "../results/2-Analysis/field/"
graph_fd <- "../graphs/"

# Files
mvt_fl  <- paste0(res_fd, "Plate_mean.tsv")

myseed <- 7

myclr <- c(Hotspot = "red", Coldspot = "blue", Control = "black")



#=================#
# Data processing #
#=================#

# Load data
mvt  <- read.delim(mvt_fl, stringsAsFactors = FALSE)

# Normalize data for each proportion
mycol <- ncol(mvt) + 1
mytag <- apply(mvt[,c(2,4)], 1, paste0, collapse = "")
for (i in unique(mytag)) {
	mvt[ mytag == i, mycol] <- mvt[ mytag == i, 8] %>% {. / max(.)}
}

# Statistics by grouping based on location and treatment
mygroups <- apply(mvt[, c(7,2)], 1, paste0, collapse = " ")
kruskal.test(mvt[, 8], mygroups)
pw <- pairwise.wilcox.test(mvt[,8], mygroups, p.adjust.method = "bonferroni")
print(pw)
write.table(pw$p.value, paste0(res_fd, "Location_treatment_comp.tsv"), sep ="\t")



#=========#
# Figures #
#=========#

if( ! dir.exists(graph_fd)) { dir.create(graph_fd, recursive = TRUE) }

pdf(paste0(graph_fd, "Fig. 4 - movement_field.pdf"), width = 6, height = 4)
# par(mar = c(3, 4, 2.1, 0) + 0.1)
par(mar = c(3, 4, 0, 0) + 0.1)
pch <- 16
pad <- 0.2

# Uncommnet to order data by type of location
mvt %<>% .[order(.[,3], decreasing = TRUE), ]
loc <- unique(mvt[,2])
myclr_vec <- myclr[mvt[,3]] %>% replace(., mvt[, 7] == "Control", myclr["Control"])
n <- vector("list", length(loc))

set.seed(myseed)
stripchart(mvt[,9] ~ mvt[,2], method = "jitter", pch = NA,vertical = TRUE, frame = FALSE, ylab = "Normalized index movement", xaxt = "n")

# Recreate the jittered x positions manually
for (i in 1:length(loc)) {
    loc_i   <- loc[i]
    mvt_tmp <- mvt[ mvt[, 2] == loc_i, ]
    spot_tp <- mvt_tmp[, 3] %>% as.character()

    set.seed(myseed + i)
    x_jitter <- jitter(rep(i, nrow(mvt_tmp)), amount = 1/25)

    y <- mvt_tmp[, 7] == "Control"
    points(x_jitter[!y] - pad, mvt_tmp[!y, 9], col = myclr[spot_tp], pch = pch)
    points(x_jitter[y] + pad,  mvt_tmp[y, 9],  col = myclr["Control"], pch = pch)

    if (i == 1) { x_diff <- ((((min(x_jitter) - pad) * 10) %>% floor(.) / 10) - 1) %>% abs() }

    n[[i]] <- c(sum(!y), sum(y)) %>% as.character()
}

# X axis
axis(1, at = 1:length(loc), labels = loc, lwd = 0, lwd.ticks = 1)
for (i in 1:length(loc)) {
    segments(i - x_diff, par("usr")[3], i + x_diff, , lwd = 2)
    text(i - pad, par("usr")[3], label = bquote(italic(.(n[[i]][1]))), xpd = TRUE, cex = 0.6, adj = c(0.5, 1.5))
    text(i + pad, par("usr")[3], label = bquote(italic(.(n[[i]][2]))), xpd = TRUE, cex = 0.6, adj = c(0.5, 1.5))
}

# legend(mean(par("usr")[1:2]), par("usr")[4], legend = names(myclr), col = myclr, pch = pch, bty = "n", horiz = TRUE, xjust = 0.5, yjust = -0.2, xpd = TRUE)

dev.off()

