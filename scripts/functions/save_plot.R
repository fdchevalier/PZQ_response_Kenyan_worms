# Title: save_plot.R
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2022-08-22
# Modified in:



#==========#
# Comments #
#==========#

# Intersect list of genomic positions with position intervals similarly to the bed intersect command.



#==========#
# Versions #
#==========#

# v0.0 - 2022-08-22: creation



#==========#
# Function #
#==========#

# Saving plots as pdf and png
## source: https://stackoverflow.com/a/26233113
save_plot <- function(filename, width = 7, height = 7, res = 300, off = FALSE) {
    if (! off) {
        # PDF device
        pdf(paste0(filename, ".pdf"), width = width, height = height)

        # PNG device
        png(paste0(filename, ".png"), width = width * res, height = height * res, res = res)
        dev.control("enable")
    } else {
        pdf_dev <- dev.list()["pdf"]
        dev.copy(which = pdf_dev)
        graphics.off()
    }
}

