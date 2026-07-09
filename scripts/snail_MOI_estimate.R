#!/usr/bin/env Rscript
# Title: snail_MOI_estimate.R
# Version: 0.0
# Author: Tim JC Anderson <tanderso@txbiomed.org>
# Created in: 2026-06-07
# Modified in:
# Licence: GPL v3



#==========#
# Comments #
#==========#

"Estimate the number of genotype per infected snail."



#==========#
# Versions #
#==========#

# v0.0 - 2026-06-07: creation



#===========#
# Functions #
#===========#

# Function to calculate MOI and genotype proportions
calc_MOI <- function(N, I){

  q <- I/N

  # estimated establishment probability
  p <- 1 - (1 - q)^(1/5)

  # MOI among infected snails
  moi <- (5*p)/q

  # exact prevalence CI
  bt <- binom.test(I, N)

  q_low  <- bt$conf.int[1]
  q_high <- bt$conf.int[2]

  # transform prevalence CI to MOI CI
  p_low  <- 1 - (1 - q_low)^(1/5)
  p_high <- 1 - (1 - q_high)^(1/5)

  moi_low  <- (5*p_low)/q_low
  moi_high <- (5*p_high)/q_high

  # probabilities for 1-5 genotypes
  probs <- dbinom(0:5, size = 5, prob = p)

  probs_inf <- probs[2:6] / q

  c(
    Prevalence = q,
    MOI = moi,
    MOI_LCL = moi_low,
    MOI_UCL = moi_high,
    Pvalue = bt$p.value,
    Pct1 = probs_inf[1]*100,
    Pct2 = probs_inf[2]*100,
    Pct3 = probs_inf[3]*100,
    Pct4 = probs_inf[4]*100,
    Pct5 = probs_inf[5]*100
  )
}



#===========#
# Variables #
#===========#

# System
options(width = 150)

# Working directory
setwd(file.path(getwd(), "scripts"))

# Folders
data_fd  <- "../data/"

# Files
fl <- paste0(data_fd, "snail_surv_inf_data.tsv")



#=================#
# Data processing #
#=================#

# Load data
df <- read.delim(fl)

dat <- data.frame(
  Village = 1:6,
  Examined = rowSums(df[, grepl("Surviving", colnames(df))]),
  Infected = rowSums(df[, grepl("Shedder", colnames(df))])
)


# Estimate MOI for each village
out <- cbind(
  dat,
  t(mapply(calc_MOI,
           N = dat$Examined,
           I = dat$Infected))
)

# Formatting
out$Prevalence <- round(out$Prevalence, 3)
out$MOI <- round(out$MOI, 2)

out$MOI_CI <- paste0(
  round(out$MOI_LCL,2),
  "–",
  round(out$MOI_UCL,2)
)

out[,c("Pct1","Pct2","Pct3","Pct4","Pct5")] <-
  round(out[,c("Pct1","Pct2","Pct3","Pct4","Pct5")],1)


# Final table

# Create formatted MOI column
out$MOI_summary <- sprintf(
  "%.2f (%.2f-%.2f)",
  out$MOI,
  out$MOI_LCL,
  out$MOI_UCL
)

final_table <- out[,c(
  "Village",
  "Examined",
  "Infected",
  "Prevalence",
  "MOI_summary",
  "Pct1",
  "Pct2",
  "Pct3",
  "Pct4",
  "Pct5"
)]

colnames(final_table) <- c(
  "Village",
  "Examined",
  "Infected",
  "Prevalence",
  "Mean_MOI_95CI",
  "%_1_genotype",
  "%_2_genotypes",
  "%_3_genotypes",
  "%_4_genotypes",
  "%_5_genotypes"
)

print(final_table, row.names = FALSE)
