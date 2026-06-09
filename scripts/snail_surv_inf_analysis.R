#!/usr/bin/env Rscript
# Title: snail_sur_inf_analysis.R
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2026-06-09
# Modified in:
# Licence: GPL v3



#==========#
# Comments #
#==========#

"Analyze snail survival and prevalence after exposure to schistosome parasites. Snails were either reared in the lab or collected in the field."



#==========#
# Versions #
#==========#

# v0.0 - 2026-06-09: creation



#==========#
# Packages #
#==========#

suppressMessages({
    library("magrittr")
    library("dplyr")
    library("tidyr")
    library("lme4")
    library("DHARMa")
    library("glmmTMB")
    library("rlang")
})



#===========#
# Variables #
#===========#

# Working directory
setwd(file.path(getwd(), "scripts"))

# Folders
data_fd  <- "../data/"
graph_fd <- "../graphs/"

# Files
fl <- paste0(data_fd, "snail_surv_inf_data.tsv")



#=================#
# Data processing #
#=================#

# Load data
df <- read.delim(fl)

df_long <- df %>%
  mutate(
    survived_lab    = Surviving.lab.snails,
    survived_field  = Surviving.field.snails,
    mortality_lab   = Lab.snail.exposed - Surviving.lab.snails,
    mortality_field = Field.snail.exposed - Surviving.field.snails,
    `non-inf_lab`   = Surviving.lab.snails - Shedder.lab,
    `non-inf_field` = Surviving.field.snails - Shedder.field,
    inf_lab         = Shedder.lab,
    inf_field       = Shedder.field
  ) %>%
  pivot_longer(
    cols      = c(survived_lab, survived_field, mortality_lab, mortality_field, `non-inf_lab`, `non-inf_field`, inf_lab, inf_field),
    names_to  = c(".value", "group"),
    names_sep = "_"
  ) %>%
  rename(died = mortality)


# Survival
cat("Proportion of surviving snails:\n")
df_long %>% {cbind(.[, c("Location", "group")], {.[,"survived"] / rowSums(.[,c("survived", "died")])} %>% round(., 2))}

## Model
cat("Survival GLM:\n")
surv_model <- glmmTMB(cbind(survived, died) ~ group + (1|Location),
               data = df_long %>% filter(survived + died > 0),
               family = betabinomial(link = "logit"))

summary(surv_model)

## Checks
sim_output <- simulateResiduals(surv_model, n = 1000)
pdf(paste0(graph_fd, "GLM_snail_survival_checks.pdf"))
plot(sim_output)                        # overall diagnostic
testDispersion(sim_output)
dev.off() %>% invisible()


# Infection
infection_long <- df_long %>% filter(`non-inf` + inf > 0)

## Model
cat("Infection GLM:\n")
infection_model <- glmer(cbind(inf, `non-inf`) ~ group + (1|Location),
               data = infection_long,
               family = binomial)
summary(infection_model)

## Checks
cat("overdispersion estimation:")
overdisp_ratio <- deviance(infection_model) / df.residual(infection_model)
print(overdisp_ratio)

sim_output <- simulateResiduals(infection_model, n = 1000)
pdf(paste0(graph_fd, "GLM_snail_infection_checks.pdf"))
plot(sim_output)                        # overall diagnostic
testDispersion(sim_output)              # confirm overdispersion is okay
testOutliers(sim_output)                # check for outliers
testQuantiles(sim_output)               # identify which quantiles fail
plotResiduals(sim_output, infection_long$group)      # by group
plotResiduals(sim_output, infection_long$Location) # by experiment
dev.off() %>% invisible()
