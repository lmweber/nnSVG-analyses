##########################
# Simulations: evaluations
# Lukas Weber, Mar 2023
##########################


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)


dir_sims <- here("outputs", "simulations")
dir_plots <- here("plots", "simulations")


# ------------
# load results
# ------------

res_list <- readRDS(here(dir_sims, "res_simulations.rds"))

length(res_list)

sim_names <- names(res_list)
sim_names

