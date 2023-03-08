######################################
# Simulations: build simulated dataset
# Lukas Weber, Mar 2023
######################################

# this script builds a simulated dataset using empirical parameters from the 
# human DLPFC dataset


library(SpatialExperiment)
library(STexampleData)
library(here)


# ---------
# load data
# ---------

spe <- Visium_humanDLPFC()

# keep only spots over tissue
spe <- spe[, colData(spe)$in_tissue == 1]
dim(spe)

