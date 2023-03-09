###########################################
# Simulations: run nnSVG - main simulations
# Lukas Weber, Mar 2023
###########################################

# interactive session on compute cluster

# qrsh -pe local 10 -l mem_free=2G,h_vmem=3G,h_fsize=200G -now n
# cd /dcs04/hicks/data/lweber/nnSVG_analyses/nnSVG-analyses
# module load conda_R/4.2.x
# R


library(SpatialExperiment)
library(here)
library(nnSVG)


dir_sims <- here("outputs", "simulations")


# ---------------------------------------
# filenames for saved simulation datasets
# ---------------------------------------

sim_names <- c(
  "sim_largeBandwidth_fullExpr", 
  "sim_largeBandwidth_medExpr", 
  "sim_largeBandwidth_lowExpr", 
  "sim_medBandwidth_fullExpr", 
  "sim_medBandwidth_medExpr", 
  "sim_medBandwidth_lowExpr", 
  "sim_smallBandwidth_fullExpr", 
  "sim_smallBandwidth_medExpr", 
  "sim_smallBandwidth_lowExpr"
)


# ---------
# run nnSVG
# ---------

# load each simulated dataset, run nnSVG, and store results

res_list <- as.list(rep(NA, length(sim_names)))
names(res_list) <- sim_names

for (s in seq_along(sim_names)) {
  
  # load simulated dataset
  spe <- readRDS(here(dir_sims, paste0("spe_", sim_names[s], ".rds")))
  
  # note: no additional filtering for simulated datasets
  
  # run nnSVG
  set.seed(123)
  spe <- nnSVG(
    spe, 
    n_threads = 10
  )
  
  print(paste0("Completed simulation: ", sim_names[s]))
  
  # store results
  res_list[[s]] <- rowData(spe)
}


# ------------
# save results
# ------------

file <- file.path(dir_sims, "res_simulations.rds")
saveRDS(res_list, file = file)

