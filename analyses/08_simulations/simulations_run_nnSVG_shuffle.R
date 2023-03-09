###############################################
# Simulations: run nnSVG - shuffled simulations
# Lukas Weber, Mar 2023
###############################################

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

sim_names_shuffle <- c(
  "sim_shuffle00", 
  "sim_shuffle01", 
  "sim_shuffle02", 
  "sim_shuffle03", 
  "sim_shuffle04", 
  "sim_shuffle05", 
  "sim_shuffle06", 
  "sim_shuffle07", 
  "sim_shuffle08", 
  "sim_shuffle09", 
  "sim_shuffle10"
)


# ---------
# run nnSVG
# ---------

# load each simulated dataset, run nnSVG, and store results

res_list <- as.list(rep(NA, length(sim_names_shuffle)))
names(res_list) <- sim_names_shuffle

for (s in seq_along(sim_names_shuffle)) {
  
  # load simulated dataset
  spe <- readRDS(here(dir_sims, paste0("spe_", sim_names_shuffle[s], ".rds")))
  
  # note: no additional filtering for simulated datasets
  
  # run nnSVG
  set.seed(123)
  spe <- nnSVG(
    spe, 
    n_threads = 10
  )
  
  print(paste0("Completed simulation: ", sim_names_shuffle[s]))
  
  # store results
  res_list[[s]] <- rowData(spe)
}


# ------------
# save results
# ------------

file <- file.path(dir_sims, "res_simulations_shuffle.rds")
saveRDS(res_list, file = file)

