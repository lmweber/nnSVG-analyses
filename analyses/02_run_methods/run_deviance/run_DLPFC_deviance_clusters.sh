#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=20G,h_fsize=100G

module load conda_R/4.1.x
/usr/bin/time -v Rscript run_DLPFC_deviance_clusters.R 2> ../../../outputs/memory/deviance/DLPFC_deviance_clusters.mem

