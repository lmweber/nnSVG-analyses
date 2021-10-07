#!/bin/bash
#$ -cwd
#$ -l mem_free=3G,h_vmem=4G,h_fsize=100G

module load conda_R/4.1.x
/usr/bin/time -v Rscript run_deviance_clusters_mOB.R 2> ../../../outputs/memory/deviance/deviance_clusters_mOB.mem

