#!/bin/bash
#$ -cwd
#$ -pe local 10
#$ -l mem_free=3G,h_vmem=4G,h_fsize=100G

module load conda_R/4.1.x
/usr/bin/time -v Rscript run_nnSVG_mOB_clusters.R 2> ../../../outputs/memory/nnSVG_mOB_clusters.mem

