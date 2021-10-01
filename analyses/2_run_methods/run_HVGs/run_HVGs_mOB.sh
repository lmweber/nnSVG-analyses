#!/bin/bash
#$ -cwd
#$ -l mem_free=2G,h_vmem=3G,h_fsize=100G

module load conda_R/4.1.x
Rscript run_HVGs_mOB.R

