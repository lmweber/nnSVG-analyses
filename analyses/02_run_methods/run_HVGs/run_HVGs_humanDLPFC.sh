#!/bin/bash
#$ -cwd
#$ -l mem_free=2G,h_vmem=3G,h_fsize=100G

module load conda_R/4.1.x
/usr/bin/time -v Rscript run_HVGs_humanDLPFC.R 2> ../../../outputs/memory/HVGs_humanDLPFC.mem

