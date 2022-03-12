#!/bin/bash
#$ -cwd
#$ -l mem_free=2G,h_vmem=3G,h_fsize=100G

module load conda_R/4.1.x
/usr/bin/time -v Rscript run_HVGs_humanDLPFC_noFilt.R 2> ../../../outputs/memory/HVGs_humanDLPFC_noFilt.mem

