#!/bin/bash
#$ -cwd
#$ -pe local 10
#$ -l mem_free=10G,h_vmem=12G,h_fsize=100G

module load conda_R/4.1.x
/usr/bin/time -v Rscript run_MoransI_humanDLPFC_noFilt.R 2> ../../../outputs/memory/MoransI_humanDLPFC_noFilt.mem

