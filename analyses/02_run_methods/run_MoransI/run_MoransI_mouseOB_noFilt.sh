#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=12G,h_fsize=100G

module load conda_R/4.1.x
/usr/bin/time -v Rscript run_MoransI_mouseOB_noFilt.R 2> ../../../outputs/memory/MoransI_mouseOB_noFilt.mem

