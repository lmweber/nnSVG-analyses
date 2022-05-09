#!/bin/bash
#$ -cwd
#$ -l mem_free=5G,h_vmem=6G,h_fsize=100G

module load conda_R/4.1.x
/usr/bin/time -v Rscript run_MoransI_mouseEmbryo_noFilt.R 2> ../../../outputs/memory/MoransI_mouseEmbryo_noFilt.mem

