#!/bin/bash
#$ -cwd
#$ -pe local 10
#$ -l bluejay,mem_free=5G,h_vmem=10G,h_fsize=100G

module load conda_R/4.1.x
/usr/bin/time -v Rscript run_SPARKX_mouseHPC_noFilt.R 2> ../../../outputs/memory/SPARKX_mouseHPC_noFilt.mem

