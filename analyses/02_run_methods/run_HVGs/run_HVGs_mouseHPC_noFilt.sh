#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=10G,h_vmem=12G,h_fsize=100G

module load conda_R/4.1.x
/usr/bin/time -v Rscript run_HVGs_mouseHPC_noFilt.R 2> ../../../outputs/memory/HVGs_mouseHPC_noFilt.mem

