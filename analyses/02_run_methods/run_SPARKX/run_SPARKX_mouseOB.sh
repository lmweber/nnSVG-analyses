#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=3G,h_vmem=4G,h_fsize=100G

module load conda_R/4.1.x
/usr/bin/time -v Rscript run_SPARKX_mouseOB.R 2> ../../../outputs/memory/SPARKX_mouseOB.mem

