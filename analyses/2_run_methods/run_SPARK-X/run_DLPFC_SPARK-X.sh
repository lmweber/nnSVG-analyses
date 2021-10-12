#!/bin/bash
#$ -cwd
#$ -pe local 4
#$ -l mem_free=3G,h_vmem=4G,h_fsize=100G

module load conda_R/4.1.x
/usr/bin/time -v Rscript run_DLPFC_SPARK-X.R 2> ../../../outputs/memory/SPARK-X/DLPFC_SPARK-X.mem

