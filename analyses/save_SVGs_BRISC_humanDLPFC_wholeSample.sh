#!/bin/bash
#$ -V
#$ -cwd
#$ -pe local 6
#$ -l mem_free=2G,h_vmem=3G,h_fsize=50G

# submit to cluster with 'qsub <filename.sh>'

module load conda_R/devel

Rscript save_SVGs_BRISC_humanDLPFC_wholeSample.R
