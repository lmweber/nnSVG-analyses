#!/bin/bash
#$ -V
#$ -cwd
#$ -pe local 6
#$ -l mem_free=1G,h_vmem=2G,h_fsize=10G

# submit to cluster with 'qsub <filename.sh>'

module load conda_R/devel

Rscript save_SVGs_BRISC_humanDLPFC_withinLayers.R
