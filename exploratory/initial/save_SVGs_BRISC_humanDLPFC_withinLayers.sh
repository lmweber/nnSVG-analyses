#!/bin/bash
#$ -V
#$ -cwd
#$ -pe local 6
#$ -l mem_free=4G,h_vmem=5G,h_fsize=100G

# submit to cluster with 'qsub <filename.sh>'

module load conda_R/devel

Rscript save_SVGs_BRISC_humanDLPFC_withinLayers.R
