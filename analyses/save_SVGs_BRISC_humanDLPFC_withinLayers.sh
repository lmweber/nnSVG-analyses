#!/bin/bash
#$ -pe local 10
#$ -l mem_free=1G,h_vmem=2G,h_fsize=10G

module load conda_R/devel

Rscript save_SVGs_BRISC_humanDLPFC_withinLayers.R
