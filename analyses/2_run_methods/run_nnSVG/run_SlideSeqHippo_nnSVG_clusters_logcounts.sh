#!/bin/bash
#$ -cwd
#$ -pe local 10
#$ -l mem_free=10G,h_vmem=12G,h_fsize=300G

module load conda_R/4.1.x
/usr/bin/time -v Rscript run_SlideSeqHippo_nnSVG_clusters_logcounts.R 2> ../../../outputs/memory/nnSVG/SlideSeqHippo_nnSVG_clusters_logcounts.mem

