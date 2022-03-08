#!/bin/bash
#$ -cwd
#$ -pe local 10
#$ -l mem_free=20G,h_vmem=24G,h_fsize=100G

module load conda_R/4.1.x
/usr/bin/time -v Rscript run_SlideSeqHippo_nnSVG_logcounts.R 2> ../../../outputs/memory/nnSVG/SlideSeqHippo_nnSVG_logcounts.mem

