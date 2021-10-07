#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=30G,h_fsize=100G

module load conda_R/4.1.x
/usr/bin/time -v Rscript run_SlideSeqHippo_HVGs.R 2> ../../../outputs/memory/HVGs/SlideSeqHippo_HVGs.mem

