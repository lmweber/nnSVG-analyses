#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=30G,h_fsize=100G

module load conda_R/4.1.x
/usr/bin/time -v Rscript run_SlideSeqHippo_deviance_onevsall.R 2> ../../../outputs/memory/deviance/SlideSeqHippo_deviance_onevsall.mem

