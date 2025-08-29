#!/bin/bash
#
#SBATCH --job-name=acmtfr_cv
#SBATCH --output=res_ACMTFR_CV.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=7-00:30:00
#SBATCH --mem-per-cpu=8000
#SBATCH --mail-user=g.r.ploeg@uva.nl
#SBATCH --mail-type=ALL

date
Rscript GOHTRANS_ACMTFR_model_selection_deltaT.R
date
