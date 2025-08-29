#!/bin/bash
#
#SBATCH --job-name=acmtfr_fit
#SBATCH --output=res_ACMTFR_fit.txt
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=7-00:30:00
#SBATCH --mem-per-cpu=8000
#SBATCH --mail-user=g.r.ploeg@uva.nl
#SBATCH --mail-type=ALL

date
Rscript GOHTRANS_ACMTFR_fitModel.R
date
