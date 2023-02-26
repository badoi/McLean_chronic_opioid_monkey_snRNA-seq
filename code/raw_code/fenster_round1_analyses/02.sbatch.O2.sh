#!/bin/bash

# https://slurm.schedmd.com/sbatch.html

#SBATCH --partition=priority        # Partition (queue)
#SBATCH --time=5:00:00              # Runtime in D-HH:MM format
#SBATCH --job-name=seurat           # Job name
#SBATCH -c 1			    # cores
#SBATCH --mem=50G                   # Memory needed per CPU or --mem
#SBATCH --output=project_%j.out     # File to which STDOUT will be written, including job ID
#SBATCH --error=project_%j.err      # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL             # Type of email notification (BEGIN, END, FAIL, ALL)

source activate r
Rscript /n/data1/cores/bcbio/PIs/kerry_ressler/fenster2020_morphine_macaques_hbc03961/01.bcbio_to_sce.R
#Rscript --vanilla -e 'rmarkdown::render("02.quality_control.Rmd")'
#Rscript --vanilla -e 'rmarkdown::render("03.clusters.v2.inhib.Rmd")'
#Rscript --vanilla -e 'rmarkdown::render("07.qc_saturation.Rmd")'
source deactivate
