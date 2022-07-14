#!/bin/bash
# Load gcc and R and jpeg header (for Hmisc)
module load gcc/8.1.0
module load R/4.0.4_gcc
module load libjpeg-turbo/2.0.2

cd ~
Rscript --vanilla naa_om_hpcc_script.R 7 7 480