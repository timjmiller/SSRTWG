#!/bin/bash
# Load gcc and R and jpeg header (for Hmisc)
module load gcc/8.1.0
module load R/4.0.4_gcc
module load libjpeg-turbo/2.0.2

echo "CPU threads: $(grep -c processor /proc/cpuinfo)"
grep 'cpu cores' /proc/cpuinfo | uniq

cd ~/SSRTWG
for sim in {7..7}
do
  for om in {1..24}
  do
    for em in {1..20}
    do
      Rscript --vanilla ~/SSRTWG/Project_0/code/naa_om_simf_fit_script_hpcc.R sim om em
    done
  done
done