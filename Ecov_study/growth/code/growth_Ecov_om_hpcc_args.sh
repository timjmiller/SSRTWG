#!/bin/bash
# Load gcc and R and jpeg header (for Hmisc)
module load gcc/8.1.0
module load R/4.0.4_gcc
module load libjpeg-turbo/2.0.2
module load openmpi/4.0.1

echo "CPU threads: $(grep -c processor /proc/cpuinfo)"
grep 'cpu cores' /proc/cpuinfo | uniq

cd ~/SSRTWG
for ((sim= $1; sim <= $2; sim++))
do
  for ((om= $3; om <= $4; om++))
  do
   for ((em= $5; em <= $6; em++))
   do
     echo $sim
     echo $om
     echo $em
     Rscript --vanilla ~/SSRTWG/Ecov_study/growth/code/M_Ecov_om_hpcc_script.R $sim $om $em &
   done
 done
done

echo "script is done"
