#!/bin/bash
# Load gcc and R and jpeg header (for Hmisc)
module load gcc/8.1.0
module load R/4.0.4_gcc
module load libjpeg-turbo/2.0.2
module load openmpi/4.0.1

# sim start = $1
# sim end = $2
# om start = $3
# om end = $4
# em start = $5
# em end = $6

echo "CPU threads: $(grep -c processor /proc/cpuinfo)"
grep 'cpu cores' /proc/cpuinfo | uniq

cd ~/SSRTWG
#Rscript --vanilla ~/SSRTWG/Project_0/code/naa_om_hpcc_script.R 7 1 1
for sim in {$1..$2}
do
  for om in {$3..$4}
  do
   for em in {$5..$6}
   do
     echo $sim
     echo $om
     echo $em
     Rscript --vanilla ~/SSRTWG/Project_0/code/naa_om_hpcc_script.R $sim $om $em &
   done
 done
done
