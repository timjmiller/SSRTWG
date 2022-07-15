#!/bin/bash
# Load gcc and R and jpeg header (for Hmisc)
module load gcc/8.1.0
module load R/4.0.4_gcc
module load libjpeg-turbo/2.0.2
module load openmpi/4.0.1

simstart = $1
simend = $2
omstart = $3
omend = $4
emstart = $5
emend = $6

echo "CPU threads: $(grep -c processor /proc/cpuinfo)"
grep 'cpu cores' /proc/cpuinfo | uniq

cd ~/SSRTWG
#Rscript --vanilla ~/SSRTWG/Project_0/code/naa_om_hpcc_script.R 7 1 1
for sim in {$simstart..$simend}
do
  for om in {$omstart..$omend}
  do
   for em in {$emstart..$emend}
   do
     echo $sim
     echo $om
     echo $em
     Rscript --vanilla ~/SSRTWG/Project_0/code/naa_om_hpcc_script.R $sim $om $em &
   done
 done
done
