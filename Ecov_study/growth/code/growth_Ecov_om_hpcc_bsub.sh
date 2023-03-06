#!/bin/bash

#$1 is first sim
#$2 is last sim
#$3 is first om
#$4 is last om
#$5 is first em
#$6 is last em
#$7 is ncores per job (20 or less being the number of ems specified)

for ((thisom= $3; thisom <= $4; thisom++)) #sim in {$simstart..$simend}
do
  for ((thisem= $5; thisem <= $6; thisem++)) #om in {$omstart..$omend}
  do
    echo $this_om
#    bsub -n $7 -q short -W 4:00 -o ~/logs/sims_${1}_to_${2}_om_${thisom}_em_${thisem}.log -R "rusage[mem=5000]" -R "span[hosts=1]" -J M_ecov_om_sim "bash ~/SSRTWG/Ecov_study/mortality/code/M_Ecov_om_hpcc_args.sh $1 $2 $thisom $thisom $thisem $thisem"
    bsub -n $7 -q short -W 4:00 -o ~/logs/short/sims_${1}_to_${2}_om_${thisom}_em_${thisem}.log -R "rusage[mem=5000]" -R "span[hosts=1]" -J ${1}_to_${2}_${thisom}_${thisem} "bash ~/SSRTWG/Ecov_study/mortality/code/M_Ecov_om_hpcc_args.sh $1 $2 $thisom $thisom $thisem $thisem"
  done
done
echo "script is done"
