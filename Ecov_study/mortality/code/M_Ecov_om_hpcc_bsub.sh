#!/bin/bash

#$1 is first sim
#$2 is last sim
#$3 is first om
#$4 is last om
#$5 is first em
#$6 is last em
#$7 is ncores per job (20 or less being the number of ems specified)

for ((this_om= $3; this_om <= $4; this_om++)) #sim in {$simstart..$simend}
do
  for ((this_em= $5; this_em <= $6; this_em++)) #om in {$omstart..$omend}
  do
    echo $this_om
    bsub -n $7 -q short -W 4:00 -o "$(~/logs/sims_$1_to_$2_om_$this_om_em_$this_em.log)" -R "rusage[mem=5000]" -R "span[hosts=1]" -J M_ecov_om_sim "bash ~/SSRTWG/Ecov_study/mortality/code/M_Ecov_om_hpcc_args.sh $1 $2 $this_om $this_om $this_em $this_em"
  done
done
echo "script is done"
