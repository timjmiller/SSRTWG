#!/bin/bash

#$1 is first sim
#$2 is last sim
#$3 is first om
#$4 is last om
#$5 is first em
#$6 is last em
#$7 is ncores per job (20 or less being the number of ems specified)

for ((this_sim= $1; this_sim <= $2; this_sim++)) #sim in {$simstart..$simend}
do
  for ((this_om= $3; this_om <= $4; this_om++)) #om in {$omstart..$omend}
  do
    echo $this_om
    bsub -n $7 -q short -W 4:00 -o ~/logs/short/om_${this_om}_sim_${this_sim}_ems_${5}_to_${6}.log -R "rusage[mem=5000]" -R "span[hosts=1]" -J naa_${5}_to_${6}_${this_om}_${this_sim} "bash ~/SSRTWG/Project_0/code/naa_om_hpcc_args.sh $this_sim $this_sim $this_om $this_om $5 $6"
  done
done
echo "script is done"
