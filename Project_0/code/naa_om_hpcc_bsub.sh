#!/bin/bash

#$1 is first sim
#$2 is last sim
#$3 is first om
#$4 is last om

for ((this_sim= $1; this_sim <= $2; this_sim++)) #sim in {$simstart..$simend}
do
  for ((this_om= $3; this_om <= $4; this_om++)) #om in {$omstart..$omend}
  do
    echo $this_om
    bsub -n 20 -q short -W 2:00 -o ~/logs/logfile.%J -R "rusage[mem=5000]" -R "span[hosts=1]" -J naa_om_sim "bash ~/SSRTWG/Project_0/code/naa_om_hpcc_args.sh $this_sim $this_sim $this_om $this_om 1 20"
  done
done
echo "script is done"
