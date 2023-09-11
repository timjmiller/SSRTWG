#!/bin/bash 
#SBATCH -n 2                                #Request cores
#SBATCH -N 1                                 #Request nodes
#SBATCH -t 12:00:00                          #Request runtime
#SBATCH -C centos7                           #Request types of nodes
#SBATCH -p sched_mit_darwin                  #Choose partition
#SBATCH --mem-per-cpu=4000                   #'000s of MB  memory per CPU
#SBATCH -o output_%j.txt                     #redirect output to output_JOBID.txt
#SBATCH -e error_%j.txt                      #redirect errors to error_JOBID.txt
#SBATCH --array=1-96
module load R/4.1.0
echo "Executing R script"
Rscript run_ems.R $SLURM_ARRAY_TASK_ID