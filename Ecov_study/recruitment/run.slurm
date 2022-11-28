#!/bin/bash 
#SBATCH -n 4 #Request 4 tasks (cores)
#SBATCH -N 1 #Request 1 node
#SBATCH -t 0-12:00                           #Request runtime of 30 minutes
#SBATCH -C centos7                           #Request only Centos7 nodes
#SBATCH -p sched_mit_darwin                  #Run on sched_engaging_default partition
#SBATCH --mem-per-cpu=4000                   #Request 4G of memory per CPU
#SBATCH -o output_%j.txt                     #redirect output to output_JOBID.txt
#SBATCH -e error_%j.txt                      #redirect errors to error_JOBID.txt
#SBATCH --mail-type=BEGIN,END                #Mail when job starts and ends
#SBATCH --mail-user=gregleebritten@gmail.com #email recipient
echo "Executing R script"
Rscript ecov_r_sim_script.r
