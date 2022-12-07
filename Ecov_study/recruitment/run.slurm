#!/bin/bash 
#SBATCH -n 4 #Request 4 tasks (cores)
#SBATCH -N 1 #Request 1 node
#SBATCH -t 12:00:00                          #Request runtime
#SBATCH -C centos7                           #Request types of nodes
#SBATCH -p sched_mit_darwin                  #Choose partition
#SBATCH --mem-per-cpu=4000                   #'000s of MB  memory per CPU
#SBATCH -o output_%j.txt                     #redirect output to output_JOBID.txt
#SBATCH -e error_%j.txt                      #redirect errors to error_JOBID.txt
#SBATCH --mail-type=BEGIN,END                #Mail when job starts and ends
#SBATCH --mail-user=gregleebritten@gmail.com #email recipient
echo "Executing R script"
Rscript ecov_r_sim_script.r
