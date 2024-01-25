blah = commandArgs(trailingOnly=TRUE)
om   = as.integer(blah[1])

library(snowfall)
library(doParallel)
library(here)

om_inputs <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "om_inputs.RDS"))
em_inputs <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "em_inputs.RDS"))
df.ems    <- readRDS(file.path(here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))
n.ems     <- nrow(df.ems)
df.oms    <- readRDS(file.path(here(),"Ecov_study", "recruitment_functions", "inputs", "df.oms.RDS"))

run_iter <- function(sim, om, em){
  cmd <- paste("Rscript --vanilla hpc_script.R", sim,om,em)
  system(cmd)
}

source(file.path(here(),"Ecov_study","recruitment_functions","code","get_failed_jobs.R"))

sfInit(parallel=TRUE, cpus=4)

XX <- as.data.frame(fail.list[[om]]$iter_em)
fail.list[[om]]$iter_em <- dplyr::filter(XX,em.fails==6)
fails <- fail.list[[om]]

if(nrow(fails$iter_em)>0){
  sfExportAll()
  trash <- sfLapply(1:nrow(fails$iter_em), function(x) run_iter(as.integer(fails$iter_em[x,1]),
                                                         as.integer(om),
                                                         as.integer(fails$iter_em[x,2])))
}

sfStop()
