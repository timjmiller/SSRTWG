blah = commandArgs(trailingOnly=TRUE)
om   = as.integer(blah[1])

library(snowfall)
library(doParallel)
library(here)

nsim <- 100

om_inputs <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "om_inputs.RDS"))
em_inputs <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "em_inputs.RDS"))
df.ems    <- readRDS(file.path(here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))
n.ems     <- nrow(df.ems)
df.oms    <- readRDS(file.path(here(),"Ecov_study", "recruitment_functions", "inputs", "df.oms.RDS"))

run_iter <- function(sim, om, em){
  cmd <- paste("Rscript --vanilla hpc_script.R", sim,om,em)
  system(cmd)
}

sfInit(parallel=TRUE, cpus=16)

for(em in 1:n.ems){
  sfExportAll()
  trash <- sfLapply(1:nsim, function(sim) run_iter(sim,om,em))
}

sfStop()
