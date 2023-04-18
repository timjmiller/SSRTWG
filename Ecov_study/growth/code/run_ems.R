library(snowfall)

### To do (no particular order):
## -Turn on random effects: NAA and Ecov (+sigma) when ready
## -Check set_* functions are working with growth options

## Make OM and EM inputs
source("growth_Ecov_om_setup.R")
source("em_setup.R")
## clear workspace otherwise gets pushed into remote sessions
## using snowfall
rm(list=ls())

run_iter <- function(sim, om, em){
  cmd <-
    paste("Rscript --vanilla growth_Ecov_om_hpcc_script.R", sim,om,em)
  system(cmd)
}

sfInit(parallel=TRUE, cpus=10)

## sfExportAll()
## run_iter(1,2,1)
## trash <- sfLapply(1:30, function(sim) run_iter(sim,2,1))

for(om in 1:2){
  for(em in 1:2){
    sfExportAll()
    trash <- sfLapply(1:20, function(sim) run_iter(sim,em,om))
  }
}
