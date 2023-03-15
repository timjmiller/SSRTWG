library(snowfall)

run_iter <- function(sim, om, em){
  cmd <-
    paste("Rscript --vanilla growth_Ecov_om_hpcc_script.R", sim,om,em)
  system(cmd)
}

sfInit(parallel=TRUE, cpus=10)
sfExportAll()

## no Ecov, growth not estimated
trash <- sfLapply(1:30, function(sim) run_iter(sim,1,2))
## no Ecov, growth estimated
trash <- sfLapply(1:30, function(sim) run_iter(sim,1,1))
