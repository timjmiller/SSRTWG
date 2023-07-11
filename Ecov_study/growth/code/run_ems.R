library(snowfall)

### To do (no particular order):
## -Turn on random effects: NAA and Ecov (+sigma) when ready
## -Check set_* functions are working with growth options

## Make OM and EM inputs
source(file.path(here(), "Ecov_study", "growth", "code", "growth_Ecov_om_setup.R"))
source(file.path(here(), "Ecov_study", "growth", "code", "em_setup.R"))

## clear workspace otherwise gets pushed into remote sessions
## using snowfall
rm(list=ls())

df.ems = readRDS(file.path(here(), "Ecov_study", "growth", "inputs", "df.ems.RDS"))
df.oms = readRDS(file.path(here(), "Ecov_study", "growth", "inputs", "df.oms.RDS"))

run_iter <- function(sim, om, em){
  cmd <-
    paste("Rscript --vanilla growth_Ecov_om_hpcc_script.R", sim, om, em)
  system(cmd)
}

sfInit(parallel=TRUE, cpus=6)

sfExportAll()
# run_iter(sim = 1, om = 2, em = 5)
# trash <- sfLapply(1:30, function(sim) run_iter(sim,2,1))

for(om in 1:nrow(df.oms)){
  for(em in 1:nrow(df.ems)){
    sfExportAll()
    trash <- sfLapply(1:30, function(sim) run_iter(sim,om,em))
  }
}

sfStop()
