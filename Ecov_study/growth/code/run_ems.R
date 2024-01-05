library(snowfall)

### To do (no particular order):
## -Turn on random effects: NAA and Ecov (+sigma) when ready
## -Check set_* functions are working with growth options

## Make OM and EM inputs
source("growth_Ecov_om_setup.R")
source("em_setup.R")

## Hack to see which fixed effects are being estimated for each EM
## pars <- lapply(em_inputs, function(x) names(fit_wham(x, do.fit=FALSE)$par))
## par0 <- Reduce(intersect, pars)
## lapply(pars, function(x) x[-which(x %in% par0)])


## clear workspace otherwise gets pushed into remote sessions
## using snowfall
rm(list=ls())

df.ems <- readRDS('../inputs/df.ems.RDS')
df.oms <- readRDS('../inputs/df.oms.RDS')

run_iter <- function(sim, om, em){
  cmd <-
    paste("Rscript --vanilla growth_Ecov_om_hpcc_script.R", sim,om,em)
  system(cmd)
}

sfInit(parallel=TRUE, cpus=18)
sfExportAll()

## run_iter(sim = 1, om = 2, em = 5)
## trash <- sfLapply(1:30, function(sim) run_iter(sim,2,1))

for(om in 1:nrow(df.oms)){
  for(em in 1:nrow(df.ems)){
    sfExportAll()
    trash <- sfLapply(1:100, function(sim) run_iter(sim,om,em))
  }
}

sfStop()
