library(snowfall)
library(here)

## Make OM and EM inputs
#source("om_setup.R")
#source("em_setup.R")

## Hack to see which fixed effects are being estimated for each EM
## pars <- lapply(em_inputs, function(x) names(fit_wham(x, do.fit=FALSE)$par))
## par0 <- Reduce(intersect, pars)
## lapply(pars, function(x) x[-which(x %in% par0)])

## clear workspace otherwise gets pushed into remote sessions
## using snowfall
rm(list=ls())

df.ems <- readRDS(file.path(here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))
df.oms <- readRDS(file.path(here(),"Ecov_study", "recruitment_functions", "inputs", "df.oms.RDS"))

run_iter <- function(sim, om, em){
  cmd <- paste("Rscript --vanilla growth_Ecov_om_hpcc_script.R", sim,om,em)
  system(cmd)
}

sfInit(parallel=TRUE, cpus=8)

sfExportAll()
# run_iter(sim = 1, om = 2, em = 5)
# trash <- sfLapply(1:30, function(sim) run_iter(sim,2,1))

for(om in 1:nrow(df.oms)){
  for(em in 1:nrow(df.ems)){
    sfExportAll()
    trash <- sfLapply(1:8, function(sim) run_iter(sim,om,em))
  }
}

sfStop()
