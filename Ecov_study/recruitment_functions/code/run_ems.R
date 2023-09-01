library(snowfall)
library(doParallel)
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
#rm(list=ls())

nsim <- 100

om_inputs <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "om_inputs.RDS"))
em_inputs <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "em_inputs.RDS"))
df.ems    <- readRDS(file.path(here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))
df.oms    <- readRDS(file.path(here(),"Ecov_study", "recruitment_functions", "inputs", "df.oms.RDS"))

run_iter <- function(sim, om, em){
  cmd <- paste("Rscript --vanilla hpc_script.R", sim,om,em)
  system(cmd)
}

x <- detectCores()      
sfInit(parallel=TRUE, cpus=x-1)

#for(om in 1:nrow(df.oms)){
for(om in 1:1){
#  for(em in 1:nrow(df.ems)){
 for(em in 2){
    sfExportAll()
    trash <- sfLapply(1:nsim, function(sim) run_iter(sim,om,em))
  }
}

sfStop()
