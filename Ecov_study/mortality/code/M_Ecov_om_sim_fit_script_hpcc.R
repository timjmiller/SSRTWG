#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

.libPaths("~/Rlib/")
library(here)
library(wham)

this_om <- as.integer(args[1])
this_em <- as.integer(args[2])
this_sim <- as.integer(args[3])

om_inputs <- readRDS(file.path(here(),"Ecov_study", "mortality", "inputs", "om_inputs.RDS"))
em_inputs <- readRDS(file.path(here(),"Ecov_study", "mortality", "inputs", "em_inputs.RDS"))
df.ems <- readRDS(file.path(here(),"Ecov_study", "mortality", "inputs", "df.ems.RDS"))
df.oms <- readRDS(file.path(here(),"Ecov_study", "mortality", "inputs", "df.oms.RDS"))
#######################################################
#need to have matching assumptions about CVs for catch and indices, too
obs_names <- c("agg_catch","agg_catch_sigma", "agg_indices", "agg_index_sigma", "catch_paa", "index_paa", 
  "Ecov_obs", "obs", "obsvec")
#######################################################

#######################################################
#I don't think we want to use the same (e.g. 1000) seeds for everything.
seeds <- readRDS(file.path(here(), "Ecov_study", "mortality", "inputs","seeds.RDS"))
#######################################################

cat(paste0("OM: ", this_om, " Sim: ", this_sim, " EM: ", this_em, "\n"))

# Set seed
om <- fit_wham(om_inputs[[this_om]], do.fit = FALSE, MakeADFun.silent = TRUE)
#seeds are different for each om
set.seed(seeds[[this_om]][this_sim])
sim_data <- om$simulate(complete=TRUE)
truth <- sim_data
#save the version for reproducibility
truth$wham_version = om$wham_version
EM_input <- em_inputs[[this_em]] # Read in the EM 
#put simulated data into the em input
EM_input$data[obs_names] = sim_data[obs_names]
#not estimating observation error in Ecov
EM_input$par$Ecov_obs_logsigma <- om_inputs[[this_om]]$par$Ecov_obs_logsigma
res <- list(truth = truth)
res$fit <- list()
#do fit withouth sdreport first
fit <- tryCatch(fit_wham(EM_input, do.sdrep=F, do.osa=F, do.retro=T, do.proj=F, MakeADFun.silent=TRUE),
  error = function(e) conditionMessage(e))
  
# Deal with issues fitting EM to non-matching OM data
# empty elements below can be used to summarize convergence information
if(!'err' %in% names(fit) & class(fit) != "character"){
  res$fit <- fit[c("wham_version", "TMB_version", "opt", "final_gradient", "rep")]
  # res$fit$wham_version <- fit$wham_version
  # res$fit$TMB_version <- fit$TMB_version
  # res$fit$opt <- fit$opt
  # res$fit$final_gradient <- fit$final_gradient
  # res$fit$rep <- fit$rep
  res$fit$mohns_rho <- tryCatch(mohns_rho(fit),
    error = function(e) conditionMessage(e))
  fit$sdrep <- tryCatch(TMB::sdreport(fit), # no bc
          error = function(e) conditionMessage(e))
  if(class(fit$sdrep) == "sdreport"){ 
    res$fit$sdrep <- list(
      "Estimate_par" = as.list(fit$sdrep, what = "Est"),
      "SE_par" = as.list(fit$sdrep, what = "Std"),
      "Estimate_rep" = as.list(fit$sdrep, what = "Est", report = TRUE),
      "SE_rep" = as.list(fit$sdrep, what = "Std", report = TRUE))
  }
}
rds.fn = file.path(here(), "Ecov_study", "mortality", "results", paste0("om", this_om, "_sim", this_sim, "_em", this_em, ".RDS"))
res_store <- readRDS(rds.fn)
res_store[[this_em]] <- res
saveRDS(res_store, file = rds.fn)
