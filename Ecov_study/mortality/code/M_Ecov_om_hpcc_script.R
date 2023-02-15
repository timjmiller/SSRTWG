#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
simi = as.integer(args[1])
omj = as.integer(args[2])
emk = as.integer(args[3])
.libPaths("~/Rlib/")
library(wham)
source(file.path(here::here(), "Ecov_study", "mortality", "code", "sim_management.R"))
verify_version()

om_inputs <- readRDS(file.path(here::here(),"Ecov_study", "mortality", "inputs", "om_inputs.RDS"))
em_inputs <- readRDS(file.path(here::here(),"Ecov_study", "mortality", "inputs", "em_inputs.RDS"))
df.ems <- readRDS(file.path(here::here(),"Ecov_study", "mortality", "inputs", "df.ems.RDS"))
df.oms <- readRDS(file.path(here::here(),"Ecov_study", "mortality", "inputs", "df.oms.RDS"))
#######################################################
#need to have matching assumptions about CVs for catch and indices, too
obs_names <- c("agg_catch","agg_catch_sigma", "agg_indices", "agg_index_sigma", "catch_paa", "index_paa", 
  "Ecov_obs", "obs", "obsvec")
#######################################################

#######################################################
#I don't think we want to use the same (e.g. 1000) seeds for everything.
seeds <- readRDS(file.path(here::here(), "Ecov_study", "mortality", "inputs","seeds.RDS"))
#######################################################

cat(paste0("OM: ", omj, " Sim: ", simi, " EM: ", emk, "\n"))
write.dir <- file.path(here::here(),"Ecov_study", "mortality", "results", paste0("om", omj))
dir.create(write.dir, recursive = T, showWarnings = FALSE)
#script.full.path <- file.path(here::here(), "Ecov_study", "mortality", "code", "M_Ecov_om_sim_fit_script_hpcc.R") 
#system(paste0("Rscript --vanilla ", script.full.path, " " , this_om, " ",  this_em, " ", this_sim, " \n"))

# Set seed
om <- fit_wham(om_inputs[[omj]], do.fit = FALSE, MakeADFun.silent = TRUE)
#seeds are different for each om
#cat(omj)
#cat(simi)
#cat(length(seeds))
set.seed(seeds[[omj]][simi])
sim_data <- om$simulate(complete=TRUE)
truth <- sim_data
#save the version for reproducibility
truth$wham_version = om$wham_version
EM_input <- em_inputs[[emk]] # Read in the EM 
#put simulated data into the em input
EM_input$data[obs_names] = sim_data[obs_names]
#not estimating observation error in Ecov
EM_input$par$Ecov_obs_logsigma[] <- om_inputs[[omj]]$par$Ecov_obs_logsigma
res <- list(truth = truth)
res$fit <- list()
#do fit withouth sdreport first
fit <- tryCatch(fit_wham(EM_input, do.sdrep=F, do.osa=F, do.retro=T, do.proj=F, MakeADFun.silent=TRUE),
  error = function(e) conditionMessage(e))
  
# Deal with issues fitting EM to non-matching OM data
# empty elements below can be used to summarize convergence information
if(!'err' %in% names(fit) & class(fit) != "character"){
  res$fit <- fit[c("wham_version", "TMB_version", "opt", "final_gradient", "rep")]
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
#Ecov_study/mortality/results/om_x/simy_em_z.RDS
rds.fn = file.path(here::here(), "Ecov_study", "mortality", "results", paste0("om", omj), paste0("sim", simi, "_em", emk, ".RDS"))
#res_store <- readRDS(rds.fn)
#res_store[[this_em]] <- res
saveRDS(res, file = rds.fn)