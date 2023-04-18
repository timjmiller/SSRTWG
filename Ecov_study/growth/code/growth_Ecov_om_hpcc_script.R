#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
simi = as.integer(args[1])
omj = as.integer(args[2])
emk = as.integer(args[3])
## .libPaths("~/Rlib/")
library(wham)
source(file.path(here::here(), "Ecov_study", "growth", "code", "sim_management.R"))
## verify_version()
## simi=5;omj=1; emk=1
om_inputs <- readRDS(file.path(here::here(),"Ecov_study", "growth", "inputs", "om_inputs.RDS"))
em_inputs <- readRDS(file.path(here::here(),"Ecov_study", "growth", "inputs", "em_inputs.RDS"))
df.ems <- readRDS(file.path(here::here(),"Ecov_study", "growth", "inputs", "df.ems.RDS"))
df.oms <- readRDS(file.path(here::here(),"Ecov_study", "growth", "inputs", "df.oms.RDS"))
#######################################################
#need to have matching assumptions about CVs for catch and indices, too

x <- data.frame(df.ems[emk,])
names(x) <- paste0('em_',names(x))
y <- data.frame(df.oms[omj,])
names(y) <- paste0('om_',names(y))
model <- cbind(im=simi, om=omj, em=emk, optimized=FALSE, sdreport=FALSE, y,x)

## only new data for growth study is the marginal lengths in
## index_pal for survey 2
obs_names <- c("agg_catch","agg_catch_sigma", "agg_indices", "agg_index_sigma", "catch_paa", "index_paa",
  "Ecov_obs", "obs", "obsvec", "index_pal")
#######################################################

#######################################################
#I don't think we want to use the same (e.g. 1000) seeds for everything.
seeds <- readRDS(file.path(here::here(), "Ecov_study", "growth", "inputs","seeds.RDS"))
#######################################################

cat(paste0("START OM: ", omj, " Sim: ", simi, " EM: ", emk, "\n"))
write.dir <- file.path(here::here(),"Ecov_study", "growth", "results", paste0("om", omj))
dir.create(write.dir, recursive = T, showWarnings = FALSE)
#script.full.path <- file.path(here::here(), "Ecov_study", "growth", "code", "M_Ecov_om_sim_fit_script_hpcc.R")
#system(paste0("Rscript --vanilla ", script.full.path, " " , this_om, " ",  this_em, " ", this_sim, " \n"))

# Set seed
om <- fit_wham(om_inputs[[omj]], do.fit = FALSE, MakeADFun.silent = TRUE)
##seeds are different for each om
#cat(omj)
#cat(simi)
#cat(length(seeds))
set.seed(seeds[[omj]][simi])
sim_data <- om$simulate(complete=TRUE)
truth <- sim_data
#save the version for reproducibility
truth$wham_version = om$wham_version
EM_input <- em_inputs[[emk]] # Read in the EM
## EM_input$map$growth_a
## EM_input$map$SD_par
#put simulated data into the em input
EM_input$data[obs_names] = sim_data[obs_names]
#not estimating observation error in Ecov
EM_input$par$Ecov_obs_logsigma[] <- om_inputs[[omj]]$par$Ecov_obs_logsigma
## the fixed effects used to generate truth

ompars <- data.frame(par=names(om$par), value=om$par) |> dplyr::filter(par!='F_devs')
ompars$par2 <- sapply(unique(ompars$par), function(x) {
  y <- which(ompars$par==x)
  if(length(y)==1) return(x)
  x <- paste(x, 1:length(y), sep='_')
  return(x)
}) %>% unlist
res <- list(truth = truth, model=model, ompars=ompars)
res$fit <- list()

#do fit withouth sdreport first
fit <- tryCatch(fit_wham(EM_input, do.sdrep=F, do.osa=F, do.retro=F, do.proj=F, MakeADFun.silent=TRUE),
  error = function(e) conditionMessage(e))

# Deal with issues fitting EM to non-matching OM data
# empty elements below can be used to summarize convergence information
if(!'err' %in% names(fit) & class(fit) != "character"){
  res$model$optimized <- TRUE
  res$fit <- fit[c("wham_version", "TMB_version", "opt", "final_gradient", "rep")]
  ## res$fit$mohns_rho <- tryCatch(mohns_rho(fit),
  ##                               error = function(e) conditionMessage(e))
  ## fit$sdrep <- tryCatch(TMB::sdreport(fit), # no bc
  ##         error = function(e) conditionMessage(e))
  ## if(class(fit$sdrep) == "sdreport"){
  ##   res$model$sdreport <- TRUE
  ##   res$fit$sdrep <- list(
  ##     "Estimate_par" = as.list(fit$sdrep, what = "Est"),
  ##     "SE_par" = as.list(fit$sdrep, what = "Std"),
  ##     "Estimate_rep" = as.list(fit$sdrep, what = "Est", report = TRUE),
  ##     "SE_rep" = as.list(fit$sdrep, what = "Std", report = TRUE))
  ## }
  empars <- data.frame(par=names(res$fit$opt$par), value=res$fit$opt$par)%>%
    dplyr::filter(!grepl(x=par,'F_devs|log_NAA'))
  empars$par2 <- sapply(unique(empars$par), function(x) {
    y <- which(empars$par==x)
    if(length(y)==1) return(x)
    x <- paste(x, 1:length(y), sep='_')
    return(x)
  }) %>% unlist
  res$empars <- empars
}
#Ecov_study/growth/results/om_x/simy_em_z.RDS
rds.fn = file.path(here::here(), "Ecov_study", "growth", "results", paste0("om", omj), paste0("sim", simi, "_em", emk, ".RDS"))
#res_store <- readRDS(rds.fn)
#res_store[[this_em]] <- res
saveRDS(res, file = rds.fn)
cat(paste0("END OM: ", omj, " Sim: ", simi, " EM: ", emk, "\n"))
