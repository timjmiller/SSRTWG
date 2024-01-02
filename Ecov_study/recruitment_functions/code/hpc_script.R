#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
simi = as.integer(args[1])
omj  = as.integer(args[2])
emk  = as.integer(args[3])
#simi <- 1; omj <- 1; emk <- 2

set.seed(simi)

library(wham)
library(here)
source(file.path(here::here(), "Ecov_study", "recruitment_functions", "code", "sim_management.R"))
#verify_version()

om_inputs <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "om_inputs.RDS"))
em_inputs <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "em_inputs.RDS"))
df.ems    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))
df.oms    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.oms.RDS"))

#extract info
x        <- data.frame(df.ems[emk,])
names(x) <- paste0('em_',names(x))
y        <- data.frame(df.oms[omj,])
names(y) <- paste0('om_',names(y))
model    <- cbind(im=simi, om=omj, em=emk, optimized=FALSE, sdreport=FALSE, y,x)

obs_names <- c("agg_catch","agg_catch_sigma", "agg_indices", "agg_index_sigma", "catch_paa", "index_paa", 
               "Ecov_obs", "obs", "obsvec")

#######################################################
#seeds <- readRDS(file.path(here::here(), "Ecov_study", "recruitment_functions", "inputs","seeds.RDS"))
#######################################################
cat(paste0("START OM: ", omj, " Sim: ", simi, " EM: ", emk, "\n"))

write.dir <- file.path(here::here(),"Ecov_study", "recruitment_functions", "results", paste0("om", omj))
dir.create(write.dir, recursive = T, showWarnings = FALSE)

om                 <- fit_wham(om_inputs[[omj]], do.fit = FALSE, MakeADFun.silent = TRUE)
sim_data           <- om$simulate(complete=TRUE)
truth              <- sim_data
truth$wham_version <- om$wham_version #save the version for reproducibility
EM_input           <- em_inputs[[emk]] # Read in the EM

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

## Build test object to test for initial 0 gradients?
test <- fit_wham(EM_input, do.fit=FALSE, do.sdrep=F, do.osa=F, do.retro=F, do.proj=F, MakeADFun.silent=TRUE)
if(any(test$gr()==0)) warning((paste0("Initial gradients 0 in OM: ", omj, " Sim: ", simi, " EM: ", emk, "\n")))


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
  fit$sdrep <- tryCatch(TMB::sdreport(fit), # no bc
           error = function(e) conditionMessage(e))
  if(class(fit$sdrep) == "sdreport"){
    res$model$sdreport <- TRUE
    res$fit$sdrep <- list(
       "Estimate_par" = as.list(fit$sdrep, what = "Est"),
       "SE_par" = as.list(fit$sdrep, what = "Std"),
       "Estimate_rep" = as.list(fit$sdrep, what = "Est", report = TRUE),
       "SE_rep" = as.list(fit$sdrep, what = "Std", report = TRUE))
    
    ##-PEELS--#######################
    
    fit$peels <- retro(fit, ran = unique(names(fit$env$par[fit$env$random])), do.sdrep = TRUE,
                       n.peels=n.yrs.peel, save.input = TRUE, MakeADFun.silent = TRUE,
                       retro.silent = TRUE, n.newton=1) #n.newton = 0 (2.25 min); n.newton=1 (3.1 min; 4.531353 mins with sdrep); n.newton=2 (3.9 min)
    # project from last peel ====
    tmp.proj.cont.ecov <- project_wham(
      model=fit$peels[[n.yrs.peel]],
      proj.opts = list(n.yrs = n.yrs.proj,  
                       use.last.F = FALSE, use.avg.F = FALSE, use.FXSPR = FALSE, use.FMSY = FALSE, 
                       proj.F = tail(truth$F, n.yrs.peel),# use Truth so can compare projections with Trut  
                       #proj.catch = NULL, 
                       avg.yrs = tail(fit$peels[[n.yrs.proj]]$input$years, yrs.avg), # this will cause selectivity to be averaged over last 5 years
                       cont.ecov = TRUE, use.last.ecov = FALSE, avg.ecov.yrs = NULL, proj.ecov = NULL, 
                       cont.Mre = NULL, 
                       #avg.rec.yrs = NULL,  #(only relevant if recr are fixed effects)
                       percentFXSPR = 100, percentFMSY = 100 #,
                       #proj_F_opt = rep(3,n.yrs.proj)
      ),
      n.newton = 3,
      do.sdrep = TRUE,
      MakeADFun.silent = TRUE,
      save.sdrep = TRUE
    )  
    
    tmp.proj.avg.ecov <- project_wham(
      model=fit$peels[[n.yrs.peel]],
      proj.opts = list(n.yrs = n.yrs.proj,  
                       use.last.F = FALSE, use.avg.F = FALSE, use.FXSPR = FALSE, use.FMSY = FALSE, 
                       proj.F = tail(truth$F, n.yrs.peel),  # use Truth so can compare projections with Truth
                       #proj.catch = NULL, 
                       avg.yrs = tail(fit$peels[[n.yrs.proj]]$input$years, yrs.avg), # this will cause selectivity to be averaged over last 5 years
                       cont.ecov = FALSE, use.last.ecov = FALSE, 
                       avg.ecov.yrs = tail(fit$peels[[n.yrs.proj]]$input$years, yrs.avg), 
                       proj.ecov = NULL, 
                       cont.Mre = NULL, 
                       #avg.rec.yrs = NULL,  #(only relevant if recr are fixed effects)
                       percentFXSPR = 100, percentFMSY = 100 #,
                       #proj_F_opt = rep(3,n.yrs.proj)
      ),
      n.newton = 3,
      do.sdrep = TRUE,
      MakeADFun.silent = TRUE,
      save.sdrep = TRUE
    )  
    
    tmp.proj.use.ecov <- project_wham(
      model=fit$peels[[n.yrs.peel]],
      proj.opts = list(n.yrs = n.yrs.proj,  
                       use.last.F = FALSE, use.avg.F = FALSE, use.FXSPR = FALSE, use.FMSY = FALSE, 
                       proj.F = tail(truth$F, n.yrs.peel),  # use Truth so can compare projections with Truth
                       #proj.catch = NULL, 
                       avg.yrs = tail(fit$peels[[n.yrs.proj]]$input$years, yrs.avg), # this will cause selectivity to be averaged over same years as ecov
                       cont.ecov = FALSE, use.last.ecov = FALSE, 
                       avg.ecov.yrs = NULL, 
                       proj.ecov = tail(fit$input$data$Ecov_obs, n.yrs.proj), 
                       cont.Mre = NULL, 
                       #avg.rec.yrs = NULL,  #(only relevant if recr are fixed effects)
                       percentFXSPR = 100, percentFMSY = 100 #,
                       #proj_F_opt = rep(3,n.yrs.proj)
      ),
      n.newton = 3,
      do.sdrep = TRUE,
      MakeADFun.silent = TRUE,
      save.sdrep = TRUE
    )  

        #res$peels <- fit$peels
    res$peels <- lapply(fit$peels, function(x) list(rep=x$rep,
                                                    sdrep= list(
                                                      "Estimate_par" = as.list(x$sdrep, what = "Est"),
                                                      "SE_par" = as.list(x$sdrep, what = "Std"),
                                                      "Estimate_rep" = as.list(x$sdrep, what = "Est", report = TRUE),
                                                      "SE_rep" = as.list(x$sdrep, what = "Std", report = TRUE) ) )
    )
    
    res$proj$cont.ecov$rep <- tmp.proj.cont.ecov$rep
    res$proj$cont.ecov$sdrep <- list(
      "Estimate_par" = as.list(tmp.proj.cont.ecov$sdrep, what = "Est"),
      "SE_par" = as.list(tmp.proj.cont.ecov$sdrep, what = "Std"),
      "Estimate_rep" = as.list(tmp.proj.cont.ecov$sdrep, what = "Est", report = TRUE),
      "SE_rep" = as.list(tmp.proj.cont.ecov$sdrep, what = "Std", report = TRUE)
    )
    
    res$proj$avg.ecov$rep <- tmp.proj.avg.ecov$rep
    res$proj$avg.ecov$sdrep <- list(
      "Estimate_par" = as.list(tmp.proj.avg.ecov$sdrep, what = "Est"),
      "SE_par" = as.list(tmp.proj.avg.ecov$sdrep, what = "Std"),
      "Estimate_rep" = as.list(tmp.proj.avg.ecov$sdrep, what = "Est", report = TRUE),
      "SE_rep" = as.list(tmp.proj.avg.ecov$sdrep, what = "Std", report = TRUE)
    )
    
    res$proj$use.ecov$rep <- tmp.proj.use.ecov$rep
    res$proj$use.ecov$sdrep <- list(
      "Estimate_par" = as.list(tmp.proj.use.ecov$sdrep, what = "Est"),
      "SE_par" = as.list(tmp.proj.use.ecov$sdrep, what = "Std"),
      "Estimate_rep" = as.list(tmp.proj.use.ecov$sdrep, what = "Est", report = TRUE),
      "SE_rep" = as.list(tmp.proj.use.ecov$sdrep, what = "Std", report = TRUE)
    )
  }
  
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

rds.fn = file.path(here::here(), "Ecov_study", "recruitment_functions", "results", paste0("om", omj), paste0("sim", simi, "_em", emk, ".RDS"))
saveRDS(res, file = rds.fn)
cat(paste0("END OM: ", omj, " Sim: ", simi, " EM: ", emk, "\n"))
