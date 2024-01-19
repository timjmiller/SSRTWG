library(here)
library(tidyverse)
library(rpart)

res.dir <- file.path(here::here(), "Ecov_study","recruitment_functions", "results")
plot.dir <- 'plots'
plot.suffix <- '_beta_fix' 

df.oms          <- readRDS(file.path(here::here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS"))
df.ems    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))
no.conv.runs <- readRDS(file.path(here(),'Ecov_study','recruitment_functions','plots_lizruns', "no.conv.runs.RDS") )

conv.runs <- no.conv.runs %>%
  rename(Sim=sim) %>%
  mutate(ok.run =(bad.opt+ bad.conv+ bad.sdrep+ bad.grad+ bad.se.big) ) %>%
  replace_na(list(ok.run=1)) %>%
  relocate(OM, EM, Sim, ok.run, bad.opt, bad.conv, bad.sdrep, bad.grad, bad.se.big) %>%
  filter(ok.run==0)  # drop unconverged runs

n_oms <- nrow(df.oms)
n_ems <- nrow(df.ems)
n_sims <- 50
nyears <- 40

n.conv.runs <- nrow(conv.runs)
ten.pct <- round(n.conv.runs/10, 0)

# get parameters, only loop through converged runs
par_mat_yearly <- matrix(NA, nrow=nrow(conv.runs)*nyears, ncol=(ncol(df.oms)+ 16) )
colnames(par_mat_yearly) <- c(names(df.oms), 'OM', 'Sim', 'EM', 'Year', 'FXSPR','Y_FXSPR', 'SSB_FXSPR', 'SPR0', 'SPR_FXSPR', 'SR_a', 'SR_b', 'FMSY', 'SSB_MSY', 'YPR_MSY', 'SPR_MSY',  'R_MSY'  )

true_mat_yearly <- matrix(NA, nrow=nrow(conv.runs)*nyears, ncol=16 )
colnames(true_mat_yearly) <- c( 'OM', 'Sim', 'EM', 'Year', 'FXSPR','Y_FXSPR', 'SSB_FXSPR', 'SPR0', 'SPR_FXSPR', 'SR_a', 'SR_b', 'FMSY', 'SSB_MSY', 'YPR_MSY', 'SPR_MSY',  'R_MSY'  )


par_mat <- matrix(NA, nrow=nrow(conv.runs), ncol=(ncol(df.oms)+ 27) )
colnames(par_mat) <- c(names(df.oms), 'OM', 'Sim', 'EM',  'mean_rec1', 'mean_rec2', 'FXSPR_static', 'Y_FXSPR_static', 'SSB_FXSPR_static', 'SPR0_static', 'SPR_FXSPR_static',  'q1', 'q2', 'NAA_sigma', 'NAA_rho','selpars_f1a', 'selpars_f1b', 'selpars_i1a' , 'selpars_i1b','selpars_i2a' , 'selpars_i2b', 'catch_paa_par', 'index_paa_pars1', 'index_paa_pars2', 'Ecov_beta', 'Ecov_process_pars1', 'Ecov_process_pars2', 'Ecov_process_pars3'  )

true_mat <- matrix(NA, nrow=nrow(conv.runs), ncol=27 )
colnames(true_mat) <- c( 'OM', 'Sim', 'EM',  'mean_rec1', 'mean_rec2', 'FXSPR_static', 'Y_FXSPR_static', 'SSB_FXSPR_static', 'SPR0_static', 'SPR_FXSPR_static',  'q1', 'q2', 'NAA_sigma', 'NAA_rho','selpars_f1a', 'selpars_f1b', 'selpars_i1a' , 'selpars_i1b','selpars_i2a' , 'selpars_i2b', 'catch_paa_par', 'index_paa_pars1', 'index_paa_pars2', 'Ecov_beta', 'Ecov_process_pars1', 'Ecov_process_pars2', 'Ecov_process_pars3'  )



t1 <- Sys.time()
k <- 1

for(irun in 1:n.conv.runs) {    #   1:n.conv.runs
  dat <- try(readRDS(file.path(res.dir, paste0("om", conv.runs$OM[irun], '/','sim',conv.runs$Sim[irun],'_','em', conv.runs$EM[irun],'.RDS') ) ) )
  if( k==1)  year1 <- dat$truth$year1_model
  if(class(dat)!='try-error'){
    
    # grab annual parameter estimates
    # if(length(re.recr[[iom]][[jem]])==2000) {
    par_mat_yearly[((k-1)*nyears+1):(k*nyears),13:21] <- cbind(rep(conv.runs$OM[irun], nyears), rep(conv.runs$Sim[irun], nyears), rep(conv.runs$EM[irun], nyears),   seq(year1, (year1+nyears-1)),  exp(dat$fit$rep$log_FXSPR), exp(dat$fit$rep$log_Y_FXSPR), exp(dat$fit$rep$log_SSB_FXSPR), exp(dat$fit$rep$log_SPR0), exp(dat$fit$rep$log_SPR_FXSPR) )
    
    if(conv.runs$EM[irun]>2) par_mat_yearly[((k-1)*nyears+1):(k*nyears),22:28] <- c(exp(dat$fit$rep$log_SR_a),  exp(dat$fit$rep$log_SR_b), exp(dat$fit$rep$log_FMSY), exp(dat$fit$rep$log_SSB_MSY), exp(dat$fit$rep$log_YPR_MSY), exp(dat$fit$rep$log_SPR_MSY), exp(dat$fit$rep$log_R_MSY) )
    
   true_mat_yearly[((k-1)*nyears+1):(k*nyears),1:16] <- cbind(rep(conv.runs$OM[irun], nyears), rep(conv.runs$Sim[irun], nyears), rep(conv.runs$EM[irun], nyears),   seq(year1, (year1+nyears-1)), exp(dat$truth$log_FXSPR), exp(dat$truth$log_Y_FXSPR), exp(dat$truth$log_SSB_FXSPR), exp(dat$truth$log_SPR0), exp(dat$truth$log_SPR_FXSPR), exp(dat$truth$log_SR_a),  exp(dat$truth$log_SR_b), exp(dat$truth$log_FMSY), exp(dat$truth$log_SSB_MSY), exp(dat$truth$log_YPR_MSY), exp(dat$truth$log_SPR_MSY), exp(dat$truth$log_R_MSY) )
    
    # grab parameters that aren't annual values
    par_mat[k, 13:39] <- c(conv.runs$OM[irun], conv.runs$Sim[irun], conv.runs$EM[irun], exp(dat$fit$rep$mean_rec_pars[1]), ifelse(conv.runs$EM[irun]>2, exp(dat$fit$rep$mean_rec_pars[2]), NA), exp(dat$fit$rep$log_FXSPR_static), exp(dat$fit$rep$log_Y_FXSPR_static), exp(dat$fit$rep$log_SSB_FXSPR_static), exp(dat$fit$rep$log_SPR0_static), exp(dat$fit$rep$log_SPR_FXSPR_static), dat$fit$rep$q[nyears, 1], dat$fit$rep$q[nyears, 2], exp(dat$fit$sdrep$Estimate_par$log_NAA_sigma), (-1 + 2/(1 + exp(-2*dat$fit$sdrep$Estimate_par$trans_NAA_rho[2])) ), (-1 + 2/(1 + exp(-2*dat$fit$sdrep$Estimate_par$logit_selpars[1,11])) ), (-1 + 2/(1 + exp(-2*dat$fit$sdrep$Estimate_par$logit_selpars[1,12])) ), (-1 + 2/(1 + exp(-2*dat$fit$sdrep$Estimate_par$logit_selpars[2,11])) ), (-1 + 2/(1 + exp(-2*dat$fit$sdrep$Estimate_par$logit_selpars[2,12])) ), (-1 + 2/(1 + exp(-2*dat$fit$sdrep$Estimate_par$logit_selpars[3,11])) ), (-1 + 2/(1 + exp(-2*dat$fit$sdrep$Estimate_par$logit_selpars[3,12])) ),  dat$fit$sdrep$Estimate_par$catch_paa_pars[1,1], dat$fit$sdrep$Estimate_par$index_paa_pars[1,1], dat$fit$sdrep$Estimate_par$index_paa_pars[2,1], dat$fit$sdrep$Estimate_par$Ecov_beta[[1]], exp(dat$fit$sdrep$Estimate_par$Ecov_process_pars[1,1]), exp(dat$fit$sdrep$Estimate_par$Ecov_process_pars[2,1]), (-1 + 2/(1 + exp(-dat$fit$sdrep$Estimate_par$Ecov_process_pars[3,1])) ) )
    
   
    true_mat[k, 1:27] <- c(conv.runs$OM[irun], conv.runs$Sim[irun], conv.runs$EM[irun], exp(dat$truth$mean_rec_pars[1]), exp(dat$truth$mean_rec_pars[2]), exp(dat$truth$log_FXSPR_static), exp(dat$truth$log_Y_FXSPR_static), exp(dat$truth$log_SSB_FXSPR_static), exp(dat$truth$log_SPR0_static), exp(dat$truth$log_SPR_FXSPR_static), dat$truth$q[nyears,1], dat$truth$q[nyears,2], exp(dat$truth$log_NAA_sigma), (-1 + 2/(1 + exp(-2*dat$truth$trans_NAA_rho[2])) ), (-1 + 2/(1 + exp(-2*dat$truth$logit_selpars[1,11])) ), (-1 + 2/(1 + exp(-2*dat$truth$logit_selpars[1,12])) ), (-1 + 2/(1 + exp(-2*dat$truth$logit_selpars[2,11])) ), (-1 + 2/(1 + exp(-2*dat$truth$logit_selpars[2,12])) ), (-1 + 2/(1 + exp(-2*dat$truth$logit_selpars[3,11])) ), (-1 + 2/(1 + exp(-2*dat$truth$logit_selpars[3,12])) ), dat$truth$catch_paa_pars[1,1], dat$truth$index_paa_pars[1,1], dat$truth$index_paa_pars[2,1], dat$truth$Ecov_beta[[1]], exp(dat$truth$Ecov_process_pars[1,1]), exp(dat$truth$Ecov_process_pars[2,1] ), (-1 + 2/(1 + exp(-dat$truth$Ecov_process_pars[3,1])) )  )
  }  #end class(dat) check
  if( k %% ten.pct==0)  print(paste0('k = ', k, ' out of ', n.conv.runs, ' at ', Sys.time() ) )
  k=k+1
  
} #iom loop
t2 <- Sys.time()
t2-t1  # Time difference of 2.479192 mins  ; session memory 1.1 GB

as_tibble(true_mat)
as_tibble(par_mat)

saveRDS(par_mat_yearly, file.path(res.dir, paste0("par_mat_yearly", plot.suffix, ".RDS") ) )
saveRDS(par_mat, file.path(res.dir, paste0("par_mat", plot.suffix, ".RDS") ) )
saveRDS(true_mat_yearly, file.path(res.dir, paste0("true_mat_yearly", plot.suffix, ".RDS") ) )
saveRDS(true_mat, file.path(res.dir, paste0("true_mat", plot.suffix, ".RDS") ) )



