library(here)
library(tidyverse)
library(rpart)

res.path <- 'E:/results_beta_fix'  # directory where simulation runs are (beta unstandardized)
res.dir <- 'results_beta_fix'   # 'results'     'results_beta_fix'   # results folder where AIC dataframes are
plot.dir <- 'plots'
plot.suffix <- '_beta_fix' 

## specify bad.grad.label and bad.se.value (these are the thresholds set in convergence_summaries.R to determine convergence) 
bad.grad.value <- 1E-6 #tim used 1E-6 (the abs of the exponent will be used for output suffix; ex: filename_grad_6.png)
bad.grad.label <- as.numeric(strsplit(as.character(bad.grad.value), split="-")[[1]][2])
bad.se.value <- 100 #tim used 100 (this value will be used for output suffix; ex: filename_se_100.png)

df.oms          <- readRDS(file.path(here::here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS"))
df.ems    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))
conv.runs <- readRDS(file.path(here::here(),'Ecov_study','recruitment_functions',res.dir, paste0("conv.runs_grad_", bad.grad.label, "_SE_", bad.se.value, plot.suffix, ".RDS")  ) )

conv.runs <- conv.runs %>%
  rename(Sim=sim) %>%
  mutate(ok.run =(bad.opt+ bad.conv+ bad.sdrep+ bad.grad+ bad.se.big) ) %>%
  replace_na(list(ok.run=1)) %>%
  relocate(OM, EM, Sim, ok.run, bad.opt, bad.conv, bad.sdrep, bad.grad, bad.se.big) %>%
  filter(ok.run==0)  # drop unconverged runs (there shouldn't be any at this point)

n_oms <- nrow(df.oms)
n_ems <- nrow(df.ems)
n_sims <- 100
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
  dat <- try(readRDS(file.path(res.path, paste0("om", conv.runs$OM[irun], '/','sim',conv.runs$Sim[irun],'_','em', conv.runs$EM[irun],'.RDS') ) ) )
  if( k==1)  year1 <- dat$truth$year1_model
  if(class(dat)!='try-error'){
    
    # grab annual parameter estimates
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
t2-t1  # Time difference of 2.479192 mins  ; session memory 1.1 GB (from liz runs, RDS w/o peels)
       # Time difference of 1.173523 hours (greg's runs w/ peels; 1.8 GB session memory)


df.oms2 <- as_tibble(cbind(OM=seq(1,256), df.oms)) 
df.ems2 <- cbind(EM=seq(1,6), df.ems)
colnames(df.ems2) <- c("EM", "EM_ecov_how", "EM_r_mod")
em_tib <- as_tibble(df.ems2) %>%
  mutate(SR=ifelse(EM_r_mod==2, 'Mean', 'BH')) %>%
  mutate(EM_mod = paste0(SR, "_", EM_ecov_how))

colnames(true_mat)[4:27] <- paste0(rep('True.'), colnames(true_mat)[4:27] )
colnames(true_mat_yearly)[5:16] <- paste0(rep('True.'), colnames(true_mat_yearly)[5:16] )

true_mat_tib <- as_tibble(true_mat)
true_mat_yearly_tib <- as_tibble(true_mat_yearly)

par_mat_tib <- as_tibble(par_mat) %>%
  select(-c(colnames(df.oms))) %>%
  left_join(df.oms2) %>%
  left_join(true_mat_tib) %>%
  mutate(True.NAA_sigmaC = True.NAA_sigma* sqrt(1-True.NAA_rho^2),
         True.Ecov_process_pars2C = True.Ecov_process_pars2* sqrt(1-True.Ecov_process_pars3^2) )  %>% #calculate 'marginal' sigma, given cor.
         mutate(RE.mean_rec1 = (mean_rec1-True.mean_rec1)/True.mean_rec1, 
         RE.mean_rec2 = ifelse(is.na(mean_rec2)==FALSE, (mean_rec2-True.mean_rec2)/True.mean_rec2, NA), 
         RE.FXSPR_static = (FXSPR_static-True.FXSPR_static)/True.FXSPR_static,
         RE.Y_FXSPR_static = (Y_FXSPR_static-True.Y_FXSPR_static)/True.Y_FXSPR_static,
         RE.SSB_FXSPR_static = (SSB_FXSPR_static-True.SSB_FXSPR_static)/True.SSB_FXSPR_static,
         RE.SPR0_static = (SPR0_static-True.SPR0_static)/True.SPR0_static,
         RE.SPR_FXSPR_static = (SPR_FXSPR_static-True.SPR_FXSPR_static)/True.SPR_FXSPR_static,
         RE.q1 = (q1-True.q1)/True.q1, RE.q2 = (q2-True.q2)/True.q2, 
         RE.NAA_sigma = (NAA_sigma-True.NAA_sigmaC)/True.NAA_sigmaC,
         RE.NAA_rho = (NAA_rho - True.NAA_rho)/True.NAA_rho, 
         RE.selpars_f1a = (selpars_f1a-True.selpars_f1a)/True.selpars_f1a,
         RE.selpars_f1b = (selpars_f1a-True.selpars_f1b)/True.selpars_f1b,
         RE.selpars_i1a = (selpars_i1a-True.selpars_i1a)/True.selpars_i1a,
         RE.selpars_i1b = (selpars_i1b-True.selpars_i1b)/True.selpars_i1b,
         RE.selpars_i2a = (selpars_i2a-True.selpars_i2a)/True.selpars_i2a,
         RE.selpars_i2b = (selpars_i2b-True.selpars_i2b)/True.selpars_i2b,
         RE.catch_paa_par = (catch_paa_par-True.catch_paa_par)/True.catch_paa_par,
         RE.index_paa_par1 = (index_paa_pars1-True.index_paa_pars1)/True.index_paa_pars1,
         RE.index_paa_par2 = (index_paa_pars2-True.index_paa_pars2)/True.index_paa_pars2,
         RE.Ecov_beta = (Ecov_beta-True.Ecov_beta)/True.Ecov_beta,
         RE.Ecov_process_pars1= (Ecov_process_pars1-True.Ecov_process_pars1)/True.Ecov_process_pars1,
         RE.Ecov_process_pars2= (Ecov_process_pars2-True.Ecov_process_pars2C)/True.Ecov_process_pars2C,
         RE.Ecov_process_pars3= (Ecov_process_pars3-True.Ecov_process_pars3)/True.Ecov_process_pars3
         )
  

par_mat_yearly_tib <- as_tibble(par_mat_yearly) %>%
  select(-c(colnames(df.oms))) %>%
  left_join(df.oms2) %>%
  left_join(true_mat_yearly_tib) %>%
  mutate(RE.FXSPR = (FXSPR-True.FXSPR)/True.FXSPR,
         RE.Y_FXSPR = (Y_FXSPR-True.Y_FXSPR)/True.Y_FXSPR,
         RE.SSB_FXSPR = (SSB_FXSPR-True.SSB_FXSPR)/True.SSB_FXSPR,
         RE.SPR0 = (SPR0-True.SPR0)/True.SPR0,
         RE.SPR_FXSPR = (SPR_FXSPR-True.SPR_FXSPR)/True.SPR_FXSPR,
         RE.SR_a = (SR_a-True.SR_a)/True.SR_a, RE.SR_b = (SR_b-True.SR_b)/True.SR_b
     #!!!!!!!!!! not calculating the MSY bias yet, as the True values are often NA or NaN  !!!!!
     #!!!!!!!!!! will likely have to make this calculation with my own function !!!!!!!
  )
         
saveRDS(par_mat_yearly_tib, file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0("par_mat_yearly", plot.suffix, ".RDS") ) )
saveRDS(par_mat_tib, file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0("par_mat", plot.suffix, ".RDS") ) )
saveRDS(true_mat_yearly_tib, file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0("true_mat_yearly", plot.suffix, ".RDS") ) )
saveRDS(true_mat_tib, file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0("true_mat", plot.suffix, ".RDS") ) )



