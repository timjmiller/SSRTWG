library(here)
library(tidyverse)
library(rpart)
library(rpart.plot)

res.path <- 'E:/results_beta_fix'  # directory where simulation runs are (beta unstandardized)
res.dir <- 'results_beta_fix'   # 'results'     'results_beta_fix'   # results folder where AIC dataframes are
plot.dir <- 'plots_beta_fix'    # 'plots_lizruns'  'plots_beta_fix'   
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

# look at summary of RE distributions ====  
RE.par.sum <- rbind(mean_rec1=summary(par_mat_tib$RE.mean_rec1)[1:6],
                    mean_rec2=summary(par_mat_tib$RE.mean_rec2)[1:6],
                    FXSPR_static=summary(par_mat_tib$RE.FXSPR_static)[1:6],
                    Y_FXSPR_static=summary(par_mat_tib$RE.Y_FXSPR_static)[1:6],
                    SSB_FXSPR_static=summary(par_mat_tib$RE.SSB_FXSPR_static)[1:6],
                    SPR0_static=summary(par_mat_tib$RE.SPR0_static)[1:6],
                    SPR_FXSPR_static=summary(par_mat_tib$RE.SPR_FXSPR_static)[1:6],
                    q1=summary(par_mat_tib$RE.q1)[1:6],
                    q2=summary(par_mat_tib$RE.q2)[1:6],
                    NAA_sigma=summary(par_mat_tib$RE.NAA_sigma)[1:6],
                    NAA_rho=summary(par_mat_tib$RE.NAA_rho)[1:6],
                    selpars_f1a=summary(par_mat_tib$RE.selpars_f1a)[1:6],
                    selpars_f1b=summary(par_mat_tib$RE.selpars_f1b)[1:6],
                    selpars_i1a=summary(par_mat_tib$RE.selpars_i1a)[1:6],
                    selpars_i1b=summary(par_mat_tib$RE.selpars_i1b)[1:6],
                    selpars_i2a=summary(par_mat_tib$RE.selpars_i2a)[1:6],
                    selpars_i2b=summary(par_mat_tib$RE.selpars_i2b)[1:6],
                    catch_paa_par=summary(par_mat_tib$RE.catch_paa_par)[1:6],
                    index_paa_par1=summary(par_mat_tib$RE.index_paa_par1)[1:6],
                    index_paa_par2=summary(par_mat_tib$RE.index_paa_par2)[1:6],
                    Ecov_beta=summary(par_mat_tib$RE.Ecov_beta)[1:6],
                    Ecov_process_pars1=summary(par_mat_tib$RE.Ecov_process_pars1)[1:6],
                    Ecov_process_pars2=summary(par_mat_tib$RE.Ecov_process_pars2)[1:6],
                    Ecov_process_pars3=summary(par_mat_tib$RE.Ecov_process_pars3)[1:6]
                    )
RE.par.sum
# strange that some of the selectivity parameters have -Inf or NaN ; checking with tim if he used a different transform
# okayyyyyyy.... need to fix the selpars, given below from tim:     ====
# selpars are transformed like that with k = 1, but the bounds are defined in data$selpars_lower, data$selpars_upper:
#   lower + (upper-lower)/(1+exp(-x))
write.csv(RE.par.sum,file=file.path(here(),'Ecov_study','recruitment_functions','tables', paste0("RelativeError.param.summary.table_grad_", bad.grad.label, "_SE_", bad.se.value, plot.suffix, ".csv")  ), row.names = FALSE )


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

unique(par_mat_tib$mean_rec2[par_mat_tib$EM<3])
unique(true_mat_tib$True.mean_rec1)


# analysis for factors influencing RE ====
par_mat_tib$R_sig <- factor(par_mat_tib$R_sig,labels = c("Rsig_0.1","Rsig_1.0"))
par_mat_tib$Ecov_effect <- factor(par_mat_tib$Ecov_effect,labels=c("Ecov_L","Ecov_H"))
par_mat_tib$Ecov_how    <- factor(par_mat_tib$Ecov_how,labels=c("Ecov_0", "Ecov_1","Ecov_2","Ecov_4"))
par_mat_tib$NAA_cor     <- factor(par_mat_tib$NAA_cor,labels=c("Rcor_L","Rcor_H"))
par_mat_tib$Fhist       <- factor(par_mat_tib$Fhist,labels=c("H-MSY","MSY") ) 
par_mat_tib$Ecov_re_cor <- factor(par_mat_tib$Ecov_re_cor,labels=c("EcovCor_L","EcovCor_H"))
par_mat_tib$obs_error   <- factor(par_mat_tib$obs_error,labels=c("ObsErr_L","ObsErr_H"))

# regression tree for correct form (or SR or Ecov) ======================================
rf_RE.SRa   <- rpart(RE.mean_rec1   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=par_mat_tib[par_mat_tib$EM>2,], control=rpart.control(cp=0.01)) # *** only EM>2 ====
imp.var <- rf_RE.SRa$frame[rf_RE.SRa$frame$var != '<leaf>',]
nodes_RE.SRa <- unique(imp.var[,1])
nodes_RE.SRa
# nothing   

rf_RE.SRb   <- rpart(RE.mean_rec2   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=par_mat_tib[par_mat_tib$EM>2,], control=rpart.control(cp=0.01)) # *** only EM>2 ====
imp.var <- rf_RE.SRb$frame[rf_RE.SRb$frame$var != '<leaf>',]
nodes_RE.SRb <- unique(imp.var[,1])
nodes_RE.SRb
# nothing 

rf_RE.NAA_sigma   <- rpart(RE.NAA_sigma   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=par_mat_tib, control=rpart.control(cp=0.01)) 
imp.var <- rf_RE.NAA_sigma$frame[rf_RE.NAA_sigma$frame$var != '<leaf>',]
nodes_RE.NAA_sigma <- unique(imp.var[,1])
nodes_RE.NAA_sigma
# "NAA_cor"   


rf_RE.NAA_rho   <- rpart(RE.NAA_rho   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=par_mat_tib, control=rpart.control(cp=0.01))
imp.var <- rf_RE.NAA_rho$frame[rf_RE.NAA_rho$frame$var != '<leaf>',]
nodes_RE.NAA_rho <- unique(imp.var[,1])
nodes_RE.NAA_rho
# "Fhist"   "R_sig"   "NAA_cor"   

rf_RE.Ecov_process_pars1   <- rpart(RE.Ecov_process_pars1   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=par_mat_tib, control=rpart.control(cp=0.01))
imp.var <- rf_RE.Ecov_process_pars1$frame[rf_RE.Ecov_process_pars1$frame$var != '<leaf>',]
nodes_RE.Ecov_process_pars1 <- unique(imp.var[,1])
nodes_RE.Ecov_process_pars1
# nothing

rf_RE.Ecov_process_pars2   <- rpart(RE.Ecov_process_pars2   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=par_mat_tib, control=rpart.control(cp=0.01))
imp.var <- rf_RE.Ecov_process_pars2$frame[rf_RE.Ecov_process_pars2$frame$var != '<leaf>',]
nodes_RE.Ecov_process_pars2 <- unique(imp.var[,1])
nodes_RE.Ecov_process_pars2
# "Ecov_re_cor"

rf_RE.Ecov_process_pars3   <- rpart(RE.Ecov_process_pars3   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=par_mat_tib, control=rpart.control(cp=0.01))
imp.var <- rf_RE.Ecov_process_pars3$frame[rf_RE.Ecov_process_pars3$frame$var != '<leaf>',]
nodes_RE.Ecov_process_pars3 <- unique(imp.var[,1])
nodes_RE.Ecov_process_pars3
# nothing


# plot the regression trees that showed some factor sensitivity ====
pdf(file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('reg_tree_NAA_RE', plot.suffix, '.pdf') ),
    height=4,width=7)
par(mfrow=c(1,2), oma=c(5,0,0,0))
prp(rf_RE.NAA_sigma,yesno=FALSE,type=4,clip.right.labs=TRUE)
mtext('a) RE(R_sigma)',adj=0,line=2.5, cex=0.9)
prp(rf_RE.NAA_rho,yesno=FALSE,type=4,clip.right.labs=TRUE)
mtext('b) RE(R_rho)',adj=0,line=2.5, cex=0.9)
# title(sub=paste0(100*round(aic.best.pct.conv,2), "% of runs converged (", aic.best.n.bad, " out of ", aic.best.n, " failed)" ),   adj=0.5, outer=TRUE)

dev.off()

# summarize the RE that we want to look at ======

par_mean_rec1_sum <- par_mat_tib %>%
  left_join(em_tib) %>%
  select(RE.mean_rec1, NAA_cor, obs_error, R_sig, Ecov_re_cor, Fhist, EM, EM_ecov_how, EM_r_mod, Ecov_how, recruit_mod )  %>%
  filter(EM>2) %>%
  mutate(mod.match=ifelse(EM_ecov_how==as.numeric(substr(Ecov_how, 6,6)) & EM_r_mod==recruit_mod, 1, 0)) %>%
  group_by( R_sig,  Fhist,  Ecov_how, EM, mod.match ) %>%
  summarise(mean.RE.mean_rec1=mean(RE.mean_rec1), var.RE.mean_rec1=var(RE.mean_rec1), 
            median.RE.mean_rec1=median(RE.mean_rec1), across(RE.mean_rec1,                                                                                                       list(p2.5=~quantile(.,probs=0.0275 ),
                                                       p97.5=~quantile(.,probs=0.975)))
            )%>%
  ungroup() %>%
  left_join(em_tib)

#  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par_mean_rec2_sum <- par_mat_tib %>%
  left_join(em_tib) %>%
  select(RE.mean_rec2, NAA_cor, obs_error, R_sig, Ecov_re_cor, Fhist, EM, EM_ecov_how, EM_r_mod, Ecov_how, recruit_mod )  %>%
  filter(EM>2) %>%
  mutate(mod.match=ifelse(EM_ecov_how==as.numeric(substr(Ecov_how, 6,6)) & EM_r_mod==recruit_mod, 1, 0)) %>%
  group_by( R_sig,  Fhist,  Ecov_how, EM, mod.match ) %>%
  summarise(mean.RE.mean_rec2=mean(RE.mean_rec2, na.rm=T), var.daic=var(RE.mean_rec2, na.rm=T), ## calculating for all EM even though only has meaning for EM>2
            median.RE.mean_rec2=median(RE.mean_rec2, na.rm=T), across(RE.mean_rec2, 
                                                                      list(p2.5=~quantile(.,probs=0.0275 , na.rm=T),
                                                                           p97.5=~quantile(.,probs=0.975, na.rm=T))) 
  )%>%
  ungroup() %>%
  left_join(em_tib)

#  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par_R_sigma_sum <- par_mat_tib %>%
  left_join(em_tib) %>%
  select(RE.NAA_sigma, NAA_cor, obs_error, R_sig, Ecov_re_cor, Fhist, EM, EM_ecov_how, EM_r_mod, Ecov_how, recruit_mod )  %>%
  mutate(mod.match=ifelse(EM_ecov_how==as.numeric(substr(Ecov_how, 6,6)) & EM_r_mod==recruit_mod, 1, 0)) %>%
  group_by( R_sig,  Fhist, NAA_cor, Ecov_how, EM, mod.match ) %>%
  summarise(mean.RE.NAA_sigma=mean(RE.NAA_sigma), var.RE.NAA_sigma=var(RE.NAA_sigma), 
            median.RE.NAA_sigma=median(RE.NAA_sigma), across(RE.NAA_sigma,                                                                                                       list(p2.5=~quantile(.,probs=0.0275 ),                                                                                p97.5=~quantile(.,probs=0.975))) 
  )%>%
  ungroup() %>%
  left_join(em_tib)

#  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par_R_cor_sum <- par_mat_tib %>%
  left_join(em_tib) %>%
  select(RE.NAA_rho, NAA_cor, obs_error, R_sig, Ecov_re_cor, Fhist, EM, EM_ecov_how, EM_r_mod, Ecov_how, recruit_mod )  %>%
  mutate(mod.match=ifelse(EM_ecov_how==as.numeric(substr(Ecov_how, 6,6)) & EM_r_mod==recruit_mod, 1, 0)) %>%
  group_by( R_sig,  Fhist, NAA_cor, Ecov_how, EM, mod.match ) %>%
  summarise(mean.RE.NAA_rho=mean(RE.NAA_rho), var.RE.NAA_rho=var(RE.NAA_rho), 
            median.RE.NAA_rho=median(RE.NAA_rho), across(RE.NAA_rho,                                                                                                       list(p2.5=~quantile(.,probs=0.0275 ),
                                                       p97.5=~quantile(.,probs=0.975))) 
  )%>%
  ungroup() %>%
  left_join(em_tib)


#  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par_Ecov1_sum <- par_mat_tib %>%
  left_join(em_tib) %>%
  select(RE.Ecov_process_pars1, Ecov_re_cor, obs_error, R_sig, Ecov_re_cor, Fhist, EM, EM_ecov_how, EM_r_mod, Ecov_how, recruit_mod )  %>%
  mutate(mod.match=ifelse(EM_ecov_how==as.numeric(substr(Ecov_how, 6,6)) & EM_r_mod==recruit_mod, 1, 0)) %>%
  group_by(   Fhist, Ecov_re_cor, Ecov_how, EM, mod.match ) %>%
  summarise( mean.RE.Ecov_process_pars1=mean(RE.Ecov_process_pars1), var.RE.Ecov_process_pars1=var(RE.Ecov_process_pars1), 
             median.RE.Ecov_process_pars1=median(RE.Ecov_process_pars1), across(RE.Ecov_process_pars1,                                                                                                       list(p2.5=~quantile(.,probs=0.0275 ) ,
                                                                                    p97.5=~quantile(.,probs=0.975)))
  )%>%
  ungroup() %>%
  left_join(em_tib)


#  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par_Ecov2_sum <- par_mat_tib %>%
  left_join(em_tib) %>%
  select(RE.Ecov_process_pars2, Ecov_re_cor, obs_error, R_sig, Ecov_re_cor, Fhist, EM, EM_ecov_how, EM_r_mod, Ecov_how, recruit_mod )  %>%
  mutate(mod.match=ifelse(EM_ecov_how==as.numeric(substr(Ecov_how, 6,6)) & EM_r_mod==recruit_mod, 1, 0)) %>%
  group_by(   Fhist, Ecov_re_cor, Ecov_how, EM, mod.match ) %>%
  summarise( mean.RE.Ecov_process_pars2=mean(RE.Ecov_process_pars2), var.RE.Ecov_process_pars=var(RE.Ecov_process_pars2), 
             median.RE.Ecov_process_pars2=median(RE.Ecov_process_pars2), across(RE.Ecov_process_pars2,                                                                                                       list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                   p97.5=~quantile(.,probs=0.975)))
  )%>%
  ungroup() %>%
  left_join(em_tib)


#  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par_Ecov3_sum <- par_mat_tib %>%
  left_join(em_tib) %>%
  select(RE.Ecov_process_pars3, Ecov_re_cor, obs_error, R_sig, Ecov_re_cor, Fhist, EM, EM_ecov_how, EM_r_mod, Ecov_how, recruit_mod )  %>%
  mutate(mod.match=ifelse(EM_ecov_how==as.numeric(substr(Ecov_how, 6,6)) & EM_r_mod==recruit_mod, 1, 0)) %>%
  group_by(  Fhist, Ecov_re_cor, Ecov_how, EM, mod.match ) %>%
  summarise(mean.RE.Ecov_process_pars3=mean(RE.Ecov_process_pars3), var.RE.Ecov_process_pars3=var(RE.Ecov_process_pars3), 
            median.RE.Ecov_process_pars3=median(RE.Ecov_process_pars3), across(RE.Ecov_process_pars3,                                                                                                       list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                  p97.5=~quantile(.,probs=0.975)))
  )%>%
  ungroup() %>%
  left_join(em_tib)
#################################################
# RE Plots  ====

RE.mean_rec1.plot <- ggplot(par_mean_rec1_sum, aes(x=EM_mod, y=median.RE.mean_rec1, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist + R_sig   ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.mean_rec1_p2.5), ymax=(RE.mean_rec1_p97.5), width=.5) )+
  ylab('RE(Recr_par_a) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-1,5)) +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in Beverton-Holt parameter a')
ggsave(RE.mean_rec1.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.mean_rec1.plot', plot.suffix, '.png') ),  height=7, width=12)

# same plot, with narrower ylim
RE.mean_rec1.plot.ylim <- ggplot(par_mean_rec1_sum, aes(x=EM_mod, y=median.RE.mean_rec1, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist + R_sig   ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.mean_rec1_p2.5), ymax=(RE.mean_rec1_p97.5), width=.5) )+
  ylab('RE(Recr_par_a) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-1,1)) +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in Beverton-Holt parameter a')
ggsave(RE.mean_rec1.plot.ylim, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.mean_rec1.plot.ylim', plot.suffix, '.png') ),  height=7, width=12)


RE.mean_rec2.plot <- ggplot(par_mean_rec2_sum, aes(x=EM_mod, y=median.RE.mean_rec2, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist + R_sig   ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.mean_rec2_p2.5), ymax=(RE.mean_rec2_p97.5), width=.5) )+
  ylab('RE(Recr_par_b) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-1,5)) +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in Beverton-Holt parameter b')
ggsave(RE.mean_rec2.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.mean_rec2.plot', plot.suffix, '.png') ),  height=7, width=12)

# same plot, with narrower ylim
RE.mean_rec2.plot.ylim <- ggplot(par_mean_rec2_sum, aes(x=EM_mod, y=median.RE.mean_rec2, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist + R_sig   ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.mean_rec2_p2.5), ymax=(RE.mean_rec2_p97.5), width=.5) )+
  ylab('RE(Recr_par_b) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-1,1)) +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in Beverton-Holt parameter b')
ggsave(RE.mean_rec2.plot.ylim, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.mean_rec2.plot.ylim', plot.suffix, '.png') ),  height=7, width=12)


RE.R_sigma.plot <- ggplot(par_R_sigma_sum, aes(x=EM_mod, y=median.RE.NAA_sigma, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist + R_sig + NAA_cor  ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.NAA_sigma_p2.5), ymax=(RE.NAA_sigma_p97.5), width=.5) )+
  ylab('RE(R_sigma) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-1,3)) +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in R_sigma')
ggsave(RE.R_sigma.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.R_sigma.plot', plot.suffix, '.png') ),  height=8, width=14)

# same plot, with narrower ylim
RE.R_sigma.plot.ylim <- ggplot(par_R_sigma_sum, aes(x=EM_mod, y=median.RE.NAA_sigma, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist + R_sig + NAA_cor  ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.NAA_sigma_p2.5), ymax=(RE.NAA_sigma_p97.5), width=.5) )+
  ylab('RE(R_sigma) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-1,1)) +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in R_sigma')
ggsave(RE.R_sigma.plot.ylim, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.R_sigma.plot.ylim', plot.suffix, '.png') ),  height=8, width=14)



RE.R_cor.plot <- ggplot(par_R_cor_sum, aes(x=EM_mod, y=median.RE.NAA_rho, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist + R_sig + NAA_cor  ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.NAA_rho_p2.5), ymax=(RE.NAA_rho_p97.5), width=.5) )+
  ylab('RE(R_sigma) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-1,3)) +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in R_rho')
ggsave(RE.R_cor.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.R_cor.plot', plot.suffix, '.png') ),  height=8, width=14)



# mean ecov is very precisely estimated; probably because  (?):
# unique(df.oms$Ecov_obs_sig)  is a fixed input (not estimated)
# [1] 0.1

RE.Ecov.mean.plot <- ggplot(par_Ecov1_sum, aes(x=EM_mod, y=median.RE.Ecov_process_pars1, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist +  Ecov_re_cor  ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.Ecov_process_pars1_p2.5), ymax=(RE.Ecov_process_pars1_p97.5), width=.5) )+
  ylab('RE(Ecov mean) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-0.25,0.25)) +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in Ecov mean')
ggsave(RE.Ecov.mean.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.Ecov.mean.plot', plot.suffix, '.png') ),  height=8, width=14)




RE.Ecov.sigma.plot <- ggplot(par_Ecov2_sum, aes(x=EM_mod, y=median.RE.Ecov_process_pars2, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist +  Ecov_re_cor  ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.Ecov_process_pars2_p2.5), ymax=(RE.Ecov_process_pars2_p97.5), width=.5) )+
  ylab('RE(Ecov_sigma) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-1,1)) +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in Ecov process sigma')
ggsave(RE.Ecov.sigma.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.Ecov.sigma.plot', plot.suffix, '.png') ),  height=8, width=14)



RE.Ecov.rho.plot <- ggplot(par_Ecov3_sum, aes(x=EM_mod, y=median.RE.Ecov_process_pars3, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist +  Ecov_re_cor  ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.Ecov_process_pars3_p2.5), ymax=(RE.Ecov_process_pars3_p97.5), width=.5) )+
  ylab('RE(Ecov_cor) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-5,5)) +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in Ecov process autocorrelation')
ggsave(RE.Ecov.rho.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.Ecov.rho.plot', plot.suffix, '.png') ),  height=8, width=14)


# same plot, with narrower ylim
RE.Ecov.rho.plot.ylim <- ggplot(par_Ecov3_sum, aes(x=EM_mod, y=median.RE.Ecov_process_pars3, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist +  Ecov_re_cor  ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.Ecov_process_pars3_p2.5), ymax=(RE.Ecov_process_pars3_p97.5), width=.5) )+
  ylab('RE(Ecov_cor) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-1,1)) +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in Ecov process autocorrelation')
ggsave(RE.Ecov.rho.plot.ylim, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.Ecov.rho.plot.ylim', plot.suffix, '.png') ),  height=8, width=14)

