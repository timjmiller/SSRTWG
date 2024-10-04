library(here)
library(tidyverse)
library(rpart)

res.path    <- file.path(here::here(),"Ecov_study", "recruitment_functions", "results")  # directory where simulation 
res.dir     <- 'results'   # 'results'     'results_beta_fix'   # results folder where AIC dataframes are
plot.dir    <- 'plots'    # 'plots_lizruns'  'plots_beta_fix'  
table.dir   <- 'tables'   # 'results'     'results_beta_fix'
plot.suffix <- ''      # '_beta_fix'   '' 

## specify bad.grad.label and bad.se.value (these are the thresholds set in convergence_summaries.R to determine convergence) 
bad.grad.value <- 1E-6 #tim used 1E-6 (the abs of the exponent will be used for output suffix; ex: filename_grad_6.png)
bad.grad.label <- as.numeric(strsplit(as.character(bad.grad.value), split="-")[[1]][2])
bad.se.value   <- 100 #tim used 100 (this value will be used for output suffix; ex: filename_se_100.png)

df.oms    <- readRDS(file.path(here::here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS"))
df.ems    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))
conv.runs <- readRDS(file.path(here::here(),'Ecov_study','recruitment_functions',res.dir, paste0("conv.runs_grad_", bad.grad.label, "_SE_", bad.se.value, plot.suffix, ".RDS")  ) )

n_oms <- nrow(df.oms)
n_ems <- nrow(df.ems)
n_sims <- 100
nyears <- 40

##--ssb relative error--########
#re.ssb <- lapply(1:n_oms,function(om){  #produces list of length n_oms with list of length n_ems, each with matrix nyearsXn_sims ***except sometimes the matrix is a list instead of a matrix!
#  lapply(1:n_ems, function(em){
#    sapply(1:n_sims, function(sim){
#      print(paste0("om ",om," em ",em," sim ",sim," re.ssb"))
#      dat <- try(readRDS(file.path(res.path, paste0("om", om, '/','sim',sim,'_','em',em,'.RDS') ) ) )
#      if(class(dat)!='try-error'){
#        (dat$fit$rep$SSB - dat$truth$SSB)/ dat$truth$SSB   
#      }
#    })
#  })
#})
#saveRDS(re.ssb, file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0('re.ssb', plot.suffix, '.RDS') ) )
re.ssb <- readRDS(file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0('re.ssb', plot.suffix, '.RDS') ) )

##--ssb rms error--########
#rmse.ssb <- lapply(1:n_oms,function(om){  #produces list of length n_oms with list of length n_ems, each with matrix nyearsXn_sims ***except sometimes the matrix is a list instead of a matrix!
#  lapply(1:n_ems, function(em){
#    sapply(1:n_sims, function(sim){
#      print(paste0("om ",om," em ",em," sim ",sim," rmse.ssb"))
#      dat <- try(readRDS(file.path(res.path, paste0("om", om, '/','sim',sim,'_','em',em,'.RDS') ) ) )
#      if(class(dat)!='try-error'){
#        sqrt( (dat$fit$rep$SSB - dat$truth$SSB)^2)   
#      }
#    })
#  })
#})
#saveRDS(rmse.ssb, file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0('rmse.ssb', plot.suffix, '.RDS') ) )
rmse.ssb <- readRDS(file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0('rmse.ssb', plot.suffix, '.RDS') ) )

nyears<-dim(re.ssb[[1]][[1]])[1]

ssb.df            <- as.data.frame(matrix(NA, nrow=n_oms*n_ems*nyears*n_sims, ncol=(13 )   ))# OM, EM, Sim, Year, RE, df.om colnames
colnames(ssb.df)  <- c('OM', 'EM', 'Sim', 'Year', 'RE','RMSE','ssb_cv','ecov_slope','opt','conv','sdrep','max_grad','SE_par_max')

re.ssb.sims    <- matrix(NA, nrow=n_sims*nyears, ncol=1)
rmse.ssb.sims  <- matrix(NA, nrow=n_sims*nyears, ncol=1)

  k <- 1
  kdf <- 1
  for(iom in 1:n_oms) {
    for (jem in 1:n_ems){
      kk <- 1
      ssb_cvs=ecov_slopes <- numeric(n_sims)
      conv <- matrix(NA,ncol=5,nrow=n_sims)
    for (ksim in 1:n_sims) {
        print(paste0("om ",iom," em ",jem," sim ",ksim))
        dat <- try(readRDS(file.path(res.path, paste0("om", iom, '/','sim',ksim,'_','em',jem,'.RDS') ) ) )
        conv[ksim,] <- as.numeric(convergence_fn(dat)[c(1,2,3,4,5)])
        if(class(dat)!='try-error'){
          if(length(re.ssb[[iom]][[jem]])==n_sims*nyears)     re.ssb.sims[((ksim-1)*nyears+1):(ksim*nyears),1]    <- re.ssb[[iom]][[jem]][,ksim]
          if(length(rmse.ssb[[iom]][[jem]])==n_sims*nyears)   rmse.ssb.sims[((ksim-1)*nyears+1):(ksim*nyears),1]  <- rmse.ssb[[iom]][[jem]][,ksim]

          #-ssb-############
          if(length(re.ssb[[iom]][[jem]])==n_sims) {
            if(length(re.ssb[[iom]][[jem]][[ksim]])>0) {
              re.ssb.sims[((ksim-1)*nyears+1):(ksim*nyears),1] <- re.ssb[[iom]][[jem]][[ksim]]
            }}
          if(length(rmse.ssb[[iom]][[jem]])==n_sims) {
            if(length(rmse.ssb[[iom]][[jem]][[ksim]])>0) {
              rmse.ssb.sims[((ksim-1)*nyears+1):(ksim*nyears),1] <- rmse.ssb[[iom]][[jem]][[ksim]]
            }}
          ssb_cvs[ksim]     <- sd(dat$truth$SSB)/mean(dat$truth$SSB)
          ecov_slopes[ksim] <- summary(lm(dat$truth$Ecov_x ~ seq(1,40)))$coefficients[2,1]
        }
        kk  <- kk+1
        kdf <- kdf+1
        
      } 
      ssb.df[((k-1)*nyears*n_sims+1):(n_sims*nyears*k),1:4] <- cbind(rep(iom, nyears*n_sims), rep(jem, nyears*n_sims), rep(seq(1,n_sims),each=nyears), rep(seq(1,nyears), n_sims))
      ssb.df[((k-1)*nyears*n_sims+1):(n_sims*nyears*k),5 ]  <- re.ssb.sims
      ssb.df[((k-1)*nyears*n_sims+1):(n_sims*nyears*k),6 ]  <- rmse.ssb.sims
      ssb.df[((k-1)*nyears*n_sims+1):(n_sims*nyears*k),7 ]  <- rep(ssb_cvs,each=nyears)
      ssb.df[((k-1)*nyears*n_sims+1):(n_sims*nyears*k),8 ]  <- rep(ecov_slopes,each=nyears)
      ssb.df[((k-1)*nyears*n_sims+1):(n_sims*nyears*k),9:13 ]  <- as.data.frame(conv[rep(1:100,each=40),])

      k <- k+1
      
      re.ssb.sims    <- matrix(NA, nrow=n_sims*nyears, ncol=1)
      rmse.ssb.sims  <- matrix(NA, nrow=n_sims*nyears, ncol=1)
    } #jem loop
  } #iom loop
  
  df.oms$OM <- 1:576
  df.ems$EM <- 1:5

ssb.df <- left_join(x=ssb.df,y=df.oms,by='OM') %>%
  left_join(x=., y=df.ems,by="EM") %>%
  mutate(obs_error=factor(obs_error,levels=c("L","H")),
         R_sig    =as.factor(R_sig),
         Fhist    =factor(Fhist,levels=c("MSY","L-H","H-MSY")),
         NAA_cor  =as.factor(NAA_cor), 
         Ecov_re_cor = as.factor(Ecov_re_cor), 
         Ecov_effect = as.factor(Ecov_effect), 
         Ecov_how    = as.factor(Ecov_how), 
         ssb_cv      = factor(case_when(ssb_cv < mean(ssb_cv) - sd(ssb_cv) ~ 'L',
                                        ssb_cv > mean(ssb_cv) + sd(ssb_cv) ~ "H",
                                        TRUE ~ 'M')), 
         ecov_slope  = factor(case_when(ecov_slope > mean(ecov_slope) + sd(ecov_slope) ~ "H",
                                        ecov_slope < mean(ecov_slope) - sd(ecov_slope) ~ 'L',
                                        TRUE ~ 'M'))) %>%
  mutate(ssb_cv     = factor(ssb_cv,levels=c("L","M","H")),
         ecov_slope = factor(ecov_slope,levels=c("L","M","H")))


colnames(ssb.df)[19:20] <- c("EM_ecov_how","EM_recruit_mod")

##-SAVE--#########################
saveRDS(ssb.df,  file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "error.ssb",  plot.suffix, ".df.RDS") ) )

