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

##--F relative error--########
re.fbar <- lapply(1:n_oms,function(om){  #produces list of length n_oms with list of length n_ems, each with matrix nyearsXn_sims
  lapply(1:n_ems, function(em){
    sapply(1:n_sims, function(sim){
      print(paste0("om ",om," em ",em," sim ",sim," re.fbar"))
      dat <- try(readRDS(file.path(res.path, paste0("om", om, '/','sim',sim,'_','em',em,'.RDS') ) ) )#,
      # silent=TRUE)   )
      if(class(dat)!='try-error'){
        (dat$fit$rep$Fbar - dat$truth$Fbar)/ dat$truth$Fbar   
      }
    })
  })
})
saveRDS(re.fbar, file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0('re.fbar', plot.suffix, '.RDS') ))

##--F rms error--########
rmse.fbar <- lapply(1:n_oms,function(om){  #produces list of length n_oms with list of length n_ems, each with matrix nyearsXn_sims
  lapply(1:n_ems, function(em){
    sapply(1:n_sims, function(sim){
      print(paste0("om ",om," em ",em," sim ",sim," rmse.fbar"))
      dat <- try(readRDS(file.path(res.path, paste0("om", om, '/','sim',sim,'_','em',em,'.RDS') ) ) )#, silent=TRUE)   )
      if(class(dat)!='try-error'){
        sqrt( (dat$fit$rep$Fbar - dat$truth$Fbar)^2)   
      }
    })
  })
})
saveRDS(rmse.fbar, file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0('re.fbar', plot.suffix, '.RDS') ))

nyears<-dim(re.fbar[[1]][[1]])[1]

fbar.df           <- matrix(NA, nrow=n_oms*n_ems*nyears*n_sims, ncol=(6 )   )# OM, EM, Sim, Year, RE, df.om colnames
colnames(fbar.df) <- c('OM', 'EM', 'Sim', 'Year', 'RE','RMSE')

re.fbar.sims   <- matrix(NA, nrow=n_sims*nyears, ncol=1)
rmse.fbar.sims <- matrix(NA, nrow=n_sims*nyears, ncol=1)

k <- 1
kdf <- 1
for(iom in 1:n_oms) {
  for (jem in 1:n_ems){
    kk <- 1
    for (ksim in 1:n_sims) {
      print(paste0("om ",iom," em ",jem," sim ",ksim))
      dat <- try(readRDS(file.path(res.path, paste0("om", iom, '/','sim',ksim,'_','em',jem,'.RDS') ) ) )
      if(class(dat)!='try-error'){
        if(length(re.fbar[[iom]][[jem]])==n_sims*nyears)    re.fbar.sims[((ksim-1)*nyears+1):(ksim*nyears),1]   <- re.fbar[[iom]][[jem]][,ksim]
        if(length(rmse.fbar[[iom]][[jem]])==n_sims*nyears)  rmse.fbar.sims[((ksim-1)*nyears+1):(ksim*nyears),1] <- rmse.fbar[[iom]][[jem]][,ksim]

        #-F-##############
        if(length(re.fbar[[iom]][[jem]])==n_sims) {
          if(length(re.fbar[[iom]][[jem]][[ksim]])>0) {
            re.fbar.sims[((ksim-1)*nyears+1):(ksim*nyears),1] <- re.fbar[[iom]][[jem]][[ksim]]
          }}
        if(length(rmse.fbar[[iom]][[jem]])==n_sims) {
          if(length(rmse.fbar[[iom]][[jem]][[ksim]])>0) {
            rmse.fbar.sims[((ksim-1)*nyears+1):(ksim*nyears),1] <- rmse.fbar[[iom]][[jem]][[ksim]]
          }}
        
      }
      kk  <- kk+1
      kdf <- kdf+1
      
    } 
    fbar.df[((k-1)*nyears*n_sims+1):(n_sims*nyears*k),1:4] <- cbind(rep(iom, nyears*n_sims), rep(jem, nyears*n_sims), rep(seq(1,n_sims),each=nyears), rep(seq(1,nyears), n_sims))
    fbar.df[((k-1)*nyears*n_sims+1):(n_sims*nyears*k),5 ]  <- re.fbar.sims
    fbar.df[((k-1)*nyears*n_sims+1):(n_sims*nyears*k),6 ]  <- rmse.fbar.sims
    
    k <- k+1
    
    re.fbar.sims   <- matrix(NA, nrow=n_sims*nyears, ncol=1)
    rmse.fbar.sims <- matrix(NA, nrow=n_sims*nyears, ncol=1)
  } #jem loop
} #iom loop

##-SAVE--#########################
saveRDS(fbar.df, file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "error.fbar", plot.suffix, ".df.RDS") ) )
