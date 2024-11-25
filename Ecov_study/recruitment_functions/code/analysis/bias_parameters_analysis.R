library(here)
library(tidyverse)

res.path    <- file.path(here::here(),"Ecov_study", "recruitment_functions", "results")  # directory where simulation 
res.dir     <- 'results'   # 'results'     'results_beta_fix'   # results folder where AIC dataframes are
plot.dir    <- 'plots'    # 'plots_lizruns'  'plots_beta_fix'  
table.dir   <- 'tables'   # 'results'     'results_beta_fix'
plot.suffix <- ''      # '_beta_fix'   '' 

RE_par   <- as.data.frame(readRDS(file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0("RE_par.RDS") ) ))
RMSE_par <- as.data.frame(readRDS(file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0("RMSE_par.RDS") ) ))

df.oms    <- readRDS(file.path(here::here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS"))
colnames(df.oms)[1] <- 'OM'

##################################

pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','RE_par_hist.pdf'),height=12,width=12)
par(mfrow=c(5,5),mar=c(2,2,2,2),oma=c(2,2,2,2))
for(i in 4:ncol(RE_par)){
  hist(RE_par[,i],xlab='',ylab='',main='')
  mtext(colnames(RE_par)[i])
}
dev.off()


pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','RMSE_par_hist.pdf'),height=12,width=12)
par(mfrow=c(5,5),mar=c(2,2,2,2),oma=c(2,2,2,2))
for(i in 4:ncol(RMSE_par)){
  hist(RMSE_par[,i],xlab='',ylab='',main='')
  mtext(colnames(RMSE_par)[i])
}
dev.off()


#################################
for(i in 1:576){
  RMSE_par$OM[RMSE_par$OM==i] = RE_par$OM[RE_par$OM==i] <- paste0('om_',i)
}

RE_par   <- left_join(RE_par,df.oms,by=c('OM'))
RMSE_par <- left_join(RMSE_par,df.oms,by=c('OM'))


RE_par$obs_error=  RMSE_par$obs_error   <- factor(RE_par$obs_error,levels=c("L","H"))
RE_par$R_sig=      RMSE_par$R_sig       <- as.factor(RE_par$R_sig)
RE_par$Fhist=      RMSE_par$Fhist       <- factor(RE_par$Fhist,levels=c("MSY","L-H","H-MSY"))
RE_par$NAA_cor=    RMSE_par$NAA_cor     <- as.factor(RE_par$NAA_cor) 
RE_par$Ecov_re_cor=RMSE_par$Ecov_re_cor <- as.factor(RE_par$Ecov_re_cor) 
RE_par$Ecov_effect=RMSE_par$Ecov_effect <- as.factor(RE_par$Ecov_effect) 
RE_par$Ecov_how=   RMSE_par$Ecov_how    <- as.factor(RE_par$Ecov_how) 
RE_par$ssb_cv=     RMSE_par$ssb_cv      <- factor(case_when(RE_par$ssb_cv < mean(RE_par$ssb_cv) - sd(RE_par$ssb_cv) ~ 'L',
                                    RE_par$ssb_cv > mean(RE_par$ssb_cv) + sd(RE_par$ssb_cv) ~ "H",
                                    TRUE ~ 'M')) 
RE_par$ssb_cv=     RMSE_par$ssb_cv      <- factor(RE_par$ssb_cv,levels=c("L","M","H"))
RE_par$ecov_slope= RMSE_par$ecov_slope  <- as.factor(case_when(RE_par$ecov_slope > mean(RE_par$ecov_slope) + sd(RE_par$ecov_slope) ~ "H",
                                       RE_par$ecov_slope < mean(RE_par$ecov_slope) - sd(RE_par$ecov_slope) ~ 'L',
                                       TRUE ~ 'M'))
RE_par$ecov_slope= RMSE_par$ecov_slope  <- factor(RE_par$ecov_slope,levels=c("L","M","H"))





##--FITS--########################################


rf_re_NAA_rho   <- rpart(NAA_rho ~ obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                          data=RE_par, control=rpart.control(cp=0.001))





lm_re_NAA_rho <-   lm(NAA_rho ~ obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                          data=RE_par)
