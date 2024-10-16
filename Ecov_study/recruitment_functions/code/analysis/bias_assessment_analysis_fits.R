library(here)
library(tidyverse)
library(rpart)
library(rpart.plot)
library(nlme)

res.path    <- file.path(here::here(),"Ecov_study", "recruitment_functions", "results")  # directory where simulation 
res.dir     <- 'results'   # 'results'     'results_beta_fix'   # results folder where AIC dataframes are
plot.dir    <- 'plots'    # 'plots_lizruns'  'plots_beta_fix'  
table.dir   <- 'tables'   # 'results'     'results_beta_fix'
plot.suffix <- ''      # '_beta_fix'   '' 

df.oms    <- readRDS(file.path(here::here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS"))
df.ems    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))

bad.grad.value <- 1E-6 #tim used 1E-6 (the abs of the exponent will be used for output suffix; ex: filename_grad_6.png)
bad.grad.label <- as.numeric(strsplit(as.character(bad.grad.value), split="-")[[1]][2])
bad.se.value   <- 100 #tim used 100 (this value will be used for output suffix; ex: filename_se_100.png)

recr.df <- as.data.frame(readRDS( file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "error.recr", plot.suffix, ".df.RDS") ) )) %>%
  filter(opt==1 & conv==0 & sdrep==0 & max_grad<bad.grad.value)

ssb.df  <- as.data.frame(readRDS( file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "error.ssb",  plot.suffix, ".df.RDS") )  )) %>%
  filter(opt==1 & conv==0 & sdrep==0 & max_grad<bad.grad.value)

fbar.df <- as.data.frame(readRDS( file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "error.fbar", plot.suffix, ".df.RDS") )  )) %>%
  filter(opt==1 & conv==0 & sdrep==0 & max_grad<bad.grad.value)


####################################################################################
#  regression trees to find nodes for bias  ====
####################################################################################
# recr ====
print('recr rf')
rf_recr_re_all   <- rpart(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=recr.df, control=rpart.control(cp=0.001))
rf_recr_re_ten   <- rpart(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=recr.df[recr.df$Year %in% 31:40,], control=rpart.control(cp=0.01))
rf_recr_re_last  <- rpart(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=recr.df[recr.df$Year==40,], control=rpart.control(cp=0.01))
rf_recr_rmse_all <- rpart(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how,
                        data=recr.df, control=rpart.control(cp=0.001))
rf_recr_rmse_ten  <- rpart(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how,
                          data=recr.df[recr.df$Year %in% 31:40,], control=rpart.control(cp=0.01))
rf_recr_rmse_last <- rpart(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=recr.df[recr.df$Year==40,], control=rpart.control(cp=0.01))

# ssb ====
print('ssb rf')
rf_ssb_re_all   <- rpart(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                        data=ssb.df, control=rpart.control(cp=0.01))
rf_ssb_re_ten   <- rpart(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=ssb.df[ssb.df$Year %in% 31:40,], control=rpart.control(cp=0.01))
rf_ssb_re_last <- rpart(RE   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how, 
                          data=ssb.df[ssb.df$Year==40,], control=rpart.control(cp=0.01))
rf_ssb_rmse_all <- rpart(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                        data=ssb.df, control=rpart.control(cp=0.01))
rf_ssb_rmse_ten <- rpart(RMSE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how, 
                          data=ssb.df[ssb.df$Year %in% 31:40,], control=rpart.control(cp=0.01))
rf_ssb_rmse_last <- rpart(RMSE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how, 
                          data=ssb.df[ssb.df$Year==40,], control=rpart.control(cp=0.01))

# fbar ====
print('fbar rf')
rf_fbar_re_all   <- rpart(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                        data=fbar.df, control=rpart.control(cp=0.01))
rf_fbar_re_ten   <- rpart(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=fbar.df[fbar.df$Year %in% 31:40,], control=rpart.control(cp=0.01))
rf_fbar_re_last <- rpart(RE   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how, 
                          data=fbar.df[fbar.df$Year==40,], control=rpart.control(cp=0.01))
rf_fbar_rmse_all <- rpart(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                        data=fbar.df, control=rpart.control(cp=0.01))
rf_fbar_rmse_ten <- rpart(RMSE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how, 
                          data=fbar.df[fbar.df$Year %in% 31:40,], control=rpart.control(cp=0.01))
rf_fbar_rmse_last <- rpart(RMSE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how, 
                          data=fbar.df[fbar.df$Year==40,], control=rpart.control(cp=0.01))


FITS <- list(
             rf_recr_re_all=rf_recr_re_all,
             rf_recr_re_ten=rf_recr_re_ten,
             rf_recr_re_last=rf_recr_re_last,
             rf_recr_rmse_all=rf_recr_rmse_all,
             rf_recr_rmse_ten=rf_recr_rmse_ten,
             rf_recr_rmse_last=rf_recr_rmse_last,

             rf_ssb_re_all=rf_ssb_re_all,
             rf_ssb_re_ten=rf_ssb_re_ten,
             rf_ssb_re_last=rf_ssb_re_last,
             rf_ssb_rmse_all=rf_ssb_rmse_all,
             rf_ssb_rmse_ten=rf_ssb_rmse_ten,
             rf_ssb_rmse_last=rf_ssb_rmse_last,

             rf_fbar_re_all=rf_fbar_re_all,
             rf_fbar_re_ten=rf_fbar_re_ten,
             rf_fbar_re_last=rf_fbar_re_last,
             rf_fbar_rmse_all=rf_fbar_rmse_all,
             rf_fbar_rmse_ten=rf_fbar_rmse_ten,
             rf_fbar_rmse_last=rf_fbar_rmse_last)

print('saving random forests')
saveRDS(FITS,file=file.path(here::here(),'Ecov_study','recruitment_functions','results','FITS_error.rds'))
print('saved random forests')


##########################################
## linear models #########################
##########################################

#recr re
print('fitting recr lm lme')
lm_recr_re_all <- lm(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=recr.df)
lm_recr_re_ten <- lm(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=recr.df[recr.df$Year %in% 31:40,])
lm_recr_re_last <- lm(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=recr.df[recr.df$Year==40,])
lme_recr_re_all  <- try(lme(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          random = ~ 1|Sim, data=recr.df))                          
lme_recr_re_ten  <- try(lme(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          random = ~ 1|Sim, data=recr.df[recr.df$Year %in% 31:40,]))                          
lme_recr_re_last <- try(lme(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          random = ~ 1|Sim, data=recr.df[recr.df$Year==40,]))                          
#recr rmse
lm_recr_rmse_all <- lm(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=recr.df)
lm_recr_rmse_ten <- lm(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=recr.df[recr.df$Year %in% 31:40,])
lm_recr_rmse_last <- lm(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=recr.df[recr.df$Year==40,])
lme_recr_rmse_all  <- try(lme(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          random = ~ 1|Sim, data=recr.df))                          
lme_recr_rmse_ten  <- try(lme(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          random = ~ 1|Sim, data=recr.df[recr.df$Year %in% 31:40,]))                          
lme_recr_rmse_last <- try(lme(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          random = ~ 1|Sim, data=recr.df[recr.df$Year==40,]))                          


#ssb re
print('fitting ssb lm lme')
lm_ssb_re_all <- lm(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=recr.df)
lm_ssb_re_ten <- lm(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=recr.df[recr.df$Year %in% 31:40,])
lm_ssb_re_last <- lm(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=recr.df[recr.df$Year==40,])
lme_ssb_re_all  <- try(lme(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          random = ~ 1|Sim, data=recr.df))                          
lme_ssb_re_ten  <- try(lme(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          random = ~ 1|Sim, data=recr.df[recr.df$Year %in% 31:40,]))                          
lme_recr_re_last <- try(lme(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          random = ~ 1|Sim, data=recr.df[recr.df$Year==40,]))                          
#ssb rmse
lm_ssb_rmse_all <- lm(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=ssb.df)
lm_ssb_rmse_ten <- lm(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=ssb.df[ssb.df$Year %in% 31:40,])
lm_ssb_rmse_last <- lm(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=ssb.df[ssb.df$Year==40,])
lme_ssb_rmse_all  <- try(lme(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          random = ~ 1|Sim, data=ssb.df))                          
lme_ssb_rmse_ten  <- try(lme(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          random = ~ 1|Sim, data=ssb.df[ssb.df$Year %in% 31:40,]))                          
lme_ssb_rmse_last <- try(lme(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          random = ~ 1|Sim, data=ssb.df[ssb.df$Year==40,]))                          


#fbar re
print('fitting fbar lm lme')
lm_fbar_re_all <- lm(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=recr.df)
lm_fbar_re_ten <- lm(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=recr.df[recr.df$Year %in% 31:40,])
lm_fbar_re_last <- lm(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=recr.df[recr.df$Year==40,])
lme_fbar_re_all  <- lme(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          random = ~ 1|Sim, data=recr.df)                          
lme_fbar_re_ten  <- lme(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          random = ~ 1|Sim, data=recr.df[recr.df$Year %in% 31:40,])                          
lme_fbar_re_last <- lme(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          random = ~ 1|Sim, data=recr.df[recr.df$Year==40,])                          
#ssb rmse
lm_fbar_rmse_all <- lm(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=fbar.df)
lm_fbar_rmse_ten <- lm(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=fbar.df[fbar.df$Year %in% 31:40,])
lm_fbar_rmse_last <- lm(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          data=fbar.df[fbar.df$Year==40,])
lme_fbar_rmse_all  <- try(lme(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          random = ~ 1|Sim, data=fbar.df))                          
lme_fbar_rmse_ten  <- try(lme(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          random = ~ 1|Sim, data=fbar.df[fbar.df$Year %in% 31:40,]))                          
lme_fbar_rmse_last <- try(lme(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM_ecov_how, 
                          random = ~ 1|Sim, data=fbar.df[fbar.df$Year==40,]))                          


FITS <- c(FITS,
        list(lm_recr_re_all=lm_recr_re_all,
             lm_recr_re_ten=lm_recr_re_ten,
             lm_recr_re_last=lm_recr_re_last,
             lm_recr_rmse_all=lm_recr_rmse_all,
             lm_recr_rmse_ten=lm_recr_rmse_ten,
             lm_recr_rmse_last=lm_recr_rmse_last,
             lme_recr_re_all=lme_recr_re_all,
             lme_recr_re_ten=lme_recr_re_ten,
             lme_recr_re_last=lme_recr_re_last,
             lme_recr_rmse_all=lme_recr_rmse_all,
             lme_recr_rmse_ten=lme_recr_rmse_ten,
             lme_recr_rmse_last=lme_recr_rmse_last,

             lm_ssb_re_all=lm_ssb_re_all,
             lm_ssb_re_ten=lm_ssb_re_ten,
             lm_ssb_re_last=lm_ssb_re_last,
             lm_ssb_rmse_all=lm_ssb_rmse_all,
             lm_ssb_rmse_ten=lm_ssb_rmse_ten,
             lm_ssb_rmse_last=lm_ssb_rmse_last,
             lme_ssb_re_all=lme_ssb_re_all,
             lme_ssb_re_ten=lme_ssb_re_ten,
             lme_ssb_re_last=lme_ssb_re_last,
             lme_ssb_rmse_all=lme_ssb_rmse_all,
             lme_ssb_rmse_ten=lme_ssb_rmse_ten,
             lme_ssb_rmse_last=lme_ssb_rmse_last,

             lm_fbar_re_all=lm_fbar_re_all,
             lm_fbar_re_ten=lm_fbar_re_ten,
             lm_fbar_re_last=lm_fbar_re_last,
             lm_fbar_rmse_all=lm_fbar_rmse_all,
             lm_fbar_rmse_ten=lm_fbar_rmse_ten,
             lm_fbar_rmse_last=lm_fbar_rmse_last,
             lme_fbar_re_all=lme_fbar_re_all,
             lme_fbar_re_ten=lme_fbar_re_ten,
             lme_fbar_re_last=lme_fbar_re_last,
             lme_fbar_rmse_all=lme_fbar_rmse_all,
             lme_fbar_rmse_ten=lme_fbar_rmse_ten,
             lme_fbar_rmse_last=lme_fbar_rmse_last)

print('saving all fits')
saveRDS(FITS,file=file.path(here::here(),'Ecov_study','recruitment_functions','results','FITS_error.rds'))
print('saved all fits')
print('end')

