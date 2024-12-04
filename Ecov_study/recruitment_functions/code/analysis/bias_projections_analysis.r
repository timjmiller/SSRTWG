library(here)
library(tidyverse)
library(wham)
library(rpart)
library(rpart.plot)
source(file.path(here::here(),'Ecov_study','recruitment_functions','code','analysis','functions_for_analysis.R'))

df.recr.proj <- as.data.frame(readRDS(file.path(here::here(), 'Ecov_study','recruitment_functions','results','df.recr.proj.RDS')))
df.ssb.proj <- as.data.frame(readRDS(file.path(here::here(), 'Ecov_study','recruitment_functions','results','df.ssb.proj.RDS')))
df.catch.proj <- as.data.frame(readRDS(file.path(here::here(), 'Ecov_study','recruitment_functions','results','df.catch.proj.RDS')))

colnames(df.recr.proj) = colnames(df.ssb.proj) = colnames(df.catch.proj) <- c("OM","EM","Sim","Year",
    "re.cont.ecov",  "re.avg.ecov",  "re.use.ecov",
    "rmse.cont.ecov","rmse.avg.ecov","rmse.use.ecov")

df.oms    <- readRDS(file.path(here::here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS"))
colnames(df.oms)[1] <- 'OM'

for(i in 1:576){
  df.recr.proj$OM[df.recr.proj$OM==i] = 
  df.ssb.proj$OM[df.ssb.proj$OM==i] =
  df.catch.proj$OM[df.catch.proj$OM==i] <- paste0('om_',i)
}

df.recr.proj  <- left_join(df.recr.proj, df.oms, by=c("OM"))
df.ssb.proj   <- left_join(df.ssb.proj, df.oms, by=c("OM"))
df.catch.proj <- left_join(df.catch.proj, df.oms, by=c("OM"))


####################################
## COVARIATES ######################
####################################
df.recr.proj$obs_error=  df.ssb.proj$obs_error=  df.catch.proj$obs_error      <- factor(df.ssb.proj$obs_error,levels=c("L","H"))
df.recr.proj$R_sig=      df.ssb.proj$R_sig=      df.catch.proj$R_sig          <- as.factor(df.ssb.proj$R_sig)
df.recr.proj$Fhist=      df.ssb.proj$Fhist=      df.catch.proj$RMSE_par$Fhist <- factor(df.ssb.proj$Fhist,levels=c("MSY","L-H","H-MSY"))
df.recr.proj$NAA_cor=    df.ssb.proj$NAA_cor=    df.catch.proj$NAA_cor        <- as.factor(df.ssb.proj$NAA_cor) 
df.recr.proj$Ecov_re_cor=df.ssb.proj$Ecov_re_cor=df.catch.proj$Ecov_re_cor    <- as.factor(df.ssb.proj$Ecov_re_cor) 
df.recr.proj$Ecov_effect=df.ssb.proj$Ecov_effect=df.catch.proj$Ecov_effect    <- as.factor(df.ssb.proj$Ecov_effect) 
df.recr.proj$Ecov_how=   df.ssb.proj$Ecov_how=   df.catch.proj$Ecov_how       <- as.factor(df.ssb.proj$Ecov_how) 
df.recr.proj$ssb_cv=     df.ssb.proj$ssb_cv=     df.catch.proj$ssb_cv         <- factor(case_when(df.ssb.proj$ssb_cv < mean(df.ssb.proj$ssb_cv) - sd(df.ssb.proj$ssb_cv) ~ 'L',
                                    df.ssb.proj$ssb_cv > mean(df.ssb.proj$ssb_cv) + sd(df.ssb.proj$ssb_cv) ~ "H",
                                    TRUE ~ 'M')) 
df.recr.proj$ssb_cv=     df.ssb.proj$ssb_cv=      df.catch.proj$ssb_cv        <- factor(df.ssb.proj$ssb_cv,levels=c("L","M","H"))
df.recr.proj$ecov_slope= df.ssb.proj$ecov_slope=  df.catch.proj$ecov_slope    <- as.factor(case_when(df.ssb.proj$ecov_slope > mean(df.ssb.proj$ecov_slope) + sd(df.ssb.proj$ecov_slope) ~ "H",
                                       df.ssb.proj$ecov_slope < mean(df.ssb.proj$ecov_slope) - sd(df.ssb.proj$ecov_slope) ~ 'L',
                                       TRUE ~ 'M'))
df.recr.proj$ecov_slope= df.ssb.proj$ecov_slope= df.catch.proj$ecov_slope     <- factor(df.ssb.proj$ecov_slope,levels=c("L","M","H"))




labels <- c(expression(italic('mean')),
            expression(sigma['obs']~'= H'),
            expression(sigma['r']~'= 0.3'),
            expression(sigma['r']~'= 0.5'),
            expression(sigma['r']~'= 1.0'),
            expression(italic('F')~'= L-H'),
            expression(italic('F')~'= H-MSY'),
            expression(rho['NAA']~' = H'), 
            expression(rho['E'['cov']]~' = H'), 
            expression(beta['E'['cov']]~' = H'), 
            expression(italic('f(E'['cov']*')')~'= 1'), 
            expression(italic('f(E'['cov']*')')~'= 2'), 
            expression(Delta*'E'['cov']*' = M'), 
            expression(Delta*'E'['cov']*' = H'), 
            expression(italic('CV'['SSB']~'= M')),
            expression(italic('CV'['SSB']~'= H')))
