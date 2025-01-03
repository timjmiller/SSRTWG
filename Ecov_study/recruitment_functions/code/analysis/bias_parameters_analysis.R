library(here)
library(tidyverse)
library(rpart)
library(rpart.plot)
source(file.path(here::here(),'Ecov_study','recruitment_functions','code','analysis','functions_for_analysis.R'))

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
for(i in c(6:10,13:ncol(RE_par))){
  hist(RE_par[,i],xlab='',ylab='',main='')
  mtext(colnames(RE_par)[i])
}
dev.off()


pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','RMSE_par_hist.pdf'),height=12,width=12)
par(mfrow=c(5,5),mar=c(2,2,2,2),oma=c(2,2,2,2))
for(i in c(6:10,13:ncol(RMSE_par))){
  hist(RMSE_par[,i],xlab='',ylab='',main='')
  mtext(colnames(RMSE_par)[i])
}
dev.off()


#################################
for(i in 1:576){
  RMSE_par$OM[RMSE_par$OM==i] = RE_par$OM[RE_par$OM==i] <- paste0('om_',i)
}

RE_par   <- left_join(RE_par,df.oms,by=c('OM'))
RE_par$Ecov_beta[!is.finite(RE_par$Ecov_beta) | RE_par$Ecov_beta=='NaN'] <- NA

RMSE_par <- left_join(RMSE_par,df.oms,by=c('OM'))

##--COVARIATES--########
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
RE_par$ecov_slope= RMSE_par$ecov_slope  <- as.factor(case_when(RMSE_par$ecov_slope > mean(RMSE_par$ecov_slope) + sd(RMSE_par$ecov_slope) ~ "+",
                                       RMSE_par$ecov_slope < mean(RMSE_par$ecov_slope) - sd(RMSE_par$ecov_slope) ~ '-',
                                       TRUE ~ '0'))
RE_par$ecov_slope= RMSE_par$ecov_slope  <- factor(RMSE_par$ecov_slope,levels=c("-","0","+"))

vars <- c("obs_error","R_sig","Fhist","NAA_cor","Ecov_re_cor","Ecov_effect","Ecov_how","ecov_slope","ssb_cv")

labels <- c(expression(sigma['obs']~'= L'),
            expression(sigma['obs']~'= H'),
            expression(sigma['r']~'= 0.1'),
            expression(sigma['r']~'= 0.3'),
            expression(sigma['r']~'= 0.5'),
            expression(sigma['r']~'= 1.0'),
            expression(italic('F')~'= MSY'),
            expression(italic('F')~'= L-H'),
            expression(italic('F')~'= H-MSY'),
            expression(rho['NAA']~' = L'), 
            expression(rho['NAA']~' = H'), 
            expression(rho['E'['cov']]~' = L'), 
            expression(rho['E'['cov']]~' = H'), 
            expression(beta['E'['cov']]~' = L'), 
            expression(beta['E'['cov']]~' = H'), 
            expression(italic('f(E'['cov']*')')~'= 0'), 
            expression(italic('f(E'['cov']*')')~'= 1'), 
            expression(italic('f(E'['cov']*')')~'= 2'), 
            expression(Delta*'E'['cov']*' = -'), 
            expression(Delta*'E'['cov']*' = 0'), 
            expression(Delta*'E'['cov']*' = +'), 
            expression(italic('CV'['SSB']~'= L')),
            expression(italic('CV'['SSB']~'= M')),
            expression(italic('CV'['SSB']~'= H'))
)

pdf(file=file.path(here::here(), 'Ecov_study','recruitment_functions','plots','par_marg.pdf'),height=6,width=4.5)
par(mfrow=c(3,1),mar=c(1,2,1,0),oma=c(6,4,2,2),cex.axis=0.6,cex.lab=0.6)
ylims <- c(0,1)
dd(RMSE_par,vars=vars,labels=rep(NA,length(labels)),yvar="Ecov_beta",ylims=ylims)
  mtext(expression("RE"),side=2,line=2.5)
abline(h=0,lty=2)
#mtext('a) Recruitment',adj=0)


dev.off()



