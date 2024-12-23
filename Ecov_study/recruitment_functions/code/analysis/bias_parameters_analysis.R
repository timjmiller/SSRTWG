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
##--FITS--########################################
#pars: mean_rec1, mean_rec2, Y_FXSPR_static, SSB_FXSPR_static, NAA_rho, Ecov_beta, 

##--RANDOM FORESTS--####
rf_re_mean_rec1   <- rpart(mean_rec1 ~ obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                          data=RE_par, control=rpart.control(cp=0.001))
rf_re_mean_rec2   <- rpart(mean_rec2 ~ obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                          data=RE_par, control=rpart.control(cp=0.001))
rf_re_Y_FXSPR_static   <- rpart(Y_FXSPR_static ~ obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                          data=RE_par, control=rpart.control(cp=0.001))
rf_re_SSB_FXSPR_static <- rpart(SSB_FXSPR_static ~ obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                          data=RE_par, control=rpart.control(cp=0.001))
rf_re_NAA_rho        <- rpart(NAA_rho ~ obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                          data=RE_par, control=rpart.control(cp=0.001))
rf_re_Ecov_beta        <- rpart(Ecov_beta ~ obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                          data=RE_par, control=rpart.control(cp=0.001))

pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','re_trees_pars.pdf'),height=10,width=10)
par(mfrow=c(2,3))
prp(rf_re_mean_rec1,yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('a)',adj=0, line=2, cex=0.9)
prp(rf_re_mean_rec2,yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('a)',adj=0, line=2, cex=0.9)
prp(rf_re_Y_FXSPR_static,yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('a)',adj=0, line=2, cex=0.9)
prp(rf_re_SSB_FXSPR_static,yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('a)',adj=0, line=2, cex=0.9)
prp(rf_re_NAA_rho,yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('a)',adj=0, line=2, cex=0.9)
prp(rf_re_Ecov_beta,yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('a)',adj=0, line=2, cex=0.9)
dev.off()

##--LINEAR MODELS--#######

lm_re_mean_rec1        <- lm(mean_rec1 ~ obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                          data=RE_par)
lm_re_mean_rec2        <- lm(mean_rec2 ~ obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                          data=RE_par)
lm_re_Y_FXSPR_static   <- lm(Y_FXSPR_static ~ obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                          data=RE_par)
lm_re_SSB_FXSPR_static <- lm(SSB_FXSPR_static ~ obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                          data=RE_par)
lm_re_NAA_rho          <- lm(NAA_rho ~ obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                          data=RE_par)
lm_re_Ecov_beta        <- lm(Ecov_beta ~ obs_error + R_sig + Fhist +  NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                          data=RE_par)


pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','re_lm_pars.pdf'),height=10,width=10)
par(mfrow=c(2,3),mar=c(1,2,1,2),oma=c(6,4,2,2),cex.axis=0.9)
ylims <- c(-10000,10000)
plotlm(lm_re_mean_rec1,add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
#plotlm(FITS[['lme_recr_re_all']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),mean=mean(recr.df$RE),nlme_try=TRUE)#,int=TRUE)
  mtext('a) RE mean_rec1',adj=0.0)
  mtext(expression('RE or'~Delta*'RE'),side=2,line=2.5)

ylims <- c(-3,3)
plotlm(lm_re_mean_rec2,add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
#plotlm(FITS[['lme_recr_re_all']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),mean=mean(recr.df$RE),nlme_try=TRUE)#,int=TRUE)
  mtext('a) RE mean_rec2',adj=0.0)
  mtext(expression('RE or'~Delta*'RE'),side=2,line=2.5)

ylims <- c(-0.01,0.01)
plotlm(lm_re_Y_FXSPR_static,add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
#plotlm(FITS[['lme_recr_re_all']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),mean=mean(recr.df$RE),nlme_try=TRUE)#,int=TRUE)
  mtext('a) RE Y_FXSPR_static',adj=0.0)
  mtext(expression('RE or'~Delta*'RE'),side=2,line=2.5)

ylims <- c(-0.01,0.01)
plotlm(lm_re_SSB_FXSPR_static,add=FALSE,ylim=ylims,labels=labels,nlme_try=FALSE)#,int=TRUE)
#plotlm(FITS[['lme_recr_re_all']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),mean=mean(recr.df$RE),nlme_try=TRUE)#,int=TRUE)
  mtext('a) RE SSB_FXSPR_static',adj=0.0)
  mtext(expression('RE or'~Delta*'RE'),side=2,line=2.5)

ylims <- c(-5,5)
plotlm(lm_re_NAA_rho,add=FALSE,ylim=ylims,labels=labels,nlme_try=FALSE)#,int=TRUE)
#plotlm(FITS[['lme_recr_re_all']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),mean=mean(recr.df$RE),nlme_try=TRUE)#,int=TRUE)
  mtext('a) RE NAA_rho',adj=0.0)
  mtext(expression('RE or'~Delta*'RE'),side=2,line=2.5)

ylims <- c(-1,1)
plotlm(lm_re_Ecov_beta,add=FALSE,ylim=ylims,labels=labels,nlme_try=FALSE)#,int=TRUE)
#plotlm(FITS[['lme_recr_re_all']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),mean=mean(recr.df$RE),nlme_try=TRUE)#,int=TRUE)
  mtext('a) RE Ecov_beta',adj=0.0)
  mtext(expression('RE or'~Delta*'RE'),side=2,line=2.5)

dev.off()
