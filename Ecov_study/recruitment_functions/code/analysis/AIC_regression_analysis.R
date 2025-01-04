#library(rpart)
#library(rpart.plot)
#library(nlme)
#library(lme4)
library(tidyverse)
library(here)

source(file.path(here::here(),'Ecov_study','recruitment_functions','code','analysis','functions_for_analysis.R'))

res.dir     <- 'results'   # 'results'     'results_beta_fix'   # results folder where AIC dataframes are
plot.dir    <- 'plots'    # 'plots_lizruns'  'plots_beta_fix'   
plot.suffix <- ''      # '_beta_fix'   '' 

#########################################
##--AIC_best--###########################
#########################################
bad.grad.value <- 1E-6 #tim used 1E-6 (the abs of the exponent will be used for output suffix; ex: filename_grad_6.png)
bad.grad.label <- as.numeric(strsplit(as.character(bad.grad.value), split="-")[[1]][2])
bad.se.value   <- 100 #tim used 100 (this value will be used for output suffix; ex: filename_se_100.png)

AIC_best             <- readRDS( file.path(here(),'Ecov_study','recruitment_functions',res.dir, paste0("AIC_best_conv.runs_grad_", bad.grad.label, "_SE_", bad.se.value, plot.suffix, ".RDS")  ) )  # EM with lowest AIC per OM-Sim (that converged)
AIC_best$obs_error   <- factor(AIC_best$obs_error,levels=c("L","H"))
AIC_best$R_sig       <- as.factor(AIC_best$R_sig)
AIC_best$Fhist       <- factor(AIC_best$Fhist,levels=c("MSY","L-H","H-MSY"))
AIC_best$NAA_cor     <- as.factor(AIC_best$NAA_cor) 
AIC_best$Ecov_re_cor <- as.factor(AIC_best$Ecov_re_cor) 
AIC_best$Ecov_effect <- as.factor(AIC_best$Ecov_effect) 
AIC_best$Ecov_how    <- as.factor(AIC_best$Ecov_how) 
AIC_best$ssb_cv      <- factor(case_when(AIC_best$ssb_cv < mean(AIC_best$ssb_cv) - sd(AIC_best$ssb_cv) ~ 'L',
                                    AIC_best$ssb_cv > mean(AIC_best$ssb_cv) + sd(AIC_best$ssb_cv) ~ "H",
                                    TRUE ~ 'M')) 
AIC_best$ssb_cv <- factor(AIC_best$ssb_cv,levels=c("L","M","H"))
AIC_best$ecov_slope  <- as.factor(case_when(AIC_best$ecov_slope > mean(AIC_best$ecov_slope) + sd(AIC_best$ecov_slope) ~ "+",
                                       AIC_best$ecov_slope < mean(AIC_best$ecov_slope) - sd(AIC_best$ecov_slope) ~ '-',
                                       TRUE ~ '0'))
AIC_best$ecov_slope <- factor(AIC_best$ecov_slope,levels=c("-","0","+"))



##--MAKE PLOTS--###########################
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

pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','AIC_marg.pdf'),height=5.75,width=5)
par(mfrow=c(3,1),mar=c(1,2,1,0),oma=c(6,4,2,2),cex.axis=0.7,cex.lab=0.7)
ylims <- c(0,1)
dd(AIC_best,vars=vars,labels=rep(NA,length(labels)),yvar="correct_SR",ylims=ylims,mean=TRUE)
  mtext(expression(italic(p(correct))),side=2,line=2.5)
  mtext(expression('a) SR Y/N'),adj=0)
dd(AIC_best,vars=vars,labels=rep(NA,length(labels)),yvar="correct_ecov",ylims=ylims,mean=TRUE)
  mtext(expression(italic(p(correct))),side=2,line=2.5)
  mtext(expression('b) E'['cov']~'Y/N'),adj=0)
dd(AIC_best,vars=vars,labels=labels,yvar="correct_form",ylims=ylims,mean=TRUE)
  mtext(expression(italic(p(correct))),side=2,line=2.5)
  mtext(expression('c) SR & E'['cov']~'Y/N'),adj=0)
dev.off()

