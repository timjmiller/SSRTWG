#library(rpart)
#library(rpart.plot)
#library(nlme)
library(tidyverse)
library(here)

source(file.path(here::here(),'Ecov_study','recruitment_functions','code','analysis','functions_for_analysis.R'))

res.dir     <- 'results'   # 'results'     'results_beta_fix'   # results folder where AIC dataframes are
plot.dir    <- 'plots'    # 'plots_lizruns'  'plots_beta_fix'   
plot.suffix <- ''      # '_beta_fix'   '' 

bad.grad.value <- 1E-6 #tim used 1E-6 (the abs of the exponent will be used for output suffix; ex: filename_grad_6.png)
bad.grad.label <- as.numeric(strsplit(as.character(bad.grad.value), split="-")[[1]][2])
bad.se.value   <- 100 #tim used 100 (this value will be used for output suffix; ex: filename_se_100.png)

AIC_weight             <- readRDS( file.path(here(),'Ecov_study','recruitment_functions',res.dir, paste0("AIC_weight_conv.runs_grad_", bad.grad.label, "_SE_", bad.se.value, plot.suffix, ".RDS")  ) ) %>% # EM with lowest AIC per OM-Sim (that converged)
    filter(AIC_rank==2)   #find the second best model
AIC_weight$obs_error   <- factor(AIC_weight$obs_error,levels=c("L","H"))
AIC_weight$R_sig       <- as.factor(AIC_weight$R_sig)
AIC_weight$Fhist       <- factor(AIC_weight$Fhist,levels=c("MSY","L-H","H-MSY"))
AIC_weight$NAA_cor     <- as.factor(AIC_weight$NAA_cor) 
AIC_weight$Ecov_re_cor <- as.factor(AIC_weight$Ecov_re_cor) 
AIC_weight$Ecov_effect <- as.factor(AIC_weight$Ecov_effect) 
AIC_weight$Ecov_how    <- as.factor(AIC_weight$Ecov_how) 
AIC_weight$ssb_cv      <- factor(case_when(AIC_weight$ssb_cv < mean(AIC_weight$ssb_cv) - sd(AIC_weight$ssb_cv) ~ 'L',
                                    AIC_weight$ssb_cv > mean(AIC_weight$ssb_cv) + sd(AIC_weight$ssb_cv) ~ "H",
                                    TRUE ~ 'M')) 
AIC_weight$ssb_cv <- factor(AIC_weight$ssb_cv,levels=c("L","M","H"))
AIC_weight$ecov_slope  <- as.factor(case_when(AIC_weight$ecov_slope > mean(AIC_weight$ecov_slope) + sd(AIC_weight$ecov_slope) ~ "+",
                                       AIC_weight$ecov_slope < mean(AIC_weight$ecov_slope) - sd(AIC_weight$ecov_slope) ~ '-',
                                       TRUE ~ '0'))
AIC_weight$ecov_slope <- factor(AIC_weight$ecov_slope,levels=c("-","0","+"))

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

pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','dAIC_marg.pdf'),height=3.5,width=5)
par(mfrow=c(1,1),mar=c(1,2,1,0),oma=c(6,4,2,2),cex.axis=0.7,cex.lab=0.7)
ylims <- c(0,1.5)
dd(AIC_weight,vars=vars,labels=labels,yvar="dAIC",ylims=ylims,mean=TRUE)
  mtext(expression(Delta*'AIC'),side=2,line=2)
  #mtext(expression('a) Constant Ecov'),adj=0)
dev.off()



##################################################################
##--Decision tree
##################################################################

tree_daic <- rpart(dAIC ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                data=AIC_weight)

pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','tree_dAIC.pdf'),height=4,width=4)
par(mfrow=c(1,1), oma=c(5,0,0,0))
prp(tree_daic,yesno=FALSE,type=4,clip.right.labs=TRUE)
  mtext('a)',adj=0,line=-1, cex=0.9)
dev.off()


