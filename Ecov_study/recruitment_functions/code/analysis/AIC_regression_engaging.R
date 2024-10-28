library(rpart)
library(rpart.plot)
library(nlme)
library(lme4)
library(tidyverse)
library(here)

source(file.path(here::here(),'Ecov_study','recruitment_functions','code','analysis','functions_for_analysis.R'))

res.dir     <- 'results'   # 'results'     'results_beta_fix'   # results folder where AIC dataframes are
plot.dir    <- 'plots'    # 'plots_lizruns'  'plots_beta_fix'   
plot.suffix <- ''      # '_beta_fix'   '' 

####################################
##--AIC--###########################
####################################
# AIC <- readRDS(file.path(here::here(),'Ecov_study','recruitment_functions','results','AIC.rds'))
# AIC$obs_error   <- as.factor(AIC$obs_error)
# AIC$R_sig       <- as.factor(AIC$R_sig)
# AIC$Fhist       <- as.factor(AIC$Fhist)
# AIC$NAA_cor     <- as.factor(AIC$NAA_cor) 
# AIC$Ecov_re_cor <- as.factor(AIC$Ecov_re_cor) 
# AIC$Ecov_effect <- as.factor(AIC$Ecov_effect) 
# AIC$Ecov_how    <- as.factor(AIC$Ecov_how) 
# AIC$ssb_cv      <- factor(case_when(AIC$ssb_cv < mean(AIC$ssb_cv) - sd(AIC$ssb_cv) ~ 'L',
#                                        AIC$ssb_cv > mean(AIC$ssb_cv) + sd(AIC$ssb_cv) ~ "H",
#                                        TRUE ~ 'M')) 
# AIC$ssb_cv <- factor(AIC$ssb_cv,levels=c("L","M","H"))
# AIC$ecov_slope  <- as.factor(case_when(AIC$ecov_slope > mean(AIC$ecov_slope) + sd(AIC$ecov_slope) ~ "H",
#                                        AIC$ecov_slope < mean(AIC$ecov_slope) - sd(AIC$ecov_slope) ~ 'L',
#                                        TRUE ~ 'M'))
# AIC$ecov_slope <- factor(AIC$ecov_slope,levels=c("L","M","H"))
# 
# ##--FIT MODELS--######################
# glm_ecov <- glm(correct_ecov ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
#                 data=AIC, family='binomial')
# glm_SR   <- glm(correct_SR   ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
#                 data=AIC, family='binomial')
# glm_form <- glm(correct_form ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
#                 data=AIC, family='binomial')
# 
# glmm_ecov <- glmer(correct_ecov ~  (1|sim) + obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
#                    data=AIC, family='binomial')
# glmm_SR   <- glmer(correct_SR   ~  (1|sim) + obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
#                    data=AIC, family='binomial')
# glmm_form <- glmer(correct_form ~  (1|sim) + obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
#                    data=AIC, family='binomial')
# 
# FITS_AIC <- list(glm_ecov,glm_SR,glm_form,
#                  glmm_ecov,glmm_SR,glmm_form)
# saveRDS(FITS_AIC,file=file.path(here::here(),'Ecov_study','recruitment_functions','results','FITS_AIC.rds'))


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
AIC_best$ecov_slope  <- as.factor(case_when(AIC_best$ecov_slope > mean(AIC_best$ecov_slope) + sd(AIC_best$ecov_slope) ~ "H",
                                       AIC_best$ecov_slope < mean(AIC_best$ecov_slope) - sd(AIC_best$ecov_slope) ~ 'L',
                                       TRUE ~ 'M'))
AIC_best$ecov_slope <- factor(AIC_best$ecov_slope,levels=c("L","M","H"))

##--FIT MODELS--######################
glm_ecov <- glm(correct_ecov ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                data=AIC_best, family='binomial')
glm_SR   <- glm(correct_SR   ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                data=AIC_best, family='binomial')
glm_form <- glm(correct_form ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                data=AIC_best, family='binomial')

glmm_ecov <- glmer(correct_ecov ~  (1|sim) + obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                   data=AIC_best, family='binomial')
glmm_SR   <- glmer(correct_SR   ~  (1|sim) + obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                   data=AIC_best, family='binomial')
glmm_form <- glmer(correct_form ~  (1|sim) + obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                   data=AIC_best, family='binomial')

FITS_AIC <- list(glm_ecov,glm_SR,glm_form,
                 glmm_ecov,glmm_SR,glmm_form)
saveRDS(FITS_AIC,file=file.path(here::here(),'Ecov_study','recruitment_functions','results','FITS_AIC.rds'))





##--MAKE PLOTS--###########################
labels <- c(expression(italic('intercept')),
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

pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','effects.pdf'),height=4.5,width=4.5)
par(mfrow=c(3,1),mar=c(1,2,1,2),oma=c(6,2,2,2),cex.axis=0.9)
ylims <- c(-5,5)
plotlm(glm_SR,add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),int=TRUE)
plotlm(glmm_SR,add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),int=TRUE)
  mtext('a) SR Y/N',adj=0.0)
plotlm(glm_ecov,add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),int=TRUE)
plotlm(glmm_ecov,add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),int=TRUE)
  mtext(expression('b) E'['cov']~'Y/N'),adj=0.0)
plotlm(glm_form,add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),int=TRUE)
plotlm(glmm_form,add=TRUE,ylim=ylims,labels=labels,int=TRUE)
  mtext(expression('c) SR & E'['cov']~'Y/N'),adj=0.0)
mtext(outer=TRUE,expression(Delta*'log['~italic('p/(1-p)')*']'),side=2)
dev.off()


##################################################################
##--Decision tree
##################################################################

tree_ecov <- rpart(correct_ecov ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                data=AIC_best)
tree_SR   <- rpart(correct_SR   ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                data=AIC_best)
tree_form <- rpart(correct_form ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                data=AIC_best)

pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','trees_AIC.pdf'),height=5,width=8)
par(mfrow=c(1,3), oma=c(5,0,0,0))
prp(tree_SR,yesno=FALSE,type=4,clip.right.labs=TRUE)
  mtext('a) SRR Y/N',adj=0,line=-1, cex=0.9)
prp(tree_ecov,yesno=FALSE,type=4,clip.right.labs=TRUE)
  mtext(expression('b) E'['cov']*'Y/N'),adj=0,line=-1, cex=0.9)
  #title(sub=paste0(100*round(aic.best.pct.conv,2), "% of runs converged (", aic.best.n.bad, " out of ", aic.best.n, " failed)" ),   adj=0.5, outer=TRUE)
prp(tree_form,yesno=FALSE,type=4,clip.right.labs=TRUE)
  mtext(expression('c) SRR and'~'E'['cov']~'Y/N'),adj=0,line=-1, cex=0.9)
dev.off()

