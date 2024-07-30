library(glmnet)
library(rpart)
library(rpart.plot)
library(lme4)
library(tidyverse)

source(file.path(here::here(),'Ecov_study','recruitment_functions','code','analysis','functions_for_analysis.r'))

##--data processing--###################
AIC <- readRDS(file.path(here::here(),'Ecov_study','recruitment_functions','results','AIC.rds'))
AIC$obs_error   <- as.factor(AIC$obs_error)
AIC$R_sig       <- as.factor(AIC$R_sig)
AIC$Fhist       <- as.factor(AIC$Fhist)
AIC$NAA_cor     <- as.factor(AIC$NAA_cor) 
AIC$Ecov_re_cor <- as.factor(AIC$Ecov_re_cor) 
AIC$Ecov_effect <- as.factor(AIC$Ecov_effect) 
AIC$Ecov_how    <- as.factor(AIC$Ecov_how) 
AIC$ssb_cv      <- factor(case_when(AIC$ssb_cv < mean(AIC$ssb_cv) - sd(AIC$ssb_cv) ~ 'L',
                                       AIC$ssb_cv > mean(AIC$ssb_cv) + sd(AIC$ssb_cv) ~ "H",
                                       TRUE ~ 'M')) 
AIC$ssb_cv <- factor(AIC$ssb_cv,levels=c("L","M","H"))
AIC$ecov_slope  <- as.factor(case_when(AIC$ecov_slope > mean(AIC$ecov_slope) + sd(AIC$ecov_slope) ~ "H",
                                       AIC$ecov_slope < mean(AIC$ecov_slope) - sd(AIC$ecov_slope) ~ 'L',
                                       TRUE ~ 'M'))
AIC$ecov_slope <- factor(AIC$ecov_slope,levels=c("L","M","H"))

##--FIT MODELS--######################
glm_ecov <- glm(correct_ecov ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                data=AIC, family='binomial')
glm_SR   <- glm(correct_SR   ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                data=AIC, family='binomial')
glm_form <- glm(correct_form ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                data=AIC, family='binomial')

glmm_ecov <- glmer(correct_ecov ~  (1|sim) + obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                   data=AIC, family='binomial')
glmm_SR   <- glmer(correct_SR   ~  (1|sim) + obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                   data=AIC, family='binomial')
glmm_form <- glmer(correct_form ~  (1|sim) + obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                   data=AIC, family='binomial')

FITS_AIC <- list(glm_ecov,glm_SR,glm_form,
                 glmm_ecov,glmm_SR,glmm_form)
saveRDS(FITS_AIC,file=file.path(here::here(),'Ecov_study','recruitment_functions','results','FITS_AIC.rds'))



##--MAKE PLOTS--###########################
labels <- c(expression(italic('intercept')),
            expression(sigma['obs']~'= L'),
            expression(sigma['r']~'= 0.3'),
            expression(sigma['r']~'= 0.5'),
            expression(sigma['r']~'= 1.0'),
            expression(italic('F')~'= L-H'),
            expression(italic('F')~'= MSY'),
            expression(rho['NAA']), 
            expression(rho['E'['cov']]), 
            expression(beta['E'['cov']]), 
            expression(italic('f(E'['cov']*')')~'= 1'), 
            expression(italic('f(E'['cov']*')')~'= 2'), 
            expression(Delta*'E'['cov']*' = M'), 
            expression(Delta*'E'['cov']*' = H'), 
            expression(italic('CV'['SSB']~'= M')),
            expression(italic('CV'['SSB']~'= H')))

pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','effects.pdf'),height=5,width=5)
par(mfrow=c(3,1),mar=c(1,2,1,2),oma=c(6,2,2,2),cex.axis=0.9)
ylims <- c(-5,5)
plotlm(glm_SR,add=FALSE,ylim=ylims,labels=rep(NA,length(labels)))
plotlm(glmm_SR,add=TRUE,ylim=ylims,labels=rep(NA,length(labels)))
  mtext('a) SR Y/N',adj=0.0)
plotlm(glm_ecov,add=FALSE,ylim=ylims,labels=rep(NA,length(labels)))
plotlm(glmm_ecov,add=TRUE,ylim=ylims,labels=rep(NA,length(labels)))
  mtext(expression('b) E'['cov']~'Y/N'),adj=0.0)
plotlm(glm_form,add=FALSE,ylim=ylims,labels=rep(NA,length(labels)))
plotlm(glmm_form,add=TRUE,ylim=ylims,labels=labels)
  mtext(expression('c) SR & E'['cov']~'Y/N'),adj=0.0)
mtext(outer=TRUE,expression(Delta*'log['~italic('p/(1-p)')*']'),side=2)
dev.off()

##################################################################
##--Sparse regression
##################################################################
#setup data
X$y <- AIC$correct_ecov
f  <- as.formula(y ~ .*.)
XX <- model.matrix(f, X)[,-1]

y <- AIC$correct_ecov
fit <- glmnet(x=XX,y=y,family='binomial')
plot(fit)

cvfit <- cv.glmnet(x=X, y=y, family='binomial')

plot(cvfit)




##################################################################
##--Decision tree
##################################################################

tree_ecov <- rpart(correct_ecov ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                data=AIC)
tree_SR   <- rpart(correct_SR   ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                data=AIC)
tree_form <- rpart(correct_form ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                data=AIC)

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

