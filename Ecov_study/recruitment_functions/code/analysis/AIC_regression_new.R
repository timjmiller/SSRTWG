library(glmnet)
library(rpart)
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
AIC$ssb_cv <- relevel(AIC$ssb_cv,ref="L")
AIC$ecov_slope  <- as.factor(case_when(AIC$ecov_slope > mean(AIC$ecov_slope) + sd(AIC$ecov_slope) ~ "H",
                                       AIC$ecov_slope < mean(AIC$ecov_slope) - sd(AIC$ecov_slope) ~ 'L',
                                       TRUE ~ 'M'))
AIC$ecov_slope <- relevel(AIC$ecov_slope,ref='L')

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


##--MAKE PLOTS--###########################
labels <- c(expression(italic('int')),
            expression(sigma['obs']~'= L'),
            expression(sigma['r']~'= L'),
            expression(sigma['r']~'= L'),
            expression(sigma['r']~'= L'),
            expression(sigma['r']~'= L'),
            expression(italic('F'['hist'])~'= L-H'),
            expression(italic('F'['hist'])~'= MSY'),
            expression(rho['NAA']), 
            expression(rho['ecov']), 
            expression(beta['ecov']), 
            expression(italic('ecov'['how'])~'= 1'), 
            expression(italic('ecov'['how'])~'= 2'), 
            expression(Delta*'ecov'), 
            expression(italic('CV'['ecov'])))

par(mfrow=c(3,1),mar=c(1,2,0,2),oma=c(4,2,2,2))
ylims <- c(-4,4)
plotlm(glm_ecov,add=FALSE,ylim=ylims,labels=rep(NA,12))
plotlm(glmm_ecov,add=TRUE,ylim=ylims,labels=rep(NA,12))

plotlm(glm_SR,add=FALSE,ylim=ylims,labels=rep(NA,12))
plotlm(glmm_SR,add=TRUE,ylim=ylims,labels=rep(NA,12))

plotlm(glm_form,add=FALSE,ylim=ylims,labels=rep(NA,12))
plotlm(glmm_form,add=TRUE,ylim=ylims,labels=labels)


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
##--Sparse regression
##################################################################




rf_SR <- rpart(correct_SR ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=AIC_best, control=rpart.control(cp=0.01))









hist(AIC$ecov_slope)
plot(AIC$Ecov_re_cor,AIC$ecov_slope)


