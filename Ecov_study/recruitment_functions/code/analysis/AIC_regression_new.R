library(glmnet)
library(rpart)

AIC <- readRDS(file.path(here::here(),'Ecov_study','recruitment_functions','results','AIC.rds'))

X <- data.frame(obs_error=as.factor(AIC$obs_error), 
                R_sig=AIC$R_sig,
                Fhist=as.factor(AIC$Fhist),
                NAA_cor=as.factor(AIC$NAA_cor), 
                Ecov_re_cor=as.factor(AIC$Ecov_re_cor), 
                Ecov_effect=as.factor(AIC$Ecov_effect), 
                Ecov_how=as.factor(AIC$Ecov_how), 
                ssb_cv=AIC$ssb_cv, 
                ecov_slope=AIC$ecov_slope)



fit_glm_ecov <- glm(AIC$correct_ecov ~ ., data=X, family='binomial')
summary(fit_glm_ecov)

fit_glm_SR <- glm(AIC$correct_SR ~ ., data=X, family='binomial')
summary(fit_glm_SR)

fit_glm_form <- glm(AIC$correct_form ~ ., data=X, family='binomial')
summary(fit_glm_form)


fit_glmnet <- glmnet(y=AIC$correct_ecov, x=X,family='binomial',alpha=0)
plot(fit_glmnet)
summary(fit_glm)



rf_SR <- rpart(correct_SR   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=AIC_best, control=rpart.control(cp=0.01))









hist(AIC$ecov_slope)



