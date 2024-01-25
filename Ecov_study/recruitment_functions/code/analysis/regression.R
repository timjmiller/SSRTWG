library("here")
library("lme4")

AIC <- readRDS(file.path(here(),'Ecov_study','recruitment_functions','results','AIC.rds'))

fit_SR_glm   <- glm(correct_SR ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect, data=AIC, family='binomial')
fit_form_glm <- glm(correct_form ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect, data=AIC, family='binomial')

fit_SR_glmer   <- glmer(correct_SR ~ (1|sim) + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect, data=AIC, family='binomial')
fit_form_glmer <- glmer(correct_form ~ (1|sim) + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect, data=AIC, family='binomial')


plotlm <- function(fit,add=FALSE,ylims){
  coef <- summary(fit)$coefficients[,1]
  ses  <- summary(fit)$coefficients[,2]
  n    <- length(coef)
  if(add==FALSE){
    #plot(1:n,coef,ylim=c(min(coef-2*ses),max(coef+2*ses)),pch=19,xaxt='n')
    plot(1:n,coef,pch='-',xaxt='n',ylim=ylims,xlim=c(1,n+0.5))
    segments(x0=1:n,x1=1:n,y0=coef-2*ses,y1=coef+2*ses)
  }
  if(add==TRUE){
    points(1:n+0.25,coef,ylim=c(min(coef-2*ses),max(coef+2*ses)),pch='-',col='red')
    segments(x0=1:n+0.25,x1=1:n+0.25,y0=coef-2*ses,y1=coef+2*ses,col='red')
  }
  abline(h=0,lty=2)
  axis(side=1,at=1:n,labels=row.names(summary(fit)$coefficients),las=2)
}

#predict(fit_SR_glm,type='response')

pdf(file.path(here(),'Ecov_study','recruitment_functions','plots','effect_size_glm_glmm.pdf'),
    height=4,width=7)
par(mfrow=c(1,2),mar=c(2,2,2,2),oma=c(5,2,2,2),cex.axis=1)
plotlm(fit_SR_glm,ylims=c(-3,2))
legend('bottomright',bty='n',legend=c('GLM','GLMM'),pch='-',lty=1,col=c('black','red'),cex=0.6)
plotlm(fit_SR_glmer,add=TRUE)

plotlm(fit_form_glm,ylims=c(-3,2))
plotlm(fit_form_glmer,add=TRUE)
mtext(outer=TRUE,side=2,expression('Effect Size [log odds scale]'))
dev.off()

############################################################
############################################################

library(rpart.plot)

AIC$R_sig <- factor(AIC$R_sig,labels = c("L","H"))
AIC$Ecov_effect <- factor(AIC$Ecov_effect,labels=c("L","H"))
#AIC$Ecov_how    <- factor(AIC$Ecov_how,labels=c("1","2","4"))
AIC$Ecov_how    <- factor(AIC$Ecov_how,labels=c("0","1"))
AIC$NAA_cor     <- factor(AIC$NAA_cor,labels=c("L","H"))
AIC$Fhist       <- factor(AIC$Fhist)
AIC$Ecov_re_cor <- factor(AIC$Ecov_re_cor,labels=c("L","H"))

#rf_SR   <- rpart(correct_SR   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect, data=AIC)
#rf_form <- rpart(correct_form ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect, data=AIC)
rf_SR   <- rpart(correct_SR   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how, data=AIC, control=rpart.control(cp=0.01))
rf_form <- rpart(correct_form ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how, data=AIC, control=rpart.control(cp=0.01))

pdf(file.path(here(),'Ecov_study','recruitment_functions','plots','reg_tree.pdf'),
    height=4,width=7)
par(mfrow=c(1,2))
prp(rf_SR,yesno=FALSE,type=4,clip.right.labs=TRUE)
mtext('a)',adj=0,line=1.5)
prp(rf_form,yesno=FALSE,type=4,clip.right.labs=TRUE)
mtext('b)',adj=0,line=1.5)
dev.off()


######################################################################
######################################################################

pdf(file.path(here(),'Ecov_study','recruitment_functions','plots','raw_boxplots_form.pdf'),
    height=5,width=7)
AIC$correct_form <- factor(AIC$correct_form,labels=c("N","Y"))
par(mfrow=c(2,3),mar=c(2,2,2,3),oma=c(2,2,2,2))
plot(AIC$correct_form ~ AIC$Ecov_effect); mtext("Ecov_effect",adj=0)
plot(AIC$correct_form ~ AIC$R_sig); mtext("R_sig",adj=0)
plot(AIC$correct_form ~ AIC$NAA_cor); mtext("NAA_cor",adj=0)
plot(AIC$correct_form ~ AIC$Fhist); mtext("Fhist",adj=0)
plot(AIC$correct_form ~ AIC$Ecov_re_cor); mtext("Ecov_re_cor",adj=0)
plot(AIC$correct_form ~ AIC$Ecov_how); mtext("Ecov_how",adj=0)
#mtext(outer=TRUE,pro)
dev.off()


pdf(file.path(here(),'Ecov_study','recruitment_functions','plots','raw_boxplots_SR.pdf'),
    height=5,width=7)
AIC$correct_SR <- factor(AIC$correct_SR,labels=c("N","Y"))
par(mfrow=c(2,3),mar=c(2,2,2,3),oma=c(2,2,2,2))
plot(AIC$correct_SR ~ AIC$Ecov_effect); mtext("Ecov_effect",adj=0)
plot(AIC$correct_SR ~ AIC$R_sig); mtext("R_sig",adj=0)
plot(AIC$correct_SR ~ AIC$NAA_cor); mtext("NAA_cor",adj=0)
plot(AIC$correct_SR ~ AIC$Fhist); mtext("Fhist",adj=0)
plot(AIC$correct_SR ~ AIC$Ecov_re_cor); mtext("Ecov_re_cor",adj=0)
plot(AIC$correct_SR ~ AIC$Ecov_how); mtext("Ecov_how",adj=0)
#mtext(outer=TRUE,side=4,"Proportion",las=0)
dev.off()

