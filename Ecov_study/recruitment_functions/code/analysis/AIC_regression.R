library("here")
#library("lme4")  #i have a package incompatability (liz)
library('tidyverse')



# read in converged/ non-converged runs corresponding to these thresholds ====
## specify bad.grad.label and bad.se.value (these are the thresholds set in convergence_summaries.R to determine convergence) 
bad.grad.value <- 1E-6 #tim used 1E-6 (the abs of the exponent will be used for output suffix; ex: filename_grad_6.png)
bad.grad.label <- as.numeric(strsplit(as.character(bad.grad.value), split="-")[[1]][2])
bad.se.value <- 100 #tim used 100 (this value will be used for output suffix; ex: filename_se_100.png)

res.dir <- 'results'       # results folder where AIC dataframes are
plot.dir <- 'plots'
plot.suffix <- '_beta_fix'


## get converged runs with these thresholds
AIC_best <- readRDS( file.path(here(),'Ecov_study','recruitment_functions',res.dir, paste0("AIC_best_conv.runs_grad_", bad.grad.label, "_SE_", bad.se.value, ".RDS")  ) )  # EM with lowest AIC per OM-Sim (that converged)
AIC_weight <- readRDS(file.path(here(),'Ecov_study','recruitment_functions',res.dir, paste0("AIC_weight_conv.runs_grad_", bad.grad.label, "_SE_", bad.se.value, ".RDS")  ))  # info for all EMs per OM-Sim

AIC_best_all <- readRDS(file.path(here(),'Ecov_study','recruitment_functions',res.dir,'AIC.rds')) #includes non=converged runs
aic.best.n <- nrow(AIC_best_all)
aic.best.n.conv <- nrow(AIC_best)
aic.best.n.bad <- aic.best.n - aic.best.n.conv
aic.best.pct.conv <- aic.best.n.conv/aic.best.n


AIC_weight_all <- readRDS(file.path(here(),'Ecov_study','recruitment_functions',res.dir,'AIC_weight.rds')) #includes non=converged runs

aic.weight.n <- nrow(AIC_weight_all)
aic.weight.n.conv <- nrow(AIC_weight)
aic.weight.n.bad <- aic.weight.n - aic.weight.n.conv
aic.weight.pct.conv <- aic.weight.n.conv/aic.weight.n

df.oms          <- readRDS(file.path(here::here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS"))
df.ems    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))


em_tib <- as_tibble(df.ems) %>%
  mutate(SR=ifelse(r_mod==2, 'Mean', 'BH')) %>%
  mutate(EM_mod = paste0(SR, "_", ecov_how), EM=seq(1,6))

############################################################
# Regression trees =====================
# use these to identify factors to use as facets in plots
############################################################

library(rpart.plot)

# analysis for best model selected ====
AIC_best$R_sig <- factor(AIC_best$R_sig,labels = c("L","H"))
AIC_best$Ecov_effect <- factor(AIC_best$Ecov_effect,labels=c("L","H"))
AIC_best$Ecov_how    <- factor(AIC_best$Ecov_how,labels=c("0", "1","2","4"))
AIC_best$NAA_cor     <- factor(AIC_best$NAA_cor,labels=c("L","H"))
AIC_best$Fhist       <- factor(AIC_best$Fhist,labels=c("H-MSY","MSY") ) 
AIC_best$Ecov_re_cor <- factor(AIC_best$Ecov_re_cor,labels=c("L","H"))

# regression tree for correct form (or SR or Ecov) ======================================
rf_SR   <- rpart(correct_SR   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how , data=AIC_best, control=rpart.control(cp=0.01))
imp.var <- rf_SR$frame[rf_SR$frame$var != '<leaf>',]
nodes_SR <- unique(imp.var[,1])
# "Fhist"   "R_sig"   "NAA_cor"

rf_ecov   <- rpart(correct_ecov   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how , data=AIC_best, control=rpart.control(cp=0.01))
imp.var <- rf_ecov$frame[rf_ecov$frame$var != '<leaf>',]
nodes_ecov <- unique(imp.var[,1])
# "Ecov_how"


rf_form <- rpart(correct_form ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how , data=AIC_best, control=rpart.control(cp=0.01))
imp.var <- rf_form$frame[rf_form$frame$var != '<leaf>',]
nodes_form <- unique(imp.var[,1])
# "R_sig"       "Fhist"       "Ecov_how"    "Ecov_effect"


pdf(file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('reg_tree_AIC_best', plot.suffix, '.pdf') ),
    height=4,width=7)
par(mfrow=c(1,3), oma=c(5,0,0,0))
prp(rf_SR,yesno=FALSE,type=4,clip.right.labs=TRUE)
mtext('a) correct SRR',adj=0,line=2.5, cex=0.9)
prp(rf_ecov,yesno=FALSE,type=4,clip.right.labs=TRUE)
mtext('b) correct Ecov_how',adj=0,line=2.5, cex=0.9)
title(sub=paste0(100*round(aic.best.pct.conv,2), "% of runs converged (", aic.best.n.bad, " out of ", aic.best.n, " failed)" ),   adj=0.5, outer=TRUE)
prp(rf_form,yesno=FALSE,type=4,clip.right.labs=TRUE)
mtext('c) correct SRR and Ecov_how',adj=0,line=2.5, cex=0.9)

dev.off()


# regression tree for dAIC ======================================
# analysis for best model selected ====
AIC_weight$R_sig <- factor(AIC_weight$R_sig,labels = c("L","H"))
AIC_weight$Ecov_effect <- factor(AIC_weight$Ecov_effect,labels=c("L","H"))
AIC_weight$Ecov_how    <- factor(AIC_weight$Ecov_how,labels=c("0", "1","2","4"))
AIC_weight$NAA_cor     <- factor(AIC_weight$NAA_cor,labels=c("L","H"))
AIC_weight$Fhist       <- factor(AIC_weight$Fhist,labels=c("H-MSY","MSY") ) 
AIC_weight$Ecov_re_cor <- factor(AIC_weight$Ecov_re_cor,labels=c("L","H"))


rf_dAIC <- rpart(dAIC ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how , data=AIC_weight, control=rpart.control(cp=0.01))
imp.var <- rf_dAIC$frame[rf_dAIC$frame$var != '<leaf>',]
nodes_dAIC <- unique(imp.var[,1])
# "R_sig"   "Fhist"   "NAA_cor"

# # this one didn't seem useful
rf_AIC_rank <- rpart(AIC_rank ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how , data=AIC_weight, control=rpart.control(cp=0.01))
imp.var <- rf_AIC_rank$frame[rf_AIC_rank$frame$var != '<leaf>',]
nodes_AIC_rank <- unique(imp.var[,1])
#"R_sig"



rf_Model_prob <- rpart(Model_prob ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how , data=AIC_weight, control=rpart.control(cp=0.01))
imp.var <- rf_Model_prob$frame[rf_Model_prob$frame$var != '<leaf>',]
nodes_Model_prob <- unique(imp.var[,1])
#"R_sig"


pdf(file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('reg_tree_AIC_prob', plot.suffix, '.pdf') ),
    height=4,width=7)
par(mfrow=c(1,3), oma=c(5,0,0,0))
prp(rf_dAIC,yesno=FALSE,type=4,clip.right.labs=TRUE)
mtext('a) dAIC',adj=0,line=2.5, cex=0.9)
prp(rf_AIC_rank,yesno=FALSE,type=4,clip.right.labs=TRUE)
mtext('b) AIC rank',adj=0,line=2.5, cex=0.9)
title(sub=paste0(100*round(aic.best.pct.conv,2), "% of runs converged (", aic.best.n.bad, " out of ", aic.best.n, " failed)" ),   adj=0.5, outer=TRUE)
prp(rf_Model_prob,yesno=FALSE,type=4,clip.right.labs=TRUE)
mtext('c) Model probability',adj=0,line=2.5, cex=0.9)

dev.off()

most.nodes <- nodes_dAIC
if(length(nodes_AIC_rank)>length(most.nodes) ) most.nodes<- nodes_AIC_rank
if(length(nodes_Model_prob)>length(most.nodes) ) most.nodes<- nodes_Model_prob

nodes.expr <- sapply(most.nodes, as.expression)
sapply(nodes.expr, eval)

# summarize mean model probability ====
AIC_w_mean <- as_tibble(AIC_weight) %>%
  mutate(mod.match=ifelse(EM_ecov_how==Ecov_how & EM_r_mod==recruit_mod, 1, 0)) %>%
  #filter(AIC_rank==1 | mod.match==1) %>%
  relocate(sim, OM, EM, mod.match, AIC_rank , Model_prob, dAIC, Ecov_how, EM_ecov_how, 
           recruit_mod, EM_r_mod) %>%
  #group_by( sapply(nodes.expr, eval), Ecov_how, EM, mod.match ) %>%
  group_by( R_sig,  Fhist, NAA_cor, Ecov_how, EM, mod.match ) %>%
    summarise(mean.prob=mean(Model_prob), var.prob=var(Model_prob), median.prob=median(Model_prob), across(Model_prob, 
                                      list(p2.5=~quantile(.,probs=0.0275 ),
                                           p97.5=~quantile(.,probs=0.975))),
                mean.daic=mean(dAIC), var.daic=var(dAIC), median.daic=median(dAIC), across(dAIC, 
                                                                  list(p2.5=~quantile(.,probs=0.0275 ),
                                                                       p97.5=~quantile(.,probs=0.975))) )%>%
  ungroup() %>%
  left_join(em_tib)


# AIC_w_merge <- as_tibble(AIC_weight) %>%
#   left_join(AIC_w_mean) 
# 
# AIC_w_merge$R_sig <- factor(AIC_w_merge$R_sig,labels = c("Rsig_0.1","Rsig_1.0"))
# AIC_w_merge$Ecov_effect <- factor(AIC_w_merge$Ecov_effect,labels=c("Ecov_L","Ecov_H"))
# AIC_w_merge$Ecov_how    <- factor(AIC_w_merge$Ecov_how,labels=c("0", "1","2","4"))
# AIC_w_merge$NAA_cor     <- factor(AIC_w_merge$NAA_cor,labels=c("L","H"))
# AIC_w_merge$Fhist       <- factor(AIC_w_merge$Fhist,labels=c("H-MSY","MSY") ) 
# AIC_w_merge$Ecov_re_cor <- factor(AIC_w_merge$Ecov_re_cor,labels=c("L","H"))

AIC_w_mean$R_sig <- factor(AIC_w_mean$R_sig,labels = c("Rsig_0.1","Rsig_1.0"))
# AIC_w_mean$Ecov_effect <- factor(AIC_w_mean$Ecov_effect,labels=c("Ecov_L", "Ecov_H"))
AIC_w_mean$Ecov_how    <- factor(AIC_w_mean$Ecov_how,labels=c("Ecov_0", "Ecov_1","Ecov_2","Ecov_4"))
AIC_w_mean$NAA_cor     <- factor(AIC_w_mean$NAA_cor,labels=c("Rcor_L","Rcor_H"))
AIC_w_mean$Fhist       <- factor(AIC_w_mean$Fhist,labels=c("H-MSY","MSY") )
#AIC_w_mean$Ecov_re_cor <- factor(AIC_w_mean$Ecov_re_cor,labels=c("EcovCor_L","EcovCor_H"))



mod.prob.plot <- ggplot(AIC_w_mean, aes(x=EM_mod, y=mean.prob, col=as.factor(mod.match) )) +
  facet_grid(R_sig  + Ecov_how~ Fhist  +NAA_cor  ) +
  geom_point(size=2) +
  geom_hline(yintercept=0.5, col='grey55') +
  # geom_errorbar(aes(ymin=(mean.prob-2*sqrt(var.prob)), ymax=(mean.prob+2*sqrt(var.prob)), width=.4) )+
  # ylab('Model Probability (Mean +/- 2 std)') +
  geom_errorbar(aes(ymin=(Model_prob_p2.5), ymax=(Model_prob_p97.5), width=.4) )+
  ylab('Model Probability (2.5 - 50 - 97.5 percentiles)') +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL))
ggsave(mod.prob.plot, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('model_probability', plot.suffix, '.png') ),  height=7, width=12)


dAIC.plot <- ggplot(AIC_w_mean, aes(x=EM_mod, y=mean.daic, col=as.factor(mod.match) )) +
  facet_grid(R_sig+Ecov_how ~ Fhist +NAA_cor    ) +
  geom_point(size=2) +
  geom_hline(yintercept=2, col='grey55') +
  geom_errorbar(aes(ymin=(mean.daic-2*sqrt(var.daic)), ymax=(mean.daic+2*sqrt(var.daic)), width=.4) )+
#  geom_errorbar(aes(ymin=(dAIC_p2.5), ymax=(dAIC_p97.5), width=.4) )+
  # ylab('Model Probability (2.5 - 50 - 97.5 percentiles)') +
  ylab('Model Probability (Mean +/- 2 std)') +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL))
ggsave(dAIC.plot, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir, paste0('model_dAIC', plot.suffix, '.png') ),  height=7, width=12)



dAIC.plot.ylim <- ggplot(AIC_w_mean, aes(x=EM_mod, y=mean.daic, col=as.factor(mod.match) )) +
  facet_grid(R_sig+Ecov_how ~ Fhist+NAA_cor  ) +
  geom_point(size=2) +
  geom_hline(yintercept=2, col='grey55') +
  geom_errorbar(aes(ymin=(mean.daic-2*sqrt(var.daic)), ymax=(mean.daic+2*sqrt(var.daic)), width=.4) )+
  #  geom_errorbar(aes(ymin=(dAIC_p2.5), ymax=(dAIC_p97.5), width=.4) )+
  # ylab('Model Probability (2.5 - 50 - 97.5 percentiles)') +
  ylab('Model Probability (Mean +/- 2 std)') +
  coord_cartesian(ylim=c(0,10) ) +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL))
ggsave(dAIC.plot.ylim, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('model_dAIC.ylim', plot.suffix, '.png') ),  height=7, width=12)


#===========================================================
#  fit glm and glmer ====
#===========================================================


fit_SR_glm   <- glm(correct_SR ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect, data=AIC, family='binomial')
fit_form_glm <- glm(correct_form ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect, data=AIC, family='binomial')

# Error with glmer:  (liz)
# Error in initializePtr() : 
#   function 'cholmod_factor_ldetA' not provided by package 'Matrix'
# answer: https://stackoverflow.com/questions/77481539/error-in-initializeptr-function-cholmod-factor-ldeta-not-provided-by-pack
# Matrix < 1.6-2 and Matrix >= 1.6-2 are binary incompatible. When you change between them, you must re-install from sources packages that link Matrix and therefore depend on the Matrix ABI:
# i have: > packageVersion("Matrix")
# [1] ‘1.6.1.1’
# packageVersion("lme4")
# [1] ‘1.1.35.1’
# Note that binaries in all repositories will be rebuilt automatically once lme4 > 1.1-35.1 is released. 

fit_SR_glmer   <- glmer(correct_SR ~ (1|sim) + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect, data=AIC, family='binomial')
fit_form_glmer <- glmer(correct_form ~ (1|sim) + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect, data=AIC, family='binomial')
#===========================================================

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

pdf(file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('effect_size_glm_glmm', plot.suffix, '.pdf') ),
    height=4,width=7)
par(mfrow=c(1,2),mar=c(2,2,2,2),oma=c(5,2,2,2),cex.axis=1)
plotlm(fit_SR_glm,ylims=c(-3,2))
legend('bottomright',bty='n',legend=c('GLM','GLMM'),pch='-',lty=1,col=c('black','red'),cex=0.6)
plotlm(fit_SR_glmer,add=TRUE)

plotlm(fit_form_glm,ylims=c(-3,2))
plotlm(fit_form_glmer,add=TRUE)
mtext(outer=TRUE,side=2,expression('Effect Size [log odds scale]'))
dev.off()

######################################################################
# Raw box plots ====
######################################################################

pdf(file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('raw_boxplots_form', plot.suffix,'.pdf') ),
    height=5,width=7)
AIC_best$correct_form <- factor(AIC_best$correct_form,labels=c("N","Y"))
par(mfrow=c(2,3),mar=c(2,2,2,3),oma=c(2,2,5,2))
plot(AIC_best$correct_form ~ AIC_best$Ecov_effect); mtext("Ecov_effect",adj=0)
title('Correct SRR and Ecov', outer=TRUE)
plot(AIC_best$correct_form ~ AIC_best$R_sig); mtext("R_sig",adj=0)
plot(AIC_best$correct_form ~ AIC_best$NAA_cor); mtext("NAA_cor",adj=0)
plot(AIC_best$correct_form ~ AIC_best$Fhist); mtext("Fhist",adj=0, xaxt="n")
plot(AIC_best$correct_form ~ AIC_best$Ecov_re_cor); mtext("Ecov_re_cor",adj=0)
plot(AIC_best$correct_form ~ AIC_best$Ecov_how); mtext("Ecov_how",adj=0)
#mtext(outer=TRUE,pro)
dev.off()


pdf(file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('raw_boxplots_SR', plot.suffix, '.pdf') ),
    height=5,width=7)
AIC_best$correct_SR <- factor(AIC_best$correct_SR,labels=c("N","Y"))
par(mfrow=c(2,3),mar=c(2,2,2,3),oma=c(2,2,5,2))
plot(AIC_best$correct_SR ~ AIC_best$Ecov_effect); mtext("Ecov_effect",adj=0)
title('Correct SRR', outer=TRUE)
plot(AIC_best$correct_SR ~ AIC_best$R_sig); mtext("R_sig",adj=0)
plot(AIC_best$correct_SR ~ AIC_best$NAA_cor); mtext("NAA_cor",adj=0)
plot(AIC_best$correct_SR ~ AIC_best$Fhist); mtext("Fhist",adj=0)
plot(AIC_best$correct_SR ~ AIC_best$Ecov_re_cor); mtext("Ecov_re_cor",adj=0)
plot(AIC_best$correct_SR ~ AIC_best$Ecov_how); mtext("Ecov_how",adj=0)
#mtext(outer=TRUE,side=4,"Proportion",las=0)
dev.off()

pdf(file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('raw_boxplots_Ecov', plot.suffix, '.pdf') ),
    height=5,width=7)
AIC_best$correct_SR <- factor(AIC_best$correct_SR,labels=c("N","Y"))
par(mfrow=c(2,3),mar=c(2,2,2,3),oma=c(2,2,5,2))
plot(AIC_best$correct_SR ~ AIC_best$Ecov_effect); mtext("Ecov_effect",adj=0)
title('Correct Ecov_how', outer=TRUE)
plot(AIC_best$correct_ecov ~ AIC_best$R_sig); mtext("R_sig",adj=0)
plot(AIC_best$correct_ecov ~ AIC_best$NAA_cor); mtext("NAA_cor",adj=0)
plot(AIC_best$correct_ecov ~ AIC_best$Fhist); mtext("Fhist",adj=0)
plot(AIC_best$correct_ecov ~ AIC_best$Ecov_re_cor); mtext("Ecov_re_cor",adj=0)
plot(AIC_best$correct_ecov ~ AIC_best$Ecov_how); mtext("Ecov_how",adj=0)
#mtext(outer=TRUE,side=4,"Proportion",las=0)
dev.off()

