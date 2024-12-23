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

FITS  <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "results", 'FITS_error.rds'))  # directory where simulation 

bad.grad.value <- 1E-6 #tim used 1E-6 (the abs of the exponent will be used for output suffix; ex: filename_grad_6.png)
bad.grad.label <- as.numeric(strsplit(as.character(bad.grad.value), split="-")[[1]][2])
bad.se.value   <- 100 #tim used 100 (this value will be used for output suffix; ex: filename_se_100.png)

recr.df <- as.data.frame(readRDS( file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "error.recr", plot.suffix, ".df.RDS") ) )) %>%
  filter(opt==1 & conv==0 & sdrep==0 & max_grad<bad.grad.value)

ssb.df  <- as.data.frame(readRDS( file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "error.ssb",  plot.suffix, ".df.RDS") )  )) %>%
  filter(opt==1 & conv==0 & sdrep==0 & max_grad<bad.grad.value)

fbar.df <- as.data.frame(readRDS( file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "error.fbar", plot.suffix, ".df.RDS") )  )) %>%
  filter(opt==1 & conv==0 & sdrep==0 & max_grad<bad.grad.value)

### TREES ############################

#recruitment
pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','trees_error_recruitment.pdf'),height=5,width=8)
par(mfrow=c(2,3), oma=c(5,0,3,0))

prp(FITS[['rf_recr_re_all']],yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('a) RE (all)',adj=0, line=2, cex=0.9)
prp(FITS[['rf_recr_re_ten']],yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('b) RE (ten)',adj=0,line=2, cex=0.9)
prp(FITS[['rf_recr_re_last']],yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('c) RE (last)',adj=0,line=2, cex=0.9)

prp(FITS[['rf_recr_rmse_all']],yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('d) RMSE (all)',adj=0, line=3, cex=0.9)
prp(FITS[['rf_recr_rmse_ten']],yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('e) RMSE (ten)',adj=0, line=3, cex=0.9)
prp(FITS[['rf_recr_rmse_last']],yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('f) RMSE (last)',adj=0, line=3, cex=0.9)
#title("Recruitment",outer=TRUE)
dev.off()

#ssb
pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','trees_error_ssb.pdf'),height=5,width=8)
par(mfrow=c(2,3), oma=c(5,0,3,0))

prp(FITS[['rf_ssb_re_all']],yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('a) RE (all)',adj=0, line=2, cex=0.9)
prp(FITS[['rf_ssb_re_ten']],yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('b) RE (ten)',adj=0,line=2, cex=0.9)
prp(FITS[['rf_ssb_re_last']],yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('c) RE (last)',adj=0,line=2, cex=0.9)

prp(FITS[['rf_ssb_rmse_all']],yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('d) RMSE (all)',adj=0, line=3, cex=0.9)
prp(FITS[['rf_ssb_rmse_ten']],yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('e) RMSE (ten)',adj=0, line=3, cex=0.9)
prp(FITS[['rf_ssb_rmse_last']],yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('f) RMSE (last)',adj=0, line=3, cex=0.9)
#title("Recruitment",outer=TRUE)
dev.off()


#fbar
pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','trees_error_fbar.pdf'),height=5,width=8)
par(mfrow=c(2,3), oma=c(5,0,3,0))

prp(FITS[['rf_fbar_re_all']],yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('a) RE (all)',adj=0, line=2, cex=0.9)
prp(FITS[['rf_fbar_re_ten']],yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('b) RE (ten)',adj=0,line=2, cex=0.9)
prp(FITS[['rf_fbar_re_last']],yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('c) RE (last)',adj=0,line=2, cex=0.9)

prp(FITS[['rf_fbar_rmse_all']],yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('d) RMSE (all)',adj=0, line=3, cex=0.9)
prp(FITS[['rf_fbar_rmse_ten']],yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('e) RMSE (ten)',adj=0, line=3, cex=0.9)
prp(FITS[['rf_fbar_rmse_last']],yesno=FALSE,type=4,clip.right.labs=TRUE,roundint=FALSE)
  mtext('f) RMSE (last)',adj=0, line=3, cex=0.9)
#title("Recruitment",outer=TRUE)
dev.off()


## LM LMM ##################################################################

##--MAKE PLOTS--###########################
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

##--RECRUITMENT--###################
pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','lm_error_recruitment.pdf'),height=4.5,width=10)
par(mfrow=c(2,3),mar=c(1,2,1,2),oma=c(6,4,2,2),cex.axis=0.9)
ylims <- c(-0.1,0.1)
plotlm(FITS[['lm_recr_re_all']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
plotlm(FITS[['lme_recr_re_all']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),mean=mean(recr.df$RE),nlme_try=TRUE)#,int=TRUE)
  mtext('a) RE (all)',adj=0.0)
  mtext(expression('RE or'~Delta*'RE'),side=2,line=2.5)

plotlm(FITS[['lm_recr_re_ten']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
plotlm(FITS[['lme_recr_re_ten']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),mean=mean(recr.df$RE[recr.df$Year %in% 31:40]),nlme_try=TRUE)#,int=TRUE)
  mtext('b) RE (ten)',adj=0.0)

plotlm(FITS[['lm_recr_re_last']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
plotlm(FITS[['lme_recr_re_last']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=TRUE,mean=mean(recr.df$RE[recr.df$Year == 40]))#,int=TRUE)
  mtext('c) RE (last)',adj=0.0)

ylims <- c(-15000,30000)
plotlm(FITS[['lm_recr_rmse_all']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
plotlm(FITS[['lme_recr_rmse_all']],add=TRUE,ylim=ylims,labels=labels,nlme_try=TRUE,mean=mean(recr.df$RMSE))#,int=TRUE)
  mtext('d) RMSE (all)',adj=0.0)
  mtext(expression('RMSE or'~Delta*'RMSE'),side=2,line=2.5)

plotlm(FITS[['lm_recr_rmse_ten']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
plotlm(FITS[['lme_recr_rmse_ten']],add=TRUE,ylim=ylims,labels=labels,nlme_try=TRUE,mean=mean(recr.df$RMSE[recr.df$Year %in% 31:40]))#,int=TRUE)
  mtext('e) RMSE (ten)',adj=0.0)

plotlm(FITS[['lm_recr_rmse_last']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
plotlm(FITS[['lme_recr_rmse_last']],add=TRUE,ylim=ylims,labels=labels,nlme_try=TRUE,mean=mean(recr.df$RMSE[recr.df$Year == 40]))#,int=TRUE)
  mtext('f) RMSE (last)',adj=0.0)

dev.off()



##--SSB--##################
pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','lm_error_ssb.pdf'),height=4.5,width=10)
par(mfrow=c(2,3),mar=c(1,2,1,2),oma=c(6,4,2,2),cex.axis=0.9)
ylims <- c(-0.1,0.1)
plotlm(FITS[['lm_ssb_re_all']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
plotlm(FITS[['lme_ssb_re_all']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=TRUE,mean=mean(ssb.df$RE))#,int=TRUE)
  mtext('a) RE (all)',adj=0.0)
  mtext(expression('RE or'~Delta*'RE'),side=2,line=2.5)

plotlm(FITS[['lm_ssb_re_ten']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
plotlm(FITS[['lme_ssb_re_ten']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=TRUE,mean=mean(ssb.df$RE[ssb.df$Year %in% 31:40]))#,int=TRUE)
  mtext('b) RE (ten)',adj=0.0)

plotlm(FITS[['lm_ssb_re_last']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
plotlm(FITS[['lme_ssb_re_last']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=TRUE,mean=mean(ssb.df$RE[ssb.df$Year == 40]))#,int=TRUE)
  mtext('c) RE (last)',adj=0.0)

ylims <- c(-15000,30000)
plotlm(FITS[['lm_ssb_rmse_all']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
plotlm(FITS[['lme_ssb_rmse_all']],add=TRUE,ylim=ylims,labels=labels,nlme_try=TRUE,mean=mean(ssb.df$RMSE))#,int=TRUE)
  mtext('d) RMSE (all)',adj=0.0)
  mtext(expression('RMSE or'~Delta*'RMSE'),side=2,line=2.5)

plotlm(FITS[['lm_ssb_rmse_ten']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
plotlm(FITS[['lme_ssb_rmse_ten']],add=TRUE,ylim=ylims,labels=labels,nlme_try=TRUE,mean=mean(ssb.df$RMSE[ssb.df$Year %in% 31:40]))#,int=TRUE)
  mtext('e) RMSE (ten)',adj=0.0)

plotlm(FITS[['lm_ssb_rmse_last']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
plotlm(FITS[['lme_ssb_rmse_last']],add=TRUE,ylim=ylims,labels=labels,nlme_try=TRUE,mean=mean(ssb.df$RMSE[ssb.df$Year == 40]))#,int=TRUE)
  mtext('f) RMSE (last)',adj=0.0)

dev.off()




##--FBAR--##############################
pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','lm_error_fbar.pdf'),height=4.5,width=10)
par(mfrow=c(2,3),mar=c(1,2,1,2),oma=c(6,4,2,2),cex.axis=0.9)
ylims <- c(-0.1,0.1)
plotlm(FITS[['lm_fbar_re_all']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
plotlm(FITS[['lme_fbar_re_all']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=TRUE,mean=mean(fbar.df$RE))#,int=TRUE)
  mtext('a) RE (all)',adj=0.0)
  mtext(expression('RE or'~Delta*'RE'),side=2,line=2.5)

plotlm(FITS[['lm_fbar_re_ten']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
plotlm(FITS[['lme_fbar_re_ten']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=TRUE,mean=mean(fbar.df$RE[fbar.df$Year %in% 31:40]))#,int=TRUE)
  mtext('b) RE (ten)',adj=0.0)

plotlm(FITS[['lm_fbar_re_last']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
plotlm(FITS[['lme_fbar_re_last']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=TRUE,mean=mean(fbar.df$RE[fbar.df$Year == 40]))#,int=TRUE)
  mtext('c) RE (last)',adj=0.0)

ylims <- c(-0.02,0.03)
plotlm(FITS[['lm_fbar_rmse_all']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
plotlm(FITS[['lme_fbar_rmse_all']],add=TRUE,ylim=ylims,labels=labels,nlme_try=TRUE,mean=mean(fbar.df$RMSE))#,int=TRUE)
  mtext('d) RMSE (all)',adj=0.0)
  mtext(expression('RMSE or'~Delta*'RMSE'),side=2,line=2.5)

plotlm(FITS[['lm_fbar_rmse_ten']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
plotlm(FITS[['lme_fbar_rmse_ten']],add=TRUE,ylim=ylims,labels=labels,nlme_try=TRUE,mean=mean(fbar.df$RMSE[fbar.df$Year %in% 31:40]))#,int=TRUE)
  mtext('e) RMSE (ten)',adj=0.0)

plotlm(FITS[['lm_fbar_rmse_last']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
plotlm(FITS[['lme_fbar_rmse_last']],add=TRUE,ylim=ylims,labels=labels,nlme_try=TRUE,mean=mean(fbar.df$RMSE[fbar.df$Year == 40]))#,int=TRUE)
  mtext('f) RMSE (last)',adj=0.0)

dev.off()






# ##--ONE BIG PLOT--###################
# ##--RECRUITMENT--######################
# pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','lm_error_assessment.pdf'),height=4*3,width=10)
# par(mfrow=c(6,3),mar=c(1,2,1,2),oma=c(6,4,2,2),cex.axis=0.9)
# ylims <- c(-0.1,0.1)
# plotlm(FITS[['lm_recr_re_all']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
# plotlm(FITS[['lme_recr_re_all']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),mean=mean(recr.df$RE),nlme_try=TRUE)#,int=TRUE)
#   mtext('a) RE (all)',adj=0.0)
#   mtext(expression('RE or'~Delta*'RE'),side=2,line=2.5)

# plotlm(FITS[['lm_recr_re_ten']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
# plotlm(FITS[['lme_recr_re_ten']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),mean=mean(recr.df$RE[recr.df$Year %in% 31:40]),nlme_try=TRUE)#,int=TRUE)
#   mtext('b) RE (ten)',adj=0.0)

# plotlm(FITS[['lm_recr_re_last']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
# plotlm(FITS[['lme_recr_re_last']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=TRUE,mean=mean(recr.df$RE[recr.df$Year == 40]))#,int=TRUE)
#   mtext('c) RE (last)',adj=0.0)

# ylims <- c(-15000,30000)
# plotlm(FITS[['lm_recr_rmse_all']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
# plotlm(FITS[['lme_recr_rmse_all']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=TRUE,mean=mean(recr.df$RMSE))#,int=TRUE)
#   mtext('d) RMSE (all)',adj=0.0)
#   mtext(expression('RMSE or'~Delta*'RMSE'),side=2,line=2.5)

# plotlm(FITS[['lm_recr_rmse_ten']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
# plotlm(FITS[['lme_recr_rmse_ten']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=TRUE,mean=mean(recr.df$RMSE[recr.df$Year %in% 31:40]))#,int=TRUE)
#   mtext('e) RMSE (ten)',adj=0.0)

# plotlm(FITS[['lm_recr_rmse_last']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
# plotlm(FITS[['lme_recr_rmse_last']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=TRUE,mean=mean(recr.df$RMSE[recr.df$Year == 40]))#,int=TRUE)
#   mtext('f) RMSE (last)',adj=0.0)
# #dev.off()

# ##--SSB--##################
# #pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','lm_error_ssb.pdf'),height=4.5,width=10)
# #par(mfrow=c(2,3),mar=c(1,2,1,2),oma=c(6,4,2,2),cex.axis=0.9)
# ylims <- c(-0.1,0.1)
# plotlm(FITS[['lm_ssb_re_all']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
# plotlm(FITS[['lme_ssb_re_all']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=TRUE,mean=mean(ssb.df$RE))#,int=TRUE)
#   mtext('a) RE (all)',adj=0.0)
#   mtext(expression('RE or'~Delta*'RE'),side=2,line=2.5)

# plotlm(FITS[['lm_ssb_re_ten']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
# plotlm(FITS[['lme_ssb_re_ten']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=TRUE,mean=mean(ssb.df$RE[ssb.df$Year %in% 31:40]))#,int=TRUE)
#   mtext('b) RE (ten)',adj=0.0)

# plotlm(FITS[['lm_ssb_re_last']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
# plotlm(FITS[['lme_ssb_re_last']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=TRUE,mean=mean(ssb.df$RE[ssb.df$Year == 40]))#,int=TRUE)
#   mtext('c) RE (last)',adj=0.0)

# ylims <- c(-15000,30000)
# plotlm(FITS[['lm_ssb_rmse_all']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
# plotlm(FITS[['lme_ssb_rmse_all']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=TRUE,mean=mean(ssb.df$RMSE))#,int=TRUE)
#   mtext('d) RMSE (all)',adj=0.0)
#   mtext(expression('RMSE or'~Delta*'RMSE'),side=2,line=2.5)

# plotlm(FITS[['lm_ssb_rmse_ten']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
# plotlm(FITS[['lme_ssb_rmse_ten']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=TRUE,mean=mean(ssb.df$RMSE[ssb.df$Year %in% 31:40]))#,int=TRUE)
#   mtext('e) RMSE (ten)',adj=0.0)

# plotlm(FITS[['lm_ssb_rmse_last']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
# plotlm(FITS[['lme_ssb_rmse_last']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=TRUE,mean=mean(ssb.df$RMSE[ssb.df$Year == 40]))#,int=TRUE)
#   mtext('f) RMSE (last)',adj=0.0)
# #dev.off()

# ##--FBAR--##############################
# #pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','lm_error_fbar.pdf'),height=4.5,width=10)
# #par(mfrow=c(2,3),mar=c(1,2,1,2),oma=c(6,4,2,2),cex.axis=0.9)
# ylims <- c(-0.1,0.1)
# plotlm(FITS[['lm_fbar_re_all']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
# plotlm(FITS[['lme_fbar_re_all']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=TRUE,mean=mean(fbar.df$RE))#,int=TRUE)
#   mtext('a) RE (all)',adj=0.0)
#   mtext(expression('RE or'~Delta*'RE'),side=2,line=2.5)

# plotlm(FITS[['lm_fbar_re_ten']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
# plotlm(FITS[['lme_fbar_re_ten']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=TRUE,mean=mean(fbar.df$RE[fbar.df$Year %in% 31:40]))#,int=TRUE)
#   mtext('b) RE (ten)',adj=0.0)

# plotlm(FITS[['lm_fbar_re_last']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
# plotlm(FITS[['lme_fbar_re_last']],add=TRUE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=TRUE,mean=mean(fbar.df$RE[fbar.df$Year == 40]))#,int=TRUE)
#   mtext('c) RE (last)',adj=0.0)

# ylims <- c(-0.02,0.04)
# plotlm(FITS[['lm_fbar_rmse_all']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
# plotlm(FITS[['lme_fbar_rmse_all']],add=TRUE,ylim=ylims,labels=labels,nlme_try=TRUE,mean=mean(fbar.df$RMSE))#,int=TRUE)
#   mtext('d) RMSE (all)',adj=0.0)
#   mtext(expression('RMSE or'~Delta*'RMSE'),side=2,line=2.5)

# plotlm(FITS[['lm_fbar_rmse_ten']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
# plotlm(FITS[['lme_fbar_rmse_ten']],add=TRUE,ylim=ylims,labels=labels,nlme_try=TRUE,mean=mean(fbar.df$RMSE[fbar.df$Year %in% 31:40]))#,int=TRUE)
#   mtext('e) RMSE (ten)',adj=0.0)

# plotlm(FITS[['lm_fbar_rmse_last']],add=FALSE,ylim=ylims,labels=rep(NA,length(labels)),nlme_try=FALSE)#,int=TRUE)
# plotlm(FITS[['lme_fbar_rmse_last']],add=TRUE,ylim=ylims,labels=labels,nlme_try=TRUE,mean=mean(fbar.df$RMSE[fbar.df$Year == 40]))#,int=TRUE)
#   mtext('f) RMSE (last)',adj=0.0)

# dev.off()










