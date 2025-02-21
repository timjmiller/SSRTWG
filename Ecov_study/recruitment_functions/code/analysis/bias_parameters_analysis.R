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
df.ems    <- readRDS(file.path(here::here(),"Ecov_study","recruitment_functions", "inputs", "df.ems.RDS"))
colnames(df.oms)[1] <- 'OM'

for(i in 1:576){
  RMSE_par$OM[RMSE_par$OM==i] = RE_par$OM[RE_par$OM==i] <- paste0('om_',i)
}

RE_par   <- left_join(RE_par,df.oms,by=c('OM'))

##--FILTERING--################
RE_par$Ecov_beta[!is.finite(RE_par$Ecov_beta) | RE_par$Ecov_beta=='NaN'] <- NA

RMSE_par <- left_join(RMSE_par,df.oms,by=c('OM'))

RE_par$recruit_mod_EM=RMSE_par$recruit_mod_EM <- NA
  RMSE_par$recruit_mod_EM[RMSE_par$EM==1] = RE_par$recruit_mod_EM[RE_par$EM==1] <- 2
  RMSE_par$recruit_mod_EM[RMSE_par$EM==2] = RE_par$recruit_mod_EM[RE_par$EM==2] <- 2
  RMSE_par$recruit_mod_EM[RMSE_par$EM==3] = RE_par$recruit_mod_EM[RE_par$EM==3] <- 3
  RMSE_par$recruit_mod_EM[RMSE_par$EM==4] = RE_par$recruit_mod_EM[RE_par$EM==4] <- 3
  RMSE_par$recruit_mod_EM[RMSE_par$EM==5] = RE_par$recruit_mod_EM[RE_par$EM==5] <- 3

RE_par$ecov_how_EM=RMSE_par$ecov_how_EM <- NA
  RMSE_par$ecov_how_EM[RMSE_par$EM==1] = RE_par$ecov_how_EM[RE_par$EM==1] <- 0
  RMSE_par$ecov_how_EM[RMSE_par$EM==2] = RE_par$ecov_how_EM[RE_par$EM==2] <- 1
  RMSE_par$ecov_how_EM[RMSE_par$EM==3] = RE_par$ecov_how_EM[RE_par$EM==3] <- 0
  RMSE_par$ecov_how_EM[RMSE_par$EM==4] = RE_par$ecov_how_EM[RE_par$EM==4] <- 1
  RMSE_par$ecov_how_EM[RMSE_par$EM==5] = RE_par$ecov_how_EM[RE_par$EM==5] <- 2

##-select only best fitting model
AIC_weight             <- readRDS( file.path(here(),'Ecov_study','recruitment_functions',res.dir, paste0("AIC_weight_conv.runs_grad_", bad.grad.label, "_SE_", bad.se.value, plot.suffix, ".RDS")  ) )
RE_par   <- RE_par[AIC_weight$AIC_rank==1,]
RMSE_par <- RMSE_par[AIC_weight$AIC_rank==1,]

RE_par$rec_a   <- RE_par$mean_rec1
RMSE_par$rec_a <- RMSE_par$mean_rec1
RE_par$rec_b   <- RE_par$mean_rec2
RMSE_par$rec_b <- RMSE_par$mean_rec2

RE_par$mean_rec1[RE_par$recruit_mod_EM==3 | RE_par$Ecov_how %in% c(1,2)] = RMSE_par$mean_rec1[RMSE_par$recruit_mod_EM==3 | RE_par$Ecov_how %in% c(1,2)] <- NA  #remove BH from mean_rec1
RE_par$mean_rec2[RE_par$recruit_mod_EM==3 | RE_par$Ecov_how %in% c(1,2)] = RMSE_par$mean_rec2[RMSE_par$recruit_mod_EM==3 | RE_par$Ecov_how %in% c(1,2)] <- NA  #remove BH from mean_rec2

RE_par$rec_a[RE_par$recruit_mod_EM==2 | RE_par$Ecov_how %in% c(1,2)] = RMSE_par$rec_a[RMSE_par$recruit_mod_EM==2 | RE_par$Ecov_how %in% c(1,2)] <- NA  #remove noSR from reca
RE_par$rec_b[RE_par$recruit_mod_EM==2 | RE_par$Ecov_how %in% c(1,2)] = RMSE_par$rec_b[RMSE_par$recruit_mod_EM==2 | RE_par$Ecov_how %in% c(1,2)] <- NA  #remove BH from recb


#RE_par$mean_rec1   <- RE_par$rec_a
#RE_par$mean_rec2   <- RE_par$rec_b
#RMSE_par$mean_rec1 <- RMSE_par$rec_a
#RMSE_par$mean_rec2 <- RMSE_par$rec_b

RE_par$Ecov_beta[RE_par$Ecov_how==0] = NA
RMSE_par$Ecov_beta[RMSE_par$Ecov_how==0] = NA

#RE_par$mean_rec1[RE_par$recruit_mod_EM==3] = RMSE_par$mean_rec1[RMSE_par$recruit_mod_EM==3] <- NA  #remove BH from mean_rec1
#RE_par$mean_rec2[RE_par$recruit_mod_EM==3] = RMSE_par$mean_rec2[RMSE_par$recruit_mod_EM==3] <- NA  #remove BH from mean_rec1


RE_par$Ecov_beta <- RE_par$Ecov_beta + 1


##-HISTOGRAMS-###########
pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','RE_par_hist.pdf'),height=12,width=12)
par(mfrow=c(5,5),mar=c(2,2,2,2),oma=c(2,2,2,2))
for(i in c(8:10,13:29,43:44)){
  x <- RE_par[,i]
  x <- x[x > quantile(x,p=0.01,na.rm=TRUE) & x < quantile(x,p=0.99,na.rm=TRUE)]
  hist(x,xlab='',ylab='',main='',col='white')
  abline(v=mean(x,na.rm=TRUE),lwd=2,col='red',lty=1)
  mtext(colnames(RE_par)[i])
}
dev.off()


#pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','RMSE_par_hist.pdf'),height=12,width=12)
#par(mfrow=c(5,5),mar=c(2,2,2,2),oma=c(2,2,2,2))
#for(i in c(6:10,13:29)){
#  x <- RMSE_par[,i]
#  x <- x[x > quantile(x,p=0.01,na.rm=TRUE) & x < quantile(x,p=0.99,na.rm=TRUE)]
#  hist(x,xlab='',ylab='',main='')
#  mtext(colnames(RMSE_par)[i])
#}
#dev.off()

##--MARGINAL MEDIANS--###############

##--COVARIATES--########
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
RE_par$ecov_slope= RMSE_par$ecov_slope  <- as.factor(case_when(RMSE_par$ecov_slope > mean(RMSE_par$ecov_slope) + sd(RMSE_par$ecov_slope) ~ "+",
                                       RMSE_par$ecov_slope < mean(RMSE_par$ecov_slope) - sd(RMSE_par$ecov_slope) ~ '-',
                                       TRUE ~ '0'))
RE_par$ecov_slope= RMSE_par$ecov_slope  <- factor(RMSE_par$ecov_slope,levels=c("-","0","+"))



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


pdf(file=file.path(here::here(), 'Ecov_study','recruitment_functions','plots','par_marg.pdf'),height=5,width=9)
par(mfrow=c(2,3),mar=c(1,2,1,0),oma=c(6,4,2,2),cex.axis=0.65,cex.lab=0.5)
ylims <- c(-1,2)
dd(RE_par,vars=vars,labels=rep(NA,length(labels)),yvar="rec_a",ylims=ylims,mean=TRUE)
  abline(h=0,lty=2)
  mtext(expression("RE"),side=2,line=2.5)
  mtext(expression("a) rec_a"),adj=0.0,cex=0.8)

ylims <- c(-1,22)
dd(RE_par,vars=vars,labels=rep(NA,length(labels)),yvar="rec_b",ylims=ylims,mean=TRUE)
  abline(h=0,lty=2)
  mtext(expression("b) rec_b"),adj=0.0,cex=0.8)

ylims <- c(-1,6)
dd(RE_par,vars=vars,labels=rep(NA,length(labels)),yvar="NAA_sigma",ylims=ylims,mean=TRUE)
  abline(h=0,lty=2)
  mtext(expression("c) NAA_sigma"),adj=0.0,cex=0.8)

ylims <- c(-0.2,0.1)
dd(RE_par,vars=vars,labels=labels,yvar="NAA_rho",ylims=ylims,mean=TRUE)
  abline(h=0,lty=2)
  mtext(expression("RE"),side=2,line=2.5)
  mtext(expression("d) NAA_rho"),adj=0.0,cex=0.8)

ylims <- c(-0.6,0.6)
dd(RE_par,vars=vars,labels=labels,yvar="Ecov_beta",ylims=ylims,mean=TRUE)
  abline(h=0,lty=2)
  mtext(expression("e) Ecov_beta"),adj=0.0,cex=0.8)

ylims <- c(-0.5,0.5)
dd(RE_par,vars=vars,labels=labels,yvar="Ecov_process_pars3",ylims=ylims,mean=TRUE)
  abline(h=0,lty=2)
  mtext(expression("f) Ecov_process_pars3"),adj=0.1,cex=0.8)

dev.off()













