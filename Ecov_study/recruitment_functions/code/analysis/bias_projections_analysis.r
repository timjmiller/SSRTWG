library(here)
library(tidyverse)
library(wham)
library(rpart)
library(rpart.plot)
source(file.path(here::here(),'Ecov_study','recruitment_functions','code','analysis','functions_for_analysis.R'))

df.recr.proj <- as.data.frame(readRDS(file.path(here::here(), 'Ecov_study','recruitment_functions','results','df.recr.proj.RDS')))
df.ssb.proj <- as.data.frame(readRDS(file.path(here::here(), 'Ecov_study','recruitment_functions','results','df.ssb.proj.RDS')))
df.catch.proj <- as.data.frame(readRDS(file.path(here::here(), 'Ecov_study','recruitment_functions','results','df.catch.proj.RDS')))

df.oms    <- readRDS(file.path(here::here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS"))
colnames(df.oms)[1] <- 'OM'

for(i in 1:576){
  df.recr.proj$OM[df.recr.proj$OM==i] = 
  df.ssb.proj$OM[df.ssb.proj$OM==i] =
  df.catch.proj$OM[df.catch.proj$OM==i] <- paste0('om_',i)
}

##--JOIN--###########
df.recr.proj  <- left_join(df.recr.proj, df.oms, by=c("OM")) %>%
    filter(Ecov_how!=0 & EM!=1) 
df.ssb.proj   <- left_join(df.ssb.proj, df.oms, by=c("OM")) %>%
    filter(Ecov_how!=0 & EM!=1) 
df.catch.proj <- left_join(df.catch.proj, df.oms, by=c("OM")) %>%
    filter(Ecov_how!=0 & EM!=1) 

##--COVARIATES--########
df.recr.proj$obs_error=  df.ssb.proj$obs_error=  df.catch.proj$obs_error      <- factor(df.ssb.proj$obs_error,levels=c("L","H"))
df.recr.proj$R_sig=      df.ssb.proj$R_sig=      df.catch.proj$R_sig          <- as.factor(df.ssb.proj$R_sig)
df.recr.proj$Fhist=      df.ssb.proj$Fhist=      df.catch.proj$RMSE_par$Fhist <- factor(df.ssb.proj$Fhist,levels=c("MSY","L-H","H-MSY"))
df.recr.proj$NAA_cor=    df.ssb.proj$NAA_cor=    df.catch.proj$NAA_cor        <- as.factor(df.ssb.proj$NAA_cor) 
df.recr.proj$Ecov_re_cor=df.ssb.proj$Ecov_re_cor=df.catch.proj$Ecov_re_cor    <- as.factor(df.ssb.proj$Ecov_re_cor) 
df.recr.proj$Ecov_effect=df.ssb.proj$Ecov_effect=df.catch.proj$Ecov_effect    <- as.factor(df.ssb.proj$Ecov_effect) 
df.recr.proj$Ecov_how=   df.ssb.proj$Ecov_how=   df.catch.proj$Ecov_how       <- as.factor(df.ssb.proj$Ecov_how) 
df.recr.proj$ssb_cv=     df.ssb.proj$ssb_cv=     df.catch.proj$ssb_cv         <- factor(case_when(df.ssb.proj$ssb_cv < mean(df.ssb.proj$ssb_cv) - sd(df.ssb.proj$ssb_cv) ~ 'L',
                                    df.ssb.proj$ssb_cv > mean(df.ssb.proj$ssb_cv) + sd(df.ssb.proj$ssb_cv) ~ "H",
                                    TRUE ~ 'M')) 
df.recr.proj$ssb_cv=     df.ssb.proj$ssb_cv=      df.catch.proj$ssb_cv        <- factor(df.ssb.proj$ssb_cv,levels=c("L","M","H"))
df.recr.proj$ecov_slope= df.ssb.proj$ecov_slope=  df.catch.proj$ecov_slope    <- as.factor(case_when(df.ssb.proj$ecov_slope > mean(df.ssb.proj$ecov_slope) + sd(df.ssb.proj$ecov_slope) ~ "+",
                                       df.ssb.proj$ecov_slope < mean(df.ssb.proj$ecov_slope) - sd(df.ssb.proj$ecov_slope) ~ '-',
                                       TRUE ~ '0'))
df.recr.proj$ecov_slope= df.ssb.proj$ecov_slope= df.catch.proj$ecov_slope     <- factor(df.ssb.proj$ecov_slope,levels=c("-","0","+"))

##--OUTLIERS--##########
df.recr.proj <- filter(df.recr.proj, re.cont.ecov <2 & re.avg.ecov <2 & re.use.ecov <2)
df.ssb.proj <- filter(df.ssb.proj, re.cont.ecov <2 & re.avg.ecov <2 & re.use.ecov <2)
df.catch.proj <- filter(df.catch.proj, re.cont.ecov <2 & re.avg.ecov <2 & re.use.ecov <2)

##--SUMMARIZE--#########
names <- factor(c('re.cont.ecov','re.avg.ecov','re.use.ecov','rmse.cont.ecov','rmse.avg.ecov','rmse.use.ecov'))

X <- array(NA,dim=c(10,6,3))
X[,1,1] <- pull(group_by(df.recr.proj,Year) %>% summarise(mean=mean(re.cont.ecov)),mean)
X[,2,1] <- pull(group_by(df.recr.proj,Year) %>% summarise(mean=mean(re.avg.ecov)),mean)
X[,3,1] <- pull(group_by(df.recr.proj,Year) %>% summarise(mean=mean(re.use.ecov)),mean)
X[,4,1] <- pull(group_by(df.recr.proj,Year) %>% summarise(mean=mean(rmse.cont.ecov)),mean)
X[,5,1] <- pull(group_by(df.recr.proj,Year) %>% summarise(mean=mean(rmse.avg.ecov)),mean)
X[,6,1] <- pull(group_by(df.recr.proj,Year) %>% summarise(mean=mean(rmse.use.ecov)),mean)
X[,1,2] <- pull(group_by(df.ssb.proj,Year) %>% summarise(mean=mean(re.cont.ecov)),mean)
X[,2,2] <- pull(group_by(df.ssb.proj,Year) %>% summarise(mean=mean(re.avg.ecov)),mean)
X[,3,2] <- pull(group_by(df.ssb.proj,Year) %>% summarise(mean=mean(re.use.ecov)),mean)
X[,4,2] <- pull(group_by(df.ssb.proj,Year) %>% summarise(mean=mean(rmse.cont.ecov)),mean)
X[,5,2] <- pull(group_by(df.ssb.proj,Year) %>% summarise(mean=mean(rmse.avg.ecov)),mean)
X[,6,2] <- pull(group_by(df.ssb.proj,Year) %>% summarise(mean=mean(rmse.use.ecov)),mean)
X[,1,3] <- pull(group_by(df.catch.proj,Year) %>% summarise(mean=mean(re.cont.ecov)),mean)
X[,2,3] <- pull(group_by(df.catch.proj,Year) %>% summarise(mean=mean(re.avg.ecov)),mean)
X[,3,3] <- pull(group_by(df.catch.proj,Year) %>% summarise(mean=mean(re.use.ecov)),mean)
X[,4,3] <- pull(group_by(df.catch.proj,Year) %>% summarise(mean=mean(rmse.cont.ecov)),mean)
X[,5,3] <- pull(group_by(df.catch.proj,Year) %>% summarise(mean=mean(rmse.avg.ecov)),mean)
X[,6,3] <- pull(group_by(df.catch.proj,Year) %>% summarise(mean=mean(rmse.use.ecov)),mean)


##--SUMMARY TIME SERIES PLOTS--###################
pdf(file=file.path(here::here(),'Ecov_study','recruitment_functions','plots','proj_time_series_error.pdf'),height=6,width=7)
par(mfrow=c(3,2),mar=c(2,4,2,2),oma=c(2,2,2,2),cex.axis=0.8)

plot(-999,xlim=c(0,10),ylim=c(-0.1,0.1),xlab='',ylab='')
    abline(h=0)
    for(i in 1:3) lines(X[,i,1],col='black',lty=i)
legend('topright',legend=c("Constant","Average","True"),lty=1:3,bty='n')
plot(-999,xlim=c(0,10),ylim=c(10000,30000),xlab='',ylab='')
for(i in 4:6) lines(X[,i,1],col='blue',lty=i-3)

plot(-999,xlim=c(0,10),ylim=c(-0.1,0.1),xlab='',ylab='')
    abline(h=0)
    for(i in 1:3) lines(X[,i,2],col='black',lty=i)
mtext(side=2,line=2.5,'Relative Error (RE)')
plot(-999,xlim=c(0,10),ylim=c(10000,100000),xlab='',ylab='')
for(i in 4:6) lines(X[,i,2],col='blue',lty=i-3)
mtext(side=2,line=2.5,'Root Mean Squared Error (RMSE)')

plot(-999,xlim=c(0,10),ylim=c(-0.1,0.1),xlab='',ylab='')
    abline(h=0)
    for(i in 1:3) lines(X[,i,3],col='black',lty=i)
mtext(side=1,line=2.5,'Forecast Years')
plot(-999,xlim=c(0,10),ylim=c(5000,20000),xlab='',ylab='')
for(i in 4:6) lines(X[,i,3],col='blue',lty=i-3)
mtext(side=1,line=2.5,'Forecast Years')

dev.off()


##--COVARIATE MODELS--##########################
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
            expression(italic('f(E'['cov']*')')~'= 1'), 
            expression(italic('f(E'['cov']*')')~'= 2'), 
            expression(Delta*'E'['cov']*' = -'), 
            expression(Delta*'E'['cov']*' = 0'), 
            expression(Delta*'E'['cov']*' = +'), 
            expression(italic('CV'['SSB']~'= L')),
            expression(italic('CV'['SSB']~'= M')),
            expression(italic('CV'['SSB']~'= H'))
)


pdf(file=file.path(here::here(),'Ecov_study','recruitment_functions','plots','proj_use_ecov_marg.pdf'),height=4.25,width=8)
par(mfrow=c(2,3),mar=c(1,2,1,0),oma=c(6,4,2,2),cex.axis=0.7,cex.lab=0.7)
ylims <- c(-0.5,0.5)
dd(df.recr.proj,vars=vars,labels=rep(NA,length(labels)),yvar="re.use.ecov",ylims=ylims)
  abline(h=0,lty=2)
  mtext(expression('a) Recruitment'),adj=0)
  mtext(expression('RE'),side=2,line=2.5)

ylims <- c(-0.15,0.15)
dd(df.ssb.proj,vars=vars,labels=rep(NA,length(labels)),yvar="re.use.ecov",ylims=ylims)
  abline(h=0,lty=2)
  mtext(expression('b) SSB'),adj=0)

ylims <- c(-0.1,0.1)
dd(df.catch.proj,vars=vars,labels=rep(NA,length(labels)),yvar="re.use.ecov",ylims=ylims)
  abline(h=0,lty=2)
  mtext(expression('c) Catch'),adj=0)

ylims <- c(0,2E4)
dd(df.recr.proj,vars=vars,labels=labels,yvar="rmse.use.ecov",ylims=ylims)
  mtext(expression('d) Recruitment'),adj=0)
  mtext(expression('RMSE'),side=2,line=2.5)

ylims <- c(0,6E4)
dd(df.ssb.proj,vars=vars,labels=labels,yvar="rmse.use.ecov",ylims=ylims)
  mtext(expression('e) SSB'),adj=0)

ylims <- c(0,1.1E4)
dd(df.catch.proj,vars=vars,labels=labels,yvar="rmse.use.ecov",ylims=ylims)
  mtext(expression('f) Catch'),adj=0)

dev.off()






############################################################################

pdf(file=file.path(here::here(),'Ecov_study','recruitment_functions','plots','proj_recr_marg.pdf'),height=4.25,width=8)
par(mfrow=c(2,3),mar=c(1,2,1,0),oma=c(6,4,2,2),cex.axis=0.7,cex.lab=0.7)
ylims <- c(-0.5,0.5)
dd(df.recr.proj,vars=vars,labels=rep(NA,length(labels)),yvar="re.cont.ecov",ylims=ylims)
  abline(h=0,lty=2)
  mtext(expression('RE'),side=2,line=2.5)
  mtext(expression('a) Constant Ecov'),adj=0)
dd(df.recr.proj,vars=vars,labels=rep(NA,length(labels)),yvar="re.avg.ecov",ylims=ylims)
  abline(h=0,lty=2)
  mtext(expression('b) Average Ecov'),adj=0)
dd(df.recr.proj,vars=vars,labels=rep(NA,length(labels)),yvar="re.use.ecov",ylims=ylims)
  abline(h=0,lty=2)
  mtext(expression('c) True Ecov'),adj=0)
ylims <- c(0,2E4)
dd(df.recr.proj,vars=vars,labels=labels,yvar="rmse.cont.ecov",ylims=ylims)
  mtext(expression('RMSE'),side=2,line=2.5)
dd(df.recr.proj,vars=vars,labels=labels,yvar="rmse.avg.ecov",ylims=ylims)
dd(df.recr.proj,vars=vars,labels=labels,yvar="rmse.use.ecov",ylims=ylims)
abline(h=0,lty=2)
dev.off()


pdf(file=file.path(here::here(),'Ecov_study','recruitment_functions','plots','proj_ssb_marg.pdf'),height=4.25,width=8)
par(mfrow=c(2,3),mar=c(1,2,1,0),oma=c(6,4,2,2),cex.axis=0.7,cex.lab=0.7)
ylims <- c(-0.15,0.15)
dd(df.ssb.proj,vars=vars,labels=rep(NA,length(labels)),yvar="re.cont.ecov",ylims=ylims)
  abline(h=0,lty=2)
  mtext(expression('RE'),side=2,line=2.5)
  mtext(expression('a) Constant Ecov'),adj=0)
dd(df.ssb.proj,vars=vars,labels=rep(NA,length(labels)),yvar="re.avg.ecov",ylims=ylims)
  abline(h=0,lty=2)
  mtext(expression('b) Average Ecov'),adj=0)
dd(df.ssb.proj,vars=vars,labels=rep(NA,length(labels)),yvar="re.use.ecov",ylims=ylims)
  abline(h=0,lty=2)
  mtext(expression('c) True Ecov'),adj=0)
ylims <- c(0,5E4)
dd(df.ssb.proj,vars=vars,labels=labels,yvar="rmse.cont.ecov",ylims=ylims)
  mtext(expression('RMSE'),side=2,line=2.5)
dd(df.ssb.proj,vars=vars,labels=labels,yvar="rmse.avg.ecov",ylims=ylims)
dd(df.ssb.proj,vars=vars,labels=labels,yvar="rmse.use.ecov",ylims=ylims)
abline(h=0,lty=2)
dev.off()


pdf(file=file.path(here::here(),'Ecov_study','recruitment_functions','plots','proj_catch_marg.pdf'),height=4.25,width=8)
par(mfrow=c(2,3),mar=c(1,2,1,0),oma=c(6,4,2,2),cex.axis=0.7,cex.lab=0.7)
ylims <- c(-0.1,0.1)
dd(df.catch.proj,vars=vars,labels=rep(NA,length(labels)),yvar="re.cont.ecov",ylims=ylims)
  abline(h=0,lty=2)
  mtext(expression('RE'),side=2,line=2.5)
  mtext(expression('a) Constant Ecov'),adj=0)
dd(df.catch.proj,vars=vars,labels=rep(NA,length(labels)),yvar="re.avg.ecov",ylims=ylims)
  abline(h=0,lty=2)
  mtext(expression('b) Average Ecov'),adj=0)
dd(df.catch.proj,vars=vars,labels=rep(NA,length(labels)),yvar="re.avg.ecov",ylims=ylims)
  abline(h=0,lty=2)
  mtext(expression('c) True Ecov'),adj=0)
ylims <- c(0,1.1E4)
dd(df.catch.proj,vars=vars,labels=labels,yvar="rmse.cont.ecov",ylims=ylims)
  mtext(expression('RMSE'),side=2,line=2.5)
dd(df.catch.proj,vars=vars,labels=labels,yvar="rmse.avg.ecov",ylims=ylims)
dd(df.catch.proj,vars=vars,labels=labels,yvar="rmse.avg.ecov",ylims=ylims)
dev.off()

cp <- 1E-8
maxdepth <- 2

rf.re.recr  <- rpart(re.use.ecov ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                data=df.recr.proj, control=rpart.control(cp=cp,maxdepth=maxdepth))
rf.re.ssb   <- rpart(re.use.ecov ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                data=df.ssb.proj, control=rpart.control(cp=cp,maxdepth=maxdepth))
rf.re.catch <- rpart(re.use.ecov ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                data=df.catch.proj, control=rpart.control(cp=cp,maxdepth=maxdepth))

rf.rmse.recr  <- rpart(rmse.use.ecov ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                data=df.recr.proj, control=rpart.control(cp=cp,maxdepth=maxdepth))
rf.rmse.ssb   <- rpart(rmse.use.ecov ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                data=df.ssb.proj, control=rpart.control(cp=cp,maxdepth=maxdepth))
rf.rmse.catch <- rpart(rmse.use.ecov ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                data=df.catch.proj, control=rpart.control(cp=cp,maxdepth=maxdepth))

pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','trees_proj.pdf'),height=5,width=8)
par(mfrow=c(2,3), oma=c(5,0,0,0))
prp(rf.re.recr,yesno=FALSE,type=4,clip.right.labs=TRUE)
  mtext(expression('a) Recruitment RE'),adj=0,line=2.5, cex=0.9)
prp(rf.re.ssb,yesno=FALSE,type=4,clip.right.labs=TRUE)
  mtext(expression('b) SSB RE'),adj=0,line=2.5, cex=0.9)
prp(rf.re.catch,yesno=FALSE,type=4,clip.right.labs=TRUE)
  mtext(expression('c) Catch RE'),adj=0,line=2.5, cex=0.9)

prp(rf.rmse.recr,yesno=FALSE,type=4,clip.right.labs=TRUE)
  mtext(expression('d) Recruitment RMSE'),adj=0,line=2.5, cex=0.9)
prp(rf.rmse.ssb,yesno=FALSE,type=4,clip.right.labs=TRUE)
  mtext(expression('e) SSN RMSE'),adj=0,line=2.5, cex=0.9)
prp(rf.rmse.catch,yesno=FALSE,type=4,clip.right.labs=TRUE)
  mtext(expression('f) Catch RMSE'),adj=0,line=2.5, cex=0.9)

dev.off()








