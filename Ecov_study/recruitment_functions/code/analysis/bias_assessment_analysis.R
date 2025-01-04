library(here)
library(tidyverse)
library(rpart)
library(rpart.plot)
library(nlme)

res.path    <- file.path(here::here(),"Ecov_study", "recruitment_functions", "results")  # directory where simulation 
res.dir     <- 'results'   # 'results'     'results_beta_fix'   # results folder where AIC dataframes are
plot.dir    <- 'plots'    # 'plots_lizruns'  'plots_beta_fix'  
table.dir   <- 'tables'   # 'results'     'results_beta_fix'
plot.suffix <- ''      # '_beta_fix'   '' 

df.oms    <- readRDS(file.path(here::here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS"))
df.ems    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))

bad.grad.value <- 1E-6 #tim used 1E-6 (the abs of the exponent will be used for output suffix; ex: filename_grad_6.png)
bad.grad.label <- as.numeric(strsplit(as.character(bad.grad.value), split="-")[[1]][2])
bad.se.value   <- 100 #tim used 100 (this value will be used for output suffix; ex: filename_se_100.png)

recr.df <- as.data.frame(readRDS( file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "error.recr", plot.suffix, ".df.RDS") ) )) %>%
  filter(opt==1 & conv==0 & sdrep==0 & max_grad<bad.grad.value)

ssb.df  <- as.data.frame(readRDS( file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "error.ssb",  plot.suffix, ".df.RDS") )  )) %>%
  filter(opt==1 & conv==0 & sdrep==0 & max_grad<bad.grad.value)

fbar.df <- as.data.frame(readRDS( file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "error.fbar", plot.suffix, ".df.RDS") )  )) %>%
  filter(opt==1 & conv==0 & sdrep==0 & max_grad<bad.grad.value)

recr.df$ecov_slope=ssb.df$ecov_slope=fbar.df$ecov_slope <- revalue(recr.df$ecov_slope, c("L" = "-", "M" = "0", "H" = "+"))

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

#mars <- c(6,3,3,3)

pdf(file=file.path(here::here(), 'Ecov_study','recruitment_functions','plots','assessment.marg.pdf'),height=10,width=9)
par(mfrow=c(4,3),mar=c(1,2,1,0),oma=c(6,4,2,2),cex.axis=0.8,cex.lab=0.8)
ylims <- c(-0.1,0.1)
dd(recr.df,vars=vars,labels=rep(NA,length(labels)),yvar="RE",ylims=ylims)
abline(h=0,lty=2)
mtext('a) Recruitment (all)',adj=0)
mtext('RE',side=2,line=2.5)

dd(ssb.df,vars=vars,labels=rep(NA,length(labels)),yvar="RE",ylims=ylims)
  mtext(expression(""),side=2,line=2.5)
abline(h=0,lty=2)
mtext('b) SSB (all)',adj=0)

dd(fbar.df,vars=vars,labels=rep(NA,length(labels)),yvar="RE",ylims=ylims)
  mtext(expression(""),side=2,line=2.5)
abline(h=0,lty=2)
mtext('c) F (all)',adj=0)

dd(recr.df[recr.df$Year==40,],vars=vars,labels=rep(NA,length(labels)),yvar="RE",ylims=ylims)
abline(h=0,lty=2)
mtext('d) Recruitment (terminal)',adj=0)
mtext('RE',side=2,line=2.5)

dd(ssb.df[ssb.df$Year==40,],vars=vars,labels=rep(NA,length(labels)),yvar="RE",ylims=ylims)
  mtext(expression(""),side=2,line=2.5)
abline(h=0,lty=2)
mtext('e) SSB (terminal)',adj=0)

dd(fbar.df[fbar.df$Year==40,],vars=vars,labels=rep(NA,length(labels)),yvar="RE",ylims=ylims)
  mtext(expression(""),side=2,line=2.5)
abline(h=0,lty=2)
mtext('f) F (terminal)',adj=0)



ylims <- c(0,1E4)
dd(recr.df,vars=vars,labels=rep(NA,length(labels)),yvar="RMSE",ylims=ylims)
  mtext(expression(""),side=2,line=2.5)
#mtext('a) Recruitment',adj=0)
mtext('RMSE',side=2,line=2.5)
mtext('g) Recruitment (all)',adj=0)

ylims <- c(0,1E4)
dd(ssb.df,vars=vars,labels=rep(NA,length(labels)),yvar="RMSE",ylims=ylims)
  mtext(expression(""),side=2,line=2.5)
#mtext('a) Recruitment',adj=0)
mtext('h) SSB (all)',adj=0)

ylims <- c(0,0.1)
dd(fbar.df,vars=vars,labels=rep(NA,length(labels)),yvar="RMSE",ylims=ylims)
  mtext(expression(""),side=2,line=2.5)
#mtext('a) Recruitment',adj=0)
mtext('i) F (all)',adj=0)

ylims <- c(0,1E4)
dd(recr.df[recr.df$Year==40,],vars=vars,labels=labels,yvar="RMSE",ylims=ylims)
  mtext(expression(""),side=2,line=2.5)
#mtext('a) Recruitment',adj=0)
mtext('RMSE',side=2,line=2.5)
mtext('j) Recruitment (terminal)',adj=0)

ylims <- c(0,1E4)
dd(ssb.df[ssb.df$Year==40,],vars=vars,labels=labels,yvar="RMSE",ylims=ylims)
  mtext(expression(""),side=2,line=2.5)
#mtext('a) Recruitment',adj=0)
mtext('k) SSB (terminal)',adj=0)

ylims <- c(0,0.1)
dd(fbar.df[fbar.df$Year==40,],vars=vars,labels=labels,yvar="RMSE",ylims=ylims)
  mtext(expression(""),side=2,line=2.5)
#mtext('a) Recruitment',adj=0)
mtext('l) F (terminal)',adj=0)

dev.off()



