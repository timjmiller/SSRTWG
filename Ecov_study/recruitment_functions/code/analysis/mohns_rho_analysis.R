library(tidyverse)
source(file.path(here::here(),'Ecov_study','recruitment_functions','code','analysis','functions_for_analysis.R'))

mrho.df <- as.data.frame(readRDS(file.path(here::here(), 'Ecov_study','recruitment_functions','results','mrho.df.RDS')))
df.oms    <- readRDS(file.path(here::here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS"))
colnames(df.oms)[1] <- 'OM'

for(i in 1:576){
  mrho.df$OM[mrho.df$OM==i] <- paste0('om_',i)
}

mrho.df  <- left_join(mrho.df, df.oms, by=c("OM"))


mrho.df$obs_error      <- factor(mrho.df$obs_error,levels=c("L","H"))
mrho.df$Fhist          <- factor(mrho.df$Fhist,levels=c("MSY","L-H","H-MSY"))
mrho.df$ssb_cv <- factor(case_when(mrho.df$ssb_cv < mean(mrho.df$ssb_cv) - sd(mrho.df$ssb_cv) ~ "L",
                                   mrho.df$ssb_cv > mean(mrho.df$ssb_cv) + sd(mrho.df$ssb_cv) ~ "H",
                                   TRUE ~ "M")) 
mrho.df$ssb_cv <- factor(mrho.df$ssb_cv,levels=c("L","M","H"))
mrho.df$ecov_slope <- as.factor(case_when(mrho.df$ecov_slope < mean(mrho.df$ecov_slope) - sd(mrho.df$ecov_slope) ~ "-",
                                          mrho.df$ecov_slope > mean(mrho.df$ecov_slope) + sd(mrho.df$ecov_slope) ~ '+',
                                          TRUE ~ '0'))
mrho.df$ecov_slope <- factor(mrho.df$ecov_slope,levels=c("-","0","+"))

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

pdf(file=file.path(here::here(), 'Ecov_study','recruitment_functions','plots','mrho.all.pdf'),height=6,width=4.5)
par(mfrow=c(3,1),mar=c(1,2,1,0),oma=c(6,4,2,2),cex.axis=0.6,cex.lab=0.6)
dd(mrho.df,vars=vars,labels=rep(NA,length(labels)),yvar="n1",ylims=c(-0.05,0.05))
  mtext(expression("Mohn's"~rho),side=2,line=2.5)
abline(h=0,lty=2)
mtext('a) Recruitment',adj=0)

dd(mrho.df,vars=vars,labels=rep(NA,length(labels)),yvar="ssb",ylims=c(-0.05,0.05))
  mtext(expression("Mohn's"~rho),side=2,line=2.5)
abline(h=0,lty=2)
mtext('b) Spawning Stock Biomass',adj=0)

dd(mrho.df,vars=vars,labels=labels,yvar="fbar",ylims=c(-0.05,0.05))
  mtext(expression("Mohn's"~rho),side=2,line=2.5)
abline(h=0,lty=2)
mtext('c) F',adj=0)

dev.off()

