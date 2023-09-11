library("here")
df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.oms.RDS"))
use.df.ems <- df.ems[1:20,]
df.oms$NAA_sig[which(is.na(df.oms$NAA_sig))] <- 0
om_ind <- 1:24


all_naa_sd_log_SSB = lapply(1:NROW(df.oms), function(y){
  res = sapply(1:100, function(x){
    print(paste0("sim_",x,".RDS"))
    sim = readRDS(file.path(here::here(),"Project_0", "results", "naa_om", paste0("om_", y), paste0("sim_",x,".RDS")))
    res = sapply(sim, function(x) {
      out <- NA
      if(length(x$truth)) out <- sd(log(x$truth$SSB))
      return(out)
    })
    return(res[which(!is.na(res))][1])
  })
  return(res)
})


#############################################################
## MAKE DF ##################################################
#############################################################
all_naa_relssb <- readRDS(file = file.path(here(),"Project_0","results", "all_naa_relssb_results.RDS"))
all_naa_relR   <- readRDS(file = file.path(here(),"Project_0","results", "all_naa_relR_results.RDS"))
all_naa_relF   <- readRDS(file = file.path(here(),"Project_0","results", "all_naa_relF_results.RDS"))
make_df_tab <- function(all_naa_relssb, all_naa_relR, all_naa_relF, est_ind = 1:length(df.ems), om_ind = NULL){
  if(!is.null(om_ind)) all_naa_relssb <- all_naa_relssb[om_ind] 
  if(!is.null(om_ind)) all_naa_relR   <- all_naa_relR[om_ind] 
  if(!is.null(om_ind)) all_naa_relF   <- all_naa_relF[om_ind] 
  df_out <- cbind.data.frame(OM = numeric(),
                             R_sig = numeric(), NAA_sig = numeric(), Fhist = character(), obs_error = character(), sd_log_SSB = numeric(),
                             relSSB = numeric(),relR=numeric(),relF=numeric())
  for(i in 1:length(all_naa_relssb)){
    tmpSSB <- lapply(all_naa_relssb[[i]], function(k){ return(colMeans(k))})
    tmpR   <- lapply(all_naa_relR[[i]],   function(k){ return(colMeans(k))})
    tmpF   <- lapply(all_naa_relF[[i]],   function(k){ return(colMeans(k))})
    df_i <- cbind.data.frame(OM = i, R_sig = df.oms$R_sig[om_ind[i]], NAA_sig = df.oms$NAA_sig[om_ind[i]], 
                             Fhist = df.oms$Fhist[om_ind[i]], obs_error = df.oms$obs_error[om_ind[i]],
                             sd_log_SSB = rep(all_naa_sd_log_SSB[[om_ind[i]]],20),
                             SR=rep(df.ems[1:20,1],100),M_est=rep(df.ems[1:20,2],100),re=rep(df.ems[1:20,3],100),
                             relSSB = unlist(tmpSSB), relR=unlist(tmpR), relF=unlist(tmpF))
    df_out <- rbind(df_out, df_i)
  }
  return(df_out)
}

df = make_df_tab(all_naa_relssb,all_naa_relR,all_naa_relF, om_ind=om_ind)



###############
sub_ssb <- filter(df, relSSB>1E6)
  tab_ssb <- apply(sub_ssb[,c(1,2,3,4,5,7,8,9)],2,table)
sub_F <- filter(df, relF<1E-6)
  tab_F <- apply(sub_F[,c(1,2,3,4,5,7,8,9)],2,table)
sub_r <- filter(df, relR>1E6)
  tab_r <- apply(sub_r[,c(1,2,3,4,5,7,8,9)],2,table)

par(mfrow=c(1,3))
hist(log10(df$relSSB),breaks=2000,xlim=c(-1,15),main='',xlab=expression('log'['10']*'(SSBfit/SSBtruth)'))
hist(log10(df$relR),breaks=1000,xlim=c(-1,15),main='',xlab=expression('log'['10']*'(Rfit/Rtruth)'))
hist(log10(df$relF),breaks=2000,xlim=c(-10,3),main='',xlab=expression('log'['10']*'(Ffit/Ftruth)'))

pdf('~/dropbox/tryyy.pdf',height=10,width=8)
par(mfrow=c(6,3),mar=c(2,2,2,2),oma=c(2,2,2,2),cex.axis=0.9)
for(i in 2:7){
  barplot(tab_ssb[[i]])
    mtext(adj=0,names(tab_ssb[i]))
  barplot(tab_r[[i]])
  barplot(tab_F[[i]])
}
dev.off()


####################
df_fit <- filter(df, between(relSSB,0,2) & between(relR,0,2) & between(relF,0,2))

fit_ssb <- lm(relSSB ~ R_sig + NAA_sig + as.factor(Fhist) + as.factor(obs_error) + sd_log_SSB + as.factor(SR) + as.factor(M_est) + as.factor(re), 
              data=df_fit)
summary(fit_ssb)

fit_r <- lm(relR ~ R_sig + NAA_sig + as.factor(Fhist) + as.factor(obs_error) + sd_log_SSB + as.factor(SR) + as.factor(M_est) + as.factor(re), 
              data=df_fit)
summary(fit_r)

fit_f <- lm(relF ~ R_sig + NAA_sig + as.factor(Fhist) + as.factor(obs_error) + sd_log_SSB + as.factor(SR) + as.factor(M_est) + as.factor(re), 
            data=df_fit)
summary(fit_f)


fit_ssb2 <- lm(relSSB ~ (R_sig + NAA_sig + as.factor(Fhist) + as.factor(obs_error) + sd_log_SSB + as.factor(SR) + as.factor(M_est) + as.factor(re))^2, 
              data=df_fit)
summary(fit_ssb2)

fit_r2 <- lm(relR ~ R_sig + NAA_sig + as.factor(Fhist) + as.factor(obs_error) + sd_log_SSB + as.factor(SR) + as.factor(M_est) + as.factor(re), 
            data=df_fit)
summary(fit_r2)

fit_f <- lm(relF ~ R_sig + NAA_sig + as.factor(Fhist) + as.factor(obs_error) + sd_log_SSB + as.factor(SR) + as.factor(M_est) + as.factor(re), 
            data=df_fit)
summary(fit_f)




plotlm <- function(fit){
  coef <- summary(fit)$coefficients[2:12,1]
  ses  <- summary(fit)$coefficients[2:12,2]
  plot(1:11,coef,ylim=c(min(coef-2*ses),max(coef+2*ses)),pch=19,xaxt='n')
  segments(x0=1:11,x1=1:11,y0=coef-2*ses,y1=coef+2*ses)
  abline(h=0,lty=2)
  axis(side=1,at=1:11,labels=row.names(summary(fit)$coefficients[2:12,]),las=2)
}


par(mfrow=c(1,3),mar=c(2,2,2,2),oma=c(8,2,2,2))
plotlm(fit_ssb)
mtext(side=2,line=2.5,expression('Coefficient'))
plotlm(fit_r)
plotlm(fit_f)

