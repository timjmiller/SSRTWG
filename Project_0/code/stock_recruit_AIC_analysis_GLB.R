library("here")
df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))

#first NAA operating models
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.oms.RDS"))
use.df.ems <- df.ems[1:20,]

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
#om_ind <- which(!is.na(df.oms$NAA_sig))
om_ind <- 1:24
#SR_rec_Mest = which(use.df.ems$re_config == "rec+1" & use.df.ems$M_est == TRUE)
SR_rec_Mest = which(use.df.ems$re_config %in% c("rec", "rec+1") & use.df.ems$M_est == TRUE)
all_naa_aic <- readRDS(file = file.path(here(),"Project_0","results", "all_naa_aic_results.RDS"))
df.oms$NAA_sig[which(is.na(df.oms$NAA_sig))] <- 0
#make a data frame with columns defining operating model characteristics and SD of true (log) SSB 
# and 1/0 for which SR assumption had lowest AIC

make_df <- function(all_aic, est_ind = SR_rec_Mest, om_ind = NULL){
  if(!is.null(om_ind)) all_aic <- all_aic[om_ind] 
  df_out <- cbind.data.frame(OM = numeric(),
    R_sig = numeric(), NAA_sig = numeric(), Fhist = character(), obs_error = character(), sd_log_SSB = numeric(),
    est_SR = numeric(), est_R = numeric())
  for(i in 1:length(all_aic)){
    #est_ind has No SR assumption first
    if(df.oms$NAA_sig[om_ind[i]] == 0) {  
        est_ind <- which(use.df.ems$re_config %in% c("rec") & use.df.ems$M_est == TRUE)
    } else est_ind <- which(use.df.ems$re_config %in% c("rec+1") & use.df.ems$M_est == TRUE)
    tmp = t(apply(all_aic[[i]][rev(est_ind),],2, function(x) {
      if(any(!is.na(x))) {
        return(x == min(x,na.rm=T))
      } else return(rep(NA, length(x)))
    })) #tmp is rows of T and F for which model had lowest AIC
    df_i <- cbind.data.frame(OM = i, R_sig = df.oms$R_sig[om_ind[i]], NAA_sig = df.oms$NAA_sig[om_ind[i]], 
        Fhist = df.oms$Fhist[om_ind[i]], obs_error = df.oms$obs_error[om_ind[i]],
      sd_log_SSB = all_naa_sd_log_SSB[[om_ind[i]]], est_SR = as.numeric(tmp[,1]), est_R = as.numeric(tmp[,2]))
    df_out <- rbind(df_out, df_i)
  }
  return(df_out)
}
temp = make_df(all_naa_aic, SR_rec_Mest, om_ind)
temp <- subset(temp, !is.na(est_SR) & !is.na(est_R))
fit0 <- glm(cbind(est_SR,est_R) ~  1, family = binomial, data = temp)
fit1 <-glm(cbind(est_SR,est_R) ~ factor(R_sig), family = binomial, data = temp)
fit2 <-glm(cbind(est_SR,est_R) ~ factor(NAA_sig), family = binomial, data = temp)
fit3 <-glm(cbind(est_SR,est_R) ~ Fhist, family = binomial, data = temp)
fit4 <-glm(cbind(est_SR,est_R) ~ obs_error, family = binomial, data = temp)
fit5 <-glm(cbind(est_SR,est_R) ~ sd_log_SSB, family = binomial, data = temp)
AIC(fit0,fit1,fit2,fit3,fit4,fit5)
fit6 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + factor(NAA_sig), family = binomial, data = temp)
fit7 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + Fhist, family = binomial, data = temp)
fit8 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + obs_error, family = binomial, data = temp)
fit9 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + sd_log_SSB, family = binomial, data = temp)
AIC(fit1,fit6,fit7,fit8,fit9)
fit10 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + sd_log_SSB + factor(NAA_sig), family = binomial, data = temp)
fit11 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + sd_log_SSB + Fhist, family = binomial, data = temp)
fit12 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + sd_log_SSB + obs_error, family = binomial, data = temp)
AIC(fit9,fit10,fit11,fit12)
fit13 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + sd_log_SSB + factor(NAA_sig) + obs_error, family = binomial, data = temp)
fit14 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + sd_log_SSB + factor(NAA_sig) + Fhist, family = binomial, data = temp)
AIC(fit10,fit13,fit14)
fit15 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + sd_log_SSB + factor(NAA_sig) + obs_error + Fhist, family = binomial, data = temp)
AIC(fit13,fit15)
fit16 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) * sd_log_SSB + factor(NAA_sig) + obs_error, family = binomial, data = temp)
fit17 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + factor(NAA_sig) * sd_log_SSB + obs_error, family = binomial, data = temp)
fit18 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + factor(NAA_sig) + obs_error * sd_log_SSB, family = binomial, data = temp)
AIC(fit13,fit16,fit17,fit18)
fit19 <-glm(cbind(est_SR,est_R) ~ (factor(R_sig) + factor(NAA_sig)) * sd_log_SSB + obs_error, family = binomial, data = temp)
fit20 <-glm(cbind(est_SR,est_R) ~ (factor(R_sig) + obs_error) * sd_log_SSB + factor(NAA_sig), family = binomial, data = temp)
AIC(fit16,fit19,fit20)

summary(fit16)
library(statmod)
qres <- qresiduals(fit16)
qqnorm(qres)
qqline(qres)
summary(qres)

df.oms[om_ind,]
qres <- qresiduals(fit0)
qqnorm(qres)
qqline(qres)

x <- cbind(temp, pred= predict(fit12, type = "response"))
library(ggplot2)
plt <- ggplot(x, aes(x = sd_log_SSB, y = pred, color = obs_error)) + 
    geom_line() + facet_grid(R_sig ~ NAA_sig, labeller = label_both)
plt <- plt + xlab("SD(log(SSB))") + ylab("P(B-H assumption best)") + 
  labs(color = "Observation\n Error") +
  scale_color_discrete(labels = c("High", "Low")) +
  theme_bw()
  #scale_linetype_discrete(labels = c("2.5 Fmsy -> Fmsy", "Fmsy"))
ggsave(here("Project_0", "paper", "pred_SR_best_NAA_oms.png"), plt)


ssbfit <- lm(sd_log_SSB ~ factor(R_sig) + factor(NAA_sig) + obs_error + Fhist, data=temp)



##--DECISION TREE--#########################
library(rpart)
library(rpart.plot)
#fit <- rpart(est_SR ~ NAA_sig + Fhist + obs_error + sd_log_SSB, data=temp, control=rpart.control(cp=0.006))
fit <- rpart(est_SR ~ NAA_sig + Fhist + obs_error, data=temp, control=rpart.control(cp=0.001))
rpart.plot(fit)









#############################################################
##-M operating models-#######################################
#############################################################
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.M.oms.RDS"))
use.df.ems = df.ems[5:24,] #only use 
use.df.ems$re_config[c(5:8,17:20)] = paste0("M_",c("iid","ar1")[match(use.df.ems$M_re_cor[c(5:8,17:20)], c("iid","ar1_y"))])
use.df.ems = use.df.ems[c(5:8,17:20,1:4,9:16),]

all_M_sd_log_SSB = lapply(1:NROW(df.oms), function(y){
  res = sapply(1:100, function(x){
    print(paste0("sim_",x,".RDS"))
    sim = readRDS(file.path(here::here(),"Project_0", "results", "M_om", paste0("om_", y), paste0("sim_",x,".RDS")))
    res = sapply(sim, function(x) {
      out <- NA
      if(length(x$truth)) out <- sd(log(x$truth$SSB))
      return(out)
    })
    return(res[which(!is.na(res))][1])
  })
  return(res)
})
all_M_sd_log_SSB[[1]]
sim = readRDS(file.path(here::here(),"Project_0", "results", "M_om", paste0("om_", 8), paste0("sim_",23,".RDS")))
om_ind <- 1:NROW(df.oms)
#SR_rec_Mest = which(use.df.ems$re_config == "rec+1" & use.df.ems$M_est == TRUE)
SR_rec_Mest = which(use.df.ems$re_config %in% c("M_iid", "M_ar1") & use.df.ems$M_est == TRUE)
all_M_aic <- readRDS(file = file.path(here(),"Project_0","results", "all_M_aic_results.RDS"))
#make a data frame with columns defining operating model characteristics and SD of true (log) SSB 
# and 1/0 for which SR assumption had lowest AIC
make_df <- function(all_aic, est_ind = SR_rec_Mest, om_ind = NULL){
  if(!is.null(om_ind)) all_aic <- all_aic[om_ind]
  df_out <- cbind.data.frame(OM = numeric(),
    R_sig = numeric(), NAA_sig = numeric(), Fhist = character(), obs_error = character(), sd_log_SSB = numeric(),
    est_SR = numeric(), est_R = numeric())
  for(i in 1:length(all_aic)){
    #est_ind has No SR assumption first
    if(df.oms$NAA_sig[om_ind[i]] == 0) {  #if there are no NAA random effects
        est_ind <- which(use.df.ems$re_config %in% c("rec") & use.df.ems$M_est == TRUE) 
    } else est_ind <- which(use.df.ems$re_config %in% c("rec+1") & use.df.ems$M_est == TRUE)
    tmp = t(apply(all_aic[[i]][rev(est_ind),],2, function(x) {
      if(any(!is.na(x))) {
        return(x == min(x,na.rm=T))
      } else return(rep(NA, length(x)))
    })) #tmp is rows of T and F for which model had lowest AIC
    df_i <- cbind.data.frame(OM = i, R_sig = df.oms$R_sig[om_ind[i]], NAA_sig = df.oms$NAA_sig[om_ind[i]], 
        Fhist = df.oms$Fhist[om_ind[i]], obs_error = df.oms$obs_error[om_ind[i]],
      sd_log_SSB = all_naa_sd_log_SSB[[om_ind[i]]], est_SR = as.numeric(tmp[,1]), est_R = as.numeric(tmp[,2]))
    df_out <- rbind(df_out, df_i)
  }
  return(df_out)
}
temp = make_df(all_naa_aic, SR_rec_Mest, om_ind)
