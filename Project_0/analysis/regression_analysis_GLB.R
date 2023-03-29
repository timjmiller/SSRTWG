library("here")
df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))
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
#SR_rec_Mest = which(use.df.ems$re_config %in% c("rec", "rec+1") & use.df.ems$M_est == TRUE)
all_naa_relssb <- readRDS(file = file.path(here(),"Project_0","results", "all_naa_relssb_results.RDS"))
df.oms$NAA_sig[which(is.na(df.oms$NAA_sig))] <- 0
#make a data frame with columns defining operating model characteristics and SD of true (log) SSB 
# and 1/0 for which SR assumption had lowest AIC
make_df_ssb <- function(all_naa_relssb, est_ind = 1:length(df.ems), om_ind = NULL){
  if(!is.null(om_ind)) all_naa_relssb <- all_naa_relssb[om_ind] 
  df_out <- cbind.data.frame(OM = numeric(),
                             R_sig = numeric(), NAA_sig = numeric(), Fhist = character(), obs_error = character(), sd_log_SSB = numeric(),
                             relSSB = numeric())
  for(i in 1:length(all_naa_relssb)){
    tmp <- lapply(all_naa_relssb[[i]], function(k){
      return(colMeans(k))
    })
    df_i <- cbind.data.frame(OM = i, R_sig = df.oms$R_sig[om_ind[i]], NAA_sig = df.oms$NAA_sig[om_ind[i]], 
                             Fhist = df.oms$Fhist[om_ind[i]], obs_error = df.oms$obs_error[om_ind[i]],
                             sd_log_SSB = rep(all_naa_sd_log_SSB[[om_ind[i]]],20), relSSB = unlist(tmp),
                             SR=rep(df.ems[1:20,1],100),M_est=rep(df.ems[1:20,2],100),re=rep(df.ems[1:20,3],100))
    df_out <- rbind(df_out, df_i)
  }
  return(df_out)
}

temp = make_df_ssb(all_naa_relssb,om_ind=om_ind)
temp <- temp[temp$relSSB < 1.5,]
fit <- lm(relSSB ~ R_sig + NAA_sig + as.factor(Fhist) + as.factor(obs_error) + sd_log_SSB + as.factor(SR) + as.factor(M_est) + as.factor(re), 
          data=temp)
summary(fit)




