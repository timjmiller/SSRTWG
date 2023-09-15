library("here")
library("tidyr")
df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))

#first NAA operating models
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.oms.RDS"))
use.df.ems <- df.ems[1:20,]

make_results <- FALSE
if(make_results) {
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
  saveRDS(all_naa_sd_log_SSB, file = here("Project_0","results", "all_naa_sd_log_SSB.RDS"))

  df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.M.oms.RDS"))
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
  saveRDS(all_M_sd_log_SSB, file = here("Project_0","results", "all_M_sd_log_SSB.RDS"))

  df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.Sel.oms.RDS"))
  all_Sel_sd_log_SSB = lapply(1:NROW(df.oms), function(y){
    res = sapply(1:100, function(x){
      print(paste0("sim_",x,".RDS"))
      sim = readRDS(file.path(here::here(),"Project_0", "results", "Sel_om", paste0("om_", y), paste0("sim_",x,".RDS")))
      res = sapply(sim, function(x) {
        out <- NA
        if(length(x$truth)) out <- sd(log(x$truth$SSB))
        return(out)
      })
      return(res[which(!is.na(res))][1])
    })
    return(res)
  })
  saveRDS(all_Sel_sd_log_SSB, file = here("Project_0","results", "all_Sel_sd_log_SSB.RDS"))

  df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.q.oms.RDS"))
  all_q_sd_log_SSB = lapply(1:NROW(df.oms), function(y){
    res = sapply(1:100, function(x){
      print(paste0("sim_",x,".RDS"))
      sim = readRDS(file.path(here::here(),"Project_0", "results", "q_om", paste0("om_", y), paste0("sim_",x,".RDS")))
      res = sapply(sim, function(x) {
        out <- NA
        if(length(x$truth)) out <- sd(log(x$truth$SSB))
        return(out)
      })
      return(res[which(!is.na(res))][1])
    })
    return(res)
  })
  saveRDS(all_q_sd_log_SSB, file = here("Project_0","results", "all_q_sd_log_SSB.RDS"))
}

coefficient_fn <- function(om = 1, df){
  fit <- unname(summary(glm(cbind(est_SR,est_R) ~ log(sd_log_SSB), family = binomial, data = subset(df, OM == om)))$coefficients[2,])
  return(c(Estimate = fit[1], SE = fit[2], lo = fit[1] + qnorm(0.025)*fit[2], hi = fit[1] + qnorm(0.975)*fit[2]))
}


all_naa_sd_log_SSB <- readRDS(file = here("Project_0","results", "all_naa_sd_log_SSB.RDS"))
all_naa_sd_log_SSB <- t(matrix(unname(unlist(all_naa_sd_log_SSB)), nrow = 100))
colnames(all_naa_sd_log_SSB) <- paste0("sim", 1:100)
df <- cbind(df.oms, all_naa_sd_log_SSB)
df <- df %>% pivot_longer(cols = paste0("sim", 1:100), names_to = "sim", values_to = "sd_log_SSB")
#df$aic_BH_ME <- sapply(all)

#om_ind <- which(!is.na(df.oms$NAA_sig))
om_ind <- 1:24
#SR_rec_Mest = which(use.df.ems$re_config == "rec+1" & use.df.ems$M_est == TRUE)
SR_rec_Mest = which(use.df.ems$re_config %in% c("rec", "rec+1") & use.df.ems$M_est == TRUE)
all_naa_aic <- readRDS(file = file.path(here(),"Project_0","results", "all_naa_aic_results.RDS"))
df.oms$NAA_sig[which(is.na(df.oms$NAA_sig))] <- 0


aic_fn <- function(x, all_aic, df.oms, df.ems, rec_mod, M_est) {
  if(!is.null(df.oms$NAA_sig)){
    if(is.na(df.oms$NAA_sig[x])) re_config <- "rec"
    else re_config <- "rec+1"  
    em_ind <- which(df.ems$SR_model == rec_mod & df.ems$M_est == M_est & df.ems$re_config == re_config)
  }
  if(!is.null(df.oms$M_sig)){
    if(df.oms$M_cor[x] == 0) re_config <- "iid"
    else re_config <- "ar1_y"  
    em_ind <- which(df.ems$SR_model == rec_mod & df.ems$M_est == M_est & df.ems$M_re_cor == re_config)
  }
  if(!is.null(df.oms$Sel_sig)){
    if(df.oms$Sel_cor[x] == 0) re_config <- "iid"
    else re_config <- "ar1_y"  
    em_ind <- which(df.ems$SR_model == rec_mod & df.ems$M_est == M_est & df.ems$sel_re_cor == re_config)
  }
  if(!is.null(df.oms$q_sig)){
    if(df.oms$q_cor[x] == 0) re_config <- "iid"
    else re_config <- "ar1"  
    em_ind <- which(df.ems$SR_model == rec_mod & df.ems$M_est == M_est & df.ems$q_re_cor == re_config)
  }
  return(all_aic[[x]][em_ind,]) #single EM/row
}
coefficient_fn <- function(om = 1, df, col = "BH_best_ME"){
  fit <- unname(summary(glm(get(col) ~ log(sd_log_SSB), family = binomial, data = subset(df, OM == om)))$coefficients[2,])
  return(c(Estimate = fit[1], SE = fit[2], lo = fit[1] + qnorm(0.025)*fit[2], hi = fit[1] + qnorm(0.975)*fit[2]))
}

df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))
# em_ind <- c(5:20,25:28)
# use.df.ems <- df.ems[em_ind,]
types <- c("naa", "M", "Sel", "q")
types.plt <- c("R, R+S","R+M", "R+Sel", "R+q")
est_conds <- cbind(c(FALSE,TRUE))
em_inds <- cbind(1:20, 5:24, c(5:20,25:28), c(5:20, 29:32))
for(i in 1:4) {
  all_aic <- readRDS(file = here("Project_0","results", paste0("all_", types[i], "_aic_results.RDS")))
  df.oms = readRDS(here("Project_0","inputs", paste0("df.", ifelse(i==1, "", paste0(types[i],".")),"oms.RDS")))
  use.df.ems <- df.ems[em_inds[,i],]
  df.oms$OM <- 1:NROW(df.oms)
  sd_log_SSB <- readRDS(file = here("Project_0","results", paste0("all_", types[i], "_sd_log_SSB.RDS")))
  sd_log_SSB <- t(matrix(unname(unlist(sd_log_SSB)), nrow = 100))
#  colnames(sd_log_SSB) <- paste0("sim", 1:100)
  df <- cbind(df.oms, sd_log_SSB)
  df <- df %>% pivot_longer(cols = paste0(1:100), names_to = "sim", values_to = "sd_log_SSB")
  df$aic_BH_ME <- c(sapply(1:length(all_aic), aic_fn, all_aic = all_aic, df.oms = df.oms, df.ems = use.df.ems, 3, TRUE))
  df$aic_R_ME <- c(sapply(1:length(all_aic), aic_fn, all_aic = all_aic, df.oms = df.oms, df.ems = use.df.ems, 2, TRUE))
  df$BH_best_ME <- as.integer(df$aic_BH_ME < df$aic_R_ME)
  df$aic_BH_MF <- c(sapply(1:length(all_aic), aic_fn, all_aic = all_aic, df.oms = df.oms, df.ems = use.df.ems, 3, FALSE))
  df$aic_R_MF <- c(sapply(1:length(all_aic), aic_fn, all_aic = all_aic, df.oms = df.oms, df.ems = use.df.ems, 2, FALSE))
  df$BH_best_MF <- as.integer(df$aic_BH_MF < df$aic_R_MF)
  if(i==1) {
    df$NAA_sig[which(is.na(df$NAA_sig))] <- 0
    df <- df %>% mutate(NAA_sig = recode(NAA_sig,
        "0" = "sigma['2+'] == 0",
        "0.25" = "sigma['2+'] == 0.25",
        "0.5" = "sigma['2+'] == 0.5")) %>%
      mutate(R_sig = recode(R_sig,
        "0.5" = "sigma[R] == 0.5",
        "1.5" = "sigma[R] == 1.5")) %>% as.data.frame
  }
  if(i ==2){
    df <- df %>% mutate(M_sig = recode(M_sig,
        "0.1" = "sigma['M'] == 0.1",
        "0.5" = "sigma['M'] == 0.5")) %>% 
      mutate(M_cor = recode(M_cor,
        "0" = "rho['M'] == 0",
        "0.9" = "rho['M'] == 0.9"))
  }
  if(i == 3){
    df <- df %>% mutate(Sel_sig = recode(Sel_sig,
        "0.1" = "sigma['Sel'] == 0.1",
        "0.5" = "sigma['Sel'] == 0.5")) %>% 
      mutate(Sel_cor = recode(Sel_cor,
        "0" = "rho['Sel'] == 0",
        "0.9" = "rho['Sel'] == 0.9"))
  }
  if(i == 4){
    df <- df %>% mutate(q_sig = recode(q_sig,
        "0.1" = "sigma['q'] == 0.1",
        "0.5" = "sigma['q'] == 0.5")) %>%
      mutate(q_cor = recode(q_cor,
        "0" = "rho['q'] == 0",
        "0.9" = "rho['q'] == 0.9"))  
  }
  df <- df %>% as.data.frame
  facs <- c("Fhist","NAA_sig", "R_sig", "M_sig", "M_cor", "Sel_sig", "Sel_cor", "q_sig", "q_cor", "obs_error")
  df[names(df) %in% facs] <- lapply(df[names(df) %in% facs], as.factor)

  for(j in c("MF","ME")){
    fit <- glm(get(paste0("BH_best_", j)) ~ factor(OM) * log(sd_log_SSB), family = binomial, data = df)
    #x <- t(sapply(1:24, coefficient_fn, df =df, col = "BH_best_MF"))
    df[[paste0("pred_",j)]] <- NA
    df[[paste0("pred_",j)]][(!is.na(df[[paste0("BH_best_", j)]]))] <- predict(fit, type = "response")
    # fit <- glm(BH_best_ME ~ factor(OM) * log(sd_log_SSB), family = binomial, data = df)
    # #x <- t(sapply(1:24, coefficient_fn, df =df, col = "BH_best_MF"))
    # df$pred_ME <- NA
    # df$pred_ME[(!is.na(df$BH_best_ME))] <- predict(fit, type = "response")
    plt <- ggplot(df, aes(x = sd_log_SSB, y = get(paste0("pred_",j)), colour = Fhist:obs_error)) + 
      scale_color_viridis_d() +
      theme_bw() +
      geom_rug(position = "jitter", sides = "b", alpha = 0.5) +
      geom_line(linewidth = 2, alpha = 0.5) + 
      xlab("SD(log(SSB))") + ylab("P(B-H assumption best)")
    if(i == 1) plt <- plt + facet_grid(R_sig ~ NAA_sig, labeller = label_parsed)
    if(i == 2) plt <- plt + facet_grid(M_sig ~ M_cor, labeller = label_parsed)
    if(i == 3) plt <- plt + facet_grid(Sel_sig ~ Sel_cor, labeller = label_parsed)
    if(i == 4) plt <- plt + facet_grid(q_sig ~ q_cor, labeller = label_parsed)
    plt
    ggsave(here("Project_0", "paper", paste0(types[i], "_om_", j, "_pred_BH_best.png")), plt, height = 12, width = 20, units = "in")
  }

}

x <- t(sapply(1:24, coefficient_fn, df =df, col = "BH_best_MF"))
x <- t(sapply(1:24, coefficient_fn, df =df, col = "BH_best_ME"))

# make_df <- function(all_aic, est_ind = SR_rec_Mest, om_ind = NULL){
#   if(!is.null(om_ind)) all_aic <- all_aic[om_ind]
#   df_out <- cbind.data.frame(OM = numeric(),
#     R_sig = numeric(), NAA_sig = numeric(), Fhist = character(), obs_error = character(), sd_log_SSB = numeric(),
#     est_SR = numeric(), est_R = numeric())
#   for(i in 1:length(all_aic)){
#     #est_ind has No SR assumption first
#     if(df.oms$NAA_sig[om_ind[i]] == 0) {
#         est_ind <- which(use.df.ems$re_config %in% c("rec") & use.df.ems$M_est == TRUE)
#     } else est_ind <- which(use.df.ems$re_config %in% c("rec+1") & use.df.ems$M_est == TRUE)
#     tmp = t(apply(all_aic[[i]][rev(est_ind),],2, function(x) {
#       if(any(!is.na(x))) {
#         return(x == min(x,na.rm=T))
#       } else return(rep(NA, length(x)))
#     })) #tmp is rows of T and F for which model had lowest AIC
#     df_i <- cbind.data.frame(OM = i, R_sig = df.oms$R_sig[om_ind[i]], NAA_sig = df.oms$NAA_sig[om_ind[i]], 
#         Fhist = df.oms$Fhist[om_ind[i]], obs_error = df.oms$obs_error[om_ind[i]],
#       sd_log_SSB = all_naa_sd_log_SSB[[om_ind[i]]], est_SR = as.numeric(tmp[,1]), est_R = as.numeric(tmp[,2]))
#     df_out <- rbind(df_out, df_i)
#   }
#   return(df_out)
# }
# temp = make_df(all_naa_aic, SR_rec_Mest, om_ind)
# fit <- glm(cbind(BH_best_MF) ~ factor(OM) * log(sd_log_SSB), family = binomial, data = temp)
# x <- t(sapply(1:24, coefficient_fn, df =temp))
# temp <- cbind(temp, pred= NA)
# temp$pred[(!is.na(temp$est_R) & !is.na(temp$est_SR))] <- predict(fit, type = "response")
# facs <- c("Fhist","NAA_sig", "R_sig", "obs_error")
# temp[facs] <- lapply(temp[facs], as.factor)
# library(ggplot2)
# plt <- ggplot(temp, aes(x = sd_log_SSB, y = pred, colour = Fhist:obs_error)) + 
#   scale_color_viridis_d() +
#   facet_grid(R_sig ~ NAA_sig, labeller = label_both) +
#   theme_bw() +
#   geom_rug(position = "jitter", sides = "b", alpha = 0.5) +
#   geom_line(linewidth = 2, alpha = 0.5) + 
#   xlab("SD(log(SSB))") + ylab("P(B-H assumption best)") 
# plt
#   #scale_linetype_discrete(labels = c("2.5 Fmsy -> Fmsy", "Fmsy"))
# ggsave(here("Project_0", "paper", "pred_SR_best_NAA_oms.png"), plt, height = 12, width = 20, units = "in")

# fit_om1 <- glm(cbind(est_SR,est_R) ~ log(sd_log_SSB), family = binomial, data = subset(temp, OM == 1))
# temp <- subset(temp, !is.na(est_SR) & !is.na(est_R))
# fit0 <- glm(cbind(est_SR,est_R) ~  1, family = binomial, data = temp)
# fit1 <-glm(cbind(est_SR,est_R) ~ factor(R_sig), family = binomial, data = temp)
# fit2 <-glm(cbind(est_SR,est_R) ~ factor(NAA_sig), family = binomial, data = temp)
# fit3 <-glm(cbind(est_SR,est_R) ~ Fhist, family = binomial, data = temp)
# fit4 <-glm(cbind(est_SR,est_R) ~ obs_error, family = binomial, data = temp)
# fit5 <-glm(cbind(est_SR,est_R) ~ sd_log_SSB, family = binomial, data = temp)
# AIC(fit0,fit1,fit2,fit3,fit4,fit5)
# fit6 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + factor(NAA_sig), family = binomial, data = temp)
# fit7 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + Fhist, family = binomial, data = temp)
# fit8 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + obs_error, family = binomial, data = temp)
# fit9 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + sd_log_SSB, family = binomial, data = temp)
# AIC(fit1,fit6,fit7,fit8,fit9)
# fit10 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + sd_log_SSB + factor(NAA_sig), family = binomial, data = temp)
# fit11 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + sd_log_SSB + Fhist, family = binomial, data = temp)
# fit12 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + sd_log_SSB + obs_error, family = binomial, data = temp)
# AIC(fit9,fit10,fit11,fit12)
# fit13 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + sd_log_SSB + factor(NAA_sig) + obs_error, family = binomial, data = temp)
# fit14 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + sd_log_SSB + factor(NAA_sig) + Fhist, family = binomial, data = temp)
# AIC(fit10,fit13,fit14)
# fit15 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + sd_log_SSB + factor(NAA_sig) + obs_error + Fhist, family = binomial, data = temp)
# AIC(fit13,fit15)
# fit16 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) * sd_log_SSB + factor(NAA_sig) + obs_error, family = binomial, data = temp)
# fit17 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + factor(NAA_sig) * sd_log_SSB + obs_error, family = binomial, data = temp)
# fit18 <-glm(cbind(est_SR,est_R) ~ factor(R_sig) + factor(NAA_sig) + obs_error * sd_log_SSB, family = binomial, data = temp)
# AIC(fit13,fit16,fit17,fit18)
# fit19 <-glm(cbind(est_SR,est_R) ~ (factor(R_sig) + factor(NAA_sig)) * sd_log_SSB + obs_error, family = binomial, data = temp)
# fit20 <-glm(cbind(est_SR,est_R) ~ (factor(R_sig) + obs_error) * sd_log_SSB + factor(NAA_sig), family = binomial, data = temp)
# AIC(fit16,fit19,fit20)

# summary(fit16)
# library(statmod)
# qres <- qresiduals(fit16)
# qqnorm(qres)
# qqline(qres)
# summary(qres)

# df.oms[om_ind,]
# qres <- qresiduals(fit0)
# qqnorm(qres)
# qqline(qres)

# x <- cbind(temp, pred= predict(fit12, type = "response"))
# library(ggplot2)
# plt <- ggplot(x, aes(x = sd_log_SSB, y = pred, color = obs_error)) + 
#     geom_line() + facet_grid(R_sig ~ NAA_sig, labeller = label_both)
# plt <- plt + xlab("SD(log(SSB))") + ylab("P(B-H assumption best)") + 
#   labs(color = "Observation\n Error") +
#   scale_color_discrete(labels = c("High", "Low")) +
#   theme_bw()
#   #scale_linetype_discrete(labels = c("2.5 Fmsy -> Fmsy", "Fmsy"))
# ggsave(here("Project_0", "paper", "pred_SR_best_NAA_oms.png"), plt)



#M operating models
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.M.oms.RDS"))
use.df.ems = df.ems[5:24,]
use.df.ems$re_config[c(5:8,17:20)] = paste0("M_",c("iid","ar1")[match(use.df.ems$M_re_cor[c(5:8,17:20)], c("iid","ar1_y"))])
use.df.ems = use.df.ems[c(5:8,17:20,1:4,9:16),]

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
