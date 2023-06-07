relres_fn <- function(df.oms, om_type = "naa", res = "ssb"){
  all_out <- lapply(1:NROW(df.oms), function(y){
    res <- lapply(1:100, function(x){
      print(paste0("om_", y, ", sim_",x))
      sim = readRDS(file.path(here::here(),"Project_0", "results", paste0(om_type,"_om"), paste0("om_", y), paste0("sim_",x,".RDS")))
      relres <- lapply(sim,function(z) {
        out <- matrix(NA, 40,3)
        if(length(z)) if(length(z$fit)) {
          if(res == "ssb") out[,1] <- z$fit$rep$SSB/z$truth$SSB
          if(res == "R") out[,1] <- z$fit$rep$NAA[,1]/z$truth$NAA[,1]
          if(res == "F") out[,1] <- z$fit$rep$F[,1]/z$truth$F[,1]
        }
        if(!is.null(z$fit$sdrep)){
          #ind <- which(!is.na(z$fit$sdrep$SE_rep$log_SSB))
          if(res == "ssb") {
            out[,2] <- z$fit$sdrep$SE_rep$log_SSB #se
            ci <- matrix(z$fit$sdrep$Estimate_rep$log_SSB - log(z$truth$SSB),40,2) + qnorm(0.975) * cbind(-out[,2],out[,2])
          }
          if(res == "R") {
            out[,2] <- z$fit$sdrep$SE_rep$log_NAA_rep[,1] #se
            ci <- matrix(z$fit$sdrep$Estimate_rep$log_NAA_rep[,1] - log(z$truth$NAA[,1]), 40,2) + qnorm(0.975) * cbind(-out[,2],out[,2])
          }
          if(res == "F") {
            # print(z$fit$sdrep$SE_rep$log_F[,1])
            # print(length(z$fit$sdrep$SE_rep$log_F[,1]))
            # print(dim(out))
            out[,2] <- z$fit$sdrep$SE_rep$log_F[,1] #se
            ci <- matrix(z$fit$sdrep$Estimate_rep$log_F[,1] - log(z$truth$F[,1]), 40,2) + qnorm(0.975) * cbind(-out[,2],out[,2])
          }
          out[,3] <- 0 >= ci[,1] & 0 <= ci[,2]
        }
        return(out)
      })
      return(relres)
    })
    return(res)
  })
  return(all_out)
}


relssb_fn <- function(df.oms, om_type = "naa"){
  relssb = lapply(1:NROW(df.oms), function(y){
    res = lapply(1:100, function(x){
      print(paste0("om_", y, ", sim_",x))
      sim = readRDS(file.path(here::here(),"Project_0", "results", paste0(om_type,"_om"), paste0("om_", y), paste0("sim_",x,".RDS")))
      relSSB = lapply(sim,function(z) {
        out <- matrix(NA, 40,3)
        if(length(z)) if(length(z$fit)) out[,1] = z$fit$rep$SSB/z$truth$SSB
        if(!is.null(z$fit$sdrep)){
          #ind <- which(!is.na(z$fit$sdrep$SE_rep$log_SSB))
          out[,2] <- z$fit$sdrep$SE_rep$log_SSB #se
            ci = matrix(z$fit$sdrep$Estimate_rep$log_SSB - log(z$truth$SSB),40,2) + qnorm(0.975) * cbind(-out[,2],out[,2])
          out[,3] <- 0 >= ci[,1] & 0 <= ci[,2]
        }
        return(out)
      })
      return(relSSB)
    })
    return(res)
  })
  return(relssb)
}

relR_fn <- function(df.oms, om_type = "naa", res =){
  relR = lapply(1:NROW(df.oms), function(y){
    res = lapply(1:100, function(x){
      print(paste0("om_", y, ", sim_",x))
      sim = readRDS(file.path(here::here(),"Project_0", "results", paste0(om_type,"_om"), paste0("om_", y), paste0("sim_",x,".RDS")))
      relR = sapply(sim,function(z) {
        out <- rep(NA, 40)
        if(length(z)) if(length(z$fit)) out = z$fit$rep$NAA[,1]/z$truth$NAA[,1]
        if(!is.null(z$fit$sdrep)){
          #ind <- which(!is.na(z$fit$sdrep$SE_rep$log_SSB))
          out[,2] <- c(z$fit$sdrep$SE_rep$log_N1_pars[1], z$fit$sdrep$SE_rep$log_NAA[,1]) #se
            ci = matrix(z$fit$sdrep$Estimate_rep$log_SSB - log(z$truth$SSB),40,2) + qnorm(0.975) * cbind(-out[,2],out[,2])
          out[,3] <- 0 >= ci[,1] & 0 <= ci[,2]
        }
        return(out)
      })
      return(relR)
    })
    return(res)
  })
  return(relR)
}

relF_fn <- function(df.oms, om_type = "naa"){
  relF = lapply(1:NROW(df.oms), function(y){
    res = lapply(1:100, function(x){
      print(paste0("om_", y, ", sim_",x))
      sim = readRDS(file.path(here::here(),"Project_0", "results", paste0(om_type,"_om"), paste0("om_", y), paste0("sim_",x,".RDS")))
      relF = sapply(sim,function(y) {
        out <- rep(NA, 40)
        if(length(y)) if(length(y$fit)) out = y$fit$rep$F[,1]/y$truth$F[,1]
        return(out)
      })
      return(relF)
    })
    return(res)
  })
  return(relF)
}
make_results <- FALSE
if(make_results) {
  df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.oms.RDS"))
  all_naa_relssb <- relssb_fn(df.oms,"naa")
  saveRDS(all_naa_relssb, file = here("Project_0","results", "all_naa_relssb_results.RDS"))
  df.M.oms = readRDS(file.path(here(),"Project_0","inputs", "df.M.oms.RDS"))
  all_M_relssb <- relssb_fn(df.M.oms,"M")
  saveRDS(all_M_relssb, file = here("Project_0","results", "all_M_relssb_results.RDS"))
  df.Sel.oms = readRDS(file.path(here(),"Project_0","inputs", "df.Sel.oms.RDS"))
  all_Sel_relssb <- relssb_fn(df.Sel.oms,"Sel")
  saveRDS(all_Sel_relssb, file = here("Project_0","results", "all_Sel_relssb_results.RDS"))
  df.q.oms = readRDS(file.path(here(),"Project_0","inputs", "df.q.oms.RDS"))
  all_q_relssb <- relssb_fn(df.q.oms,"q")
  saveRDS(all_q_relssb, file = here("Project_0","results", "all_q_relssb_results.RDS"))

  df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.oms.RDS"))
  all_naa_relR <- relres_fn(df.oms,om_type = "naa", res = "R")
  saveRDS(all_naa_relR, file = here("Project_0","results", "all_naa_relR_results.RDS"))
  df.M.oms = readRDS(file.path(here(),"Project_0","inputs", "df.M.oms.RDS"))
  all_M_relR <- relres_fn(df.M.oms,"M", res = "R")
  saveRDS(all_M_relR, file = here("Project_0","results", "all_M_relR_results.RDS"))
  df.Sel.oms = readRDS(file.path(here(),"Project_0","inputs", "df.Sel.oms.RDS"))
  all_Sel_relR <- relres_fn(df.Sel.oms,"Sel", res = "R")
  saveRDS(all_Sel_relR, file = here("Project_0","results", "all_Sel_relR_results.RDS"))
  df.q.oms = readRDS(file.path(here(),"Project_0","inputs", "df.q.oms.RDS"))
  all_q_relR <- relres_fn(df.q.oms,"q", res = "R")
  saveRDS(all_q_relR, file = here("Project_0","results", "all_q_relR_results.RDS"))

  df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.oms.RDS"))
  all_naa_relF <- relres_fn(df.oms,"naa", res = "F")
  saveRDS(all_naa_relF, file = here("Project_0","results", "all_naa_relF_results.RDS"))
  df.M.oms = readRDS(file.path(here(),"Project_0","inputs", "df.M.oms.RDS"))
  all_M_relF <- relres_fn(df.M.oms,"M", res = "F")
  saveRDS(all_M_relF, file = here("Project_0","results", "all_M_relF_results.RDS"))
  df.Sel.oms = readRDS(file.path(here(),"Project_0","inputs", "df.Sel.oms.RDS"))
  all_Sel_relF <- relres_fn(df.Sel.oms,"Sel", res = "F")
  saveRDS(all_Sel_relF, file = here("Project_0","results", "all_Sel_relF_results.RDS"))
  df.q.oms = readRDS(file.path(here(),"Project_0","inputs", "df.q.oms.RDS"))
  all_q_relF <- relres_fn(df.q.oms,"q", res = "F")
  saveRDS(all_q_relF, file = here("Project_0","results", "all_q_relF_results.RDS"))
}


custom_boxplot_stat <- function(x){#, n.yrs = 1, n.sim = 100) {
  x <- x[which(!is.na(x))]
  n <- length(x)
  bnds95 <- qbinom(c(0.025,0.975), n, 0.5)/n # 95% CI bounds for median (sims x years)
  bnds80 <- qbinom(c(0.1,0.9), n, 0.5)/n # 80% CI bounds for median (sims x years)
  r <- quantile(x, probs = c(bnds95[1], bnds80[1], 0.5, bnds80[2], bnds95[2]))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
get_df_ems <- function(type = "naa"){
  df.ems = readRDS(here::here("Project_0","inputs", "df.ems.RDS"))
  df.ems$pe <- c("R","R+S","R+M","R+Sel","R+q")[match(df.ems$re_config, c("rec","rec+1","M_re", "sel_re","q_re"))]
  em_ind <- 1:20
  if(type != "naa"){
    if(type == "M") {
      em_ind <- 5:24
      ind <- which(df.ems$re_config =="M_re")
      df.ems$pe[ind] <- paste0(df.ems$pe[ind], "(", df.ems$M_re_cor[ind], ")") 
    }
    if(type == "Sel") {
      em_ind <- c(5:20,25:28)
      ind <- which(df.ems$re_config =="sel_re")
      df.ems$pe[ind] <- paste0(df.ems$pe[ind], "(", df.ems$sel_re_cor[ind], ")") 
    }
    if(type == "q") {
      em_ind <- c(5:20,29:32)
      ind <- which(df.ems$re_config =="q_re")
      df.ems$pe[ind] <- paste0(df.ems$pe[ind], "(", df.ems$q_re_cor[ind], ")") 
    }
  }
  df.ems <- df.ems[em_ind,]
  df.ems$meanR <- "Mean R"
  df.ems$meanR[df.ems$SR_model ==3] <- "B-H"
  return(df.ems)
}
naa_om_df.ems <- get_df_ems()
M_om_df.ems <- get_df_ems("M")
Sel_om_df.ems <- get_df_ems("Sel")
q_om_df.ems <- get_df_ems("q")

make_plot_df <- function(type = "naa", res = all_naa_om_relssb, is_SE = FALSE) {
  library(reshape2)
  res <- melt(res)
  names(res) <- c("year", "column", "value", "em", "sim","om")
  res <- res %>% mutate(column = recode(column,
      "1" = "relerror",
      "2" = "cv",
      "3" = "in_ci",
    )) %>% as.data.frame   
  df <- res %>% pivot_wider(names_from = column, values_from = value) %>% as.data.frame
  print(dim(df))
  if(is_SE) df <- filter(df, !is.na(cv))
  print(dim(df))
  df$relerror = df$relerror - 1

  df <- df %>% group_by(om, em, year) %>%
    reframe(stats = custom_boxplot_stat(relerror)) %>% as.data.frame
  df$type <- c("ymin", "lower", "middle", "upper", "ymax")
  library(tidyr)
  df.ems <- get_df_ems(type)
  if(type == "naa") {
    df.oms <- readRDS(here::here("Project_0", "inputs", "df.oms.RDS"))
  } else {
    df.oms <- readRDS(here::here("Project_0", "inputs", paste0("df.", type, ".oms.RDS")))
  }
  print(head(df.oms))
  df <- df %>% pivot_wider(names_from = type, values_from = stats) %>% as.data.frame
  df$em_pe = df.ems$pe#[df$em]
  df <- cbind(df, df.ems[df$em,])
  df <- cbind(df, df.oms[df$om,])
  if(type == "naa") {
    df$NAA_sig[which(is.na(df$NAA_sig))] <- 0
    df <- df %>% mutate(NAA_sig = recode(NAA_sig,
        "0" = "sigma['2+'] == 0",
        "0.25" = "sigma['2+'] == 0.25",
        "0.5" = "sigma['2+'] == 0.5")) %>%
      mutate(R_sig = recode(R_sig,
        "0.5" = "sigma[R] == 0.5",
        "1.5" = "sigma[R] == 1.5"))
  }
  if(type == "M") {
    df <- df %>% mutate(M_sig = recode(M_sig,
        "0.1" = "sigma['M'] == 0.1",
        "0.5" = "sigma['M'] == 0.5")) %>% 
      mutate(M_cor = recode(M_cor,
        "0" = "rho['M'] == 0",
        "0.9" = "rho['M'] == 0.9"))
  }
  if(type == "Sel") {
    df <- df %>% mutate(Sel_sig = recode(Sel_sig,
        "0.1" = "sigma['Sel'] == 0.1",
        "0.5" = "sigma['Sel'] == 0.5")) %>% 
      mutate(Sel_cor = recode(Sel_cor,
        "0" = "rho['Sel'] == 0",
        "0.9" = "rho['Sel'] == 0.9"))
  }
  if(type == "q") {
    df <- df %>% mutate(q_sig = recode(q_sig,
        "0.1" = "sigma['q'] == 0.1",
        "0.5" = "sigma['q'] == 0.5")) %>%
      mutate(q_cor = recode(q_cor,
        "0" = "rho['q'] == 0",
        "0.9" = "rho['q'] == 0.9"))  
  }
  df <- df %>% mutate(obs_error = recode(obs_error,
      "L" = "Low obs error (indices, age comp)",
      "H" = "High obs error (indices, age comp)"))
  df <- df %>% mutate(Fhist = recode(Fhist,
      "H-MSY" = "F history: High->FMSY",
      "MSY" = "F history: FMSY"))
  df <- df %>% as.data.frame
  facs <- c("Fhist","NAA_sig", "R_sig", "M_sig", "M_cor", "Sel_sig", "Sel_cor", "q_sig", "q_cor", "obs_error")
  df[names(df) %in% facs] <- lapply(df[names(df) %in% facs], as.factor)
  return(df)
}

library(here)
library(dplyr)
library(tidyr)
all_naa_om_relssb <- readRDS(file = here("Project_0","results", "all_naa_relssb_results.RDS"))
library(reshape2)
temp <- melt(all_naa_om_relssb)
names(temp) <- c("year", "column", "value", "em", "sim","om")
temp <- temp %>% mutate(column = recode(column,
    "1" = "relbias",
    "2" = "cv",
    "3" = "in_ci",
  )) %>% as.data.frame
temp <- temp %>% pivot_wider(names_from = column, values_from = value) %>% as.data.frame

estM <- c(FALSE,FALSE,TRUE,TRUE)
estSR<- c("Mean R", "B-H", "Mean R", "B-H")
##############################################################################################################
#naa oms

#SSB
plot.df <- make_plot_df(is_SE = TRUE)
for(i in 1:4){
  plt <- ggplot(filter(plot.df, M_est== estM[i] & meanR == estSR[i] )) + scale_fill_viridis_d() + 
      geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
      facet_grid(R_sig + NAA_sig ~ Fhist + obs_error, labeller = labeller(R_sig = label_parsed, NAA_sig = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
      geom_line(aes(x = year, y = middle, colour = pe), position = position_dodge(0.1), linewidth = 1) + #geom_point(position = position_dodge(0.1), size = 1) + 
      geom_ribbon(aes(x = year, ymin=ymin,ymax=ymax, fill = pe), alpha=0.3, position = position_dodge(0.1)) +
      theme_bw() + coord_cartesian(ylim = c(-0.4, 0.4)) + ylab("Relative Bias of SSB") + xlab("Year") +
      ggtitle(paste0("EM: ", ifelse(estM[i], "M estimated", "M = 0.2"), " and ", ifelse(estSR[i]=="B-H", "BH assumed", "no SRR"))) + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM Process Error", fill = "EM Process Error")
  plt
  ggsave(here("Project_0", "paper", paste0("naa_om_", ifelse(estSR[i]=="B-H","BH","R"), "_M", ifelse(estM[i], "E","F"), "_relbias_ssb.png")), plt, width = 20, height = 12, units = "in")
}

#F
all_naa_om_relF <- readRDS(file = here("Project_0","results", "all_naa_relF_results.RDS"))
plot.df <- make_plot_df("naa", all_naa_om_relF, is_SE = TRUE)
for(i in 1:4){
  plt <- ggplot(filter(plot.df, M_est== estM[i] & meanR == estSR[i] )) + scale_fill_viridis_d() + 
      geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
      facet_grid(R_sig + NAA_sig ~ Fhist + obs_error, labeller = labeller(R_sig = label_parsed, NAA_sig = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
      geom_line(aes(x = year, y = middle, colour = pe), position = position_dodge(0.1), linewidth = 1) + #geom_point(position = position_dodge(0.1), size = 1) + 
      geom_ribbon(aes(x = year, ymin=ymin,ymax=ymax, fill = pe), alpha=0.3, position = position_dodge(0.1)) +
      theme_bw() + coord_cartesian(ylim = c(-0.4, 0.4)) + ylab("Relative Bias of F") + xlab("Year") +
      ggtitle(paste0("EM: ", ifelse(estM[i], "M estimated", "M = 0.2"), " and ", ifelse(estSR[i]=="B-H", "BH assumed", "no SRR"))) + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM Process Error", fill = "EM Process Error")
  plt
  ggsave(here("Project_0", "paper", paste0("naa_om_", ifelse(estSR[i]=="B-H","BH","R"), "_M", ifelse(estM[i], "E","F"), "_relbias_F.png")), plt, width = 20, height = 12, units = "in")
}

#R
all_naa_om_relR <- readRDS(file = here("Project_0","results", "all_naa_relR_results.RDS"))
plot.df <- make_plot_df("naa", all_naa_om_relR, is_SE = TRUE)
for(i in 1:4){
  plt <- ggplot(filter(plot.df, M_est== estM[i] & meanR == estSR[i] )) + scale_fill_viridis_d() + 
      geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
      facet_grid(R_sig + NAA_sig ~ Fhist + obs_error, labeller = labeller(R_sig = label_parsed, NAA_sig = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
      geom_line(aes(x = year, y = middle, colour = pe), position = position_dodge(0.1), linewidth = 1) + #geom_point(position = position_dodge(0.1), size = 1) + 
      geom_ribbon(aes(x = year, ymin=ymin,ymax=ymax, fill = pe), alpha=0.3, position = position_dodge(0.1)) +
      theme_bw() + coord_cartesian(ylim = c(-0.4, 0.4)) + ylab("Relative Bias of Recruitment") + xlab("Year") +
      ggtitle(paste0("EM: ", ifelse(estM[i], "M estimated", "M = 0.2"), " and ", ifelse(estSR[i]=="B-H", "BH assumed", "no SRR"))) + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM Process Error", fill = "EM Process Error")
  plt
  ggsave(here("Project_0", "paper", paste0("naa_om_", ifelse(estSR[i]=="B-H","BH","R"), "_M", ifelse(estM[i], "E","F"), "_relbias_R.png")), plt, width = 20, height = 12, units = "in")
}

##############################################################################################################
#M oms
all_M_relssb <- readRDS(file = here("Project_0","results", "all_M_relssb_results.RDS"))
plot.df <- make_plot_df("M", all_M_relssb, is_SE = TRUE)

for(i in 1:4){
  plt <- ggplot(filter(plot.df, M_est== estM[i] & meanR == estSR[i] )) + scale_fill_viridis_d() + 
      geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
      facet_grid(M_sig + M_cor ~ Fhist + obs_error, labeller = labeller(M_sig = label_parsed, M_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
      geom_line(aes(x = year, y = middle, colour = pe), position = position_dodge(0.1), linewidth = 1) + #geom_point(position = position_dodge(0.1), size = 1) + 
      geom_ribbon(aes(x = year, ymin=ymin,ymax=ymax, fill = pe), alpha=0.3, position = position_dodge(0.1)) +
      theme_bw() + coord_cartesian(ylim = c(-0.4, 0.4)) + ylab("Relative Bias of SSB") + xlab("Year") +
      ggtitle(paste0("EM: ", ifelse(estM[i], "M estimated", "M = 0.2"), " and ", ifelse(estSR[i]=="B-H", "BH assumed", "no SRR"))) + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM Process Error", fill = "EM Process Error")
  plt
  ggsave(here("Project_0", "paper", paste0("M_om_", ifelse(estSR[i]=="B-H","BH","R"), "_M", ifelse(estM[i], "E","F"), "_relbias_ssb.png")), plt, width = 20, height = 12, units = "in")
}

#F
all_M_om_relF <- readRDS(file = here("Project_0","results", "all_M_relF_results.RDS"))
plot.df <- make_plot_df("M", all_M_om_relF, is_SE = TRUE)
for(i in 1:4){
  plt <- ggplot(filter(plot.df, M_est== estM[i] & meanR == estSR[i] )) + scale_fill_viridis_d() + 
      geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
      facet_grid(M_sig + M_cor ~ Fhist + obs_error, labeller = labeller(M_sig = label_parsed, M_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
      geom_line(aes(x = year, y = middle, colour = pe), position = position_dodge(0.1), linewidth = 1) + #geom_point(position = position_dodge(0.1), size = 1) + 
      geom_ribbon(aes(x = year, ymin=ymin,ymax=ymax, fill = pe), alpha=0.3, position = position_dodge(0.1)) +
      theme_bw() + coord_cartesian(ylim = c(-0.4, 0.4)) + ylab("Relative Bias of F") + xlab("Year") +
      ggtitle(paste0("EM: ", ifelse(estM[i], "M estimated", "M = 0.2"), " and ", ifelse(estSR[i]=="B-H", "BH assumed", "no SRR"))) + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM Process Error", fill = "EM Process Error")
  plt
  ggsave(here("Project_0", "paper", paste0("M_om_", ifelse(estSR[i]=="B-H","BH","R"), "_M", ifelse(estM[i], "E","F"), "_relbias_F.png")), plt, width = 20, height = 12, units = "in")
}

#R
all_M_om_relR <- readRDS(file = here("Project_0","results", "all_M_relR_results.RDS"))
plot.df <- make_plot_df("M", all_M_om_relR, is_SE = TRUE)
for(i in 1:4){
  plt <- ggplot(filter(plot.df, M_est== estM[i] & meanR == estSR[i] )) + scale_fill_viridis_d() + 
      geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
      facet_grid(M_sig + M_cor ~ Fhist + obs_error, labeller = labeller(M_sig = label_parsed, M_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
      geom_line(aes(x = year, y = middle, colour = pe), position = position_dodge(0.1), linewidth = 1) + #geom_point(position = position_dodge(0.1), size = 1) + 
      geom_ribbon(aes(x = year, ymin=ymin,ymax=ymax, fill = pe), alpha=0.3, position = position_dodge(0.1)) +
      theme_bw() + coord_cartesian(ylim = c(-0.4, 0.4)) + ylab("Relative Bias of Recruitment") + xlab("Year") +
      ggtitle(paste0("EM: ", ifelse(estM[i], "M estimated", "M = 0.2"), " and ", ifelse(estSR[i]=="B-H", "BH assumed", "no SRR"))) + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM Process Error", fill = "EM Process Error")
  plt
  ggsave(here("Project_0", "paper", paste0("M_om_", ifelse(estSR[i]=="B-H","BH","R"), "_M", ifelse(estM[i], "E","F"), "_relbias_R.png")), plt, width = 20, height = 12, units = "in")
}

##############################################################################################################
#Sel oms

all_Sel_relssb <- readRDS(file = here("Project_0","results", "all_Sel_relssb_results.RDS"))
plot.df <- make_plot_df("Sel", all_Sel_relssb, is_SE = TRUE)
for(i in 1:4){
  plt <- ggplot(filter(plot.df, M_est== estM[i] & meanR == estSR[i] )) + scale_fill_viridis_d() + 
      geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
      facet_grid(Sel_sig + Sel_cor ~ Fhist + obs_error, labeller = labeller(Sel_sig = label_parsed, Sel_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
      geom_line(aes(x = year, y = middle, colour = pe), position = position_dodge(0.1), linewidth = 1) + #geom_point(position = position_dodge(0.1), size = 1) + 
      geom_ribbon(aes(x = year, ymin=ymin,ymax=ymax, fill = pe), alpha=0.3, position = position_dodge(0.1)) +
      theme_bw() + coord_cartesian(ylim = c(-0.4, 0.4)) + ylab("Relative Bias of SSB") + xlab("Year") +
      ggtitle(paste0("EM: ", ifelse(estM[i], "M estimated", "M = 0.2"), " and ", ifelse(estSR[i]=="B-H", "BH assumed", "no SRR"))) + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM Process Error", fill = "EM Process Error")
  plt
  ggsave(here("Project_0", "paper", paste0("Sel_om_", ifelse(estSR[i]=="B-H","BH","R"), "_M", ifelse(estM[i], "E","F"), "_relbias_ssb.png")), plt, width = 20, height = 12, units = "in")
}

#F
all_Sel_om_relF <- readRDS(file = here("Project_0","results", "all_Sel_relF_results.RDS"))
plot.df <- make_plot_df("Sel", all_Sel_om_relF, is_SE = TRUE)
for(i in 1:4){
  plt <- ggplot(filter(plot.df, M_est== estM[i] & meanR == estSR[i] )) + scale_fill_viridis_d() + 
      geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
      facet_grid(Sel_sig + Sel_cor ~ Fhist + obs_error, labeller = labeller(Sel_sig = label_parsed, Sel_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
      geom_line(aes(x = year, y = middle, colour = pe), position = position_dodge(0.1), linewidth = 1) + #geom_point(position = position_dodge(0.1), size = 1) + 
      geom_ribbon(aes(x = year, ymin=ymin,ymax=ymax, fill = pe), alpha=0.3, position = position_dodge(0.1)) +
      theme_bw() + coord_cartesian(ylim = c(-0.4, 0.4)) + ylab("Relative Bias of F") + xlab("Year") +
      ggtitle(paste0("EM: ", ifelse(estM[i], "M estimated", "M = 0.2"), " and ", ifelse(estSR[i]=="B-H", "BH assumed", "no SRR"))) + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM Process Error", fill = "EM Process Error")
  plt
  ggsave(here("Project_0", "paper", paste0("Sel_om_", ifelse(estSR[i]=="B-H","BH","R"), "_M", ifelse(estM[i], "E","F"), "_relbias_F.png")), plt, width = 20, height = 12, units = "in")
}

#R
all_Sel_om_relR <- readRDS(file = here("Project_0","results", "all_Sel_relR_results.RDS"))
plot.df <- make_plot_df("Sel", all_Sel_om_relR, is_SE = TRUE)
for(i in 1:4){
  plt <- ggplot(filter(plot.df, M_est== estM[i] & meanR == estSR[i] )) + scale_fill_viridis_d() + 
      geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
      facet_grid(Sel_sig + Sel_cor ~ Fhist + obs_error, labeller = labeller(Sel_sig = label_parsed, Sel_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
      geom_line(aes(x = year, y = middle, colour = pe), position = position_dodge(0.1), linewidth = 1) + #geom_point(position = position_dodge(0.1), size = 1) + 
      geom_ribbon(aes(x = year, ymin=ymin,ymax=ymax, fill = pe), alpha=0.3, position = position_dodge(0.1)) +
      theme_bw() + coord_cartesian(ylim = c(-0.4, 0.4)) + ylab("Relative Bias of Recruitment") + xlab("Year") +
      ggtitle(paste0("EM: ", ifelse(estM[i], "M estimated", "M = 0.2"), " and ", ifelse(estSR[i]=="B-H", "BH assumed", "no SRR"))) + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM Process Error", fill = "EM Process Error")
  plt
  ggsave(here("Project_0", "paper", paste0("Sel_om_", ifelse(estSR[i]=="B-H","BH","R"), "_M", ifelse(estM[i], "E","F"), "_relbias_R.png")), plt, width = 20, height = 12, units = "in")
}

##############################################################################################################
#q oms
all_q_relssb <- readRDS(file = here("Project_0","results", "all_q_relssb_results.RDS"))
plot.df <- make_plot_df("q", all_q_relssb, is_SE = TRUE)
for(i in 1:4){
  plt <- ggplot(filter(plot.df, M_est== estM[i] & meanR == estSR[i] )) + scale_fill_viridis_d() + 
      geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
      facet_grid(q_sig + q_cor ~ Fhist + obs_error, labeller = labeller(q_sig = label_parsed, q_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
      geom_line(aes(x = year, y = middle, colour = pe), position = position_dodge(0.1), linewidth = 1) + #geom_point(position = position_dodge(0.1), size = 1) + 
      geom_ribbon(aes(x = year, ymin=ymin,ymax=ymax, fill = pe), alpha=0.3, position = position_dodge(0.1)) +
      theme_bw() + coord_cartesian(ylim = c(-0.4, 0.4)) + ylab("Relative Bias of SSB") + xlab("Year") +
      ggtitle(paste0("EM: ", ifelse(estM[i], "M estimated", "M = 0.2"), " and ", ifelse(estSR[i]=="B-H", "BH assumed", "no SRR"))) + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM Process Error", fill = "EM Process Error")
  plt
  ggsave(here("Project_0", "paper", paste0("q_om_", ifelse(estSR[i]=="B-H","BH","R"), "_M", ifelse(estM[i], "E","F"), "_relbias_ssb.png")), plt, width = 20, height = 12, units = "in")
}

#F
all_q_om_relF <- readRDS(file = here("Project_0","results", "all_q_relF_results.RDS"))
plot.df <- make_plot_df("q", all_q_om_relF, is_SE = TRUE)
for(i in 1:4){
  plt <- ggplot(filter(plot.df, M_est== estM[i] & meanR == estSR[i] )) + scale_fill_viridis_d() + 
      geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
      facet_grid(q_sig + q_cor ~ Fhist + obs_error, labeller = labeller(q_sig = label_parsed, q_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
      geom_line(aes(x = year, y = middle, colour = pe), position = position_dodge(0.1), linewidth = 1) + #geom_point(position = position_dodge(0.1), size = 1) + 
      geom_ribbon(aes(x = year, ymin=ymin,ymax=ymax, fill = pe), alpha=0.3, position = position_dodge(0.1)) +
      theme_bw() + coord_cartesian(ylim = c(-0.4, 0.4)) + ylab("Relative Bias of F") + xlab("Year") +
      ggtitle(paste0("EM: ", ifelse(estM[i], "M estimated", "M = 0.2"), " and ", ifelse(estSR[i]=="B-H", "BH assumed", "no SRR"))) + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM Process Error", fill = "EM Process Error")
  plt
  ggsave(here("Project_0", "paper", paste0("q_om_", ifelse(estSR[i]=="B-H","BH","R"), "_M", ifelse(estM[i], "E","F"), "_relbias_F.png")), plt, width = 20, height = 12, units = "in")
}

#R
all_q_om_relR <- readRDS(file = here("Project_0","results", "all_q_relR_results.RDS"))
plot.df <- make_plot_df("q", all_q_om_relR, is_SE = TRUE)
for(i in 1:4){
  plt <- ggplot(filter(plot.df, M_est== estM[i] & meanR == estSR[i] )) + scale_fill_viridis_d() + 
      geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
      facet_grid(q_sig + q_cor ~ Fhist + obs_error, labeller = labeller(q_sig = label_parsed, q_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
      geom_line(aes(x = year, y = middle, colour = pe), position = position_dodge(0.1), linewidth = 1) + #geom_point(position = position_dodge(0.1), size = 1) + 
      geom_ribbon(aes(x = year, ymin=ymin,ymax=ymax, fill = pe), alpha=0.3, position = position_dodge(0.1)) +
      theme_bw() + coord_cartesian(ylim = c(-0.4, 0.4)) + ylab("Relative Bias of Recruitment") + xlab("Year") +
      ggtitle(paste0("EM: ", ifelse(estM[i], "M estimated", "M = 0.2"), " and ", ifelse(estSR[i]=="B-H", "BH assumed", "no SRR"))) + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM Process Error", fill = "EM Process Error")
  plt
  ggsave(here("Project_0", "paper", paste0("q_om_", ifelse(estSR[i]=="B-H","BH","R"), "_M", ifelse(estM[i], "E","F"), "_relbias_R.png")), plt, width = 20, height = 12, units = "in")
}

