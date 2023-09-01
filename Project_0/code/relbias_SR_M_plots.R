x = readRDS(file.path(here::here(),"Project_0", "results", paste0("naa","_om"), paste0("om_", 1), paste0("sim_",1,".RDS")))
get_SR_estimates_fn <- function(om_type = "naa"){
  df.ems = readRDS(here::here("Project_0","inputs", "df.ems.RDS"))
  df.oms = readRDS(here("Project_0","inputs", paste0("df.", ifelse(om_type == "naa", "", paste0(om_type,".")),"oms.RDS")))
  if(om_type == "naa") {
    em_ind <- 1:20
  }
  if(om_type == "M") {
    em_ind <- 5:24
  }
  if(om_type == "Sel") {
    em_ind <- c(5:20,25:28)
  }
  if(om_type == "q") {
    em_ind <- c(5:20,28:32)
  }
  df.ems <- df.ems[em_ind,]

  all_out <- lapply(1:NROW(df.oms), function(y){
    if(om_type == "naa") {
      if(is.na(df.oms$NAA_sig[y])) em_ind <- 1:4
      else em_ind <- 5:8
      print("naa")
      print(em_ind)
    }
    if(om_type == "M") {
      if(df.oms$M_cor[y] == 0) em_ind <- which(df.ems$re_config == "M_re" & df.ems$M_re_cor == "iid")
      else em_ind <- which(df.ems$re_config == "M_re" & df.ems$M_re_cor == "ar1_y")
    }
    if(om_type == "Sel") {
      if(df.oms$Sel_cor[y] == 0) em_ind <- which(df.ems$re_config == "sel_re" & df.ems$sel_re_cor == "iid")
      else em_ind <- which(df.ems$re_config == "sel_re" & df.ems$sel_re_cor == "ar1_y")
    }
    if(om_type == "q") {
      if(df.oms$q_cor[y] == 0) em_ind <- which(df.ems$re_config == "q_re" & df.ems$q_re_cor == "iid")
      else em_ind <- which(df.ems$re_config == "q_re" & df.ems$M_re_cor == "ar1")
    }
    em_ind <- em_ind[which(df.ems$SR_model[em_ind] == 3)]
    res <- lapply(1:100, function(x){
      print(paste0("om_", y, ", sim_",x))
      sim = readRDS(file.path(here::here(),"Project_0", "results", paste0(om_type,"_om"), paste0("om_", y), paste0("sim_",x,".RDS")))
      relres <- lapply(sim[em_ind],function(z) {
        out <- matrix(NA, 2,3)
        if(length(z)) if(length(z$fit)) {
          # print(names(z))
          # print(names(z$fit))
          # print(names(z$fit$rep))
          # print(exp(z$fit$rep$log_SR_a[1])/exp(z$truth$log_SR_a[1]))
          out[1,1] <- exp(z$fit$rep$log_SR_a[1])/exp(z$truth$log_SR_a[1])
          out[2,1] <- exp(z$fit$rep$log_SR_b[1])/exp(z$truth$log_SR_b[1])
        }
        if(!is.null(z$fit$sdrep)){
          #ind <- which(!is.na(z$fit$sdrep$SE_rep$log_SSB))
            out[1,2] <- z$fit$sdrep$SE_rep$log_SR_a[1] #se
            out[2,2] <- z$fit$sdrep$SE_rep$log_SR_b[1] #se
            ci <- matrix(c(z$fit$rep$log_SR_a[1],z$fit$rep$log_SR_b[1]) - c(z$truth$log_SR_a[1],z$truth$log_SR_b[1]),2,2) + 
              qnorm(0.975) * cbind(-out[,2],out[,2])
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
df.ems <- readRDS(here::here("Project_0","inputs", "df.ems.RDS"))
temp <- get_SR_estimates_fn("naa")
temp <- get_SR_estimates_fn(df.oms,"naa", df.ems)

make_results <- FALSE
if(make_results) {
  naa_relSR_results <- get_SR_estimates_fn("naa")
  saveRDS(naa_relSR_results, file = here("Project_0","results", "naa_relSR_results.RDS"))
  M_relSR_results <- get_SR_estimates_fn("M")
  saveRDS(M_relSR_results, file = here("Project_0","results", "M_relSR_results.RDS"))
  Sel_relSR_results <- get_SR_estimates_fn("Sel")
  saveRDS(Sel_relSR_results, file = here("Project_0","results", "Sel_relSR_results.RDS"))
  q_relSR_results <- get_SR_estimates_fn("q")
  saveRDS(q_relSR_results, file = here("Project_0","results", "q_relSR_results.RDS"))
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
make_plot_df <- function(type = "naa", res = all_naa_om_relssb, is_SE = FALSE) {
  library(reshape2)
  res <- melt(res)
  names(res) <- c("par", "column", "value", "em", "sim","om")
  res <- res %>% mutate(par = recode(par,
      "1" = "a",
      "2" = "b",
    )) %>% as.data.frame   
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

all_naa_mohns_rho <-  readRDS(file = here("Project_0","results", "all_naa_mohns_rho_results.RDS"))
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
df <- make_plot_df(type = "naa", res = all_naa_mohns_rho)

# naa_om_df.ems <- get_df_ems()
# M_om_df.ems <- get_df_ems("M")
# Sel_om_df.ems <- get_df_ems("Sel")
# q_om_df.ems <- get_df_ems("q")
