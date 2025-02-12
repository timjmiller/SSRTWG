library(here)
library(ggplot2)
library(ggh4x)
library(dplyr)
library(tidyr)
library(ggpattern)
df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))


all_naa_aic <- readRDS(file = file.path(here(),"Project_0","results", "all_naa_aic_results.RDS"))
aic_fn <- function(all, est_ind, om_ind = NULL){
  if(!is.null(om_ind)) all <- all[om_ind]
  out <- sapply(all, function(res){
    tmp = apply(res[est_ind,],2, function(x) {
      if(any(!is.na(x))) {
        return(x == min(x,na.rm=T))
      } else return(rep(NA, length(x)))
    })
    return(apply(tmp,1,sum,na.rm=T))
  })  
  return(out)
}

make_plot_df <- function(aic_res, df.ems, df.oms, em_ind = 1:20, SR_est =FALSE, M_est = FALSE){
  use.df.ems <- df.ems[em_ind,]
  SR_M_ind = which(use.df.ems$SR_model == ifelse(SR_est,3,2) & use.df.ems$M_est == M_est)
  aic_res <- t(aic_fn(aic_res, SR_M_ind))
  df <- cbind(df.oms[-1], aic_res)
  print(df)
  df <- df %>% pivot_longer(cols = as.character(1:5), names_to = "EM", values_to = "n") %>% as.data.frame
  this.df.ems <- use.df.ems[SR_M_ind,]
  print(SR_M_ind)
  print(this.df.ems)
  print(df$EM)
  df <- cbind(this.df.ems[as.integer(df$EM),],df)
  print(head(df))
  df$correct_EM_PE <- "No"
  if(max(em_ind)==20) {
    df <- df %>% mutate(EM_process_error = recode(EM,
        "1" = "R",
        "2" = "R+S",
        "3" = "R+M (iid)",
        "4" = "R+Sel (iid)",
        "5" = "R+q (iid)"
      ))
    df$NAA_sig[which(is.na(df$NAA_sig))] <- 0
    df$correct_EM_PE[df$NAA_sig == 0 & df$EM_process_error == "R"] <- "Yes"
    df$correct_EM_PE[df$NAA_sig >  0 & df$EM_process_error == "R+S"] <- "Yes"
    df <- df %>% 
      mutate(NAA_sig = recode(NAA_sig,
        "0" = "sigma['2+'] == 0",
        "0.25" = "sigma['2+'] == 0.25",
        "0.5" = "sigma['2+'] == 0.5")) %>% 
      mutate(R_sig = recode(R_sig,
        "0.5" = "sigma[italic(R)] == 0.5",
        "1.5" = "sigma[italic(R)] == 1.5")) %>% as.data.frame
  }
  if(max(em_ind)==24) {
    df <- df %>% 
      mutate(EM_process_error = recode(EM,
        "1" = "R+S",
        "2" = "R+M (iid)",
        "3" = "R+Sel (iid)",
        "4" = "R+q (iid)",
        "5" = "R+M (AR1)"
      )) 
    df$correct_EM_PE[df$M_cor == 0 & df$EM_process_error == "R+M (iid)"] <- "Yes"
    df$correct_EM_PE[df$M_cor >  0 & df$EM_process_error == "R+M (AR1)"] <- "Yes"
    df <- df %>% 
      mutate(M_sig = recode(M_sig,
        "0.1" = "sigma[italic(M)] == 0.1",
        "0.5" = "sigma[italic(M)] == 0.5")) %>% 
      mutate(M_cor = recode(M_cor,
        "0" = "rho[italic(M)] == 0",
        "0.9" = "rho[italic(M)] == 0.9")) %>% as.data.frame
  }
  if(max(em_ind)==28) {
    df <- df %>% 
      mutate(EM_process_error = recode(EM,
        "1" = "R+S",
        "2" = "R+M (iid)",
        "3" = "R+Sel (iid)",
        "4" = "R+q (iid)",
        "5" = "R+Sel (AR1)"
      ))    
    df$correct_EM_PE[df$Sel_cor == 0 & df$EM_process_error == "R+Sel (iid)"] <- "Yes"
    df$correct_EM_PE[df$Sel_cor >  0 & df$EM_process_error == "R+Sel (AR1)"] <- "Yes"
    df <- df %>% 
      mutate(Sel_sig = recode(Sel_sig,
        "0.1" = "sigma[Sel] == 0.1",
        "0.5" = "sigma[Sel] == 0.5")) %>% 
      mutate(Sel_cor = recode(Sel_cor,
        "0" = "rho[Sel] == 0",
        "0.9" = "rho[Sel] == 0.9")) %>% as.data.frame
  }
  if(max(em_ind)==32) {
    df <- df %>% 
      mutate(EM_process_error = recode(EM,
        "1" = "R+S",
        "2" = "R+M (iid)",
        "3" = "R+Sel (iid)",
        "4" = "R+q (iid)",
        "5" = "R+q (AR1)"
      ))
    df$correct_EM_PE[df$q_cor == 0 & df$EM_process_error == "R+q (iid)"] <- "Yes"
    df$correct_EM_PE[df$q_cor >  0 & df$EM_process_error == "R+q (AR1)"] <- "Yes"
    df <- df %>% 
      mutate(q_sig = recode(q_sig,
        "0.1" = "sigma[italic(q)] == 0.1",
        "0.5" = "sigma[italic(q)] == 0.5")) %>%
      mutate(q_cor = recode(q_cor,
        "0" = "rho[italic(q)] == 0",
        "0.9" = "rho[italic(q)] == 0.9"))  %>% as.data.frame
  }
  df <- df %>% mutate(Fhist = recode(Fhist,
        "H-MSY"  = "2.5*italic(F)[MSY] %->% italic(F)[MSY]",
        "MSY" = "italic(F)[MSY]"))
  df <- df %>% mutate(obs_error = recode(obs_error,
      "L" = "Low observation error",
      "H" = "High observation error"))

  return(df)
}
 # make_plot_df(all_aic, df.ems, df.oms, em_ind = em_inds[,i], SR_est = est_conds[j,1], M_est = est_conds[j,2])


df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))
# em_ind <- c(5:20,25:28)
# use.df.ems <- df.ems[em_ind,]
types <- c("naa", "M", "Sel", "q")
types.plt <- c("R, R+S","R+M", "R+Sel", "R+q")
est_conds <- cbind(c(FALSE,FALSE,TRUE,TRUE), c(FALSE,TRUE,FALSE,TRUE))
em_inds <- cbind(1:20, 5:24, c(5:20,25:28), c(5:20, 29:32))
est_config <- c("M known, No SR", "M estimated, No SR", "M known, SR", "M estimated, SR")
EM_process_error <- c("R","R+S","R+M (iid)","R+M (AR1)","R+Sel (iid)","R+Sel (AR1)","R+q (iid)","R+q (AR1)")

out <- list()
for(i in 1:4){
  all_aic <- readRDS(file = here("Project_0","results", paste0("all_", types[i], "_aic_results.RDS")))
  df.oms = readRDS(here("Project_0","inputs", paste0("df.", ifelse(i==1, "", paste0(types[i],".")),"oms.RDS")))
  if(i>1) re_nms <- paste0(types[i],c("_sig","_cor"))
  if(i ==1) re_nms <- c("R_sig","NAA_sig")
  for(j in 1:4){
    print(c(i,j))
    if(j == 1) out[[i]] <- cbind(make_plot_df(all_aic, df.ems, df.oms, em_ind = em_inds[,i], SR_est = est_conds[j,1], M_est = est_conds[j,2]), est_config = est_config[j])
    else out[[i]] <- rbind(out[[i]], cbind(make_plot_df(all_aic, df.ems, df.oms, em_ind = em_inds[,i], SR_est = est_conds[j,1], M_est = est_conds[j,2]), est_config = est_config[j]))
    print(head(out[[i]]))

  }
  out[[i]]$EM_process_error <- factor(out[[i]]$EM_process_error, levels = EM_process_error)
  out[[i]]$est_config <- factor(out[[i]]$est_config, levels = est_config)
}
names(out) <- types
theme_set(theme_bw())
theme_update(strip.text = element_text(size = rel(2)), strip.placement = "outside", strip.background = element_rect(), #fill = "transparent"), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)))

R_S_aic_plt <- ggplot(out$naa, aes(x = est_config, y = n)) + 
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
  facet_nested(obs_error+ Fhist ~ R_sig + NAA_sig, 
    labeller = labeller(R_sig = label_parsed, NAA_sig = label_parsed, Fhist = label_parsed)) +
  #geom_col(position = "fill", mapping = aes(fill = EM)) + 
  coord_cartesian(ylim = c(0, 1)) + ylab("Probability lowest AIC") + xlab("EM M and Stock Recruit Assumptions") +
  geom_col_pattern(mapping = aes(fill = EM_process_error, pattern = correct_EM_PE), position ="fill",
   pattern_color = "black",
   pattern_fill = "black",
   pattern_angle = 45,
   pattern_spacing = 0.2,
   pattern_density = 0.05
   # pattern_spacing = 0.025,
   # pattern_key_scale_factor = 0.6
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0, pattern = "none"), order = 1), pattern = "none") +
  scale_pattern_manual(values = c(Yes = "stripe", No = "none")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(fill = "EM Process Error")
R_S_aic_plt

R_M_aic_plt <- ggplot(out$M, aes(x = est_config, y = n)) + 
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
  facet_nested(obs_error+ Fhist ~ M_sig + M_cor, 
    labeller = labeller(M_sig = label_parsed, M_cor = label_parsed, Fhist = label_parsed)) +
  #geom_col(position = "fill", mapping = aes(fill = EM)) + 
  coord_cartesian(ylim = c(0, 1)) + ylab("Probability lowest AIC") + xlab("EM M and Stock Recruit Assumptions") +
  geom_col_pattern(mapping = aes(fill = EM_process_error, pattern = correct_EM_PE), position ="fill",
   pattern_color = "black",
   pattern_fill = "black",
   pattern_angle = 45,
   pattern_spacing = 0.2,
   pattern_density = 0.05) +
  guides(fill = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0, pattern = "none"), order = 1), pattern = "none") +
  scale_pattern_manual(values = c(Yes = "stripe", No = "none")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(fill = "EM Process Error")
R_M_aic_plt

R_Sel_aic_plt <- ggplot(out$Sel, aes(x = est_config, y = n)) + 
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
  facet_nested(obs_error+ Fhist ~ Sel_sig + Sel_cor, 
    labeller = labeller(Sel_sig = label_parsed, Sel_cor = label_parsed, Fhist = label_parsed)) +
  #geom_col(position = "fill", mapping = aes(fill = EM)) + 
  coord_cartesian(ylim = c(0, 1)) + ylab("Probability lowest AIC") + xlab("EM M and Stock Recruit Assumptions") +
  geom_col_pattern(mapping = aes(fill = EM_process_error, pattern = correct_EM_PE), position ="fill",
   pattern_color = "black",
   pattern_fill = "black",
   pattern_angle = 45,
   pattern_spacing = 0.2,
   pattern_density = 0.05) +
  guides(fill = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0, pattern = "none"), order = 1), pattern = "none") +
  scale_pattern_manual(values = c(Yes = "stripe", No = "none")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(fill = "EM Process Error")
R_Sel_aic_plt

R_q_aic_plt <- ggplot(out$q, aes(x = est_config, y = n)) + 
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
  facet_nested(obs_error+ Fhist ~ q_sig + q_cor, 
    labeller = labeller(q_sig = label_parsed, q_cor = label_parsed, Fhist = label_parsed)) +
  #geom_col(position = "fill", mapping = aes(fill = EM)) + 
  coord_cartesian(ylim = c(0, 1)) + ylab("Probability lowest AIC") + xlab("EM M and Stock Recruit Assumptions") +
  geom_col_pattern(mapping = aes(fill = EM_process_error, pattern = correct_EM_PE), position ="fill",
   pattern_color = "black",
   pattern_fill = "black",
   pattern_angle = 45,
   pattern_spacing = 0.2,
   pattern_density = 0.05) +
  guides(fill = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0, pattern = "none"), order = 1), pattern = "none") +
  scale_pattern_manual(values = c(Yes = "stripe", No = "none")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(fill = "EM Process Error")
R_q_aic_plt

#all in one figure

library(patchwork)
theme_set(theme_bw())
theme_update(strip.placement = "outside", strip.background = element_rect(), title = element_text(size = rel(1.5)))
R_Sel_aic_plt_alt <- R_Sel_aic_plt + 
# geom_col(position = "fill", show.legend = TRUE)
  geom_col_pattern(mapping = aes(fill = EM_process_error, pattern = correct_EM_PE), position ="fill", show.legend = TRUE,
   pattern_color = "black",
   pattern_fill = "black",
   pattern_angle = 45,
   pattern_spacing = 0.2,
   pattern_density = 0.05)
cairo_pdf(here("Project_0","manuscript", "pe_aic_plots.pdf"), width = 30*2/3, height = 20*2/3)
design <- c(area(1,1,1,1), area(2,1,2,1), area(1,2,1,2), area(2,2,2,2))
(R_S_aic_plt +  xlab("") + labs(title = "A: R, R+S OMs") + theme(axis.title.x = element_blank(), strip.text.y=element_blank(), axis.text.x=element_blank(), legend.position="none")) +# , strip.text=element_text(margin=margin(b = 5)), plot.margin = margin(b = 1, t = 0))) + 
(R_M_aic_plt + labs(title = "B: R+M OMs") + theme(legend.position="none", strip.text.y=element_blank())) + #axis.title.x = element_blank(), strip.text=element_text(margin=margin(b = 5)), axis.text.x=element_blank(), plot.margin = margin(b = 1, t = 0))) + 
(R_Sel_aic_plt_alt + xlab("") + labs(title = "C: R+Sel OMs") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())) +# , strip.text=element_text(margin=margin(b = 5)), plot.margin = margin(b = 1, t = 0))) + 
(R_q_aic_plt + labs(title = "D: R+q OMs") + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position = "none")) + 
  plot_layout(design = design, axis_titles = "collect")
dev.off()

##########################################################################################################################################
# Stock-recruit AIC
##########################################################################################################################################


coefficient_fn <- function(om = 1, df){
  fit <- unname(summary(glm(cbind(est_SR,est_R) ~ log(sd_log_SSB), family = binomial, data = subset(df, OM == om)))$coefficients[2,])
  return(c(Estimate = fit[1], SE = fit[2], lo = fit[1] + qnorm(0.025)*fit[2], hi = fit[1] + qnorm(0.975)*fit[2]))
}
coefficient_fn <- function(om = 1, df, col = "BH_best_ME"){
  fit <- unname(summary(glm(get(col) ~ log(sd_log_SSB), family = binomial, data = subset(df, OM == om)))$coefficients[2,])
  return(c(Estimate = fit[1], SE = fit[2], lo = fit[1] + qnorm(0.025)*fit[2], hi = fit[1] + qnorm(0.975)*fit[2]))
}


df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.oms.RDS"))
all_naa_sd_log_SSB <- readRDS(file = here("Project_0","results", "all_naa_sd_log_SSB.RDS"))
all_naa_sd_log_SSB <- t(matrix(unname(unlist(all_naa_sd_log_SSB)), nrow = 100))
colnames(all_naa_sd_log_SSB) <- paste0("sim", 1:100)
df <- cbind(df.oms, all_naa_sd_log_SSB)
df <- df %>% pivot_longer(cols = paste0("sim", 1:100), names_to = "sim", values_to = "sd_log_SSB")

om_ind <- 1:24
#SR_rec_Mest = which(use.df.ems$re_config == "rec+1" & use.df.ems$M_est == TRUE)
use.df.ems <- df.ems[1:20,]
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

df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))
# em_ind <- c(5:20,25:28)
# use.df.ems <- df.ems[em_ind,]
types <- c("naa", "M", "Sel", "q")
types.plt <- c("R, R+S","R+M", "R+Sel", "R+q")
est_conds <- cbind(c(FALSE,TRUE))
em_inds <- cbind(1:20, 5:24, c(5:20,25:28), c(5:20, 29:32))
pred_dfs <- list()

for(i in 1:4) {
  all_aic <- readRDS(file = here("Project_0","results", paste0("all_", types[i], "_aic_results.RDS")))
  df.oms = readRDS(here("Project_0","inputs", paste0("df.", ifelse(i==1, "", paste0(types[i],".")),"oms.RDS")))
  use.df.ems <- df.ems[em_inds[,i],]
  df.oms$OM <- 1:NROW(df.oms)
  sd_log_SSB <- readRDS(file = here("Project_0","results", paste0("all_", types[i], "_sd_log_SSB.RDS")))
  sd_log_SSB <- t(matrix(unname(unlist(sd_log_SSB)), nrow = 100))
#  colnames(sd_log_SSB) <- paste0("sim", 1:100)
  df <- cbind(df.oms, sd_log_SSB)
  pred_df <- cbind(df.oms, t(sapply(1:NROW(df.oms), \(x) seq(min(sd_log_SSB[x,]), max(sd_log_SSB[x,]),length.out = 100))))
  df <- df %>% pivot_longer(cols = paste0(1:100), names_to = "sim", values_to = "sd_log_SSB")
  pred_df <- pred_df %>% pivot_longer(cols = paste0(1:100), names_to = NULL, values_to = "sd_log_SSB")
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
        "0.5" = "sigma[italic(R)] == 0.5",
        "1.5" = "sigma[italic(R)] == 1.5")) %>% as.data.frame
  }
  if(i ==2){
    df <- df %>% mutate(M_sig = recode(M_sig,
        "0.1" = "sigma[italic(M)] == 0.1",
        "0.5" = "sigma[italic(M)] == 0.5")) %>% 
      mutate(M_cor = recode(M_cor,
        "0" = "rho[italic(M)] == 0",
        "0.9" = "rho[italic(M)] == 0.9"))
  }
  if(i == 3){
    df <- df %>% mutate(Sel_sig = recode(Sel_sig,
        "0.1" = "sigma[Sel] == 0.1",
        "0.5" = "sigma[Sel] == 0.5")) %>% 
      mutate(Sel_cor = recode(Sel_cor,
        "0" = "rho[Sel] == 0",
        "0.9" = "rho[Sel] == 0.9"))
  }
  if(i == 4){
    df <- df %>% mutate(q_sig = recode(q_sig,
        "0.1" = "sigma[italic(q)] == 0.1",
        "0.5" = "sigma[italic(q)] == 0.5")) %>%
      mutate(q_cor = recode(q_cor,
        "0" = "rho[italic(q)] == 0",
        "0.9" = "rho[italic(q)] == 0.9"))  
  }
  df <- df %>% mutate(Fhist = recode(Fhist,
        "H-MSY"  = "2.5*italic(F)[MSY] %->% italic(F)[MSY]",
        "MSY" = "italic(F)[MSY]"))
  df <- df %>% mutate(obs_error = recode(obs_error,
      "L" = "Low observation error",
      "H" = "High observation error"))
  df <- df %>% as.data.frame
  facs <- c("OM", "Fhist","NAA_sig", "R_sig", "M_sig", "M_cor", "Sel_sig", "Sel_cor", "q_sig", "q_cor", "obs_error")
  fac_names <- names(df)[names(df) %in% facs]
  df[fac_names] <- lapply(df[fac_names], as.factor)
  pred_df[fac_names] <- df[fac_names]
  pred_df <- pred_df %>% as.data.frame
  pred_df <- rbind(cbind(pred_df, type = "MF"), cbind(pred_df, type = "ME")) %>% as.data.frame

  fit_MF <- glm(BH_best_MF ~ factor(OM) * log(sd_log_SSB), family = binomial, data = df)
  fit_ME <- glm(BH_best_ME ~ factor(OM) * log(sd_log_SSB), family = binomial, data = df)
  #x <- t(sapply(1:24, coefficient_fn, df =df, col = "BH_best_MF"))
  df$pred <- predict(fit_MF, newdata = df, type = "response")
  predict_MF <- predict(fit_MF, newdata = subset(pred_df, type == "ME"), type = "link", se.fit = TRUE)
  predict_ME <- predict(fit_ME, newdata = subset(pred_df, type == "ME"), type = "link", se.fit = TRUE)
  pred_df$pred_link <- c(predict_MF$fit,predict_ME$fit)
  pred_df$pred_sd <- c(predict_MF$se.fit,predict_ME$se.fit)
  pred_df$pred <- 1/(1 + exp(-pred_df$pred_link))
  pred_df$pred_lo  <- 1/(1 + exp(-(pred_df$pred_link + qnorm(0.025)*pred_df$pred_sd)))
  pred_df$pred_hi  <- 1/(1 + exp(-(pred_df$pred_link + qnorm(0.975)*pred_df$pred_sd)))
  pred_df$obs <- df$sd_log_SSB
  pred_df <- pred_df %>% mutate(type = recode(type,
      "MF" = "M known",
      "ME" = "M estimated"))
  pred_df$log_sd_log_SSB <- log(pred_df$sd_log_SSB)
  pred_df$log_obs <- log(pred_df$obs)
  pred_dfs[[i]] <- pred_df
}
    R_S_sr_aic_plt <- ggplot(pred_dfs[[1]], aes(x = log_sd_log_SSB, y = pred, colour = type, fill =  type)) + 
      facet_nested(obs_error+ Fhist ~ R_sig + NAA_sig, 
        labeller = labeller(R_sig = label_parsed, NAA_sig = label_parsed, Fhist = label_parsed)) +
      scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
      scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
      geom_line(linewidth = 2, alpha = 0.5) + labs(fill = "EM M Assumption", colour = "EM M Assumption") +
      geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), linetype = 0, alpha = 0.4) +
      geom_rug(data = subset(pred_dfs[[1]], type == "M known"), mapping = aes(x=log_obs), colour = "black", position = "jitter", sides = "b", alpha = 0.5) +
      xlab("log(SD(log(SSB)))") + ylab("Probability SR assumption with lowest AIC") + xlim(-3,1)
    R_S_sr_aic_plt

    R_M_sr_aic_plt <- ggplot(pred_dfs[[2]], aes(x = log_sd_log_SSB, y = pred, colour = type, fill =  type)) + 
      facet_nested(obs_error+ Fhist ~ M_sig + M_cor, 
        labeller = labeller(M_sig = label_parsed, M_cor = label_parsed, Fhist = label_parsed)) +
      scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
      scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
      geom_line(linewidth = 2, alpha = 0.5) + labs(fill = "EM M Assumption", colour = "EM M Assumption") +
      geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), linetype = 0, alpha = 0.4) +
      geom_rug(data = subset(pred_dfs[[2]], type == "M known"), mapping = aes(x=log_obs), colour = "black", position = "jitter", sides = "b", alpha = 0.5) +
      xlab("log(SD(log(SSB)))") + ylab("Probability SR assumption with lowest AIC") + xlim(-3,1)
    R_M_sr_aic_plt

    R_Sel_sr_aic_plt <- ggplot(pred_dfs[[3]], aes(x = log_sd_log_SSB, y = pred, colour = type, fill =  type)) + 
      facet_nested(obs_error+ Fhist ~ Sel_sig + Sel_cor, 
        labeller = labeller(Sel_sig = label_parsed, Sel_cor = label_parsed, Fhist = label_parsed)) +
      scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
      scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
      geom_line(linewidth = 2, alpha = 0.5) + labs(fill = "EM M Assumption", colour = "EM M Assumption") +
      geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), linetype = 0, alpha = 0.4) +
      geom_rug(data = subset(pred_dfs[[3]], type == "M known"), mapping = aes(x=log_obs), colour = "black", position = "jitter", sides = "b", alpha = 0.5) +
      xlab("log(SD(log(SSB)))") + ylab("Probability SR assumption with lowest AIC") + xlim(-3,1)
    R_Sel_sr_aic_plt

    R_q_sr_aic_plt <- ggplot(pred_dfs[[4]], aes(x = log_sd_log_SSB, y = pred, colour = type, fill =  type)) + 
      facet_nested(obs_error+ Fhist ~ q_sig + q_cor, 
        labeller = labeller(q_sig = label_parsed, q_cor = label_parsed, Fhist = label_parsed)) +
      scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
      scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
      geom_line(linewidth = 2, alpha = 0.5) + labs(fill = "EM M Assumption", colour = "EM M Assumption") +
      geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), linetype = 0, alpha = 0.4) +
      geom_rug(data = subset(pred_dfs[[4]], type == "M known"), mapping = aes(x=log_obs), colour = "black", position = "jitter", sides = "b", alpha = 0.5) +
      xlab("log(SD(log(SSB)))") + ylab("Probability SR assumption with lowest AIC") + xlim(-3,1)
    R_q_sr_aic_plt

theme_set(theme_bw())
theme_update(strip.placement = "outside", strip.background = element_rect(), title = element_text(size = rel(1.5)))
cairo_pdf(here("Project_0","manuscript", "sr_aic_plots.pdf"), width = 30*2/3, height = 20*2/3)
design <- c(area(1,1,1,1), area(2,1,2,1), area(1,2,1,2), area(2,2,2,2))
(R_S_sr_aic_plt +  xlab("") + labs(title = "A: R, R+S OMs") + theme(axis.title.x = element_blank(), strip.text.y=element_blank(), axis.text.x=element_blank(), legend.position="none")) +# , strip.text=element_text(margin=margin(b = 5)), plot.margin = margin(b = 1, t = 0))) + 
(R_M_sr_aic_plt + labs(title = "B: R+M OMs") + theme(legend.position="none", strip.text.y=element_blank())) + #axis.title.x = element_blank(), strip.text=element_text(margin=margin(b = 5)), axis.text.x=element_blank(), plot.margin = margin(b = 1, t = 0))) + 
(R_Sel_sr_aic_plt + xlab("") + labs(title = "C: R+Sel OMs") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())) +# , strip.text=element_text(margin=margin(b = 5)), plot.margin = margin(b = 1, t = 0))) + 
(R_q_sr_aic_plt + labs(title = "D: R+q OMs") + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position = "none")) + 
  plot_layout(design = design, axis_titles = "collect")
dev.off()
