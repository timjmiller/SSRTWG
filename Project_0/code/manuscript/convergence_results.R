library(here)
library(ggplot2)
library(ggh4x)
library(dplyr)
library(tidyr)
library(patchwork)
df.ems = readRDS(here("Project_0", "inputs", "df.ems.RDS"))

convergence_fn <- function(em){
  out <- rep(NA,5) #convergence type  1, and 2
  if(!is.character(em)) {
    if(!is.null(em$fit$opt)) {
      out[1] <- 1 
      out[2] <- em$fit$opt$conv
    }
    if(!is.null(em$fit$sdrep)) out[3] <- as.integer(sum(sapply(em$fit$sdrep$SE_par, function(g) any(is.nan(g)))))
    if(!is.null(em$fit$final_gradient)) out[4] <- max(abs(em$fit$final_gradient))
    if(!is.null(em$fit$sdrep)) {
      maxs <- sapply(em$fit$sdrep$SE_par, function(g) ifelse(any(!is.na(g)), max(g,na.rm=TRUE), NA))
      out[5] <- ifelse(any(!is.na(maxs)), max(maxs, na.rm =TRUE), NA)
    }
  }
  return(out)
}
df.oms = readRDS(here("Project_0", "inputs", "df.oms.RDS"))

#convergence information:
#1: 0/1 model finished optimizing
#2: nlminb convergence flag
#3: max gradient value
#4: number of NaNs in SEs for parameters, 0 = good invertible hessian
#5: maximum non-NaN SE estimate


conv_fn <- function(all, est_ind, om_ind = NULL){
  #est_ind is 1 ROW of df.ems.
  if(!is.null(om_ind)) all <- all[om_ind]
  out <- sapply(all, function(res){
    tmp <- sapply(res, function(x) {
      return(x[est_ind,])
    })
    x <- t(tmp)
    nconv <- c(
      sum(!is.na(x[,1])), 
      sum(!is.na(x[,2]) & x[,2] == 0), 
      sum(!is.na(x[,4]) & x[,4] < 1e-6),
      sum(!is.na(x[,3]) & x[,3] == 0),
      sum(!is.na(x[,5]) & x[,3] == 0 & x[,5] < 100)
    )
    nconv <- c(nconv, NROW(x))
    return(nconv)
  })
  return(t(out))
}


conv_res_plotting_fn <- function(conv_res, df.oms, df.ems, em_ind, om_type = "naa"){

  om_ind <- 1:NROW(df.oms)
  #em_ind <- which(df.ems$M_est == M_est & df.ems$SR_model == ifelse(SR_est, 3,2))
  for(i in em_ind){
    conv <- conv_fn(conv_res, i, om_ind)
    colnames(conv) <- c(paste0("Type", 1:5, "_n_pass"), "n_sim")
    res <- df.oms[om_ind,]
    for(k in names(df.ems)) res <- cbind(res, df.ems[[k]][i])
    names(res)[NCOL(df.oms)+1:NCOL(df.ems)] <- names(df.ems)
    res <- cbind(res, conv)
    if(i == em_ind[1]){
      df <- res
    } else{
      df <- rbind(df, res)

    }
  }

  for(i in 1:5){
    df[[paste0("Type", i, "_p_pass")]] <- df[[paste0("Type", i, "_n_pass")]]/df$n_sim
    df[[paste0("Type", i, "_ci_lo")]] <- sapply(1:NROW(df), function(x) {
      binom.test(df[[paste0("Type", i, "_n_pass")]][x], df$n_sim[x], p = df[[paste0("Type", i, "_p_pass")]][x])$conf.int[1]
    })
    df[[paste0("Type", i, "_ci_hi")]]  <- sapply(1:NROW(df), function(x) {
      binom.test(df[[paste0("Type", i, "_n_pass")]][x], df$n_sim[x], p = df[[paste0("Type", i, "_p_pass")]][x])$conf.int[2]
    })
  }
  # facs <- c("Ecov_obs_sig", "Fhist", "Ecov_re_sig","Ecov_re_cor", "obs_error", "NAA_M_re", "re_config")
  # df[facs] <- lapply(df[facs], factor)
  # df_mod <- df %>%
  df <- df %>%
    mutate(obs_error = recode(obs_error,
      "L" = "Low observation error",
      "H" = "High observation error"
    ))
  df <- df %>%
    mutate(EM_process_error = recode(re_config,
      "rec" = "R",
      "rec+1" = "R+S",
      "M_re" = "R+M",
      "q_re" = "R+q",
      "sel_re" = "R+Sel"
    ))
  ind <- which(df$M_re_cor == "iid" | df$sel_re_cor == "iid" | df$q_re_cor == "iid")
  df$EM_process_error[ind] <- paste0(df$EM_process_error[ind], " (iid)")
  ind <- which(df$M_re_cor == "ar1_y" | df$sel_re_cor == "ar1_y" | df$q_re_cor == "ar1")
  df$EM_process_error[ind] <- paste0(df$EM_process_error[ind], " (AR1)")

  df <- df %>% mutate(Fhist = recode(Fhist,
        "H-MSY"  = "2.5*italic(F)[MSY] %->% italic(F)[MSY]",
        "MSY" = "italic(F)[MSY]"))

  df$correct_EM_PE <- "No"
  if(om_type == "naa") {
    df$NAA_sig[which(is.na(df$NAA_sig))] <- 0
    df$correct_EM_PE[df$NAA_sig == 0 & df$EM_process_error == "R"] <- "Yes"
    df$correct_EM_PE[df$NAA_sig >  0 & df$EM_process_error == "R+S"] <- "Yes"
    df <- df %>% mutate(R_sig = recode(R_sig,
          "0.5" = "sigma[R] == 0.5",
          "1.5" = "sigma[R] == 1.5"))
    df <- df %>% mutate(NAA_sig = recode(NAA_sig,
        "0" = 'sigma["2+"] == 0',
        "0.25" = 'sigma["2+"] == 0.25',
        "0.5" = 'sigma["2+"] == 0.5'))
  }
  if(om_type == "M") {
    df$correct_EM_PE[df$M_cor == 0 & df$EM_process_error == "R+M (iid)"] <- "Yes"
    df$correct_EM_PE[df$M_cor >  0 & df$EM_process_error == "R+M (AR1)"] <- "Yes"
    df <- df %>% mutate(M_sig = recode(M_sig,
        "0.1" = 'sigma[italic(M)] == 0.1',
        "0.5" = 'sigma[italic(M)] == 0.5'))
    df <- df %>% mutate(M_cor = recode(M_cor,
        "0" = 'rho[italic(M)] == 0',
        "0.9" = 'rho[italic(M)] == 0.9'))
  }
  if(om_type == "Sel") {
    df$correct_EM_PE[df$Sel_cor == 0 & df$EM_process_error == "R+Sel (iid)"] <- "Yes"
    df$correct_EM_PE[df$Sel_cor >  0 & df$EM_process_error == "R+Sel (AR1)"] <- "Yes"
    df <- df %>% mutate(Sel_sig = recode(Sel_sig,
        "0.1" = 'sigma[Sel] == 0.1',
        "0.5" = 'sigma[Sel] == 0.5'))
    df <- df %>% mutate(Sel_cor = recode(Sel_cor,
        "0" = 'rho[Sel] == 0',
        "0.9" = 'rho[Sel] == 0.9'))
  }
  if(om_type == "q") {
    df$correct_EM_PE[df$q_cor == 0 & df$EM_process_error == "R+q (iid)"] <- "Yes"
    df$correct_EM_PE[df$q_cor >  0 & df$EM_process_error == "R+q (AR1)"] <- "Yes"
    df <- df %>% mutate(q_sig = recode(q_sig,
        "0.1" = 'sigma[italic(q)] == 0.1',
        "0.5" = 'sigma[italic(q)] == 0.5'))
    df <- df %>% mutate(q_cor = recode(q_cor,
        "0" = 'rho[italic(q)] == 0',
        "0.9" = 'rho[italic(q)] == 0.9'))
  }
  df <- df %>% pivot_longer(
    cols = Type1_p_pass:Type5_ci_hi,
    names_to = c("Type", "locate"),
    names_pattern = "Type(.)_(.*)" #magic
  )

  df <- df %>% pivot_wider(names_from = locate, values_from = value) %>% as.data.frame
  df$Type <- as.numeric(df$Type)
  return(df)
}

df.oms = readRDS(here("Project_0", "inputs", "df.oms.RDS"))
df.ems = readRDS(here("Project_0", "inputs", "df.ems.RDS"))

naa_om_conv_res <- readRDS(here("Project_0", "results", "naa_om_convergence_results.RDS"))

ind <- which(df.ems$M_est == FALSE & df.ems$SR_model == 2)
ind <- ind[which(ind<21)]
all_res_mod_Mknown_NoSR <- conv_res_plotting_fn(naa_om_conv_res, df.oms, df.ems, em_ind = ind)
ind <- which(df.ems$M_est == TRUE & df.ems$SR_model == 2)
ind <- ind[which(ind<21)]
all_res_mod_Mest_NoSR <- conv_res_plotting_fn(naa_om_conv_res, df.oms, df.ems, em_ind = ind)
ind <- which(df.ems$M_est == FALSE & df.ems$SR_model == 3)
ind <- ind[which(ind<21)]
all_res_mod_Mknown_SR <- conv_res_plotting_fn(naa_om_conv_res, df.oms, df.ems, em_ind = ind)
ind <- which(df.ems$M_est == TRUE & df.ems$SR_model == 3)
ind <- ind[which(ind<21)]
all_res_mod_Mest_SR <- conv_res_plotting_fn(naa_om_conv_res, df.oms, df.ems, em_ind = ind)
est_config <- c("M known, No SR", "M estimated, No SR", "M known, SR", "M estimated, SR")
all_res_mod <- rbind(
  cbind(all_res_mod_Mknown_NoSR, est_config = "M known, No SR"),
  cbind(all_res_mod_Mest_NoSR, est_config = "M estimated, No SR"),
  cbind(all_res_mod_Mknown_SR, est_config = "M known, SR"),
  cbind(all_res_mod_Mest_SR, est_config = "M estimated, SR")
)
all_res_mod$est_config <- factor(all_res_mod$est_config, levels = est_config)
EM_process_error <- c("R","R+S","R+M (iid)","R+M (AR1)","R+Sel (iid)","R+Sel (AR1)","R+q (iid)","R+q (AR1)")
all_res_mod$EM_process_error <- factor(all_res_mod$EM_process_error, levels = EM_process_error)

theme_set(theme_bw())
theme_update(strip.text = element_text(size = rel(2)), strip.placement = "outside", strip.background = element_rect(), #fill = "transparent"), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)))

head(all_res_mod)
temp <- subset(all_res_mod, Type == 4)
R_S_om_p_conv_type_4 <- ggplot(temp, aes(x = est_config, y = p_pass)) + 
    facet_nested(obs_error+ Fhist ~ R_sig + NAA_sig, 
      labeller = labeller(R_sig = label_parsed, NAA_sig = label_parsed, Fhist = label_parsed)) +
    coord_cartesian(ylim = c(0, 1)) + ylab("Probability of convergence") + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7)) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_S_om_p_conv_type_4
cairo_pdf(here("Project_0","manuscript", "R_S_om_p_convergence_type_4.pdf"), width = 24, height = 16)
R_S_om_p_conv_type_4
dev.off()

temp <- subset(all_res_mod, Type == 3)
R_S_om_p_conv_type_3 <- ggplot(temp, aes(x = est_config, y = p_pass)) + 
    facet_nested(obs_error+ Fhist ~ R_sig + NAA_sig, 
      labeller = labeller(R_sig = label_parsed, NAA_sig = label_parsed, Fhist = label_parsed)) +
    coord_cartesian(ylim = c(0, 1)) + ylab("Probability of convergence") + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7)) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_S_om_p_conv_type_3
cairo_pdf(here("Project_0","manuscript", "R_S_om_p_convergence_type_3.pdf"), width = 24, height = 16)
R_S_om_p_conv_type_3
dev.off()

df.oms = readRDS(here("Project_0", "inputs", "df.M.oms.RDS"))
df.ems = readRDS(here("Project_0", "inputs", "df.ems.RDS"))[5:24,]
M_om_conv_res <- readRDS(here("Project_0", "results", "M_om_convergence_results.RDS"))
ind <- which(df.ems$M_est == FALSE & df.ems$SR_model == 2)
all_res_mod_Mknown_NoSR <- conv_res_plotting_fn(M_om_conv_res, df.oms, df.ems, em_ind = ind, om_type = "M")
ind <- which(df.ems$M_est == TRUE & df.ems$SR_model == 2)
all_res_mod_Mest_NoSR <- conv_res_plotting_fn(M_om_conv_res, df.oms, df.ems, em_ind = ind, om_type = "M")
ind <- which(df.ems$M_est == FALSE & df.ems$SR_model == 3)
all_res_mod_Mknown_SR <- conv_res_plotting_fn(M_om_conv_res, df.oms, df.ems, em_ind = ind, om_type = "M")
ind <- which(df.ems$M_est == TRUE & df.ems$SR_model == 3)
all_res_mod_Mest_SR <- conv_res_plotting_fn(M_om_conv_res, df.oms, df.ems, em_ind = ind, om_type = "M")
est_config <- c("M known, No SR", "M estimated, No SR", "M known, SR", "M estimated, SR")
all_res_mod <- rbind(
  cbind(all_res_mod_Mknown_NoSR, est_config = "M known, No SR"),
  cbind(all_res_mod_Mest_NoSR, est_config = "M estimated, No SR"),
  cbind(all_res_mod_Mknown_SR, est_config = "M known, SR"),
  cbind(all_res_mod_Mest_SR, est_config = "M estimated, SR")
)
all_res_mod$est_config <- factor(all_res_mod$est_config, levels = est_config)
EM_process_error <- c("R","R+S","R+M (iid)","R+M (AR1)","R+Sel (iid)","R+Sel (AR1)","R+q (iid)","R+q (AR1)")
all_res_mod$EM_process_error <- factor(all_res_mod$EM_process_error, levels = EM_process_error)

temp <- subset(all_res_mod, Type == 4)
R_M_om_p_conv_type_4 <- ggplot(temp, aes(x = est_config, y = p_pass)) + 
    facet_nested(obs_error+ Fhist ~ M_sig + M_cor, 
      labeller = labeller(M_cor = label_parsed, M_sig = label_parsed, Fhist = label_parsed)) +
    coord_cartesian(ylim = c(0, 1)) + ylab("Probability of convergence") + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7)) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_M_om_p_conv_type_4
cairo_pdf(here("Project_0","manuscript", "R_M_om_p_convergence_type_4.pdf"), width = 24, height = 16)
R_M_om_p_conv_type_4
dev.off()

temp <- subset(all_res_mod, Type == 3)
R_M_om_p_conv_type_3 <- ggplot(temp, aes(x = est_config, y = p_pass)) + 
    facet_nested(obs_error+ Fhist ~ M_sig + M_cor, 
      labeller = labeller(M_cor = label_parsed, M_sig = label_parsed, Fhist = label_parsed)) +
    coord_cartesian(ylim = c(0, 1)) + ylab("Probability of convergence") + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7)) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_M_om_p_conv_type_3
cairo_pdf(here("Project_0","manuscript", "R_M_om_p_convergence_type_3.pdf"), width = 24, height = 16)
R_M_om_p_conv_type_3
dev.off()


df.oms = readRDS(here("Project_0", "inputs", "df.Sel.oms.RDS"))
df.ems = readRDS(here("Project_0", "inputs", "df.ems.RDS"))[c(5:20, 25:28),]
Sel_om_conv_res <- readRDS(here("Project_0", "results", "Sel_om_convergence_results.RDS"))

ind <- which(df.ems$M_est == FALSE & df.ems$SR_model == 2)
all_res_mod_Mknown_NoSR <- conv_res_plotting_fn(Sel_om_conv_res, df.oms, df.ems, em_ind = ind, om_type = "Sel")
ind <- which(df.ems$M_est == TRUE & df.ems$SR_model == 2)
all_res_mod_Mest_NoSR <- conv_res_plotting_fn(Sel_om_conv_res, df.oms, df.ems, em_ind = ind, om_type = "Sel")
ind <- which(df.ems$M_est == FALSE & df.ems$SR_model == 3)
all_res_mod_Mknown_SR <- conv_res_plotting_fn(Sel_om_conv_res, df.oms, df.ems, em_ind = ind, om_type = "Sel")
ind <- which(df.ems$M_est == TRUE & df.ems$SR_model == 3)
all_res_mod_Mest_SR <- conv_res_plotting_fn(Sel_om_conv_res, df.oms, df.ems, em_ind = ind, om_type = "Sel")
est_config <- c("M known, No SR", "M estimated, No SR", "M known, SR", "M estimated, SR")
all_res_mod <- rbind(
  cbind(all_res_mod_Mknown_NoSR, est_config = "M known, No SR"),
  cbind(all_res_mod_Mest_NoSR, est_config = "M estimated, No SR"),
  cbind(all_res_mod_Mknown_SR, est_config = "M known, SR"),
  cbind(all_res_mod_Mest_SR, est_config = "M estimated, SR")
)
all_res_mod$est_config <- factor(all_res_mod$est_config, levels = est_config)
EM_process_error <- c("R","R+S","R+M (iid)","R+M (AR1)","R+Sel (iid)","R+Sel (AR1)","R+q (iid)","R+q (AR1)")
all_res_mod$EM_process_error <- factor(all_res_mod$EM_process_error, levels = EM_process_error)

temp <- subset(all_res_mod, Type == 4)
R_Sel_om_p_conv_type_4 <- ggplot(temp, aes(x = est_config, y = p_pass)) + 
    facet_nested(obs_error+ Fhist ~ Sel_sig + Sel_cor, 
      labeller = labeller(Sel_cor = label_parsed, Sel_sig = label_parsed, Fhist = label_parsed)) +
    coord_cartesian(ylim = c(0, 1)) + ylab("Probability of convergence") + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7)) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_Sel_om_p_conv_type_4
cairo_pdf(here("Project_0","manuscript", "R_Sel_om_p_convergence_type_4.pdf"), width = 24, height = 16)
R_Sel_om_p_conv_type_4
dev.off()

temp <- subset(all_res_mod, Type == 3)
R_Sel_om_p_conv_type_3 <- ggplot(temp, aes(x = est_config, y = p_pass)) + 
    facet_nested(obs_error+ Fhist ~ Sel_sig + Sel_cor, 
      labeller = labeller(Sel_cor = label_parsed, Sel_sig = label_parsed, Fhist = label_parsed)) +
    coord_cartesian(ylim = c(0, 1)) + ylab("Probability of convergence") + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7)) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_Sel_om_p_conv_type_3
cairo_pdf(here("Project_0","manuscript", "R_Sel_om_p_convergence_type_3.pdf"), width = 24, height = 16)
R_Sel_om_p_conv_type_3
dev.off()

df.oms = readRDS(here("Project_0", "inputs", "df.q.oms.RDS"))
df.ems = readRDS(here("Project_0", "inputs", "df.ems.RDS"))[c(5:20, 29:32),]
q_om_conv_res <- readRDS(here("Project_0", "results", "q_om_convergence_results.RDS"))
ind <- which(df.ems$M_est == FALSE & df.ems$SR_model == 2)
all_res_mod_Mknown_NoSR <- conv_res_plotting_fn(q_om_conv_res, df.oms, df.ems, em_ind = ind, om_type = "q")
ind <- which(df.ems$M_est == TRUE & df.ems$SR_model == 2)
all_res_mod_Mest_NoSR <- conv_res_plotting_fn(q_om_conv_res, df.oms, df.ems, em_ind = ind, om_type = "q")
ind <- which(df.ems$M_est == FALSE & df.ems$SR_model == 3)
all_res_mod_Mknown_SR <- conv_res_plotting_fn(q_om_conv_res, df.oms, df.ems, em_ind = ind, om_type = "q")
ind <- which(df.ems$M_est == TRUE & df.ems$SR_model == 3)
all_res_mod_Mest_SR <- conv_res_plotting_fn(q_om_conv_res, df.oms, df.ems, em_ind = ind, om_type = "q")
est_config <- c("M known, No SR", "M estimated, No SR", "M known, SR", "M estimated, SR")
all_res_mod <- rbind(
  cbind(all_res_mod_Mknown_NoSR, est_config = "M known, No SR"),
  cbind(all_res_mod_Mest_NoSR, est_config = "M estimated, No SR"),
  cbind(all_res_mod_Mknown_SR, est_config = "M known, SR"),
  cbind(all_res_mod_Mest_SR, est_config = "M estimated, SR")
)
all_res_mod$est_config <- factor(all_res_mod$est_config, levels = est_config)
EM_process_error <- c("R","R+S","R+M (iid)","R+M (AR1)","R+Sel (iid)","R+Sel (AR1)","R+q (iid)","R+q (AR1)")
all_res_mod$EM_process_error <- factor(all_res_mod$EM_process_error, levels = EM_process_error)

temp <- subset(all_res_mod, Type == 4)
R_q_om_p_conv_type_4 <- ggplot(temp, aes(x = est_config, y = p_pass)) + 
    facet_nested(obs_error+ Fhist ~ q_sig + q_cor, 
      labeller = labeller(q_cor = label_parsed, q_sig = label_parsed, Fhist = label_parsed)) +
    coord_cartesian(ylim = c(0, 1)) + ylab("Probability of convergence") + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7)) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_q_om_p_conv_type_4
cairo_pdf(here("Project_0","manuscript", "R_q_om_p_convergence_type_4.pdf"), width = 24, height = 16)
R_q_om_p_conv_type_4
dev.off()

temp <- subset(all_res_mod, Type == 3)
R_q_om_p_conv_type_3 <- ggplot(temp, aes(x = est_config, y = p_pass)) + 
    facet_nested(obs_error+ Fhist ~ q_sig + q_cor, 
      labeller = labeller(q_cor = label_parsed, q_sig = label_parsed, Fhist = label_parsed)) +
    coord_cartesian(ylim = c(0, 1)) + ylab("Probability of convergence") + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7)) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_q_om_p_conv_type_3
cairo_pdf(here("Project_0","manuscript", "R_q_om_p_convergence_type_3.pdf"), width = 24, height = 16)
R_q_om_p_conv_type_3
dev.off()


theme_set(theme_bw())
theme_update(strip.placement = "outside", strip.background = element_rect(), title = element_text(size = rel(1.5)))
R_Sel_om_p_conv_type_4_alt <- R_Sel_om_p_conv_type_4 + 
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21, show.legend = TRUE) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7)) 
R_Sel_om_p_conv_type_4_alt

cairo_pdf(here("Project_0","manuscript", "type_4_convergence_plots.pdf"), width = 30*2/3, height = 20*2/3)
design <- c(area(1,1,1,1), area(2,1,2,1), area(1,2,1,2), area(2,2,2,2))
(R_S_om_p_conv_type_4 +  xlab("") + labs(title = "A: R, R+S OMs") + theme(axis.title.x = element_blank(), strip.text.y=element_blank(), axis.text.x=element_blank(), legend.position="none")) +# , strip.text=element_text(margin=margin(b = 5)), plot.margin = margin(b = 1, t = 0))) + 
(R_M_om_p_conv_type_4 + labs(title = "C: R+M OMs") + theme(legend.position="none", strip.text.y=element_blank())) + #axis.title.x = element_blank(), strip.text=element_text(margin=margin(b = 5)), axis.text.x=element_blank(), plot.margin = margin(b = 1, t = 0))) + 
(R_Sel_om_p_conv_type_4_alt + xlab("") + labs(title = "B: R+Sel OMs") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())) +# , strip.text=element_text(margin=margin(b = 5)), plot.margin = margin(b = 1, t = 0))) + 
(R_q_om_p_conv_type_4 + labs(title = "D: R+q OMs") + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position = "none")) + 
  plot_layout(design = design, axis_titles = "collect")
dev.off()

R_Sel_om_p_conv_type_3_alt <- R_Sel_om_p_conv_type_3 + 
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21, show.legend = TRUE) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7)) 

cairo_pdf(here("Project_0","manuscript", "type_3_convergence_plots.pdf"), width = 30*2/3, height = 20*2/3)
design <- c(area(1,1,1,1), area(2,1,2,1), area(1,2,1,2), area(2,2,2,2))
(R_S_om_p_conv_type_3 +  xlab("") + labs(title = "A: R, R+S OMs") + theme(axis.title.x = element_blank(), strip.text.y=element_blank(), axis.text.x=element_blank(), legend.position="none")) +# , strip.text=element_text(margin=margin(b = 5)), plot.margin = margin(b = 1, t = 0))) + 
(R_M_om_p_conv_type_3 + labs(title = "C: R+M OMs") + theme(legend.position="none", strip.text.y=element_blank())) + #axis.title.x = element_blank(), strip.text=element_text(margin=margin(b = 5)), axis.text.x=element_blank(), plot.margin = margin(b = 1, t = 0))) + 
(R_Sel_om_p_conv_type_3_alt + xlab("") + labs(title = "B: R+Sel OMs") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())) +# , strip.text=element_text(margin=margin(b = 5)), plot.margin = margin(b = 1, t = 0))) + 
(R_q_om_p_conv_type_3 + labs(title = "D: R+q OMs") + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position = "none")) + 
  plot_layout(design = design, axis_titles = "collect")
dev.off()


