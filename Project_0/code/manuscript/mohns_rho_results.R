library(here)
library(ggplot2)
library(ggh4x)
library(dplyr)
library(tidyr)
library(patchwork)
library(reshape2)

median_ci_fn <- function(x, alpha = 0.05){#, n.yrs = 1, n.sim = 100) {
  x <- x[which(!is.na(x))]
  n <- length(x)
  bnds <- qbinom(c(alpha/2,1-alpha/2), n, 0.5)/n # 95% CI bounds for median
  r <- quantile(x, probs = c(bnds[1], 0.5, bnds[2]))
  names(r) <- c("lo", "middle", "hi")
  r
}
make_mohns_rho_df <- function(om_type = "naa", res = all_naa_mohns_rho) {
  df.ems <- readRDS(here::here("Project_0","inputs", "df.ems.RDS"))
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
    em_ind <- c(5:20,29:32)
  }
  df.ems <- df.ems[em_ind,]

  res <- melt(res)
  names(res) <- c("rho", "em", "sim","om")
  res$Type = c("SSB","F","R")
  res <- cbind(df.ems[res$em,], res)

  df <- res %>% group_by(om, em, Type) %>%
    reframe(stats = median_ci_fn(rho)) %>% as.data.frame
  df$type <- c("lo", "middle", "hi")
  if(om_type == "naa") {
    df.oms <- readRDS(here::here("Project_0", "inputs", "df.oms.RDS"))
  } else {
    df.oms <- readRDS(here::here("Project_0", "inputs", paste0("df.", om_type, ".oms.RDS")))
  }
  df <- df %>% pivot_wider(names_from = type, values_from = stats) %>% as.data.frame
  df <- cbind(df, df.oms[df$om,])
  df <- cbind(df, df.ems[df$em,])
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
  EM_process_error <- c("R","R+S","R+M (iid)","R+M (AR1)","R+Sel (iid)","R+Sel (AR1)","R+q (iid)","R+q (AR1)")
  df$EM_process_error <- factor(df$EM_process_error, levels = EM_process_error)

  df$correct_EM_PE <- "No"
  if(om_type == "naa") {
    df$NAA_sig[which(is.na(df$NAA_sig))] <- 0
    df$correct_EM_PE[df$NAA_sig == 0 & df$EM_process_error == "R"] <- "Yes"
    df$correct_EM_PE[df$NAA_sig >  0 & df$EM_process_error == "R+S"] <- "Yes"
    df <- df %>% mutate(NAA_sig = recode(NAA_sig,
        "0" = "sigma['2+'] == 0",
        "0.25" = "sigma['2+'] == 0.25",
        "0.5" = "sigma['2+'] == 0.5")) %>%
      mutate(R_sig = recode(R_sig,
        "0.5" = "sigma[italic(R)] == 0.5",
        "1.5" = "sigma[italic(R)] == 1.5"))
  }

  if(om_type == "M") {
    df$correct_EM_PE[df$M_cor == 0 & df$EM_process_error == "R+M (iid)"] <- "Yes"
    df$correct_EM_PE[df$M_cor >  0 & df$EM_process_error == "R+M (AR1)"] <- "Yes"
    df <- df %>% mutate(M_sig = recode(M_sig,
        "0.1" = "sigma[italic(M)] == 0.1",
        "0.5" = "sigma[italic(M)] == 0.5")) %>% 
      mutate(M_cor = recode(M_cor,
        "0" = "rho[italic(M)] == 0",
        "0.9" = "rho[italic(M)] == 0.9"))
  }
  if(om_type == "Sel") {
    df$correct_EM_PE[df$Sel_cor == 0 & df$EM_process_error == "R+Sel (iid)"] <- "Yes"
    df$correct_EM_PE[df$Sel_cor >  0 & df$EM_process_error == "R+Sel (AR1)"] <- "Yes"
    df <- df %>% mutate(Sel_sig = recode(Sel_sig,
        "0.1" = "sigma[Sel] == 0.1",
        "0.5" = "sigma[Sel] == 0.5")) %>% 
      mutate(Sel_cor = recode(Sel_cor,
        "0" = "rho[Sel] == 0",
        "0.9" = "rho[Sel] == 0.9"))
  }
  if(om_type == "q") {
    df$correct_EM_PE[df$q_cor == 0 & df$EM_process_error == "R+q (iid)"] <- "Yes"
    df$correct_EM_PE[df$q_cor >  0 & df$EM_process_error == "R+q (AR1)"] <- "Yes"
    df <- df %>% mutate(q_sig = recode(q_sig,
        "0.1" = "sigma[italic(q)] == 0.1",
        "0.5" = "sigma[italic(q)] == 0.5")) %>%
      mutate(q_cor = recode(q_cor,
        "0" = "rho[italic(q)] == 0",
        "0.9" = "rho[italic(q)] == 0.9"))  
  }
  df <- df %>%
    mutate(obs_error = recode(obs_error,
      "L" = "Low observation error",
      "H" = "High observation error"))
  df <- df %>% mutate(Fhist = recode(Fhist,
        "H-MSY"  = "2.5*italic(F)[MSY] %->% italic(F)[MSY]",
        "MSY" = "italic(F)[MSY]"))
  df <- df %>% as.data.frame
  facs <- c("om", "Fhist","NAA_sig", "R_sig", "M_sig", "M_cor", "Sel_sig", "Sel_cor", "q_sig", "q_cor", "obs_error", "correct_EM_PE")
  fac_names <- names(df)[names(df) %in% facs]
  df[fac_names] <- lapply(df[fac_names], as.factor)
  est_config <- c("M known, No SR", "M estimated, No SR", "M known, SR", "M estimated, SR")
  df$est_config <- est_config[1]
  df$est_config[df$M_est & df$SR_model == 2] <- est_config[2]
  df$est_config[!df$M_est & df$SR_model == 3] <- est_config[3]
  df$est_config[df$M_est & df$SR_model == 3] <- est_config[4]
  df$est_config <- factor(df$est_config, levels = est_config)
  return(df)
}

all_naa_mohns_rho <-  readRDS(file = here("Project_0","results", "all_naa_mohns_rho_results.RDS"))
naa_mohns_rho_df <- make_mohns_rho_df(om_type = "naa", res = all_naa_mohns_rho)
all_Sel_mohns_rho <-  readRDS(file = here("Project_0","results", "all_Sel_mohns_rho_results.RDS"))
Sel_mohns_rho_df <- make_mohns_rho_df(om_type = "Sel", res = all_Sel_mohns_rho)
all_M_mohns_rho <-  readRDS(file = here("Project_0","results", "all_M_mohns_rho_results.RDS"))
M_mohns_rho_df <- make_mohns_rho_df(om_type = "M", res = all_M_mohns_rho)
all_q_mohns_rho <-  readRDS(file = here("Project_0","results", "all_q_mohns_rho_results.RDS"))
q_mohns_rho_df <- make_mohns_rho_df(om_type = "q", res = all_q_mohns_rho)

ylims <- c(-0.4,0.4)
Values <- c("SSB","F","R")
Ylabs <- list(expression("Mohn's"~rho(SSB)), expression("Mohn's"~rho(italic(F))), expression("Mohn's"~rho(italic(R))))
for(i in 1:3)
  theme_set(theme_bw())
  theme_update(strip.text = element_text(size = rel(2)), strip.placement = "outside", strip.background = element_rect(), #fill = "transparent"), 
        axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
        legend.title = element_text(size = rel(2)))

  R_S_plt <- ggplot(filter(naa_mohns_rho_df, Type == Values[i]), aes(x = est_config, y = middle)) + 
      facet_nested(obs_error + Fhist ~ R_sig + NAA_sig, labeller = labeller(Fhist = label_parsed, R_sig = label_parsed, NAA_sig = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
      coord_cartesian(ylim = ylims) + ylab(Ylabs[[i]]) + xlab("EM M and Stock Recruit Assumptions") +
      scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
      scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
      geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
      geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
      geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7)) + 
      scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
      geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
      labs(colour = "EM process error", fill = "EM process error") +
      guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

  R_M_plt <- ggplot(filter(M_mohns_rho_df, Type == Values[i]), aes(x = est_config, y = middle)) + 
      facet_nested(obs_error + Fhist ~ M_sig + M_cor, labeller = labeller(Fhist = label_parsed, M_sig = label_parsed, M_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
      coord_cartesian(ylim = ylims) + ylab(Ylabs[[i]]) + xlab("EM M and Stock Recruit Assumptions") +
      scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
      scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
      geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
      geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
      geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7)) + 
      scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
      geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
      labs(colour = "EM process error", fill = "EM process error") +
      guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

  R_Sel_plt <- ggplot(filter(Sel_mohns_rho_df, Type == Values[i]), aes(x = est_config, y = middle)) + 
      facet_nested(obs_error + Fhist ~ Sel_sig + Sel_cor, labeller = labeller(Fhist = label_parsed, Sel_sig = label_parsed, Sel_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
      coord_cartesian(ylim = ylims) + ylab(Ylabs[[i]]) + xlab("EM M and Stock Recruit Assumptions") +
      scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
      scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
      geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
      geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
      geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7)) + 
      scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
      geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
      labs(colour = "EM process error", fill = "EM process error") +
      guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

  R_q_plt <- ggplot(filter(q_mohns_rho_df, Type == Values[i]), aes(x = est_config, y = middle)) + 
      facet_nested(obs_error + Fhist ~ q_sig + q_cor, labeller = labeller(Fhist = label_parsed, q_sig = label_parsed, q_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
      coord_cartesian(ylim = ylims) + ylab(Ylabs[[i]]) + xlab("EM M and Stock Recruit Assumptions") +
      scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
      scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
      geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
      geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
      geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7)) + 
      scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
      geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
      labs(colour = "EM process error", fill = "EM process error") +
      guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

  R_Sel_plt_alt <- R_Sel_plt + 
      geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21, show.legend = TRUE) + 
      geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7)) 
  
  theme_set(theme_bw())
  theme_update(strip.placement = "outside", strip.background = element_rect(), title = element_text(size = rel(1.5)))
  cairo_pdf(here("Project_0","manuscript", paste0("mohns_rho_", Values[i],"_plots.pdf")), width = 30*2/3, height = 20*2/3)
  design <- c(area(1,1,1,1), area(2,1,2,1), area(1,2,1,2), area(2,2,2,2))
  (R_S_plt +  xlab("") + labs(title = "A: R, R+S OMs") + theme(axis.title.x = element_blank(), strip.text.y=element_blank(), axis.text.x=element_blank(), legend.position="none")) +
  (R_M_plt + labs(title = "C: R+M OMs") + theme(legend.position="none", strip.text.y=element_blank())) + 
  (R_Sel_plt_alt + xlab("") + labs(title = "B: R+Sel OMs") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())) +
  (R_q_plt + labs(title = "D: R+q OMs") + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position = "none")) + 
    plot_layout(design = design, axis_titles = "collect")
  dev.off()
}

# Mohns rho Fbar
theme_set(theme_bw())
theme_update(strip.text = element_text(size = rel(2)), strip.placement = "outside", strip.background = element_rect(), #fill = "transparent"), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)))

R_S_plt <- ggplot(filter(naa_mohns_rho_df, Type == "F"), aes(x = est_config, y = middle)) + 
    coord_cartesian(ylim = ylims) + ylab(expression("Mohn's"~rho(italic(F)))) + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    facet_nested(obs_error + Fhist ~ R_sig + NAA_sig, labeller = labeller(Fhist = label_parsed, R_sig = label_parsed, NAA_sig = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    geom_point(position = position_dodge(0.7), mapping = aes(colour = EM_process_error, fill = EM_process_error, size = correct_EM_PE, shape = correct_EM_PE)) + 
    scale_shape_manual(values = c('No'=21, 'Yes'=23)) + 
    scale_size_manual(values = c('No'=2, 'Yes'=4)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error", shape = "EM process error\ncorrect?", size = "EM process error\ncorrect?") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), fill = guide_legend(order = 1)) +
    coord_cartesian(ylim = ylims) + ylab(expression("Mohn's"~rho(italic(F)))) + xlab("EM Configuration")
R_S_plt

R_M_plt <- ggplot(filter(M_mohns_rho_df, Type == "F"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ M_sig + M_cor, labeller = labeller(Fhist = label_parsed, M_sig = label_parsed, M_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    coord_cartesian(ylim = ylims) + ylab(expression("Mohn's"~rho(italic(F)))) + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), mapping = aes(colour = EM_process_error, fill = EM_process_error, size = correct_EM_PE, shape = correct_EM_PE)) + 
    scale_shape_manual(values = c('No'=21, 'Yes'=23)) + 
    scale_size_manual(values = c('No'=2, 'Yes'=4)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error", shape = "EM process error\ncorrect?", size = "EM process error\ncorrect?") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), fill = guide_legend(order = 1)) +
    coord_cartesian(ylim = ylims) + ylab(expression("Mohn's"~rho(italic(F)))) + xlab("EM Configuration")
R_M_plt


R_Sel_plt <- ggplot(filter(Sel_mohns_rho_df, Type == "F"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ Sel_sig + Sel_cor, labeller = labeller(Fhist = label_parsed, Sel_sig = label_parsed, Sel_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    coord_cartesian(ylim = ylims) + ylab(expression("Mohn's"~rho(italic(F)))) + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), mapping = aes(colour = EM_process_error, fill = EM_process_error, size = correct_EM_PE, shape = correct_EM_PE)) + 
    scale_shape_manual(values = c('No'=21, 'Yes'=23)) + 
    scale_size_manual(values = c('No'=2, 'Yes'=4)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error", shape = "EM process error\ncorrect?", size = "EM process error\ncorrect?") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), fill = guide_legend(order = 1)) +
    coord_cartesian(ylim = ylims) + ylab(expression("Mohn's"~rho(italic(F)))) + xlab("EM Configuration")
R_Sel_plt

R_q_plt <- ggplot(filter(q_mohns_rho_df, Type == "F"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ q_sig + q_cor, labeller = labeller(Fhist = label_parsed, q_sig = label_parsed, q_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    coord_cartesian(ylim = ylims) + ylab(expression("Mohn's"~rho(italic(F)))) + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), mapping = aes(colour = EM_process_error, fill = EM_process_error, size = correct_EM_PE, shape = correct_EM_PE)) + 
    scale_shape_manual(values = c('No'=21, 'Yes'=23)) + 
    scale_size_manual(values = c('No'=2, 'Yes'=4)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error", shape = "EM process error\ncorrect?", size = "EM process error\ncorrect?") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), fill = guide_legend(order = 1)) +
    coord_cartesian(ylim = ylims) + ylab(expression("Mohn's"~rho(italic(F)))) + xlab("EM Configuration")
R_q_plt

R_Sel_plt <- R_Sel_plt + 
  geom_point(position = position_dodge(0.7), mapping = aes(colour = EM_process_error, fill = EM_process_error, size = correct_EM_PE, shape = correct_EM_PE), show.legend = TRUE) + 
  geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) + #, show.legend = TRUE) +
  guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), fill = guide_legend(order = 1))


theme_set(theme_bw())
theme_update(strip.placement = "outside", strip.background = element_rect(), title = element_text(size = rel(1.5)))
cairo_pdf(here("Project_0","manuscript", "mohns_rho_F_plots.pdf"), width = 30*2/3, height = 20*2/3)
design <- c(area(1,1,1,1), area(2,1,2,1), area(1,2,1,2), area(2,2,2,2))
(R_S_plt +  xlab("") + labs(title = "A: R, R+S OMs") + theme(axis.title.x = element_blank(), strip.text.y=element_blank(), axis.text.x=element_blank(), legend.position="none")) +
(R_M_plt + labs(title = "C: R+M OMs") + theme(legend.position="none", strip.text.y=element_blank())) + 
(R_Sel_plt + xlab("") + labs(title = "B: R+Sel OMs") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())) +
(R_q_plt + labs(title = "D: R+q OMs") + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position = "none")) + 
  plot_layout(design = design, axis_titles = "collect")
dev.off()


# Mohns rho Recruit
theme_set(theme_bw())
theme_update(strip.text = element_text(size = rel(2)), strip.placement = "outside", strip.background = element_rect(), #fill = "transparent"), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)))

R_S_plt <- ggplot(filter(naa_mohns_rho_df, Type == "R"), aes(x = est_config, y = middle)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    facet_nested(obs_error + Fhist ~ R_sig + NAA_sig, labeller = labeller(Fhist = label_parsed, R_sig = label_parsed, NAA_sig = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    geom_point(position = position_dodge(0.7), mapping = aes(colour = EM_process_error, fill = EM_process_error, size = correct_EM_PE, shape = correct_EM_PE)) + 
    scale_shape_manual(values = c('No'=21, 'Yes'=23)) + 
    scale_size_manual(values = c('No'=2, 'Yes'=4)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error", shape = "EM process error\ncorrect?", size = "EM process error\ncorrect?") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), fill = guide_legend(order = 1)) +
    coord_cartesian(ylim = ylims) + ylab(expression("Mohn's"~rho(italic(R)))) + xlab("EM Configuration")
R_S_plt

R_M_plt <- ggplot(filter(M_mohns_rho_df, Type == "R"), aes(x = est_config, y = middle)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    facet_nested(obs_error + Fhist ~ M_sig + M_cor, labeller = labeller(Fhist = label_parsed, M_sig = label_parsed, M_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    geom_point(position = position_dodge(0.7), mapping = aes(colour = EM_process_error, fill = EM_process_error, size = correct_EM_PE, shape = correct_EM_PE)) + 
    scale_shape_manual(values = c('No'=21, 'Yes'=23)) + 
    scale_size_manual(values = c('No'=2, 'Yes'=4)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error", shape = "EM process error\ncorrect?", size = "EM process error\ncorrect?") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), fill = guide_legend(order = 1)) +
    coord_cartesian(ylim = ylims) + ylab(expression("Mohn's"~rho(italic(R)))) + xlab("EM Configuration")
R_M_plt


R_Sel_plt <- ggplot(filter(Sel_mohns_rho_df, Type == "R"), aes(x = est_config, y = middle)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    facet_nested(obs_error + Fhist ~ Sel_sig + Sel_cor, labeller = labeller(Fhist = label_parsed, Sel_sig = label_parsed, Sel_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    geom_point(position = position_dodge(0.7), mapping = aes(colour = EM_process_error, fill = EM_process_error, size = correct_EM_PE, shape = correct_EM_PE)) + 
    scale_shape_manual(values = c('No'=21, 'Yes'=23)) + 
    scale_size_manual(values = c('No'=2, 'Yes'=4)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error", shape = "EM process error\ncorrect?", size = "EM process error\ncorrect?") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), fill = guide_legend(order = 1)) +
    coord_cartesian(ylim = ylims) + ylab(expression("Mohn's"~rho(italic(R)))) + xlab("EM Configuration")
R_Sel_plt

R_q_plt <- ggplot(filter(q_mohns_rho_df, Type == "R"), aes(x = est_config, y = middle)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    facet_nested(obs_error + Fhist ~ q_sig + q_cor, labeller = labeller(Fhist = label_parsed, q_sig = label_parsed, q_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    geom_point(position = position_dodge(0.7), mapping = aes(colour = EM_process_error, fill = EM_process_error, size = correct_EM_PE, shape = correct_EM_PE)) + 
    scale_shape_manual(values = c('No'=21, 'Yes'=23)) + 
    scale_size_manual(values = c('No'=2, 'Yes'=4)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error", shape = "EM process error\ncorrect?", size = "EM process error\ncorrect?") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), fill = guide_legend(order = 1)) +
    coord_cartesian(ylim = ylims) + ylab(expression("Mohn's"~rho(italic(R)))) + xlab("EM Configuration")
R_q_plt

R_Sel_plt <- R_Sel_plt + 
  geom_point(position = position_dodge(0.7), mapping = aes(colour = EM_process_error, fill = EM_process_error, size = correct_EM_PE, shape = correct_EM_PE), show.legend = TRUE) + 
  geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) + #, show.legend = TRUE) +
  guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), fill = guide_legend(order = 1))


theme_set(theme_bw())
theme_update(strip.placement = "outside", strip.background = element_rect(), title = element_text(size = rel(1.5)))
cairo_pdf(here("Project_0","manuscript", "mohns_rho_R_plots.pdf"), width = 30*2/3, height = 20*2/3)
design <- c(area(1,1,1,1), area(2,1,2,1), area(1,2,1,2), area(2,2,2,2))
(R_S_plt +  xlab("") + labs(title = "A: R, R+S OMs") + theme(axis.title.x = element_blank(), strip.text.y=element_blank(), axis.text.x=element_blank(), legend.position="none")) +
(R_M_plt + labs(title = "C: R+M OMs") + theme(legend.position="none", strip.text.y=element_blank())) + 
(R_Sel_plt + xlab("") + labs(title = "B: R+Sel OMs") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())) +
(R_q_plt + labs(title = "D: R+q OMs") + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position = "none")) + 
  plot_layout(design = design, axis_titles = "collect")
dev.off()
