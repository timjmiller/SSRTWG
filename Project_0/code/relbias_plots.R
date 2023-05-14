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

relR_fn <- function(df.oms, om_type = "naa"){
  relR = lapply(1:NROW(df.oms), function(y){
    res = lapply(1:100, function(x){
      print(paste0("om_", y, ", sim_",x))
      sim = readRDS(file.path(here::here(),"Project_0", "results", paste0(om_type,"_om"), paste0("om_", y), paste0("sim_",x,".RDS")))
      relR = sapply(sim,function(y) {
        out <- rep(NA, 40)
        if(length(y)) if(length(y$fit)) out = y$fit$rep$NAA[,1]/y$truth$NAA[,1]
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

  # df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.oms.RDS"))
  # all_naa_relR <- relR_fn(df.oms,"naa")
  # saveRDS(all_naa_relR, file = here("Project_0","results", "all_naa_relR_results.RDS"))
  # df.M.oms = readRDS(file.path(here(),"Project_0","inputs", "df.M.oms.RDS"))
  # all_M_relR <- relR_fn(df.M.oms,"M")
  # saveRDS(all_M_relR, file = here("Project_0","results", "all_M_relR_results.RDS"))
  # df.Sel.oms = readRDS(file.path(here(),"Project_0","inputs", "df.Sel.oms.RDS"))
  # all_Sel_relR <- relR_fn(df.Sel.oms,"Sel")
  # saveRDS(all_Sel_relR, file = here("Project_0","results", "all_Sel_relR_results.RDS"))
  # df.q.oms = readRDS(file.path(here(),"Project_0","inputs", "df.q.oms.RDS"))
  # all_q_relR <- relR_fn(df.q.oms,"q")
  # saveRDS(all_q_relR, file = here("Project_0","results", "all_q_relR_results.RDS"))

  # df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.oms.RDS"))
  # all_naa_relF <- relF_fn(df.oms,"naa")
  # saveRDS(all_naa_relF, file = here("Project_0","results", "all_naa_relF_results.RDS"))
  # df.M.oms = readRDS(file.path(here(),"Project_0","inputs", "df.M.oms.RDS"))
  # all_M_relF <- relF_fn(df.M.oms,"M")
  # saveRDS(all_M_relF, file = here("Project_0","results", "all_M_relF_results.RDS"))
  # df.Sel.oms = readRDS(file.path(here(),"Project_0","inputs", "df.Sel.oms.RDS"))
  # all_Sel_relF <- relF_fn(df.Sel.oms,"Sel")
  # saveRDS(all_Sel_relF, file = here("Project_0","results", "all_Sel_relF_results.RDS"))
  # df.q.oms = readRDS(file.path(here(),"Project_0","inputs", "df.q.oms.RDS"))
  # all_q_relF <- relF_fn(df.q.oms,"q")
  # saveRDS(all_q_relF, file = here("Project_0","results", "all_q_relF_results.RDS"))
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
library(dplyr)
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
  plot.df <- res %>% pivot_wider(names_from = column, values_from = value) %>% as.data.frame
  print(dim(plot.df))
  if(is_SE) plot.df <- filter(plot.df, !is.na(cv))
  print(dim(plot.df))
  plot.df$relerror = plot.df$relerror - 1

  plot.df <- plot.df %>% group_by(om, em, year) %>%
    reframe(stats = custom_boxplot_stat(relerror)) %>% as.data.frame
  plot.df$type <- c("ymin", "lower", "middle", "upper", "ymax")
  library(tidyr)
  df.ems <- get_df_ems(type)
  if(type == "naa") {
    df.oms <- readRDS(here::here("Project_0", "inputs", "df.oms.RDS"))
  } else {
    df.oms <- readRDS(here::here("Project_0", "inputs", paste0("df.", type, ".oms.RDS")))
  }
  print(head(df.oms))
  plot.df <- plot.df %>% pivot_wider(names_from = type, values_from = stats) %>% as.data.frame
  plot.df$em_pe = df.ems$pe#[plot.df$em]
  plot.df <- cbind(plot.df, df.ems[plot.df$em,])
  plot.df <- cbind(plot.df, df.oms[plot.df$om,])
  if(type == "naa") plot.df$NAA_sig[which(is.na(plot.df$NAA_sig))] <- 0
  # facs <- c("em_pe", "Fhist", "R_sig","NAA_sig", "obs_error")
  # plot.df[facs] <- lapply(plot.df[facs], factor)
  return(plot.df)
}

plot.df <- make_plot_df(is_SE = TRUE)

temp <- filter(plot.df, em_pe == "R+M" & !M_est & meanR == "Mean R" )
unique(temp$em)

temp <- filter(plot.df, !M_est & meanR == "Mean R" & om == 1 & year == 1)

plt <- ggplot(filter(plot.df, !M_est & meanR == "Mean R" )) + scale_fill_viridis_d() + 
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    facet_grid(R_sig + NAA_sig ~ Fhist + obs_error, labeller = label_both) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    geom_line(aes(x = year, y = middle, colour = pe), position = position_dodge(0.1), linewidth = 1) + #geom_point(position = position_dodge(0.1), size = 1) + 
    geom_ribbon(aes(x = year, ymin=ymin,ymax=ymax, fill = pe), alpha=0.3, position = position_dodge(0.1)) +
    theme_bw() + coord_cartesian(ylim = c(-0.4, 0.4)) + ylab("Relative Bias of SSB") + xlab("Year") +
    ggtitle("EM: M = 0.2 and no SRR") + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM Process Error", fill = "EM Process Error")
#    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .05, position = position_dodge(0.1)) + 
plt
ggsave(here("Project_0", "paper", "naa_om_R_MF_relbias_ssb.png"), plt, width = 20, height = 12, units = "in")

plt <- ggplot(filter(plot.df, !M_est & meanR == "B-H" )) + scale_fill_viridis_d() + 
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    facet_grid(R_sig + NAA_sig ~ Fhist + obs_error, labeller = label_both) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    geom_line(aes(x = year, y = middle, colour = pe), position = position_dodge(0.1), linewidth = 1) + #geom_point(position = position_dodge(0.1), size = 1) + 
    geom_ribbon(aes(x = year, ymin=ymin,ymax=ymax, fill = pe), alpha=0.3, position = position_dodge(0.1)) +
    theme_bw() + coord_cartesian(ylim = c(-0.4, 0.4)) + ylab("Relative Bias of SSB") + xlab("Year") +
    ggtitle("EM: M = 0.2 and BH assumed") + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM Process Error", fill = "EM Process Error")
#    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .05, position = position_dodge(0.1)) + 
plt
ggsave(here("Project_0", "paper", "naa_om_SR_MF_relbias_ssb.png"), plt, width = 20, height = 12, units = "in")

plt <- ggplot(filter(plot.df, M_est & meanR == "Mean R" )) + scale_fill_viridis_d() + 
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    facet_grid(R_sig + NAA_sig ~ Fhist + obs_error, labeller = label_both) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    geom_line(aes(x = year, y = middle, colour = pe), position = position_dodge(0.1), linewidth = 1) + #geom_point(position = position_dodge(0.1), size = 1) + 
    geom_ribbon(aes(x = year, ymin=ymin,ymax=ymax, fill = pe), alpha=0.3, position = position_dodge(0.1)) +
    theme_bw() + coord_cartesian(ylim = c(-0.4, 0.4)) + ylab("Relative Bias of SSB") + xlab("Year") +
    ggtitle("EM: M estimated and no SRR") + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM Process Error", fill = "EM Process Error")
#    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .05, position = position_dodge(0.1)) + 
plt
ggsave(here("Project_0", "paper", "naa_om_R_ME_relbias_ssb.png"), plt, width = 20, height = 12, units = "in")

plt <- ggplot(filter(plot.df, M_est & meanR == "B-H" )) + scale_fill_viridis_d() + 
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    facet_grid(R_sig + NAA_sig ~ Fhist + obs_error, labeller = label_both) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    geom_line(aes(x = year, y = middle, colour = pe), position = position_dodge(0.1), linewidth = 1) + #geom_point(position = position_dodge(0.1), size = 1) + 
    geom_ribbon(aes(x = year, ymin=ymin,ymax=ymax, fill = pe), alpha=0.3, position = position_dodge(0.1)) +
    theme_bw() + coord_cartesian(ylim = c(-0.4, 0.4)) + ylab("Relative Bias of SSB") + xlab("Year") +
    ggtitle("EM: M estimated and BH assumed") + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM Process Error", fill = "EM Process Error")
#    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .05, position = position_dodge(0.1)) + 
plt
ggsave(here("Project_0", "paper", "naa_om_SR_ME_relbias_ssb.png"), plt, width = 20, height = 12, units = "in")

all_M_relssb <- readRDS(file = here("Project_0","results", "all_M_relssb_results.RDS"))
plot.df <- make_plot_df("M", all_M_relssb, is_SE = TRUE)

plt <- ggplot(filter(plot.df, !M_est & meanR == "Mean R" )) + scale_fill_viridis_d() + 
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    facet_grid(M_sig + M_cor ~ Fhist + obs_error, labeller = label_both) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    geom_line(aes(x = year, y = middle, colour = pe), position = position_dodge(0.1), linewidth = 1) + #geom_point(position = position_dodge(0.1), size = 1) + 
    geom_ribbon(aes(x = year, ymin=ymin,ymax=ymax, fill = pe), alpha=0.3, position = position_dodge(0.1)) +
    theme_bw() + coord_cartesian(ylim = c(-0.4, 0.4)) + ylab("Relative Bias of SSB") + xlab("Year") +
    ggtitle("EM: M = 0.2 and no SRR") + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM Process Error", fill = "EM Process Error")
#    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .05, position = position_dodge(0.1)) + 
plt
ggsave(here("Project_0", "paper", "M_om_R_MF_relbias_ssb.png"), plt, width = 20, height = 12, units = "in")

all_Sel_relssb <- readRDS(file = here("Project_0","results", "all_Sel_relssb_results.RDS"))
plot.df <- make_plot_df("Sel", all_Sel_relssb, is_SE = TRUE)

plt <- ggplot(filter(plot.df, !M_est & meanR == "Mean R" )) + scale_fill_viridis_d() + 
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    facet_grid(Sel_sig + Sel_cor ~ Fhist + obs_error, labeller = label_both) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    geom_line(aes(x = year, y = middle, colour = pe), position = position_dodge(0.1), linewidth = 1) + #geom_point(position = position_dodge(0.1), size = 1) + 
    geom_ribbon(aes(x = year, ymin=ymin,ymax=ymax, fill = pe), alpha=0.3, position = position_dodge(0.1)) +
    theme_bw() + coord_cartesian(ylim = c(-0.4, 0.4)) + ylab("Relative Bias of SSB") + xlab("Year") +
    ggtitle("EM: M = 0.2 and no SRR") + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM Process Error", fill = "EM Process Error")
#    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .05, position = position_dodge(0.1)) + 
plt
ggsave(here("Project_0", "paper", "Sel_om_R_MF_relbias_ssb.png"), plt, width = 20, height = 12, units = "in")

all_q_relssb <- readRDS(file = here("Project_0","results", "all_q_relssb_results.RDS"))
plot.df <- make_plot_df("q", all_q_relssb, is_SE = TRUE)

plt <- ggplot(filter(plot.df, !M_est & meanR == "Mean R" )) + scale_fill_viridis_d() + 
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    facet_grid(q_sig + q_cor ~ Fhist + obs_error, labeller = label_both) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    geom_line(aes(x = year, y = middle, colour = pe), position = position_dodge(0.1), linewidth = 1) + #geom_point(position = position_dodge(0.1), size = 1) + 
    geom_ribbon(aes(x = year, ymin=ymin,ymax=ymax, fill = pe), alpha=0.3, position = position_dodge(0.1)) +
    theme_bw() + coord_cartesian(ylim = c(-0.4, 0.4)) + ylab("Relative Bias of SSB") + xlab("Year") +
    ggtitle("EM: M = 0.2 and no SRR") + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM Process Error", fill = "EM Process Error")
#    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .05, position = position_dodge(0.1)) + 
plt
ggsave(here("Project_0", "paper", "q_om_R_MF_relbias_ssb.png"), plt, width = 20, height = 12, units = "in")
