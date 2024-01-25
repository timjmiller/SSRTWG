library(here)
library(ggplot2)
library(dplyr)
df.ems = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.ems.RDS"))
df.oms = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.oms.RDS"))
F40_rel_error <- readRDS(file.path(here::here(),"Ecov_study","mortality", "results", "F40_rel_error_results.RDS"))
SSB40_rel_error <- readRDS(file.path(here::here(),"Ecov_study","mortality", "results", "SSB40_rel_error_results.RDS"))
om_inputs = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "om_inputs.RDS"))

custom_boxplot_stat <- function(x){#, n.yrs = 1, n.sim = 100) {
  x <- x[which(!is.na(x))]
  n <- length(x)
  bnds95 <- qbinom(c(0.025,0.975), n, 0.5)/n # 95% CI bounds for median (sims x years)
  bnds80 <- qbinom(c(0.1,0.9), n, 0.5)/n # 80% CI bounds for median (sims x years)
  r <- quantile(x, probs = c(bnds95[1], bnds80[1], 0.5, bnds80[2], bnds95[2]))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

ind <- which(df.oms$NAA_M_re == "rec+1" & df.oms$Fhist=="MSY" & df.oms$Ecov_effect == 0.5 & df.oms$obs_error=="L")

temp <- readRDS(file.path(here::here(),"Ecov_study","mortality", "results", "aggregated_results", paste0("om_", 122, "_em_", 10, ".RDS")))

temp <- readRDS(file.path(here::here(),"Ecov_study","mortality", "results", "om122", paste0("sim", 1, "_em", 10, ".RDS")))

#all EM PE assumptions
plot_df_fn <- function(df.ems, df.oms, M_est = FALSE, Ecov_est = FALSE, error_res = F40_rel_error) {
  for(h in 1:3) {
    re_mod <- c("rec", "rec+1", "rec+M")[h]
    #EM:  M fixed, mean_M estimated, OM and EM RE assumption match
    em_ind <- which(df.ems$Ecov_est==Ecov_est & df.ems$M_est == M_est)
    #df.ems. <- df.ems[,]
    om_ind <- which(df.oms$NAA_M_re == re_mod)
    for(i in om_ind){
      out <- matrix(NA,length(em_ind), 7)
      for(j in em_ind){
        x <- error_res[[i]][[j]]
        x <- x[which(!is.na(x[,2])),1] #remove fits without standard error estimates
        out[which(em_ind==j),] <- c(median(x, na.rm = TRUE), sd(x, na.rm=T)/sqrt(sum(!is.na(x))), custom_boxplot_stat(x))#quantile(x, probs = c(0.025,0.975)))
      }
      colnames(out) <- c("bias_est", "bias_se", "ymin", "lower", "middle", "upper", "ymax")
      if(i == om_ind[1]) res <- cbind.data.frame(df.oms[rep(i,length(em_ind)),], df.ems[em_ind,], out)
      else res <- rbind.data.frame(res,cbind.data.frame(df.oms[rep(i,length(em_ind)),], df.ems[em_ind,], out))
    }
    print(head(res))
    if(h == 1) {
      all_res <- res
    } else {
      all_res <- rbind.data.frame(all_res, res)
    }
  }
  facs <- c("Ecov_obs_sig", "Fhist", "Ecov_re_sig","Ecov_re_cor", "obs_error", "NAA_M_re", "re_config")
  all_res[facs] <- lapply(all_res[facs], factor)
  all_res_mod <- all_res %>%
    mutate(Ecov_obs_sig = recode(Ecov_obs_sig,
      "0.1" = "sigma[ecov] == 0.1",
      "0.5" = "sigma[ecov] == 0.5"
    ))
  all_res_mod <- all_res_mod %>%
    mutate(obs_error = recode(obs_error,
      "L" = "Low obs error (indices, age comp)",
      "H" = "High obs error (indices, age comp)"
    ))
  all_res_mod <- all_res_mod %>%
    mutate(NAA_M_re = recode(NAA_M_re,
      "rec" = "OM: R",
      "rec+1" = "OM: R+S",
      "rec+M" = "OM: R+M"
    ))
  all_res_mod <- all_res_mod %>%
    mutate(re_config = recode(re_config,
      "rec" = "R",
      "rec+1" = "R+S",
      "rec+M" = "R+M"
    ))
  all_res_mod <- all_res_mod %>% mutate(Fhist = recode(Fhist,
      "H-MSY" = "2.5*F[MSY] %->% F[MSY]",
      "MSY" = "F[MSY]"))
  all_res_mod <- all_res_mod %>%
    mutate(Ecov_re_sig = recode(Ecov_re_sig,
      "0.1" = "sigma[Ecov] == 0.1",
      "0.5" = "sigma[Ecov] == 0.5"
    ))
  all_res_mod <- all_res_mod %>%
    mutate(Ecov_re_cor = recode(Ecov_re_cor,
      "0" = "rho[Ecov] == 0",
      "0.5" = "rho[Ecov] == 0.5"
    ))
  return(all_res_mod)
}

all_res <- plot_df_fn(df.ems, df.oms, Ecov_est = FALSE, M_est = FALSE, error_res = F40_rel_error)
plt <- ggplot(all_res, aes(x = Ecov_effect, y = bias_est, colour = re_config)) + scale_colour_viridis_d() + 
    ggtitle(bquote(EM:~beta[M] == log(0.2)*","~beta[Ecov]==0)) +
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor  ~ NAA_M_re + Fhist + obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
    theme_bw() + coord_cartesian(ylim = c(-0.5, 0.5)) + ylab(bquote(Median~Relative~Error~of~hat(F)[40])) + xlab(expression(beta[Ecov])) +
    theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM process error") +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .05, position = position_dodge(0.1))
plt
ggsave(here("Ecov_study","mortality", "paper", "F40_rel_error_M_fixed_Ecov_beta_fixed.png"), plt, width = 20, height = 12, units = "in")
remove(all_res)

all_res <- plot_df_fn(df.ems, df.oms, Ecov_est = FALSE, M_est = TRUE, error_res = F40_rel_error)
plt <- ggplot(all_res, aes(x = Ecov_effect, y = bias_est, colour = re_config)) + scale_colour_viridis_d() + 
    ggtitle(bquote(EM:~beta[M]~Estimated*","~beta[Ecov]==0)) +
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor  ~ NAA_M_re + Fhist + obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
    theme_bw() + coord_cartesian(ylim = c(-0.5, 0.5)) + ylab(bquote(Median~Relative~Error~of~hat(F)[40])) + xlab(expression(beta[Ecov])) +
    theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM process error") +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .05, position = position_dodge(0.1))
plt
ggsave(here("Ecov_study","mortality", "paper", "F40_rel_error_M_estimated_Ecov_beta_fixed.png"), plt, width = 20, height = 12, units = "in")
remove(all_res)

all_res <- plot_df_fn(df.ems, df.oms, Ecov_est = TRUE, M_est = FALSE, error_res = F40_rel_error)
plt <- ggplot(all_res, aes(x = Ecov_effect, y = bias_est, colour = re_config)) + scale_colour_viridis_d() + 
    ggtitle(bquote(EM:~beta[M] == log(0.2)*","~beta[Ecov]~Estimated)) +
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor  ~ NAA_M_re + Fhist + obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
    theme_bw() + coord_cartesian(ylim = c(-0.5, 0.5)) + ylab(bquote(Median~Relative~Error~of~hat(F)[40])) + xlab(expression(beta[Ecov])) +
    theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM process error") +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .05, position = position_dodge(0.1))
plt
ggsave(here("Ecov_study","mortality", "paper", "F40_rel_error_M_fixed_Ecov_beta_estimated.png"), plt, width = 20, height = 12, units = "in")
remove(all_res)

all_res <- plot_df_fn(df.ems, df.oms, Ecov_est = TRUE, M_est = TRUE, error_res = F40_rel_error)
plt <- ggplot(all_res, aes(x = Ecov_effect, y = bias_est, colour = re_config)) + scale_colour_viridis_d() + 
    ggtitle(bquote(EM:~beta[M]~Estimated*","~beta[Ecov]~Estimated)) +
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor  ~ NAA_M_re + Fhist + obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
    theme_bw() + coord_cartesian(ylim = c(-0.5, 0.5)) + ylab(bquote(Median~Relative~Error~of~hat(F)[40])) + xlab(expression(beta[Ecov])) +
    theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM process error") +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .05, position = position_dodge(0.1))
plt
ggsave(here("Ecov_study","mortality", "paper", "F40_rel_error_M_estimated_Ecov_beta_estimated.png"), plt, width = 20, height = 12, units = "in")
remove(all_res)

#######################
SSB40

all_res <- plot_df_fn(df.ems, df.oms, Ecov_est = FALSE, M_est = FALSE, error_res = SSB40_rel_error)
plt <- ggplot(all_res, aes(x = Ecov_effect, y = bias_est, colour = re_config)) + scale_colour_viridis_d() + 
    ggtitle(bquote(EM:~beta[M] == log(0.2)*","~beta[Ecov]==0)) +
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor  ~ NAA_M_re + Fhist + obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
    theme_bw() + coord_cartesian(ylim = c(-0.5, 0.5)) + ylab(bquote(Median~Relative~Error~of~widehat(SSB)[40])) + xlab(expression(beta[Ecov])) +
    theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM process error") +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .05, position = position_dodge(0.1))
plt
ggsave(here("Ecov_study","mortality", "paper", "F40_rel_error_M_fixed_Ecov_beta_fixed.png"), plt, width = 20, height = 12, units = "in")
remove(all_res)

all_res <- plot_df_fn(df.ems, df.oms, Ecov_est = FALSE, M_est = TRUE, error_res = F40_rel_error)
plt <- ggplot(all_res, aes(x = Ecov_effect, y = bias_est, colour = re_config)) + scale_colour_viridis_d() + 
    ggtitle(bquote(EM:~beta[M]~Estimated*","~beta[Ecov]==0)) +
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor  ~ NAA_M_re + Fhist + obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
    theme_bw() + coord_cartesian(ylim = c(-0.5, 0.5)) + ylab(bquote(Median~Relative~Error~of~hat(F)[40])) + xlab(expression(beta[Ecov])) +
    theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM process error") +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .05, position = position_dodge(0.1))
plt
ggsave(here("Ecov_study","mortality", "paper", "F40_rel_error_M_estimated_Ecov_beta_fixed.png"), plt, width = 20, height = 12, units = "in")
remove(all_res)

all_res <- plot_df_fn(df.ems, df.oms, Ecov_est = TRUE, M_est = FALSE, error_res = F40_rel_error)
plt <- ggplot(all_res, aes(x = Ecov_effect, y = bias_est, colour = re_config)) + scale_colour_viridis_d() + 
    ggtitle(bquote(EM:~beta[M] == log(0.2)*","~beta[Ecov]~Estimated)) +
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor  ~ NAA_M_re + Fhist + obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
    theme_bw() + coord_cartesian(ylim = c(-0.5, 0.5)) + ylab(bquote(Median~Relative~Error~of~hat(F)[40])) + xlab(expression(beta[Ecov])) +
    theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM process error") +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .05, position = position_dodge(0.1))
plt
ggsave(here("Ecov_study","mortality", "paper", "F40_rel_error_M_fixed_Ecov_beta_estimated.png"), plt, width = 20, height = 12, units = "in")
remove(all_res)

all_res <- plot_df_fn(df.ems, df.oms, Ecov_est = TRUE, M_est = TRUE, error_res = F40_rel_error)
plt <- ggplot(all_res, aes(x = Ecov_effect, y = bias_est, colour = re_config)) + scale_colour_viridis_d() + 
    ggtitle(bquote(EM:~beta[M]~Estimated*","~beta[Ecov]~Estimated)) +
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor  ~ NAA_M_re + Fhist + obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
    theme_bw() + coord_cartesian(ylim = c(-0.5, 0.5)) + ylab(bquote(Median~Relative~Error~of~hat(F)[40])) + xlab(expression(beta[Ecov])) +
    theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM process error") +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .05, position = position_dodge(0.1))
plt
ggsave(here("Ecov_study","mortality", "paper", "F40_rel_error_M_estimated_Ecov_beta_estimated.png"), plt, width = 20, height = 12, units = "in")
remove(all_res)
