library(here)
library(ggplot2)
library(dplyr)
df.ems = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.ems.RDS"))
df.oms = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.oms.RDS"))
all_beta_bias <- readRDS(file.path(here::here(),"Ecov_study","mortality", "results", "ecov_beta_bias_results.RDS"))

custom_boxplot_stat <- function(x){#, n.yrs = 1, n.sim = 100) {
  x <- x[which(!is.na(x))]
  n <- length(x)
  bnds95 <- qbinom(c(0.025,0.975), n, 0.5)/n # 95% CI bounds for median (sims x years)
  bnds80 <- qbinom(c(0.1,0.9), n, 0.5)/n # 80% CI bounds for median (sims x years)
  r <- quantile(x, probs = c(bnds95[1], bnds80[1], 0.5, bnds80[2], bnds95[2]))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
plot_df_fn <- function(df.ems, df.oms, M_est = FALSE) {
  for(i in 1:3) {
    re_mod <- c("rec", "rec+1", "rec+M")[i]
    #EM:  M fixed, Ecov_beta estimated, OM and EM RE assumption match
    df.ems. <- df.ems[df.ems$Ecov_est,]
    em_ind <- which(df.ems.$M_est == M_est & df.ems.$re_config == re_mod)
    om_ind <- which(df.oms$NAA_M_re == re_mod) #om and em match
    res <- t(sapply(all_beta_bias[om_ind], function(x) {
      x <- x[[em_ind]]
      x <- x[which(!is.na(x[,2])),1] #remove fits without standard error estimates
      #x <- x[which(x[,2] < 100),1] #remove fits with bad (really big) standard error estimates
      #print(x)

      # out <- c(mean(x,na.rm=T), sd(x, na.rm=T)/sqrt(sum(!is.na(x))))
      # out <- c(out, out[1] + qt(0.025, sum(!is.na(x))) * out[2])
      # out <- c(out, out[1] + qt(0.975, sum(!is.na(x))) * out[2])
      out <- c(median(x, na.rm = TRUE), sd(x, na.rm=T)/sqrt(sum(!is.na(x))), custom_boxplot_stat(x))#quantile(x, probs = c(0.025,0.975)))
      return(out)
    }))
    colnames(res) <- c("bias_est", "bias_se", "ymin", "lower", "middle", "upper", "ymax")
    res <- cbind.data.frame(df.oms[om_ind,], res)
    print(head(res))
    if(i == 1) {
      all_res <- res
    } else {
      all_res <- rbind.data.frame(all_res, res)
    }
  }
  facs <- c("Ecov_obs_sig", "Fhist", "Ecov_re_sig","Ecov_re_cor", "obs_error", "NAA_M_re")
  all_res[facs] <- lapply(all_res[facs], factor)
  all_res_mod <- all_res %>%
    mutate(Ecov_obs_sig = recode(Ecov_obs_sig,
      "0.1" = "SD(Ecov obs) = 0.1",
      "0.5" = "SD(Ecov obs) = 0.5"
    ))
  # all_res_mod <- all_res_mod %>%
  #   mutate(Ecov_re_sig = recode(Ecov_re_sig,
  #     "0.1" = "sigma(Ecov) = 0.1",
  #     "0.5" = "sigma(Ecov) = 0.5"
  #   ))
  # all_res_mod <- all_res_mod %>%
  #   mutate(Ecov_re_cor = recode(Ecov_re_cor,
  #     "0" = "rho(Ecov) = 0",
  #     "0.5" = "rho(Ecov) = 0.5"
  #   ))
  all_res_mod <- all_res_mod %>%
    mutate(obs_error = recode(obs_error,
      "L" = "Low obs error (indices, age comp)",
      "H" = "High obs error (indices, age comp)"
    ))
  all_res_mod <- all_res_mod %>%
    mutate(NAA_M_re = recode(NAA_M_re,
      "rec" = "R",
      "rec+1" = "R+S",
      "rec+M" = "R+M"
    ))
  all_res_mod <- all_res_mod %>% mutate(Fhist = recode(Fhist,
      "H-MSY" = "F history: High->FMSY",
      "MSY" = "F history: FMSY"))
  return(all_res_mod)
}

all_res <- plot_df_fn(df.ems, df.oms, M_est = FALSE)
plt <- ggplot(all_res, aes(x = Ecov_effect, y = bias_est, colour = Ecov_re_sig:Ecov_re_cor)) + scale_colour_viridis_d() + 
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig + obs_error ~ NAA_M_re + Fhist, labeller = labeller(obs_error = label_wrap_gen(width = 35))) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(-1, 1)) + ylab(bquote(Median~Bias~of~beta[Ecov])) + xlab(expression(beta[Ecov])) +
    ggtitle(bquote(beta[M]==log(0.2))) + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = expression(sigma[Ecov]:rho[Ecov])) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .05, position = position_dodge(0.1))
plt
ggsave(here("Ecov_study","mortality", "paper", "Ecov_beta_bias_M_fixed.png"), plt, width = 20, height = 12, units = "in")
remove(all_res)

all_res <- plot_df_fn(df.ems, df.oms, M_est = TRUE)
plt <- ggplot(all_res, aes(x = Ecov_effect, y = bias_est, colour = Ecov_re_sig:Ecov_re_cor)) + scale_colour_viridis_d() + 
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig + obs_error ~ NAA_M_re + Fhist, labeller = labeller(obs_error = label_wrap_gen(width = 35))) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(-1, 1)) + ylab(bquote(Median~Bias~of~beta[Ecov])) + xlab(expression(beta[Ecov])) +
    ggtitle(bquote(beta[M]~Estimated)) + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = expression(sigma[Ecov]:rho[Ecov])) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .05, position = position_dodge(0.1))
plt
ggsave(here("Ecov_study","mortality", "paper", "Ecov_beta_bias_M_estimated.png"), plt, width = 20, height = 12, units = "in")
remove(all_res)

#all EM PE assumptions
plot_df_fn <- function(df.ems, df.oms, M_est = FALSE) {
  for(h in 1:3) {
    re_mod <- c("rec", "rec+1", "rec+M")[h]
    #EM:  M fixed, Ecov_beta estimated, OM and EM RE assumption match
    df.ems. <- df.ems[df.ems$Ecov_est,]
    em_ind <- which(df.ems.$M_est == M_est) # all 3 EM re_mods  
    om_ind <- which(df.oms$NAA_M_re == re_mod) #om and em match
    for(i in om_ind){
      out <- matrix(NA,length(em_ind), 7)
      for(j in em_ind){
        x <- all_beta_bias[[i]][[j]]
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
      "0.1" = "SD(Ecov obs) = 0.1",
      "0.5" = "SD(Ecov obs) = 0.5"
    ))
  all_res_mod <- all_res_mod %>%
    mutate(obs_error = recode(obs_error,
      "L" = "Low obs error (indices, age comp)",
      "H" = "High obs error (indices, age comp)"
    ))
  all_res_mod <- all_res_mod %>%
    mutate(NAA_M_re = recode(NAA_M_re,
      "rec" = "R",
      "rec+1" = "R+S",
      "rec+M" = "R+M"
    ))
  all_res_mod <- all_res_mod %>%
    mutate(re_config = recode(re_config,
      "rec" = "R",
      "rec+1" = "R+S",
      "rec+M" = "R+M"
    ))
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
library(scales)
par(mfrow = c(2,4))
for(i in LETTERS[1:8]) show_col(viridis_pal(option=i, begin = 0.5)(4))
cols = c(viridis_pal(option="E", begin = 0.5, alpha = 0.5)(4),viridis_pal(option="G", begin = 0.5, alpha = 0.5)(4),viridis_pal(option="H", begin = 0.5, alpha = 0.5)(4))

all_res <- plot_df_fn(df.ems, df.oms, M_est = TRUE)
plt <- ggplot(all_res, aes(x = Ecov_effect, y = bias_est, colour = re_config)) + scale_colour_viridis_d() + 
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor  ~ NAA_M_re + Fhist + obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_wrap_gen(width = 15))) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(-0.5, 0.5)) + ylab(bquote(Median~Bias~of~beta[Ecov])) + xlab(expression(beta[Ecov])) +
    ggtitle(bquote(beta[M]~Estimated)) + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM process error") +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .05, position = position_dodge(0.1))
plt
ggsave(here("Ecov_study","mortality", "paper", "Ecov_beta_bias_all_PE_effect_M_estimated.png"), plt, width = 20, height = 12, units = "in")
remove(all_res)

all_res <- plot_df_fn(df.ems, df.oms, M_est = FALSE)

plt <- ggplot(all_res, aes(x = Ecov_effect, y = bias_est, colour = re_config)) + scale_colour_viridis_d() + 
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor  ~ NAA_M_re + Fhist + obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_wrap_gen(width = 15))) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(-0.5, 0.5)) + ylab(bquote(Median~Bias~of~beta[Ecov])) + xlab(expression(beta[Ecov])) +
    ggtitle(bquote(beta[M]==log(0.2))) + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "EM process error") +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .05, position = position_dodge(0.1))
plt
ggsave(here("Ecov_study","mortality", "paper", "Ecov_beta_bias_all_PE_effect_M_fixed.png"), plt, width = 20, height = 12, units = "in")
remove(all_res)
