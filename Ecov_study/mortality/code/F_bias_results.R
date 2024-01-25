library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
df.ems = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.ems.RDS"))
df.oms = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.oms.RDS"))
F_bias <- readRDS(file.path(here::here(),"Ecov_study","mortality", "results", "F_results.RDS"))
library(reshape2)
#res <- melt(F_bias)

custom_boxplot_stat <- function(x){#, n.yrs = 1, n.sim = 100) {
  x <- x[which(!is.na(x))]
  n <- length(x)
  bnds95 <- qbinom(c(0.025,0.975), n, 0.5)/n # 95% CI bounds for median (sims x years)
  bnds80 <- qbinom(c(0.1,0.9), n, 0.5)/n # 80% CI bounds for median (sims x years)
  r <- quantile(x, probs = c(bnds95[1], bnds80[1], 0.5, bnds80[2], bnds95[2]))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
plot_df_fn <- function(df.ems, df.oms, M_est=FALSE) {

  print(which(df.oms$NAA_M_re == "rec+1"))
  for(i in 1:3) {
    re_mod <- c("rec", "rec+1", "rec+M")[i]
    #EM:  M fixed, mean_M estimated, OM and EM RE assumption match
    em_ind <- which(df.ems$Ecov_est & df.ems$M_est == M_est & df.ems$re_config == re_mod)
#    print(em_ind)
    om_ind <- which(df.oms$NAA_M_re == re_mod) #om and em match
#    print(om_ind)
    res <- lapply(om_ind, function(x) {
      F_res <- t(sapply(1:40, function(y) {
        #print(om_ind)
        r_e <- F_bias[[x]][,em_ind,y,2]/F_bias[[x]][,em_ind,y,1]-1
        out <- c(median(r_e, na.rm = TRUE), sd(r_e, na.rm=T)/sqrt(sum(!is.na(r_e))), custom_boxplot_stat(r_e))#quantile(x, probs = c(0.025,0.975)))
        return(out)
      }))
      colnames(F_res) <- c("bias_est", "bias_se", "ymin", "lower", "middle", "upper", "ymax")
      return(F_res)
    })
    # print(res[[1]])
    # print(dim(res[[1]]))
    res <- reshape2::melt(res)
    colnames(res) <- c("year", "type", "value", "om")
    res$om <- om_ind[res$om]

    # print(head(res))
    # print(table(res[,1]))
    # print(table(res[,2]))
    # print(table(res[,4]))
    res <- cbind.data.frame(df.oms[res$om,], res)
    res<- res %>% tidyr::pivot_wider(names_from = type, values_from = value) %>% as.data.frame
     # print(head(res))
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
      "rec" = "R",
      "rec+1" = "R+S",
      "rec+M" = "R+M"
    ))
  all_res_mod <- all_res_mod %>% mutate(Fhist = recode(Fhist,
      "H-MSY" = "F history: High->FMSY",
      "MSY" = "F history: FMSY"))
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
all_res <- plot_df_fn(df.ems, df.oms)
temp <- filter(all_res, year %in% c(1,21,40))
temp$Ecov_effect <- factor(temp$Ecov_effect)
temp$year <- factor(temp$year)
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
plt <- ggplot(temp, aes(x = Ecov_effect, y = bias_est, colour = year)) + scale_colour_viridis_d(alpha = 0.7) + 
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    #geom_line(position = position_dodge(0.2), linewidth = 1) + 
    geom_point(position = position_dodge(0.6), size = 4) + 
    facet_grid(Ecov_obs_sig + Ecov_re_sig + Ecov_re_cor ~ NAA_M_re + Fhist+ obs_error , 
      labeller = labeller(obs_error = label_wrap_gen(width = 20), Ecov_obs_sig = label_parsed, Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed)) + 
    theme_bw() + 
    #coord_cartesian(ylim = c(-0.25, 0.25)) + 
    ylab(bquote(Median~relative~error*"(F)")) + xlab(expression(beta[Ecov])) +
    theme(plot.title = element_text(hjust = 0.6)) + labs(colour = "Time") +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .5, position = position_dodge(0.6), linewidth = 1)
plt
ggsave(here("Ecov_study","mortality", "paper", "F_bias_OM_EM_PE_match_M_fixed.png"), plt, width = 20, height = 12, units = "in")
remove(all_res)

all_res <- plot_df_fn(df.ems, df.oms, M_est = TRUE)
temp <- filter(all_res, year %in% c(1,21,40))
temp$Ecov_effect <- factor(temp$Ecov_effect)
temp$year <- factor(temp$year)
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
plt <- ggplot(temp, aes(x = Ecov_effect, y = bias_est, colour = year)) + scale_colour_viridis_d(alpha = 0.7) + 
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    #geom_line(position = position_dodge(0.2), linewidth = 1) + 
    geom_point(position = position_dodge(0.6), size = 4) + 
    facet_grid(Ecov_obs_sig + Ecov_re_sig + Ecov_re_cor ~ NAA_M_re + Fhist+ obs_error , 
      labeller = labeller(obs_error = label_wrap_gen(width = 20), Ecov_obs_sig = label_parsed, Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed)) + 
    theme_bw() + 
    coord_cartesian(ylim = c(-0.5, 0.5)) + 
    ylab(bquote(Median~relative~error*"(F)")) + xlab(expression(beta[Ecov])) +
    theme(plot.title = element_text(hjust = 0.6)) + labs(colour = "Time") +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .5, position = position_dodge(0.6), linewidth = 1)
plt
ggsave(here("Ecov_study","mortality", "paper", "F_bias_OM_EM_PE_match_M_estimated.png"), plt, width = 20, height = 12, units = "in")
remove(all_res)

#all EM PE assumptions
plot_df_fn <- function(df.ems, df.oms, M_est = FALSE, Ecov_est = FALSE) {
  for(h in 1:3) {
    re_mod <- c("rec", "rec+1", "rec+M")[h]
    #EM:  M fixed, mean_M estimated, OM and EM RE assumption match
    em_ind <- which(df.ems$Ecov_est==Ecov_est & df.ems$M_est == M_est) # all 3 EM re_mods  
    om_ind <- which(df.oms$NAA_M_re == re_mod) #om and em match
    # for(i in om_ind){
    #   out <- matrix(NA,length(em_ind), 7)
    #   for(j in em_ind){
    res <- lapply(om_ind, function(x) {
      F_res <- lapply(1:40, function(y) {
        #print(om_ind)
        out <- matrix(NA,length(em_ind), 7)
        for(j in em_ind){
          r_e <- F_bias[[x]][,j,y,2]/F_bias[[x]][,j,y,1]-1
          out[which(em_ind==j),] <- c(median(r_e, na.rm = TRUE), sd(r_e, na.rm=T)/sqrt(sum(!is.na(r_e))), custom_boxplot_stat(r_e))
        }
#        print(out)
        return(out)
        colnames(out) <- c("bias_est", "bias_se", "ymin", "lower", "middle", "upper", "ymax")
      })
      # print(length(F_res))
      # print(F_res[[1]])
      return(F_res)
    })
    # print(res[[1]])
    # print(dim(res[[1]]))
    res <- reshape2::melt(res)
    print(head(res))
    print(table(res[,1]))
    print(table(res[,2]))
    print(table(res[,4]))
    colnames(res) <- c("em_config","type","value","year", "om")
    res$om <- om_ind[res$om]

    res <- cbind.data.frame(df.oms[res$om,], res)
#      colnames(out) <- c("bias_est", "bias_se", "ymin", "lower", "middle", "upper", "ymax")
      # if(i == om_ind[1]) res <- cbind.data.frame(df.oms[rep(i,length(em_ind)),], df.ems[em_ind,], out)
      # else res <- rbind.data.frame(res,cbind.data.frame(df.oms[rep(i,length(em_ind)),], df.ems[em_ind,], out))
    print(head(res))
    if(h == 1) {
      all_res <- res
    } else {
      all_res <- rbind.data.frame(all_res, res)
    }
  }
  all_res$type <- c("bias_est", "bias_se", "ymin", "lower", "middle", "upper", "ymax")[all_res$type]
  all_res$em_config <- c("rec", "rec+1", "rec+M")[all_res$em_config]
  print(head(all_res))

  all_res<- all_res %>% tidyr::pivot_wider(names_from = type, values_from = value) %>% as.data.frame
  print(head(all_res))
  facs <- c("Ecov_obs_sig", "Fhist", "Ecov_re_sig","Ecov_re_cor", "obs_error", "NAA_M_re", "em_config")
  all_res[facs] <- lapply(all_res[facs], factor)
  all_res_mod <- all_res
  # all_res_mod <- all_res %>%
  #   mutate(Ecov_obs_sig = recode(Ecov_obs_sig,
  #     "0.1" = "sigma[ecov] == 0.1",
  #     "0.5" = "sigma[ecov] == 0.5"
  #   ))
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
    mutate(em_config = recode(em_config,
      "rec" = "EM: R",
      "rec+1" = "EM: R+S",
      "rec+M" = "EM: R+M"
    ))
  all_res_mod <- all_res_mod %>% mutate(Fhist = recode(Fhist,
      "H-MSY" = "2.5*F[MSY] %->% F[MSY]",
      "MSY" = "F[MSY]"))
  all_res_mod <- all_res_mod %>%
    mutate(Ecov_effect = recode(Ecov_effect,
      "0" = "beta[Ecov] == 0",
      "0.25" = "beta[Ecov] == 0.25",
      "0.5" = "beta[Ecov] == 0.5"
    ))
  # all_res_mod <- all_res_mod %>%
  #   mutate(Ecov_re_sig = recode(Ecov_re_sig,
  #     "0.1" = "sigma[Ecov] == 0.1",
  #     "0.5" = "sigma[Ecov] == 0.5"
  #   ))
  # all_res_mod <- all_res_mod %>%
  #   mutate(Ecov_re_cor = recode(Ecov_re_cor,
  #     "0" = "rho[Ecov] == 0",
  #     "0.5" = "rho[Ecov] == 0.5"
  #   ))
  return(all_res_mod)
}
all_res <- plot_df_fn(df.ems, df.oms)
temp <- filter(all_res, year %in% c(1,21,40))
temp$Ecov_effect <- factor(temp$Ecov_effect)
temp$year <- factor(temp$year)
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
pdodge <- 0.8
plt <- ggplot(temp, aes(x = year, y = bias_est, colour = Ecov_obs_sig:Ecov_re_sig:Ecov_re_cor)) + scale_colour_viridis_d() + 
    ggtitle(bquote(beta[M]==log(0.2)*","~beta[Ecov]==0)) +
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_point(position = position_dodge(pdodge), size = 2) + 
    facet_grid(em_config+Ecov_effect  ~ NAA_M_re + Fhist + obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_effect = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
    theme_bw() +
    coord_cartesian(ylim = c(-0.3, 0.3)) + 
    ylab(bquote(Median~Relative~Error~of~hat(F))) + xlab("Time") + 
    theme(plot.title = element_text(hjust = 0.5), panel.spacing = unit(0.01, "lines")) + labs(colour = bquote(sigma[ecov]:sigma[Ecov]:rho[Ecov])) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .01, linewidth = 1, position = position_dodge(pdodge))
plt
ggsave(here("Ecov_study","mortality", "paper", "F_bias_all_PE_effect_M_fixed_beta_fixed.png"), plt, width = 20, height = 12, units = "in")
remove(all_res)

all_res <- plot_df_fn(df.ems, df.oms, Ecov_est = TRUE)
temp <- filter(all_res, year %in% c(1,21,40))
temp$Ecov_effect <- factor(temp$Ecov_effect)
temp$year <- factor(temp$year)
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
pdodge <- 0.8
plt <- ggplot(temp, aes(x = year, y = bias_est, colour = Ecov_obs_sig:Ecov_re_sig:Ecov_re_cor)) + scale_colour_viridis_d() + 
    ggtitle(bquote(beta[M]==log(0.2)*","~beta[Ecov]~Estimated)) +
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_point(position = position_dodge(pdodge), size = 2) + 
    facet_grid(em_config+Ecov_effect  ~ NAA_M_re + Fhist + obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_effect = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
    theme_bw() +
    coord_cartesian(ylim = c(-0.3, 0.3)) + 
    ylab(bquote(Median~Relative~Error~of~hat(F))) + xlab("Time") + 
    theme(plot.title = element_text(hjust = 0.5), panel.spacing = unit(0.01, "lines")) + labs(colour = bquote(sigma[ecov]:sigma[Ecov]:rho[Ecov])) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .01, linewidth = 1, position = position_dodge(pdodge))
plt
ggsave(here("Ecov_study","mortality", "paper", "F_bias_all_PE_effect_M_fixed_beta_estimated.png"), plt, width = 20, height = 12, units = "in")
remove(all_res)

all_res <- plot_df_fn(df.ems, df.oms, Ecov_est = TRUE, M_est = TRUE)
temp <- filter(all_res, year %in% c(1,21,40))
temp$Ecov_effect <- factor(temp$Ecov_effect)
temp$year <- factor(temp$year)
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
pdodge <- 0.8
plt <- ggplot(temp, aes(x = year, y = bias_est, colour = Ecov_obs_sig:Ecov_re_sig:Ecov_re_cor)) + scale_colour_viridis_d() + 
    ggtitle(bquote(beta[M]~Estimated*","~beta[Ecov]~Estimated)) +
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_point(position = position_dodge(pdodge), size = 2) + 
    facet_grid(em_config+Ecov_effect  ~ NAA_M_re + Fhist + obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_effect = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
    theme_bw() +
    coord_cartesian(ylim = c(-0.53, 0.53)) + 
    ylab(bquote(Median~Relative~Error~of~hat(F))) + xlab("Time") + 
    theme(plot.title = element_text(hjust = 0.5), panel.spacing = unit(0.01, "lines")) + labs(colour = bquote(sigma[ecov]:sigma[Ecov]:rho[Ecov])) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .01, linewidth = 1, position = position_dodge(pdodge))
plt
ggsave(here("Ecov_study","mortality", "paper", "F_bias_all_PE_effect_M_estimated_beta_estimated.png"), plt, width = 20, height = 12, units = "in")
remove(all_res)

all_res <- plot_df_fn(df.ems, df.oms, Ecov_est = FALSE, M_est = TRUE)
temp <- filter(all_res, year %in% c(1,21,40))
temp$Ecov_effect <- factor(temp$Ecov_effect)
temp$year <- factor(temp$year)
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
pdodge <- 0.8
plt <- ggplot(temp, aes(x = year, y = bias_est, colour = Ecov_obs_sig:Ecov_re_sig:Ecov_re_cor)) + scale_colour_viridis_d() + 
    ggtitle(bquote(beta[M]~Estimated*","~beta[Ecov]==0)) +
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_point(position = position_dodge(pdodge), size = 2) + 
    facet_grid(em_config+Ecov_effect  ~ NAA_M_re + Fhist + obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_effect = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
    theme_bw() +
    coord_cartesian(ylim = c(-0.53, 0.53)) + 
    ylab(bquote(Median~Relative~Error~of~hat(F))) + xlab("Time") + 
    theme(plot.title = element_text(hjust = 0.5), panel.spacing = unit(0.01, "lines")) + labs(colour = bquote(sigma[ecov]:sigma[Ecov]:rho[Ecov])) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .01, linewidth = 1, position = position_dodge(pdodge))
plt
ggsave(here("Ecov_study","mortality", "paper", "F_bias_all_PE_effect_M_estimated_beta_fixed.png"), plt, width = 20, height = 12, units = "in")
remove(all_res)
