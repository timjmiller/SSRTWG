library(here)
library(ggplot2)
library(dplyr)
df.ems = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.ems.RDS"))
df.oms = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.oms.RDS"))
aic_res <- readRDS(file.path(here(),"Ecov_study","mortality", "results", "aic_results.RDS"))

aic_fn <- function(all, est_ind, om_ind = NULL){
  if(!is.null(om_ind)) all <- all[om_ind]
  out <- sapply(all, function(res){
    tmp <- sapply(res[est_ind], function(x) return(x))
    #print(est_ind)
    #print(length(res))
    #print(tmp)
    tmp = apply(tmp,1, function(x) {
      if(any(!is.na(x))) {
        return(x == min(x,na.rm=T))
      } else return(rep(NA, length(x)))
    })
    return(apply(tmp,1,sum,na.rm=T))
  })  
  return(out)
}

aic_rank_df_fn <- function(df.ems, df.oms, M_est = FALSE) {
  for(i in 1:3) {
    
    re_mod <- c("rec", "rec+1", "rec+M")[i]
    #EM:  M fixed, OM and EM RE assumption match
    em_ind <- which(df.ems$M_est == M_est & df.ems$re_config == re_mod)
    om_ind <- which(df.oms$NAA_M_re == re_mod) #om and em match
    #Mfixed_rec <- t(aic_wt_fn(aic_res, em_ind, om_ind))
    Mfixed_rec <- t(aic_fn(aic_res, em_ind, om_ind))
    sapply(om_ind, function(i) {
      length(aic_res[[i]][[em_ind[2]]])
    })
    colnames(Mfixed_rec) <- c("effect_est", "effect_0")
    res <- cbind(df.oms[om_ind,], Mfixed_rec)
    temp <- t(sapply(1:NROW(res), function(x) {
      if(res$Ecov_effect[x] == 0) out <- c(res$effect_0[x], res$effect_est[x])
      if(res$Ecov_effect[x] > 0) out <- c(res$effect_est[x], res$effect_0[x])
      return(out)
    }))
    colnames(temp) <- c("right","wrong")
    res <- cbind(res,temp)
    if(i == 1) {
      all_res <- res
    } else {
      all_res <- rbind(all_res, res)
    }
  }
  all_res$emp_p <- all_res$right/(all_res$right + all_res$wrong)
  all_res$ci_lo <- sapply(1:NROW(all_res), function(x) {
    binom.test(all_res$right[x], all_res$right[x] + all_res$wrong[x], p = all_res$emp_p[x])$conf.int[1]
  })
  all_res$ci_hi <- sapply(1:NROW(all_res), function(x) {
    binom.test(all_res$right[x], all_res$right[x] + all_res$wrong[x], p = all_res$emp_p[x])$conf.int[2]
  })
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
      "rec" = "OM=EM: R",
      "rec+1" = "OM=EM: R+S",
      "rec+M" = "OM=EM: R+M"
    ))
  all_res_mod <- all_res_mod %>%
    mutate(Fhist = recode(Fhist,
      "H-MSY" = "High->FMSY",
      "MSY" = "FMSY"
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

all_res <- aic_rank_df_fn(df.ems, df.oms, M_est = FALSE)

plt <- ggplot(all_res, aes(x = Ecov_effect, y = emp_p, colour = Ecov_re_sig:Ecov_re_cor)) + scale_colour_viridis_d() + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig+obs_error ~ NAA_M_re+Fhist, 
      labeller = labeller(Ecov_re_sig = label_parsed, obs_error = label_wrap_gen(width = 40), Ecov_obs_sig = label_parsed)) +
    # facet_grid(Ecov_obs_sig:obs_error ~ NAA_M_re:Fhist, labeller = label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(correct Ecov effect assumption)") + xlab(expression(beta[Ecov])) +
    ggtitle(bquote(beta[M]==log(0.2))) + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = bquote(sigma[Ecov]:rho[Ecov])) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .05, position = position_dodge(0.1)) + 
    geom_hline(aes(yintercept=0.5), linewidth = 1, linetype = "dashed", colour = "red") + scale_x_continuous(breaks = c(0, 0.25, 0.5))
plt

ggsave(here("Ecov_study","mortality", "paper", "proportion_correct_effect_M_fixed.png"), plt, width = 20, height = 12, units = "in")
remove(all_res, all_res_mod)

all_res <- aic_rank_df_fn(df.ems, df.oms, M_est = TRUE)

plt <- ggplot(all_res, aes(x = Ecov_effect, y = emp_p, colour = Ecov_re_sig:Ecov_re_cor)) + scale_colour_viridis_d() + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig+obs_error ~ NAA_M_re+Fhist, 
      labeller = labeller(Ecov_re_sig = label_parsed, obs_error = label_wrap_gen(width = 40), Ecov_obs_sig = label_parsed)) +
    # facet_grid(Ecov_obs_sig:obs_error ~ NAA_M_re:Fhist, labeller = label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(correct Ecov effect assumption)") + xlab(expression(beta[Ecov])) +
    ggtitle(bquote(beta[M]~estimated)) + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = bquote(sigma[Ecov]:rho[Ecov])) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .05, position = position_dodge(0.1)) + 
    geom_hline(aes(yintercept=0.5), linewidth = 1, linetype = "dashed", colour = "red") + scale_x_continuous(breaks = c(0, 0.25, 0.5))
plt

ggsave(here("Ecov_study","mortality", "paper", "proportion_correct_effect_M_estimated.png"), plt, width = 20, height = 12, units = "in")
remove(all_res, all_res_mod)

aic_wt_fn <- function(all, est_ind, om_ind = NULL){
  if(!is.null(om_ind)) all <- all[om_ind]
  out <- sapply(all, function(res){
    tmp <- sapply(res[est_ind], function(x) return(x))
    tmp <- apply(tmp,1, function(x) {
      if(any(!is.na(x))) {
        out <- rep(0, length(x))
        out <- exp(-0.5*(x-min(x,na.rm =TRUE)))
        out[which(is.na(out))] <- 0
        out <- out/sum(out)
        return(out)
      } else return(rep(NA, length(x)))
    })
    return(apply(tmp,1,mean,na.rm=T))
  })  
  return(out)
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

for(i in 1:3) {
  
  re_mod <- c("rec", "rec+1", "rec+M")[i]
  #EM:  M fixed, OM and EM RE assumption match
  em_ind <- which(df.ems$M_est == FALSE & df.ems$re_config == re_mod)
  om_ind <- which(df.oms$NAA_M_re == re_mod) #om and em match
  Mfixed_rec <- t(aic_wt_fn(aic_res, em_ind, om_ind))
  #Mfixed_rec <- t(aic_fn(aic_res, em_ind, om_ind))
  sapply(om_ind, function(i) {
    length(aic_res[[i]][[em_ind[2]]])
  })
  colnames(Mfixed_rec) <- c("effect_est", "effect_0")
  res <- cbind(df.oms[om_ind,], Mfixed_rec)
  temp <- t(sapply(1:NROW(res), function(x) {
    if(res$Ecov_effect[x] == 0) out <- c(res$effect_0[x], res$effect_est[x])
    if(res$Ecov_effect[x] > 0) out <- c(res$effect_est[x], res$effect_0[x])
    return(out)
  }))
  colnames(temp) <- c("right","wrong")
  res <- cbind(res,temp)
  if(i == 1) {
    all_res <- res
  } else {
    all_res <- rbind(all_res, res)
  }
}
all_res$emp_p <- all_res$right/(all_res$right + all_res$wrong)
all_res$ci_lo <- sapply(1:NROW(all_res), function(x) {
  binom.test(all_res$right[x], all_res$right[x] + all_res$wrong[x], p = all_res$emp_p[x])$conf.int[1]
})
all_res$ci_hi <- sapply(1:NROW(all_res), function(x) {
  binom.test(all_res$right[x], all_res$right[x] + all_res$wrong[x], p = all_res$emp_p[x])$conf.int[2]
})


# mod <- glm(cbind(right,wrong) ~  factor(NAA_M_re), data = all_res, family = binomial)
# glm(cbind(right,wrong) ~ factor(NAA_M_re)*factor(Ecov_effect), data = all_res, family = binomial)
# stepfits <- step(mod, scope = ~ factor(Ecov_obs_sig) + factor(Ecov_re_sig) + factor(Ecov_re_cor) + factor(Ecov_effect) + factor(Fhist) + factor(obs_error) + factor(NAA_M_re))
# mod <- glm(cbind(right,wrong) ~ factor(Ecov_obs_sig) + factor(Ecov_re_sig) + factor(Ecov_re_cor) + factor(Ecov_effect) + factor(Fhist) + factor(obs_error) + factor(NAA_M_re), data = all_res, family = binomial)
# stepfits <-step(mod, scope = ~ (factor(Ecov_obs_sig) + factor(Ecov_re_sig) + factor(Ecov_re_cor) + factor(Ecov_effect) + factor(Fhist) + factor(obs_error) + factor(NAA_M_re))^2, direction = "forward")
# mod <- glm(cbind(right,wrong) ~ (factor(Ecov_obs_sig) + factor(Ecov_re_sig) + factor(Ecov_effect) + factor(Fhist) + factor(NAA_M_re) + factor(Ecov_re_cor) + factor(obs_error))^2  - factor(Ecov_obs_sig):factor(obs_error) - factor(Ecov_re_cor):factor(Ecov_re_sig) - factor(Ecov_re_cor):factor(Ecov_obs_sig), data = all_res, family = binomial)
# stepfits <-step(mod, scope = ~ (factor(Ecov_obs_sig) + factor(Ecov_re_sig) + factor(Ecov_re_cor) + factor(Ecov_effect) + factor(Fhist) + factor(obs_error) + factor(NAA_M_re))^3, direction = "forward")
# mod <- glm(cbind(right,wrong) ~ (factor(Ecov_obs_sig) + factor(Ecov_re_sig) + factor(Ecov_effect) + factor(Fhist) + factor(NAA_M_re) + factor(Ecov_re_cor) + factor(obs_error))^2  - factor(Ecov_re_cor):factor(obs_error), data = all_res, family = binomial)
# final_form <- paste0("cbind(right,wrong) ~ ", as.character(stepfits$formula)[[3]], " - factor(Ecov_re_sig):factor(Ecov_effect):factor(obs_error) - factor(Ecov_re_sig):factor(Fhist):factor(obs_error)")
# mod <- glm(as.formula(final_form), data = all_res, family = binomial)
# stepfits <-step(mod, scope = ~ (factor(Ecov_obs_sig) + factor(Ecov_re_sig) + factor(Ecov_re_cor) + factor(Ecov_effect) + factor(Fhist) + factor(obs_error) + factor(NAA_M_re))^4, direction = "forward")

# all_res$predict_p <- predict(mod, type = "response")
# library(ggplot2)
# plt <- ggplot(all_res, aes(x = factor(Ecov_obs_sig):factor(Fhist), y = predict_p, fill = factor(Ecov_re_sig):factor(Ecov_re_cor))) + 
#     geom_col(position = "dodge") + facet_grid(NAA_M_re ~ factor(Ecov_effect):factor(obs_error), labeller = label_both)
# plt

#facs <- c("Ecov_obs_sig", "Fhist", "Ecov_re_sig","Ecov_re_cor","Ecov_effect", "obs_error", "NAA_M_re")
facs <- c("Ecov_obs_sig", "Fhist", "Ecov_re_sig","Ecov_re_cor", "obs_error", "NAA_M_re")
all_res[facs] <- lapply(all_res[facs], factor)
all_res_mod <- all_res %>%
  mutate(Ecov_obs_sig = recode(Ecov_obs_sig,
    "0.1" = "Ecov obs SD = 0.1",
    "0.5" = "Ecov obs SD = 0.5"
  ))
all_res_mod <- all_res_mod %>%
  mutate(obs_error = recode(obs_error,
    "L" = "Low obs error (indices, age comp)",
    "H" = "High obs error (indices, age comp)"
  ))
all_res_mod <- all_res_mod %>%
  mutate(NAA_M_re = recode(NAA_M_re,
    "rec" = "R",
    "rec+1" = "NAA",
    "rec+M" = "R + M"
  ))

plt <- ggplot(all_res_mod, aes(x = Ecov_effect, y = emp_p, colour = Ecov_re_sig:Ecov_re_cor)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig:obs_error ~ NAA_M_re:Fhist, labeller = label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(correct Ecov effect assumption)") + xlab("Ecov effect size") +
    ggtitle("EM: M = 0.2") + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .05, position = position_dodge(0.1)) + 
    geom_hline(aes(yintercept=0.5), linewidth = 1, linetype = "dashed", colour = "red")
plt

ggsave(here("Ecov_study","mortality", "paper", "proportion_correct_effect_M_fixed.png"), plt)
remove(all_res, all_res_mod)

for(i in 1:3) {
  re_mod <- c("rec", "rec+1", "rec+M")[i]
  #EM:  M estimated, OM and EM RE assumption match
  em_ind <- which(df.ems$M_est == TRUE & df.ems$re_config == re_mod)
  om_ind <- which(df.oms$NAA_M_re == re_mod) #om and em match
  Mfixed_rec <- t(aic_fn(aic_res, em_ind, om_ind))
  colnames(Mfixed_rec) <- c("effect_est", "effect_0")
  res <- cbind(df.oms[om_ind,], Mfixed_rec)
  temp <- t(sapply(1:NROW(res), function(x) {
    if(res$Ecov_effect[x] == 0) out <- c(res$effect_0[x], res$effect_est[x])
    if(res$Ecov_effect[x] > 0) out <- c(res$effect_est[x], res$effect_0[x])
    return(out)
  }))
  colnames(temp) <- c("right","wrong")
  res <- cbind(res,temp)
  if(i == 1) {
    all_res <- res
  } else {
    all_res <- rbind(all_res, res)
  }
}
all_res$emp_p <- all_res$right/(all_res$right + all_res$wrong)
facs <- c("Ecov_obs_sig", "Fhist", "Ecov_re_sig","Ecov_re_cor","obs_error", "NAA_M_re")
all_res[facs] <- lapply(all_res[facs], factor)
all_res$emp_p <- all_res$right/(all_res$right + all_res$wrong)
all_res$ci_lo <- sapply(1:NROW(all_res), function(x) {
  binom.test(all_res$right[x], all_res$right[x] + all_res$wrong[x], p = all_res$emp_p[x])$conf.int[1]
})
all_res$ci_hi <- sapply(1:NROW(all_res), function(x) {
  binom.test(all_res$right[x], all_res$right[x] + all_res$wrong[x], p = all_res$emp_p[x])$conf.int[2]
})

all_res_mod <- all_res %>%
  mutate(Ecov_obs_sig = recode(Ecov_obs_sig,
    "0.1" = "Ecov obs SD = 0.1",
    "0.5" = "Ecov obs SD = 0.5"
  ))
all_res_mod <- all_res_mod %>%
  mutate(obs_error = recode(obs_error,
    "L" = "Low obs error (indices, age comp)",
    "H" = "High obs error (indices, age comp)"
  ))
all_res_mod <- all_res_mod %>%
  mutate(NAA_M_re = recode(NAA_M_re,
    "rec" = "R",
    "rec+1" = "NAA",
    "rec+M" = "R + M"
  ))

plt <- ggplot(all_res_mod, aes(x = Ecov_effect, y = emp_p, colour = Ecov_re_sig:Ecov_re_cor)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig:obs_error ~ NAA_M_re:Fhist, labeller = label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(correct Ecov effect assumption)") + xlab("Ecov effect size") +
    ggtitle("EM: M estimated") + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .05, position = position_dodge(0.1)) + 
    geom_hline(aes(yintercept=0.5), linewidth = 1, linetype = "dashed", colour = "red")
plt
ggsave(here("Ecov_study","mortality", "paper", "proportion_correct_effect_M_estimated.png"), plt)

make_df_fn <- function(M_est = TRUE){
  #all EM PE assumptions
  for(i in 1:3) {
    re_mod <- c("rec", "rec+1", "rec+M")[i]
    #EM:  M estimated
    em_ind <- which(df.ems$M_est == M_est) #all EM PE assumptions
    om_ind <- which(df.oms$NAA_M_re == re_mod) #om and em match
    Mfixed_rec <- aic_fn(aic_res, em_ind, om_ind)
    res <- cbind(df.oms[rep(om_ind, each = length(em_ind)),], df.ems[rep(em_ind, length(om_ind)),], n = c(Mfixed_rec))
    if(i == 1) {
      all_res <- res
    } else {
      all_res <- rbind(all_res, res)
    }
  }
  facs <- c("Ecov_obs_sig", "Fhist", "Ecov_re_sig","Ecov_re_cor","obs_error", "NAA_M_re", "re_config", "Ecov_est")
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
  all_res_mod <- all_res_mod %>% mutate(Fhist = recode(Fhist,
      "H-MSY" = "High->FMSY",
      "MSY" = "FMSY"))
  all_res_mod <- all_res_mod %>%
    mutate(beta_ecov_fixed = recode(Ecov_est,
      "TRUE" = "beta[Ecov]==0",
      "FALSE" = "beta[Ecov]~estimated"
    ))

  #1: correct effect assumption and RE will be 1, 
  #2: correct RE, wrong effect assumption
  #3: wrong RE, correct effect assumption
  #4: wrong RE, wrong effect assumption
  all_res_mod$correct <- 0 
  all_res_mod$correct[all_res_mod$re_config == all_res_mod$NAA_M_re & all_res_mod$Ecov_effect == 0 & all_res_mod$Ecov_est == FALSE] <- 1
  all_res_mod$correct[all_res_mod$re_config == all_res_mod$NAA_M_re & all_res_mod$Ecov_effect > 0 & all_res_mod$Ecov_est == TRUE] <- 1
  all_res_mod$correct[all_res_mod$re_config == all_res_mod$NAA_M_re & all_res_mod$Ecov_effect == 0 & all_res_mod$Ecov_est == TRUE] <- 2
  all_res_mod$correct[all_res_mod$re_config == all_res_mod$NAA_M_re & all_res_mod$Ecov_effect > 0 & all_res_mod$Ecov_est == FALSE] <- 2
  all_res_mod$correct[all_res_mod$re_config != all_res_mod$NAA_M_re & all_res_mod$Ecov_effect == 0 & all_res_mod$Ecov_est == FALSE] <- 3
  all_res_mod$correct[all_res_mod$re_config != all_res_mod$NAA_M_re & all_res_mod$Ecov_effect > 0 & all_res_mod$Ecov_est == TRUE] <- 3
  all_res_mod$correct[all_res_mod$re_config != all_res_mod$NAA_M_re & all_res_mod$Ecov_effect == 0 & all_res_mod$Ecov_est == TRUE] <- 4
  all_res_mod$correct[all_res_mod$re_config != all_res_mod$NAA_M_re & all_res_mod$Ecov_effect > 0 & all_res_mod$Ecov_est == FALSE] <- 4
  all_res_mod$correct <- factor(all_res_mod$correct)
  temp <- all_res_mod 
  all_res_mod$Ecov_effect <- factor(all_res_mod$Ecov_effect)
  all_res_mod <- all_res_mod %>%
    mutate(NAA_M_re = recode(NAA_M_re,
      "R" = "OM: R",
      "R+S" = "OM: R+S",
      "R+M" = "OM: R+M"
    ))

  all_res_mod <- all_res_mod %>%
    mutate(correct = recode(correct,
      "1" = "Correct Effect, Correct PE",
      "2" = "Wrong Effect, Correct PE",
      "3" = "Correct Effect, Wrong PE",
      "4" = "Wrong Effect, Wrong PE"
    ))
    return(all_res_mod)
}

all_res_mod <- make_df_fn(M_est = TRUE)
plt <- ggplot(all_res_mod, aes(x = Ecov_effect, y = n, fill = correct)) + scale_fill_viridis_d() + 
    geom_col(position = "fill") + 
    #geom_col(position = position_dodge(0.1)) + 
    facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ NAA_M_re+Fhist+obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed)) +
    #facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ NAA_M_re+Fhist+obs_error, labeller = label_wrap_gen(width = 15)) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("Proportion ranked best") + xlab(expression(beta[Ecov])) + labs(fill = "EM assumption") +
    ggtitle(bquote(beta[M]~estimated)) + theme(plot.title = element_text(hjust = 0.5))
plt
ggsave(here("Ecov_study","mortality", "paper", "proportion_correct_PE_effect_M_estimated.png"), plt, height = 12, width = 20, units = "in")
temp <- all_res_mod %>%
  mutate(beta_ecov_fixed = recode(Ecov_est,
    "TRUE" = "NO",
    "FALSE" = "YES"
  ))
temp$beta_ecov_fixed <- factor(temp$beta_ecov_fixed)
plt <- ggplot(temp, aes(x = Ecov_effect, y = n, fill = re_config:beta_ecov_fixed)) + scale_fill_viridis_d() + 
    geom_col(position = "fill") + 
    #geom_col(position = position_dodge(0.1)) + 
    facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ NAA_M_re+Fhist+obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed)) +
    #facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ NAA_M_re+Fhist+obs_error, labeller = label_wrap_gen(width = 15)) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("Proportion ranked best") + xlab(expression(beta[Ecov])) + labs(fill = bquote(atop("EM Assumptions",PE*":"*beta[Ecov]==0))) +
    ggtitle(bquote(beta[M]~estimated)) + theme(plot.title = element_text(hjust = 0.5))
plt
ggsave(here("Ecov_study","mortality", "paper", "proportion_best_AIC_M_estimated.png"), plt, height = 12, width = 20, units = "in")

all_res_mod <- make_df_fn(M_est = FALSE)
plt <- ggplot(all_res_mod, aes(x = Ecov_effect, y = n, fill = correct)) + scale_fill_viridis_d() + 
    geom_col(position = "fill") + 
    #geom_col(position = position_dodge(0.1)) + 
    facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ NAA_M_re+Fhist+obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed)) +
    #facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ NAA_M_re+Fhist+obs_error, labeller = label_wrap_gen(width = 15)) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("Proportion ranked best") + xlab(expression(beta[Ecov])) + labs(fill = "EM assumption") +
    ggtitle(bquote(beta[M]==log(0.2))) + theme(plot.title = element_text(hjust = 0.5))
plt
ggsave(here("Ecov_study","mortality", "paper", "proportion_correct_PE_effect_M_fixed.png"), plt, height = 12, width = 20, units = "in")

temp <- all_res_mod %>%
  mutate(beta_ecov_fixed = recode(Ecov_est,
    "TRUE" = "NO",
    "FALSE" = "YES"
  ))
temp$beta_ecov_fixed <- factor(temp$beta_ecov_fixed)
plt <- ggplot(temp, aes(x = Ecov_effect, y = n, fill = re_config:beta_ecov_fixed)) + scale_fill_viridis_d() + 
    geom_col(position = "fill") + 
    #geom_col(position = position_dodge(0.1)) + 
    facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ NAA_M_re+Fhist+obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed)) +
    #facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ NAA_M_re+Fhist+obs_error, labeller = label_wrap_gen(width = 15)) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("Proportion ranked best") + xlab(expression(beta[Ecov])) + labs(fill = bquote(atop("EM Assumptions",PE*":"*beta[Ecov]==0))) +
    ggtitle(bquote(beta[M]==log(0.2))) + theme(plot.title = element_text(hjust = 0.5))
plt
ggsave(here("Ecov_study","mortality", "paper", "proportion_best_AIC_M_fixed.png"), plt, height = 12, width = 20, units = "in")



x <- subset(temp, NAA_M_re == "R+M" & Ecov_effect == "0")
aggregate(x$n, x["re_config"], sum)
x <- subset(temp, NAA_M_re == "R+M" & Ecov_effect == "0.25")
aggregate(x$n, x["re_config"], sum)
x <- subset(temp, NAA_M_re == "R+M" & Ecov_effect == "0.5")
aggregate(x$n, x["re_config"], sum)

g <- ggplot_gtable(ggplot_build(plt))
stript <- grep('strip-t', g$layout$name)
#fills <- c("red","green","blue","yellow")
fills <- cols
k <- 1
for (i in stript) {
j <- grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
k <- k+1
}
grid.draw(g)
