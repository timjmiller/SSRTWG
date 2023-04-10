library(here)
library(ggplot2)
df.ems = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.ems.RDS"))


#NAA oms: ems = 1-20
df.oms = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.oms.RDS"))

aic_res = lapply(1:NROW(df.oms), function(y){
  print(y)
  res = sapply(1:100, function(x){
    aic <- sapply(1:NROW(df.ems), function(z) {
      fit = try(readRDS(file.path(here::here(),"Ecov_study","mortality", "results", paste0("om", y), paste0("sim", x, "_em", z, ".RDS"))), silent = TRUE)
      out <- NA
      if(!is.character(fit)) if(!is.null(fit$fit$opt)) out = 2*(fit$fit$opt$obj + length(fit$fit$opt$par))
      return(out)
    })
    return(aic)
  })
})

fit = try(readRDS(file.path(here::here(),"Ecov_study","mortality", "results", paste0("om", 1), paste0("sim", 1, "_em", 1, ".RDS"))), silent = TRUE)


saveRDS(aic_res, file = file.path(here(),"Ecov_study","mortality", "results", "aic_results.RDS"))
aic_res <- readRDS(file.path(here(),"Ecov_study","mortality", "results", "aic_results.RDS"))

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

for(i in 1:3) {
  
  re_mod <- c("rec", "rec+1", "rec+M")[i]
  #EM:  M fixed, OM and EM RE assumption match
  em_ind <- which(df.ems$M_est == FALSE & df.ems$re_config == re_mod)
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
remove(all_res, all_res_mod)
#bias of ecov beta
all_beta_bias = lapply(1:NROW(df.oms), function(y){
  res = sapply(1:20, function(x) {
    bias <- sapply(which(df.ems$Ecov_est), function(z) { #only the estimating models that estimate Ecov_beta
      fit = try(readRDS(file.path(here::here(),"Ecov_study","mortality", "results", paste0("om", y), paste0("sim", x, "_em", z, ".RDS"))), silent = TRUE)
      out <- NA
      if(!is.character(fit)) if(!is.null(fit$fit$opt)) out = fit$fit$opt$par["Ecov_beta"] - df.oms$Ecov_effect[y]
      return(out)
    })
    return(bias)
  })
  return(res)
})
saveRDS(all_beta_bias, file = file.path(here::here(),"Ecov_study","mortality", "results", "ecov_beta_bias_results.RDS"))
all_beta_bias <- readRDS(file.path(here::here(),"Ecov_study","mortality", "results", "ecov_beta_bias_results.RDS"))

for(i in 1:3) {
  re_mod <- c("rec", "rec+1", "rec+M")[i]
  #EM:  M fixed, Ecov_beta estimated, OM and EM RE assumption match
  df.ems. <- df.ems[df.ems$Ecov_est,]
  em_ind <- which(!df.ems.$M_est & df.ems.$re_config == re_mod)
  om_ind <- which(df.oms$NAA_M_re == re_mod) #om and em match
  res <- t(sapply(all_beta_bias[om_ind], function(x) {
    x <- x[em_ind, ,drop = F]
    out <- c(apply(x,1,mean,na.rm=T), apply(x,1,sd, na.rm=T)/sqrt(apply(x, 1, function(y) sum(!is.na(y)))))
    out <- c(out, out[1] + qt(0.025, apply(x, 1, function(y) sum(!is.na(y)))) * out[2])
    out <- c(out, out[1] + qt(0.975, apply(x, 1, function(y) sum(!is.na(y)))) * out[2])
    return(out)
  }))
  colnames(res) <- c("bias_est", "bias_se", "bias_ci_lo", "bias_ci_hi")
  res <- cbind.data.frame(df.oms[om_ind,], res)
  print(head(res))
  if(i == 1) {
    all_res <- res
  } else {
    all_res <- rbind.data.frame(all_res, res)
  }
}
#all_res$bias_z = all_res$bias_est/all_res$bias_se
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

plt <- ggplot(all_res_mod, aes(x = Ecov_effect, y = bias_est, colour = Ecov_re_sig:Ecov_re_cor)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig:obs_error ~ NAA_M_re:Fhist, labeller = label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(-5, 5)) + ylab("Bias of Ecov effect (beta)") + xlab("Ecov effect size") +
    ggtitle("EM: M fixed") + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = bias_ci_lo, ymax = bias_ci_hi), width = .05, position = position_dodge(0.1)) + 
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "red")
plt

ggsave(here("Ecov_study","mortality", "paper", "Ecov_beta_bias_M_fixed.png"), plt)
remove(all_res, all_res_mod)

for(i in 1:3) {
  re_mod <- c("rec", "rec+1", "rec+M")[i]
  #EM:  M est, Ecov_beta estimated, OM and EM RE assumption match
  df.ems. <- df.ems[df.ems$Ecov_est,]
  em_ind <- which(df.ems.$M_est & df.ems.$re_config == re_mod)
  om_ind <- which(df.oms$NAA_M_re == re_mod) #om and em match
  res <- t(sapply(all_beta_bias[om_ind], function(x) {
    x <- x[em_ind, ,drop = F]
    out <- c(apply(x,1,mean,na.rm=T), apply(x,1,sd, na.rm=T)/sqrt(apply(x, 1, function(y) sum(!is.na(y)))))
    out <- c(out, out[1] + qt(0.025, apply(x, 1, function(y) sum(!is.na(y)))) * out[2])
    out <- c(out, out[1] + qt(0.975, apply(x, 1, function(y) sum(!is.na(y)))) * out[2])
    return(out)
  }))
  colnames(res) <- c("bias_est", "bias_se", "bias_ci_lo", "bias_ci_hi")
  res <- cbind.data.frame(df.oms[om_ind,], res)
  if(i == 1) {
    all_res <- res
  } else {
    all_res <- rbind.data.frame(all_res, res)
  }
}
#all_res$bias_z = all_res$bias_est/all_res$bias_se
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

plt <- ggplot(all_res_mod, aes(x = Ecov_effect, y = bias_est, colour = Ecov_re_sig:Ecov_re_cor)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig:obs_error ~ NAA_M_re:Fhist, labeller = label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(-5, 5)) + ylab("Bias of Ecov effect (beta)") + xlab("Ecov effect size") +
    ggtitle("EM: M fixed") + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = bias_ci_lo, ymax = bias_ci_hi), width = .05, position = position_dodge(0.1)) + 
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "red")
plt

ggsave(here("Ecov_study","mortality", "paper", "Ecov_beta_bias_M_estimated.png"), plt)

#temp = try(readRDS(file.path(here::here(),"Ecov_study","mortality", "results", paste0("om", 1), paste0("sim", 1, "_em", 1, ".RDS"))), silent = TRUE)

ssb_res = lapply(1:NROW(df.oms), function(y){
  print(y)
  res <- array(NA, dim = c(40, 12, 40, 2)) #nsims, nems, nyears, 2=(true, est)
  for(i in 1:40)for(j in 1:NROW(df.ems))  {
    fit = try(readRDS(file.path(here::here(),"Ecov_study","mortality", "results", paste0("om", y), paste0("sim", i, "_em", j, ".RDS"))), silent = TRUE)
    if(!is.character(fit)) if(!is.null(fit$fit$opt)) {
      res[i,j,,1] <-fit$truth$SSB
      res[i,j,,2] <- fit$fit$rep$SSB
    }
  }
  return(res)
})
################################################################
#depletion (SSB terminal/SSB in year 1)
depletion_error <- data.frame(om = integer(), em = integer(), sim = integer(), depletion = numeric())
for(h in 1:length(ssb_res)) for(i in 1:40)for(j in 1:NROW(df.ems))  {
  true <- ssb_res[[h]][i,j,40,1]/ssb_res[[h]][i,j,1,1] #true depletion
  est <- ssb_res[[h]][i,j,40,2]/ssb_res[[h]][i,j,1,2] #estimated depletion
  depletion_error = rbind.data.frame(depletion_error, c(h, j, i, (est-true)/true))
}
colnames(depletion_error) <- c("om", "em", "sim", "depletion")

for(i in 1:3) {
  re_mod <- c("rec", "rec+1", "rec+M")[i]
  #EM:  M est, Ecov_beta estimated, OM and EM RE assumption match
  #df.ems. <- df.ems[df.ems$Ecov_est,]
  em_ind <- which(df.ems$Ecov_est & df.ems$M_est & df.ems$re_config == re_mod)
  om_ind <- which(df.oms$NAA_M_re == re_mod) #om and em match
  res <- t(sapply(om_ind, function(x) {
    x <- subset(depletion_error, om == x)
    sapply(em_ind, function(y){ #only 1 em_ind
      x <- subset(x, em == y)
      n <- sum(!is.na(x$depletion))
      out <- c(mean(x$depletion, na.rm = T), sd(x$depletion, na.rm =T)/sqrt(n))
      out <- c(out,  out[1] + qt(0.025, n) * out[2],
        out[1] + qt(0.975, n) * out[2])
      return(out)
    })
  }))
  colnames(res) <- c("bias_est", "bias_se", "bias_ci_lo", "bias_ci_hi")
  res <- cbind.data.frame(df.oms[om_ind,], res)
  if(i == 1) {
    all_res <- res
  } else {
    all_res <- rbind.data.frame(all_res, res)
  }
}
#all_res$bias_z = all_res$bias_est/all_res$bias_se
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

plt <- ggplot(all_res_mod, aes(x = Ecov_effect, y = bias_est, colour = Ecov_re_sig:Ecov_re_cor)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig:obs_error ~ NAA_M_re:Fhist, labeller = label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(-0.2, 0.2)) + ylab("Relative Bias of Depeletion") + xlab("Ecov effect size") +
    ggtitle("EM: M Estimated") + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = bias_ci_lo, ymax = bias_ci_hi), width = .05, position = position_dodge(0.1)) + 
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "red")
plt

which(df.ems$Ecov_est & df.ems$M_est & df.ems$re_config == "rec")
temp <- subset(depletion_error, om == 3 & em == 1)
temp <- subset(all_res, NAA_M_re == "rec+M" & Fhist == "H-MSY" & Ecov_obs_sig == 0.1 & obs_error == "H" & Ecov_effect == 0.25 & Ecov_re_sig == 0.5 & Ecov_re_cor == 0.5)

#########################################################

