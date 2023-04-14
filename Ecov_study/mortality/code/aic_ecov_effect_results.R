library(here)
library(ggplot2)
library(dplyr)
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

for(i in 1:3) {
  
  re_mod <- c("rec", "rec+1", "rec+M")[i]
  #EM:  M fixed, OM and EM RE assumption match
  em_ind <- which(df.ems$M_est == FALSE & df.ems$re_config == re_mod)
  om_ind <- which(df.oms$NAA_M_re == re_mod) #om and em match
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
