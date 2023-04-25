library(here)
library(ggplot2)
library(dplyr)
df.ems = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.ems.RDS"))
df.oms = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.oms.RDS"))
term_ssb_error <- readRDS(file = file.path(here(),"Ecov_study","mortality", "results", "ssb_term_bias_res.RDS"))
ssb_res <- readRDS(file = file.path(here(),"Ecov_study","mortality", "results", "ssb_results.RDS"))

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
    ggtitle("EM: M and Ecov effect Estimated") + theme(plot.title = element_text(hjust = 0.5)) + labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = bias_ci_lo, ymax = bias_ci_hi), width = .05, position = position_dodge(0.1)) + 
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "red")
plt

ggsave(here("Ecov_study","mortality", "paper", "depletion_bias_Ecov_M_estimated.png"), plt)
