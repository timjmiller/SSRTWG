library(here)
library(ggplot2)
library(dplyr)
library(ggh4x)
library(tidyr)
library(patchwork)
library(scales)
library(ggpattern)
df.ems = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.ems.RDS"))
df.oms = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.oms.RDS"))
mean_M_bias <- readRDS(file.path(here::here(),"Ecov_study","mortality", "results", "mean_M_bias_results.RDS"))

custom_boxplot_stat <- function(x){#, n.yrs = 1, n.sim = 100) {
  x <- x[which(!is.na(x))]
  n <- length(x)
  bnds95 <- qbinom(c(0.025,0.975), n, 0.5)/n # 95% CI bounds for median (sims x years)
  bnds80 <- qbinom(c(0.1,0.9), n, 0.5)/n # 80% CI bounds for median (sims x years)
  r <- quantile(x, probs = c(bnds95[1], bnds80[1], 0.5, bnds80[2], bnds95[2]))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}


#all EM PE assumptions
plot_df_fn <- function(df.ems, df.oms, Ecov_est = FALSE) {
  for(h in 1:3) {
    re_mod <- c("rec", "rec+1", "rec+M")[h]
    #EM:  M fixed, mean_M estimated, OM and EM RE assumption match
    df.ems. <- df.ems[df.ems$M_est,]
    em_ind <- which(df.ems.$Ecov_est==Ecov_est) # all 3 EM re_mods  
    om_ind <- which(df.oms$NAA_M_re == re_mod)
    for(i in om_ind){
      out <- matrix(NA,length(em_ind), 7)
      for(j in em_ind){
        x <- mean_M_bias[[i]][[j]]
        x <- x[which(!is.na(x[,2])),1] #remove fits without standard error estimates
        out[which(em_ind==j),] <- c(median(x, na.rm = TRUE), sd(x, na.rm=T)/sqrt(sum(!is.na(x))), custom_boxplot_stat(x))#quantile(x, probs = c(0.025,0.975)))
      }
      colnames(out) <- c("bias_est", "bias_se", "ymin", "lower", "middle", "upper", "ymax")
      if(i == om_ind[1]) res <- cbind.data.frame(df.oms[rep(i,length(em_ind)),], df.ems.[em_ind,], out)
      else res <- rbind.data.frame(res,cbind.data.frame(df.oms[rep(i,length(em_ind)),], df.ems.[em_ind,], out))
    }
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
      "0.1" = "sigma[italic(e)] == 0.1",
      "0.5" = "sigma[italic(e)] == 0.5"
    ))
  all_res_mod <- all_res_mod %>%
    mutate(Ecov_re_sig = recode(Ecov_re_sig,
      "0.1" = "sigma[italic(E)] == 0.1",
      "0.5" = "sigma[italic(E)] == 0.5"
    ))
  all_res_mod <- all_res_mod %>%
    mutate(Ecov_re_cor = recode(Ecov_re_cor,
      "0" = "rho[italic(E)] == 0",
      "0.5" = "rho[italic(E)] == 0.5"
    ))
  all_res_mod <- all_res_mod %>%
    mutate(beta_Ecov = recode(as.character(Ecov_est),
      "TRUE" = 'beta[Ecov]*" estimated"',
      "FALSE" = 'beta[Ecov]==0'
    ))
  all_res_mod <- all_res_mod %>%
    mutate(obs_error = recode(obs_error,
      "L" = "Low observation error",
      "H" = "High observation error"
    ))
  all_res_mod <- all_res_mod %>%
    mutate(OM_process_error = recode(NAA_M_re,
      "rec" = "R OMs",
      "rec+1" = "R+S OMs",
      "rec+M" = "R+M OMs"
    ))
  all_res_mod <- all_res_mod %>%
    mutate(EM_process_error = recode(re_config,
      "rec" = "R",
      "rec+1" = "R+S",
      "rec+M" = "R+M"
    ))
  all_res_mod <- all_res_mod %>% mutate(Fhist = recode(Fhist,
      "H-MSY"  = "2.5*italic(F)[MSY] %->% italic(F)[MSY]",
      "MSY" = "italic(F)[MSY]"))
  return(all_res_mod)
}

theme_set(theme_bw())
theme_update(strip.text = element_text(size = rel(1.5)), strip.placement = "outside", strip.background = element_rect(), #fill = "transparent"), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)), legend.title.align=0.5)

all_res <- rbind(plot_df_fn(df.ems, df.oms, Ecov_est = FALSE),plot_df_fn(df.ems, df.oms, Ecov_est = TRUE))
temp <- subset(all_res, obs_error == "Low observation error" & Fhist == "2.5*italic(F)[MSY] %->% italic(F)[MSY]")
plt <- ggplot(temp, aes(x = Ecov_effect, y = bias_est, colour = EM_process_error)) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) + 
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    scale_x_continuous(breaks = c(0,0.25,0.5)) + 
    facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor  ~ OM_process_error + beta_Ecov, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_parsed, Ecov_obs_sig = label_parsed, beta_Ecov = label_parsed)) +
    coord_cartesian(ylim = c(-1, 1)) + ylab(bquote(Median~Error~beta[M])) + xlab(expression("True "*beta[Ecov])) +
    labs(colour = "EM process error") +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .05, position = position_dodge(0.1))
plt

cairo_pdf(here("Ecov_study","mortality","manuscript", "median_M_bias.pdf"), width = 30*2/3, height = 20*2/3)
print(plt)
dev.off()

