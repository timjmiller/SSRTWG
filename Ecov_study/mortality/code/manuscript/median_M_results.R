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
median_M_bias <- readRDS(file.path(here::here(),"Ecov_study","mortality", "results", "mean_M_bias_results.RDS"))
conv_res <- readRDS(here("Ecov_study","mortality", "results", "convergence_results.RDS"))

custom_boxplot_stat <- function(x, alpha = 0.05){#, n.yrs = 1, n.sim = 100) {
  x <- x[which(!is.na(x))]
  n <- length(x)
  prob_ci <- qbinom(c(alpha/2,1-alpha/2), n, 0.5)/n # 95% CI bounds for median (sims x years)
  r <- quantile(x, probs = c(prob_ci[1], 0.5, prob_ci[2]))
  names(r) <- c("lower", "middle", "upper")
  r
}

#see convergence_results.R
conv_fn <- function(om, em, conv_res, Type = 3){
  x <- conv_res[[om]][[em]]
  if(Type == 1) ind <- which(!is.na(x[,1]))
  if(Type == 2) ind <- which(!is.na(x[,2]) & x[,2] == 0)
  if(Type == 3) ind <- which(!is.na(x[,3]) & x[,3] == 0)
  if(Type == 4) ind <- which(!is.na(x[,4]) & x[,4] < 1e-6)
  if(Type == 5) ind <- which(!is.na(x[,5]) & x[,3] == 0 & x[,5] < 10)
  return(ind)
}
# types of convergence
#convergence information:
#1: 0/1 model finished optimizing
#2: nlminb convergence flag
#3: number of NaNs in SEs for parameters, 0 = good invertible hessian
#4: max gradient value < 1e-6
#5: maximum non-NaN SE estimate < 10

#all EM PE assumptions
plot_df_fn <- function(df.ems, df.oms, Ecov_est = FALSE, conv_type = 3) {
  cnames <- c("bias_est", "bias_se", "rmse","lower", "middle", "upper", "true_beta_sd", "beta_sd_est", "beta_sd_rel_bias", "beta_sd_mean_bias", "beta_sd_median_bias", "pCI", "pCI_low", "pCI_hi", "nfit")
  for(h in 1:3) {
    re_mod <- c("rec", "rec+1", "rec+M")[h]
    #EM:  M fixed, mean_M estimated, OM and EM RE assumption match
    df.ems. <- df.ems[df.ems$M_est,]
    # em_ind <- which(df.ems.$Ecov_est==Ecov_est) # all 3 EM re_mods  
    em_ind_median_M_res <- which(df.ems.$Ecov_est == Ecov_est) # all 3 EM re_mods, M_results just has the 6 EMs where ecov_beta is estimated (this is NOT the first 6 EMS in df.ems)
    om_ind <- which(df.oms$NAA_M_re == re_mod)
    for(i in om_ind){
      out <- matrix(NA,length(em_ind_median_M_res), length(cnames))
      for(j in em_ind_median_M_res){
        em_ind_conv_res <- which(df.ems$M_est & df.ems$Ecov_est==Ecov_est & df.ems$re_config == df.ems.$re_config[j]) # conv_results has all 12 EMs, must find the one that corresponds to the em here
        # print(c(em_ind_conv_res,j))
        # stop()
        om_em_res <- median_M_bias[[i]][[j]]
        if(!is.null(conv_type)) om_em_res <- om_em_res[conv_fn(i,em_ind_conv_res,conv_res,Type = conv_type),,drop = FALSE] #make sure subset is consistent with convergence results
        rmse <- sqrt(mean(om_em_res[,1]^2, na.rm = TRUE))
        ci_pCI <- c(NA,NA)
        if(NROW(om_em_res)>0) ci_pCI <- binom.test(sum(om_em_res[,3], na.rm= TRUE), sum(!is.na(om_em_res[,3])), p = mean(om_em_res[,3], na.rm = TRUE))$conf.int[1:2]
        out[which(em_ind_median_M_res==j),] <- c(
          median(om_em_res[,1], na.rm = TRUE), 
          sd(om_em_res[,1], na.rm=T)/sqrt(sum(!is.na(om_em_res[,1]))), 
          rmse, 
          custom_boxplot_stat(om_em_res[,1]), 
          sd(om_em_res[,1], na.rm=T), #sd of the beta estimate is the same as the sd of the bias because the true value is constant.
          mean(om_em_res[,2],na.rm = TRUE),
          median(om_em_res[,2]/sd(om_em_res[,1],na.rm=T) -1, na.rm=T),
          mean(om_em_res[,2],na.rm = TRUE) - sd(om_em_res[,1], na.rm=T),
          median(om_em_res[,2],na.rm = TRUE) - sd(om_em_res[,1], na.rm=T),
          mean(om_em_res[,3], na.rm = TRUE),
          ci_pCI,
          sum(!is.na(om_em_res[,1])))
      }
      colnames(out) <- cnames
      if(i == om_ind[1]) res <- cbind.data.frame(df.oms[rep(i,length(em_ind_median_M_res)),], df.ems.[em_ind_median_M_res,], out)
      else res <- rbind.data.frame(res,cbind.data.frame(df.oms[rep(i,length(em_ind_median_M_res)),], df.ems.[em_ind_median_M_res,], out))
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
    mutate(obs_error = recode(obs_error,
      "L" = "Low observation error",
      "H" = "High observation error"
    ))
  all_res_mod <- all_res_mod %>%
    mutate(oe = recode(obs_error,
        "Low observation error" = "Low OE",
        "High observation error" = "High OE"
    ))
  all_res_mod$oe  <- factor(all_res_mod$oe, levels = c("Low OE", "High OE"))

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
  all_res_mod <- all_res_mod %>% 
    mutate(Fhist = recode(Fhist,
      "H-MSY"  = "2.5*italic(F)[MSY] %->% italic(F)[MSY]",
      "MSY" = "italic(F)[MSY]"))
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
      "TRUE" = 'beta[italic(E)]*" estimated"',
      "FALSE" = 'beta[italic(E)]==0'
    ))
  return(all_res_mod)
}

theme_set(theme_bw())
theme_update(strip.text = element_text(size = rel(1.5)), strip.placement = "outside", strip.background = element_rect(), #fill = "transparent"), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)), legend.title.align=0.5)

all_res <- rbind(plot_df_fn(df.ems, df.oms, Ecov_est = FALSE, conv_type = 3),plot_df_fn(df.ems, df.oms, Ecov_est = TRUE, conv_type = 3))
temp <- subset(all_res, obs_error == "Low observation error" & Fhist == "2.5*italic(F)[MSY] %->% italic(F)[MSY]")
plt <- ggplot(temp, aes(x = Ecov_effect, y = bias_est, colour = EM_process_error)) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) + 
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    scale_x_continuous(breaks = c(0,0.25,0.5)) + 
    facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor  ~ OM_process_error + beta_Ecov, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_parsed, Ecov_obs_sig = label_parsed, beta_Ecov = label_parsed)) +
    coord_cartesian(ylim = c(-0.5, 0.5)) + ylab(bquote(ME(hat(beta)[italic(M)]))) + xlab(expression("True "*beta[italic(E)])) +
    labs(colour = "EM process error") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = .05, position = position_dodge(0.1))
plt
cairo_pdf(here("Ecov_study","mortality","manuscript", "beta_M_bias_main.pdf"), width = 30*2/3, height = 20*2/3)
print(plt)
dev.off()

plt <- ggplot(temp, aes(x = Ecov_effect, y = beta_sd_median_bias, colour = EM_process_error)) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) + 
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    scale_x_continuous(breaks = c(0,0.25,0.5)) + 
    facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor  ~ OM_process_error + beta_Ecov, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_parsed, Ecov_obs_sig = label_parsed, beta_Ecov = label_parsed)) +
    coord_cartesian(ylim = c(-0.5, 0.5)) + ylab(bquote(ME(widehat(SE)(hat(beta)[italic(M)])))) + xlab(expression("True "*beta[italic(E)])) +
    labs(colour = "EM process error")
plt
cairo_pdf(here("Ecov_study","mortality","manuscript", "se_beta_M_bias_main.pdf"), width = 30*2/3, height = 20*2/3)
print(plt)
dev.off()

subset(df.oms, NAA_M_re == "rec+M" & Ecov_obs_sig == 0.1 & Ecov_re_sig == 0.5 & Ecov_re_cor == 0.5 & Fhist == "H-MSY" & obs_error == "L" & Ecov_effect == 0.5)$Model
# OM 69
which(df.ems$M_est & df.ems$Ecov_est & df.ems$re_config == "rec+M")
# conv_res EM 5
which(df.ems$Ecov_est[df.ems$M_est] & df.ems$re_config[df.ems$M_est] == "rec+M")
# median_M_res index 3

x <- median_M_bias[[69]][[3]][conv_fn(69,5,conv_res,Type = 3),,drop = FALSE]

mean(x[,1] < qnorm(0.975) * x[,2] & x[,1] > qnorm(0.025)*x[,2], na.rm = T)
x <- as.data.frame(x)
x[,2] <- x[,2] - mean(x[,2])
names(x) = c("est", "se", "inCI")
            plot(x[,1],x[,2])

theme_set(theme_bw())
theme_update(strip.text = element_text(size = rel(1.5)), strip.placement = "outside", strip.background = element_rect(), #fill = "transparent"), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)), legend.title.align=0.5)
plt <- ggplot(data = x, aes(x = est, y = se)) +
            geom_smooth(method = "lm", se=TRUE, color="black", formula = y ~ x) +
            geom_point() + ylab(bquote(widehat(SE)(hat(beta)[italic(M)])- SE(hat(beta)[italic(M)]))) + xlab(bquote(hat(beta)[italic(M)]-beta[italic(M)]))

plt
cairo_pdf(here("Ecov_study","mortality","manuscript", "om_69_em_5_beta_M_se_beta_M_lm_plot.pdf"), width = 30*2/3, height = 20*2/3)
print(plt)
dev.off()


plt <- ggplot(temp, aes(x = Ecov_effect, y = pCI, colour = EM_process_error)) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) + 
    geom_hline(aes(yintercept=0.95), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    scale_x_continuous(breaks = c(0,0.25,0.5)) + 
    facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor  ~ OM_process_error + beta_Ecov, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, Ecov_obs_sig = label_parsed, beta_Ecov = label_parsed)) +
    coord_cartesian(ylim = c(0.8, 1)) + ylab(bquote(CI~Coverage(hat(beta)[italic(M)]))) + xlab(expression("True "*beta[italic(E)])) +
    geom_errorbar(aes(ymin = pCI_low, ymax = pCI_hi), width = .05, position = position_dodge(0.1)) +
    labs(colour = "EM process error")
plt
cairo_pdf(here("Ecov_study","mortality","manuscript", "beta_M_CI_coverage_main.pdf"), width = 30*2/3, height = 20*2/3)
print(plt)
dev.off()

plt <- ggplot(temp, aes(x = Ecov_effect, y = rmse, colour = EM_process_error)) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor  ~ OM_process_error + beta_Ecov, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, Ecov_obs_sig = label_parsed, beta_Ecov = label_parsed)) +
    ylab(bquote(RMSE(hat(beta)[italic(M)]))) + xlab(expression("True "*beta[italic(E)])) + 
    coord_cartesian(ylim = c(0.01, 100), clip = "on") +
    scale_y_log10(breaks = c(1e-2,1e-1,1e-0,1e1,1e2), labels = parse(text = c("10^-2","10^-1","10^-0","10^1", "10^2"))) + #ylim(c(0,10)) +
    scale_x_continuous(breaks = c(0,0.25,0.5)) + labs(colour = "EM process error")
plt
cairo_pdf(here("Ecov_study","mortality","manuscript", "beta_M_rmse_main.pdf"), width = 30*2/3, height = 20*2/3)
print(plt)
dev.off()

OMs <- levels(all_res$OM_process_error)
OMs_lab <- c("Rom","RSom","RMom")
for(i in 1:length(OMs)){
  temp <- subset(all_res, OM_process_error == OMs[i])
  plt <- ggplot(temp, aes(x = Ecov_effect, y = bias_est, colour = EM_process_error)) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) + 
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    scale_x_continuous(breaks = c(0,0.25,0.5)) + 
    facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ beta_Ecov + Fhist + oe,
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, Ecov_obs_sig = label_parsed, Fhist = label_parsed, beta_Ecov = label_parsed)) +
    coord_cartesian(ylim = c(-0.5, 0.5)) + ylab(bquote(ME(hat(beta)[italic(M)]))) + xlab(expression("True "*beta[italic(E)])) +
    labs(colour = "EM process error") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = .05, position = position_dodge(0.1))
  plt
  cairo_pdf(here("Ecov_study","mortality","manuscript", paste0("beta_M_bias_",OMs_lab[i],".pdf")), width = 30*2/3, height = 20*2/3)
  print(plt)
  dev.off()

  plt <- ggplot(temp, aes(x = Ecov_effect, y = beta_sd_median_bias, colour = EM_process_error)) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) + 
      geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
      geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
      scale_x_continuous(breaks = c(0,0.25,0.5)) + 
      facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor  ~ beta_Ecov + Fhist + oe, 
        labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, Ecov_obs_sig = label_parsed, Fhist = label_parsed, beta_Ecov = label_parsed)) +
      coord_cartesian(ylim = c(-0.5, 0.5)) + ylab(bquote(ME(widehat(SE)(hat(beta)[italic(M)])))) + xlab(expression("True "*beta[italic(E)])) +
      labs(colour = "EM process error")
  plt
  cairo_pdf(here("Ecov_study","mortality","manuscript", paste0("se_beta_M_bias_",OMs_lab[i],".pdf")), width = 30*2/3, height = 20*2/3)
  print(plt)
  dev.off()

  
  plt <- ggplot(temp, aes(x = Ecov_effect, y = pCI, colour = EM_process_error)) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) + 
      geom_hline(aes(yintercept=0.95), linewidth = 2, linetype = "dashed", colour = "grey") +
      geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
      scale_x_continuous(breaks = c(0,0.25,0.5)) + 
      facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ beta_Ecov + Fhist + oe,
        labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, Ecov_obs_sig = label_parsed, beta_Ecov = label_parsed, Fhist = label_parsed)) +
      coord_cartesian(ylim = c(0.8, 1)) + ylab(bquote(CI~Coverage(hat(beta)[italic(M)]))) + xlab(expression("True "*beta[italic(E)])) +
      geom_errorbar(aes(ymin = pCI_low, ymax = pCI_hi), width = .05, position = position_dodge(0.1)) +
      labs(colour = "EM process error")
  plt
  cairo_pdf(here("Ecov_study","mortality","manuscript", paste0("beta_M_CI_coverage_", OMs_lab[i], ".pdf")), width = 30*2/3, height = 20*2/3)
  print(plt)
  dev.off()

  plt <- ggplot(temp, aes(x = Ecov_effect, y = rmse, colour = EM_process_error)) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) + 
      geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
      facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ beta_Ecov + Fhist + oe,
        labeller = labeller(Ecov_obs_sig = label_parsed, Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed,  Fhist = label_parsed, beta_Ecov = label_parsed)) +
      ylab(bquote(RMSE(hat(beta)[italic(M)]))) + xlab(expression("True "*beta[italic(E)])) + 
      scale_x_continuous(breaks = c(0,0.25,0.5)) + labs(colour = "EM process error") +  
      coord_cartesian(ylim = c(0.01, 100), clip = "on") +
      scale_y_log10(breaks = c(1e-2,1e-1,1e-0,1e1,1e2), labels = parse(text = c("10^-2","10^-1","10^-0","10^1", "10^2")))
  plt
  cairo_pdf(here("Ecov_study","mortality","manuscript", paste0("beta_M_rmse_",OMs_lab[i],".pdf")), width = 30*2/3, height = 20*2/3)
  print(plt)
  dev.off()
}
