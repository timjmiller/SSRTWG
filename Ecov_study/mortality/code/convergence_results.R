library(here)
library(ggplot2)
library(dplyr)
df.ems = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.ems.RDS"))
df.oms = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.oms.RDS"))

conv_res <- readRDS(file.path(here(),"Ecov_study","mortality", "results", "convergence_results.RDS"))
#convergence information:
#1: 0/1 model finished optimizing
#2: nlminb convergence flag
#3: number of NaNs in SEs for parameters, 0 = good invertible hessian
#4: max gradient value
#5: maximum non-NaN SE estimate

conv_fn <- function(all, est_ind, om_ind = NULL){
  if(!is.null(om_ind)) all <- all[om_ind]
  out <- sapply(all, function(res){
    tmp <- c(sapply(res[est_ind], function(x) {
      #print(est_ind)
      #print(length(res))
      #print(tmp)
      nconv <- c(
        sum(!is.na(x[,1])), 
        sum(!is.na(x[,2]) & x[,2] == 0), 
        sum(!is.na(x[,3]) & x[,3] == 0),
        sum(!is.na(x[,4]) & x[,4] < 1e-6),
        sum(!is.na(x[,5]) & x[,3] == 0 & x[,5] < 10)
      )
      nconv <- c(nconv, NROW(x))
      return(nconv)
    }))
    return(tmp)
  })  
  return(t(out))
}
#M estimated
conv_res_plotting_fn <- function(conv_res, M_est = TRUE, Ecov_est = TRUE){

  for(i in 1:3) for(j in 1:3){
    
    re_mods <- c("rec", "rec+1", "rec+M")
    #re_mod <- c("rec", "rec+1", "rec+M")[i]
    #EM:  M fixed, OM and EM RE assumption match
    om_ind <- which(df.oms$NAA_M_re == re_mods[i]) #om and em match
    em_ind <- which(df.ems$Ecov_est == Ecov_est & df.ems$M_est == M_est & df.ems$re_config == re_mods[j])
    conv <- conv_fn(conv_res, em_ind, om_ind)
    colnames(conv) <- c(paste0("Type", 1:5, "_n_pass"), "n_sim")
    res <- cbind(df.oms[om_ind,], df.ems[em_ind,], conv)
    if(i == 1 & j == 1) {
      all_res <- res
    } else {
      all_res <- rbind(all_res, res)
    }
  }
  for(i in 1:5){
    all_res[[paste0("Type", i, "_p_pass")]] <- all_res[[paste0("Type", i, "_n_pass")]]/all_res$n_sim
    all_res[[paste0("Type", i, "_ci_lo")]] <- sapply(1:NROW(all_res), function(x) {
      binom.test(all_res[[paste0("Type", i, "_n_pass")]][x], all_res$n_sim[x], p = all_res[[paste0("Type", i, "_p_pass")]][x])$conf.int[1]
    })
    all_res[[paste0("Type", i, "_ci_hi")]]  <- sapply(1:NROW(all_res), function(x) {
      binom.test(all_res[[paste0("Type", i, "_n_pass")]][x], all_res$n_sim[x], p = all_res[[paste0("Type", i, "_p_pass")]][x])$conf.int[2]
    })
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
      "L" = "Low obs error \n (indices, age comp)",
      "H" = "High obs error \n (indices, age comp)"
    ))
  all_res_mod <- all_res_mod %>%
    mutate(OM_process_error = recode(NAA_M_re,
      "rec" = "OM: R",
      "rec+1" = "OM: R+S",
      "rec+M" = "OM: R+M"
    ))
  all_res_mod <- all_res_mod %>%
    mutate(EM_process_error = recode(re_config,
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
}

all_res_mod <- conv_res_plotting_fn(conv_res, M_est = FALSE, Ecov_est = TRUE)
temp <- subset(all_res_mod, re_config == NAA_M_re)
plt <- ggplot(all_res_mod, aes(x = Ecov_effect, y = Type3_p_pass, colour = EM_process_error)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ OM_process_error+Fhist+obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(convergence: invertible Hessian)") + xlab(expression(beta[Ecov])) +
    labs(colour = "EM process error") + geom_errorbar(aes(ymin = Type3_ci_lo, ymax = Type3_ci_hi), width = .05, position = position_dodge(0.1)) + 
    geom_hline(aes(yintercept=0.5), linewidth = 1, linetype = "dashed", colour = "red") + ggtitle(bquote(beta[M]==log(0.2)*","~beta[Ecov]~estimated)) + theme(plot.title = element_text(hjust = 0.5))
plt
ggsave(here("Ecov_study","mortality", "paper", "proportion_good_hessian_ecov_effect_est_M_fixed.png"), plt, width = 20, height = 12, units = "in")

all_res_mod <- conv_res_plotting_fn(conv_res, M_est = TRUE, Ecov_est = TRUE)
plt <- ggplot(all_res_mod, aes(x = Ecov_effect, y = Type3_p_pass, colour = EM_process_error)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ OM_process_error+Fhist+obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(convergence: invertible Hessian)") + xlab(expression(beta[Ecov])) +
    labs(colour = "EM process error") + geom_errorbar(aes(ymin = Type3_ci_lo, ymax = Type3_ci_hi), width = .05, position = position_dodge(0.1)) + 
    geom_hline(aes(yintercept=0.5), linewidth = 1, linetype = "dashed", colour = "red") + ggtitle(bquote(beta[M]~estimated*","~beta[Ecov]~estimated)) + theme(plot.title = element_text(hjust = 0.5))
plt
ggsave(here("Ecov_study","mortality", "paper", "proportion_good_hessian_ecov_effect_est_M_est.png"), plt, width = 20, height = 12, units = "in")

