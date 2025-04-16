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

make_df_fn <- function(M_est = TRUE){
  #all EM PE assumptions
  for(i in 1:3) {
    re_mod <- c("rec", "rec+1", "rec+M")[i]
    #EM:  M estimated
    em_ind <- which(df.ems$M_est == M_est) #all EM PE assumptions
    print(length(em_ind))
    om_ind <- which(df.oms$NAA_M_re == re_mod) #om PE
    print(length(om_ind))
    res <- aic_fn(aic_res, em_ind, om_ind) #n_ems x n_oms
    print(res[1,])
    print(res[,1])
    res <- cbind(df.oms[rep(om_ind, each = length(em_ind)),], df.ems[rep(em_ind, length(om_ind)),], n = c(res))
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
      "H-MSY" = "2.5*F[MSY] %->% F[MSY]",
      "MSY" = "F[MSY]"))
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
  
  all_res_mod$correct_PE <- all_res_mod$correct_Ecov_effect <- "No"
  all_res_mod$correct_type <- "No"
  all_res_mod$correct_PE[which(all_res_mod$correct %in% c(1,2))] <- "Yes"
  all_res_mod$correct_Ecov_effect[which(all_res_mod$correct %in% c(1,3))] <- "Yes"
  all_res_mod$correct_type[all_res_mod$correct %in% c(1,2)] <- "PE"
  all_res_mod$correct_type[all_res_mod$correct %in% c(1,3)] <- "Ecov_effect"
  
  all_res_mod$correct_type <- factor(all_res_mod$correct_type)
  all_res_mod$Ecov_effect <- factor(all_res_mod$Ecov_effect)
  all_res_mod$correct <- factor(all_res_mod$correct)
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

# all_res_mod <- make_df_fn(M_est = TRUE)
all_res_mod <- rbind(make_df_fn(M_est = TRUE), make_df_fn(M_est = FALSE))
temp <- subset(all_res_mod, obs_error == "Low obs error (indices, age comp)" & Fhist == "2.5*F[MSY] %->% F[MSY]")
x <- subset(temp, Ecov_obs_sig == "sigma[ecov] == 0.1" & M_est== TRUE & Ecov_re_sig == "sigma[Ecov] == 0.1" & Ecov_re_cor == "rho[Ecov] == 0" & NAA_M_re == "OM: R" & Ecov_effect == "0")

plt <- ggplot(all_res_mod, aes(x = Ecov_effect, y = n, fill = correct)) + scale_fill_viridis_d() + 
    geom_col(position = "fill") + 
    #geom_col(position = position_dodge(0.1)) + 
    facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ NAA_M_re+Fhist+obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
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
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("Proportion ranked best") + xlab(expression(beta[Ecov])) + labs(fill = bquote(atop("EM Assumptions",PE*":"*beta[Ecov]==0))) +
    ggtitle(bquote(beta[M]~estimated)) + theme(plot.title = element_text(hjust = 0.5))
plt
ggsave(here("Ecov_study","mortality", "paper", "proportion_best_AIC_M_estimated.png"), plt, height = 12, width = 20, units = "in")

all_res_mod <- make_df_fn(M_est = FALSE)
plt <- ggplot(all_res_mod, aes(x = Ecov_effect, y = n, fill = correct)) + scale_fill_viridis_d() + 
    geom_col(position = "fill") + 
    #geom_col(position = position_dodge(0.1)) + 
    facet_grid(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ NAA_M_re+Fhist+obs_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
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
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, obs_error = label_wrap_gen(width = 15), Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
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
