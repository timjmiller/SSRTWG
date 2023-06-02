library(here)
library(ggplot2)
library(dplyr)
df.ems = readRDS(here("Project_0", "inputs", "df.ems.RDS"))

convergence_fn <- function(em){
  out <- rep(NA,5) #convergence type  1, and 2
  if(!is.character(em)) {
    if(!is.null(em$fit$opt)) {
      out[1] <- 1 
      out[2] <- em$fit$opt$conv
    }
    if(!is.null(em$fit$sdrep)) out[3] <- as.integer(sum(sapply(em$fit$sdrep$SE_par, function(g) any(is.nan(g)))))
    if(!is.null(em$fit$final_gradient)) out[4] <- max(abs(em$fit$final_gradient))
    if(!is.null(em$fit$sdrep)) {
      maxs <- sapply(em$fit$sdrep$SE_par, function(g) ifelse(any(!is.na(g)), max(g,na.rm=TRUE), NA))
      out[5] <- ifelse(any(!is.na(maxs)), max(maxs, na.rm =TRUE), NA)
    }
  }
  return(out)
}
df.oms = readRDS(here("Project_0", "inputs", "df.oms.RDS"))
naa_conv_res = lapply(1:NROW(df.oms), function(y){
  res = lapply(1:100, function(x){
    print(paste0("om_", y, ", sim_",x))
    sim = readRDS(here("Project_0", "results", "naa_om", paste0("om_", y), paste0("sim_",x,".RDS")))
    conv_res = t(sapply(sim, convergence_fn))
    return(conv_res)
  })
})
saveRDS(naa_conv_res, file = here("Project_0", "results", "naa_om_convergence_results.RDS"))

df.oms = readRDS(here("Project_0", "inputs", "df.M.oms.RDS"))
M_conv_res = lapply(1:NROW(df.oms), function(y){
  res = lapply(1:100, function(x){
    print(paste0("om_", y, ", sim_",x))
    sim = readRDS(here("Project_0", "results", "M_om", paste0("om_", y), paste0("sim_",x,".RDS")))
    conv_res = t(sapply(sim, convergence_fn))
    return(conv_res)
  })
})
saveRDS(M_conv_res, file = here("Project_0", "results", "M_om_convergence_results.RDS"))

df.oms = readRDS(here("Project_0", "inputs", "df.Sel.oms.RDS"))
Sel_conv_res = lapply(1:NROW(df.oms), function(y){
  res = lapply(1:100, function(x){
    print(paste0("om_", y, ", sim_",x))
    sim = readRDS(here("Project_0", "results", "Sel_om", paste0("om_", y), paste0("sim_",x,".RDS")))
    conv_res = t(sapply(sim, convergence_fn))
    return(conv_res)
  })
})
saveRDS(Sel_conv_res, file = here("Project_0", "results", "Sel_om_convergence_results.RDS"))

df.oms = readRDS(here("Project_0", "inputs", "df.q.oms.RDS"))
q_conv_res = lapply(1:NROW(df.oms), function(y){
  res = lapply(1:100, function(x){
    print(paste0("om_", y, ", sim_",x))
    sim = readRDS(here("Project_0", "results", "q_om", paste0("om_", y), paste0("sim_",x,".RDS")))
    conv_res = t(sapply(sim, convergence_fn))
    return(conv_res)
  })
})
saveRDS(q_conv_res, file = here("Project_0", "results", "q_om_convergence_results.RDS"))

#convergence information:
#1: 0/1 model finished optimizing
#2: nlminb convergence flag
#3: max gradient value
#4: number of NaNs in SEs for parameters, 0 = good invertible hessian
#5: maximum non-NaN SE estimate

conv_fn <- function(all, est_ind, om_ind = NULL){
  #est_ind is 1 ROW of df.ems.
  if(!is.null(om_ind)) all <- all[om_ind]
  out <- sapply(all, function(res){
    tmp <- sapply(res, function(x) {
      return(x[est_ind,])
    })
    x <- t(tmp)
    nconv <- c(
      sum(!is.na(x[,1])), 
      sum(!is.na(x[,2]) & x[,2] == 0), 
      sum(!is.na(x[,4]) & x[,4] < 1e-6),
      sum(!is.na(x[,3]) & x[,3] == 0),
      sum(!is.na(x[,5]) & x[,3] == 0 & x[,5] < 100)
    )
    nconv <- c(nconv, NROW(x))
    return(nconv)
  })
  return(t(out))
}
naa_om_conv_res <- readRDS(here("Project_0", "results", "naa_om_convergence_results.RDS"))
ind <- which(df.ems$M_est == FALSE & df.ems$SR_model == 2)
ind <- ind[which(ind<21)]
x <- conv_fn(naa_om_conv_res, ind[1])

conv_res_plotting_fn <- function(conv_res, df.oms, df.ems, em_ind){

  om_ind <- 1:NROW(df.oms)
  print(om_ind)
  print(em_ind)
  #em_ind <- which(df.ems$M_est == M_est & df.ems$SR_model == ifelse(SR_est, 3,2))
  for(i in em_ind){
    conv <- conv_fn(conv_res, i, om_ind)
    colnames(conv) <- c(paste0("Type", 1:5, "_n_pass"), "n_sim")
    print(dim(conv))
    print(head(conv))
    res <- df.oms[om_ind,]
    for(k in names(df.ems)) res <- cbind(res, df.ems[[k]][i])
    names(res)[NCOL(df.oms)+1:NCOL(df.ems)] <- names(df.ems)
    res <- cbind(res, conv)
    print(head(res))
    if(i == em_ind[1]){
      all_res <- res
    } else{
      all_res <- rbind(all_res, res)

    }
  }
  print(dim(all_res))
  print(head(all_res))


  for(i in 1:5){
    all_res[[paste0("Type", i, "_p_pass")]] <- all_res[[paste0("Type", i, "_n_pass")]]/all_res$n_sim
    all_res[[paste0("Type", i, "_ci_lo")]] <- sapply(1:NROW(all_res), function(x) {
      binom.test(all_res[[paste0("Type", i, "_n_pass")]][x], all_res$n_sim[x], p = all_res[[paste0("Type", i, "_p_pass")]][x])$conf.int[1]
    })
    all_res[[paste0("Type", i, "_ci_hi")]]  <- sapply(1:NROW(all_res), function(x) {
      binom.test(all_res[[paste0("Type", i, "_n_pass")]][x], all_res$n_sim[x], p = all_res[[paste0("Type", i, "_p_pass")]][x])$conf.int[2]
    })
  }
  # facs <- c("Ecov_obs_sig", "Fhist", "Ecov_re_sig","Ecov_re_cor", "obs_error", "NAA_M_re", "re_config")
  # all_res[facs] <- lapply(all_res[facs], factor)
  # all_res_mod <- all_res %>%
  all_res <- all_res %>%
    mutate(obs_error = recode(obs_error,
      "L" = "Low obs error \n (indices, age comp)",
      "H" = "High obs error \n (indices, age comp)"
    ))
  all_res <- all_res %>%
    mutate(EM_process_error = recode(re_config,
      "rec" = "R",
      "rec+1" = "R+S",
      "M_re" = "R+M",
      "q_re" = "R+q",
      "sel_re" = "R+Sel"
    ))
  ind <- which(all_res$M_re_cor == "iid" | all_res$sel_re_cor == "iid" | all_res$q_re_cor == "iid")
  all_res$EM_process_error[ind] <- paste0(all_res$EM_process_error[ind], "(iid)")
  ind <- which(all_res$M_re_cor == "ar1_y" | all_res$sel_re_cor == "ar1_y" | all_res$q_re_cor == "ar1")
  all_res$EM_process_error[ind] <- paste0(all_res$EM_process_error[ind], "(ar1)")
  all_res <- all_res %>%
    mutate(Fhist = recode(Fhist,
      "H-MSY" = "F history: High->FMSY",
      "MSY" = "F history: FMSY"
    ))
  all_res <- all_res %>% pivot_longer(
    cols = Type1_p_pass:Type5_ci_hi,
    names_to = c("Type", "locate"),
    names_pattern = "Type(.)_(.*)" #magic
  )

  all_res <- all_res %>% pivot_wider(names_from = locate, values_from = value) %>% as.data.frame
  all_res$Type <- as.numeric(all_res$Type)
  if(!is.null(all_res$NAA_sig)) all_res$NAA_sig[which(is.na(all_res$NAA_sig))] <- 0
  return(all_res)
}

df.oms = readRDS(here("Project_0", "inputs", "df.oms.RDS"))
df.ems = readRDS(here("Project_0", "inputs", "df.ems.RDS"))
naa_om_conv_res <- readRDS(here("Project_0", "results", "naa_om_convergence_results.RDS"))

dummy <- data.frame(vals = 1:5, def = c(
"finished optimizing",
"nlminb convergence flag = 0",
"max gradient value < 1e-6",
"0 NaNs in SEs for parameters",
"0 NANs and max SE estimate < 100"
))
ind <- which(df.ems$M_est == FALSE & df.ems$SR_model == 2)
ind <- ind[which(ind<21)]
all_res_mod <- conv_res_plotting_fn(naa_om_conv_res, df.oms, df.ems, em_ind = ind)
head(all_res_mod)
plt <- ggplot(all_res_mod, aes(x = Type, y = p_pass, colour = EM_process_error)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 2) + 
    facet_grid(obs_error+ Fhist ~ R_sig + NAA_sig, 
      #labeller = label_wrap_gen(width=15)) + #, labeller = label_parsed) + 
      labeller = labeller(R_sig = label_both, NAA_sig = label_both, obs_error = label_wrap_gen(width=30))) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(convergence)") + xlab("Convergence Type") +
    ggtitle("EM: No SR, M = 0.2") + theme(plot.title = element_text(hjust = 0.5)) + #labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .05, position = position_dodge(0.1)) 
plt
ggsave(here("Project_0","paper", "naa_om_p_convergence_meanR_M_fixed.png"), plt, height = 8, width = 12, units = "in")

ind <- which(df.ems$M_est == TRUE & df.ems$SR_model == 2)
ind <- ind[which(ind<21)]
all_res_mod <- conv_res_plotting_fn(naa_om_conv_res, df.oms, df.ems, em_ind = ind)
plt <- ggplot(all_res_mod, aes(x = Type, y = p_pass, colour = EM_process_error)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 2) + 
    facet_grid(obs_error+ Fhist ~ R_sig + NAA_sig, 
      #labeller = label_wrap_gen(width=15)) + #, labeller = label_parsed) + 
      labeller = labeller(R_sig = label_both, NAA_sig = label_both, obs_error = label_wrap_gen(width=30))) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(convergence)") + xlab("Convergence Type") +
    ggtitle("EM: No SR, M estimated") + theme(plot.title = element_text(hjust = 0.5)) + #labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .05, position = position_dodge(0.1)) #+
plt
ggsave(here("Project_0","paper", "naa_om_p_convergence_meanR_M_estimated.png"), plt, height = 8, width = 12, units = "in")

ind <- which(df.ems$M_est == FALSE & df.ems$SR_model == 3)
ind <- ind[which(ind<21)]
all_res_mod <- conv_res_plotting_fn(naa_om_conv_res, df.oms, df.ems, em_ind = ind)
head(all_res_mod)
plt <- ggplot(all_res_mod, aes(x = Type, y = p_pass, colour = EM_process_error)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 2) + 
    facet_grid(obs_error+ Fhist ~ R_sig + NAA_sig, 
      #labeller = label_wrap_gen(width=15)) + #, labeller = label_parsed) + 
      labeller = labeller(R_sig = label_both, NAA_sig = label_both, obs_error = label_wrap_gen(width=30))) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(convergence)") + xlab("Convergence Type") +
    ggtitle("EM: B-H, M = 0.2") + theme(plot.title = element_text(hjust = 0.5)) + #labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .05, position = position_dodge(0.1)) 
plt
ggsave(here("Project_0","paper", "naa_om_p_convergence_BH_M_fixed.png"), plt, height = 8, width = 12, units = "in")

ind <- which(df.ems$M_est == TRUE & df.ems$SR_model == 3)
ind <- ind[which(ind<21)]
all_res_mod <- conv_res_plotting_fn(naa_om_conv_res, df.oms, df.ems, em_ind = ind)
plt <- ggplot(all_res_mod, aes(x = Type, y = p_pass, colour = EM_process_error)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 2) + 
    facet_grid(obs_error+ Fhist ~ R_sig + NAA_sig, 
      #labeller = label_wrap_gen(width=15)) + #, labeller = label_parsed) + 
      labeller = labeller(R_sig = label_both, NAA_sig = label_both, obs_error = label_wrap_gen(width=30))) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(convergence)") + xlab("Convergence Type") +
    ggtitle("EM: B-H, M estimated") + theme(plot.title = element_text(hjust = 0.5)) + #labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .05, position = position_dodge(0.1)) #+
plt
ggsave(here("Project_0","paper", "naa_om_p_convergence_BH_M_estimated.png"), plt, height = 8, width = 12, units = "in")

df.oms = readRDS(here("Project_0", "inputs", "df.M.oms.RDS"))
df.ems = readRDS(here("Project_0", "inputs", "df.ems.RDS"))[5:24,]
M_om_conv_res <- readRDS(here("Project_0", "results", "M_om_convergence_results.RDS"))
ind <- which(df.ems$M_est == FALSE & df.ems$SR_model == 2)
all_res_mod <- conv_res_plotting_fn(M_om_conv_res, df.oms, df.ems, em_ind = ind)
plt <- ggplot(all_res_mod, aes(x = Type, y = p_pass, colour = EM_process_error)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 2) + 
    facet_grid(obs_error+ Fhist ~ M_sig + M_cor, 
      #labeller = label_wrap_gen(width=15)) + #, labeller = label_parsed) + 
      labeller = labeller(M_sig = label_both, M_cor = label_both, obs_error = label_wrap_gen(width=30))) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(convergence)") + xlab("Convergence Type") +
    ggtitle("EM: No SR, M = 0.2") + theme(plot.title = element_text(hjust = 0.5)) + #labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .05, position = position_dodge(0.1)) 
plt
ggsave(here("Project_0","paper", "M_om_p_convergence_meanR_M_fixed.png"), plt, height = 8, width = 12, units = "in")

ind <- which(df.ems$M_est == TRUE & df.ems$SR_model == 2)
all_res_mod <- conv_res_plotting_fn(M_om_conv_res, df.oms, df.ems, em_ind = ind)
plt <- ggplot(all_res_mod, aes(x = Type, y = p_pass, colour = EM_process_error)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 2) + 
    facet_grid(obs_error+ Fhist ~ M_sig + M_cor, 
      #labeller = label_wrap_gen(width=15)) + #, labeller = label_parsed) + 
      labeller = labeller(M_sig = label_both, M_cor = label_both, obs_error = label_wrap_gen(width=30))) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(convergence)") + xlab("Convergence Type") +
    ggtitle("EM: No SR, M estimated") + theme(plot.title = element_text(hjust = 0.5)) + #labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .05, position = position_dodge(0.1)) 
plt
ggsave(here("Project_0","paper", "M_om_p_convergence_meanR_M_estimated.png"), plt, height = 8, width = 12, units = "in")

ind <- which(df.ems$M_est == FALSE & df.ems$SR_model == 3)
all_res_mod <- conv_res_plotting_fn(M_om_conv_res, df.oms, df.ems, em_ind = ind)
plt <- ggplot(all_res_mod, aes(x = Type, y = p_pass, colour = EM_process_error)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 2) + 
    facet_grid(obs_error+ Fhist ~ M_sig + M_cor, 
      #labeller = label_wrap_gen(width=15)) + #, labeller = label_parsed) + 
      labeller = labeller(M_sig = label_both, M_cor = label_both, obs_error = label_wrap_gen(width=30))) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(convergence)") + xlab("Convergence Type") +
    ggtitle("EM: B-H, M = 0.2") + theme(plot.title = element_text(hjust = 0.5)) + #labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .05, position = position_dodge(0.1)) 
plt
ggsave(here("Project_0","paper", "M_om_p_convergence_BH_M_fixed.png"), plt, height = 8, width = 12, units = "in")

ind <- which(df.ems$M_est == TRUE & df.ems$SR_model == 3)
all_res_mod <- conv_res_plotting_fn(M_om_conv_res, df.oms, df.ems, em_ind = ind)
plt <- ggplot(all_res_mod, aes(x = Type, y = p_pass, colour = EM_process_error)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 2) + 
    facet_grid(obs_error+ Fhist ~ M_sig + M_cor, 
      labeller = labeller(M_sig = label_both, M_cor = label_both, obs_error = label_wrap_gen(width=30))) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(convergence)") + xlab("Convergence Type") +
    ggtitle("EM: B-H, M estimated") + theme(plot.title = element_text(hjust = 0.5)) + #labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .05, position = position_dodge(0.1)) 
plt
ggsave(here("Project_0","paper", "M_om_p_convergence_BH_M_estimated.png"), plt, height = 8, width = 12, units = "in")


df.oms = readRDS(here("Project_0", "inputs", "df.Sel.oms.RDS"))
df.ems = readRDS(here("Project_0", "inputs", "df.ems.RDS"))[c(5:20, 25:28),]
Sel_om_conv_res <- readRDS(here("Project_0", "results", "Sel_om_convergence_results.RDS"))
ind <- which(df.ems$M_est == FALSE & df.ems$SR_model == 2)
all_res_mod <- conv_res_plotting_fn(Sel_om_conv_res, df.oms, df.ems, em_ind = ind)
plt <- ggplot(all_res_mod, aes(x = Type, y = p_pass, colour = EM_process_error)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 2) + 
    facet_grid(obs_error+ Fhist ~ Sel_sig + Sel_cor, 
      labeller = labeller(Sel_sig = label_both, Sel_cor = label_both, obs_error = label_wrap_gen(width=30))) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(convergence)") + xlab("Convergence Type") +
    ggtitle("EM: No SR, M = 0.2") + theme(plot.title = element_text(hjust = 0.5)) + #labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .05, position = position_dodge(0.1)) 
plt
ggsave(here("Project_0","paper", "Sel_om_p_convergence_meanR_M_fixed.png"), plt, height = 8, width = 12, units = "in")

ind <- which(df.ems$M_est == TRUE & df.ems$SR_model == 2)
all_res_mod <- conv_res_plotting_fn(Sel_om_conv_res, df.oms, df.ems, em_ind = ind)
plt <- ggplot(all_res_mod, aes(x = Type, y = p_pass, colour = EM_process_error)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 2) + 
    facet_grid(obs_error+ Fhist ~ Sel_sig + Sel_cor, 
      labeller = labeller(Sel_sig = label_both, Sel_cor = label_both, obs_error = label_wrap_gen(width=30))) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(convergence)") + xlab("Convergence Type") +
    ggtitle("EM: No SR, M estimated") + theme(plot.title = element_text(hjust = 0.5)) + #labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .05, position = position_dodge(0.1)) 
plt
ggsave(here("Project_0","paper", "Sel_om_p_convergence_meanR_M_estimated.png"), plt, height = 8, width = 12, units = "in")

ind <- which(df.ems$M_est == FALSE & df.ems$SR_model == 3)
all_res_mod <- conv_res_plotting_fn(Sel_om_conv_res, df.oms, df.ems, em_ind = ind)
plt <- ggplot(all_res_mod, aes(x = Type, y = p_pass, colour = EM_process_error)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 2) + 
    facet_grid(obs_error+ Fhist ~ Sel_sig + Sel_cor, 
      labeller = labeller(Sel_sig = label_both, Sel_cor = label_both, obs_error = label_wrap_gen(width=30))) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(convergence)") + xlab("Convergence Type") +
    ggtitle("EM: B-H, M = 0.2") + theme(plot.title = element_text(hjust = 0.5)) + #labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .05, position = position_dodge(0.1)) 
plt
ggsave(here("Project_0","paper", "Sel_om_p_convergence_BH_M_fixed.png"), plt, height = 8, width = 12, units = "in")

ind <- which(df.ems$M_est == TRUE & df.ems$SR_model == 3)
all_res_mod <- conv_res_plotting_fn(Sel_om_conv_res, df.oms, df.ems, em_ind = ind)
plt <- ggplot(all_res_mod, aes(x = Type, y = p_pass, colour = EM_process_error)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 2) + 
    facet_grid(obs_error+ Fhist ~ Sel_sig + Sel_cor, 
      labeller = labeller(Sel_sig = label_both, Sel_cor = label_both, obs_error = label_wrap_gen(width=30))) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(convergence)") + xlab("Convergence Type") +
    ggtitle("EM: B-H, M estimated") + theme(plot.title = element_text(hjust = 0.5)) + #labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .05, position = position_dodge(0.1)) 
plt
ggsave(here("Project_0","paper", "Sel_om_p_convergence_BH_M_estimated.png"), plt, height = 8, width = 12, units = "in")


df.oms = readRDS(here("Project_0", "inputs", "df.q.oms.RDS"))
df.ems = readRDS(here("Project_0", "inputs", "df.ems.RDS"))[c(5:20, 29:32),]
q_om_conv_res <- readRDS(here("Project_0", "results", "q_om_convergence_results.RDS"))
ind <- which(df.ems$M_est == FALSE & df.ems$SR_model == 2)
all_res_mod <- conv_res_plotting_fn(q_om_conv_res, df.oms, df.ems, em_ind = ind)
plt <- ggplot(all_res_mod, aes(x = Type, y = p_pass, colour = EM_process_error)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 2) + 
    facet_grid(obs_error+ Fhist ~ q_sig + q_cor, 
      labeller = labeller(q_sig = label_both, q_cor = label_both, obs_error = label_wrap_gen(width=30))) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(convergence)") + xlab("Convergence Type") +
    ggtitle("EM: No SR, M = 0.2") + theme(plot.title = element_text(hjust = 0.5)) + #labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .05, position = position_dodge(0.1)) 
plt
ggsave(here("Project_0","paper", "q_om_p_convergence_meanR_M_fixed.png"), plt, height = 8, width = 12, units = "in")

ind <- which(df.ems$M_est == TRUE & df.ems$SR_model == 2)
all_res_mod <- conv_res_plotting_fn(q_om_conv_res, df.oms, df.ems, em_ind = ind)
plt <- ggplot(all_res_mod, aes(x = Type, y = p_pass, colour = EM_process_error)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 2) + 
    facet_grid(obs_error+ Fhist ~ q_sig + q_cor, 
      labeller = labeller(q_sig = label_both, q_cor = label_both, obs_error = label_wrap_gen(width=30))) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(convergence)") + xlab("Convergence Type") +
    ggtitle("EM: No SR, M estimated") + theme(plot.title = element_text(hjust = 0.5)) + #labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .05, position = position_dodge(0.1)) 
plt
ggsave(here("Project_0","paper", "q_om_p_convergence_meanR_M_estimated.png"), plt, height = 8, width = 12, units = "in")

ind <- which(df.ems$M_est == FALSE & df.ems$SR_model == 3)
all_res_mod <- conv_res_plotting_fn(q_om_conv_res, df.oms, df.ems, em_ind = ind)
plt <- ggplot(all_res_mod, aes(x = Type, y = p_pass, colour = EM_process_error)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 2) + 
    facet_grid(obs_error+ Fhist ~ q_sig + q_cor, 
      labeller = labeller(q_sig = label_both, q_cor = label_both, obs_error = label_wrap_gen(width=30))) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(convergence)") + xlab("Convergence Type") +
    ggtitle("EM: B-H, M = 0.2") + theme(plot.title = element_text(hjust = 0.5)) + #labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .05, position = position_dodge(0.1)) 
plt
ggsave(here("Project_0","paper", "q_om_p_convergence_BH_M_fixed.png"), plt, height = 8, width = 12, units = "in")

ind <- which(df.ems$M_est == TRUE & df.ems$SR_model == 3)
all_res_mod <- conv_res_plotting_fn(q_om_conv_res, df.oms, df.ems, em_ind = ind)
plt <- ggplot(all_res_mod, aes(x = Type, y = p_pass, colour = EM_process_error)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 2) + 
    facet_grid(obs_error+ Fhist ~ q_sig + q_cor, 
      labeller = labeller(q_sig = label_both, q_cor = label_both, obs_error = label_wrap_gen(width=30))) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("P(convergence)") + xlab("Convergence Type") +
    ggtitle("EM: B-H, M estimated") + theme(plot.title = element_text(hjust = 0.5)) + #labs(colour = "Ecov SD:Ecov Cor") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .05, position = position_dodge(0.1)) 
plt
ggsave(here("Project_0","paper", "q_om_p_convergence_BH_M_estimated.png"), plt, height = 8, width = 12, units = "in")
