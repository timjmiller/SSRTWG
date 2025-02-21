library(here)
library(ggplot2)
library(ggh4x)
library(dplyr)
library(tidyr)
library(patchwork)
library(reshape2)

get_SR_M_estimates_fn <- function(om_type = "naa", M_or_SR = "SR"){
  df.ems <- readRDS(here::here("Project_0","inputs", "df.ems.RDS"))
  df.oms <- readRDS(here("Project_0","inputs", paste0("df.", ifelse(om_type == "naa", "", paste0(om_type,".")),"oms.RDS")))
  if(om_type == "naa") {
    em_ind <- 1:20
  }
  if(om_type == "M") {
    em_ind <- 5:24
  }
  if(om_type == "Sel") {
    em_ind <- c(5:20,25:28)
  }
  if(om_type == "q") {
    em_ind <- c(5:20,29:32)
  }
  df.ems <- df.ems[em_ind,]
  print(df.ems)
  if(M_or_SR == "SR") em_ind <- which(df.ems$SR_model == 3)
  if(M_or_SR == "M") em_ind <- which(df.ems$M_est)
  df.ems <- df.ems[em_ind,]
  print(df.ems)
  all_out <- lapply(1:NROW(df.oms), function(y){
    res <- lapply(1:100, function(x){
      print(paste0("om_", y, ", sim_",x))
      sim = readRDS(file.path(here::here(),"Project_0", "results", paste0(om_type,"_om"), paste0("om_", y), paste0("sim_",x,".RDS")))
      relres <- lapply(sim[em_ind],function(z) {
        if(M_or_SR == "SR"){
          out <- matrix(NA, 2,3)
          if(length(z)) if(length(z$fit)) {
            out[1,1] <- exp(z$fit$rep$log_SR_a[1])/exp(z$truth$log_SR_a[1])
            out[2,1] <- exp(z$fit$rep$log_SR_b[1])/exp(z$truth$log_SR_b[1])
          }
          if(!is.null(z$fit$sdrep)){
            out[1,2] <- z$fit$sdrep$SE_rep$log_SR_a[1] #se
            out[2,2] <- z$fit$sdrep$SE_rep$log_SR_b[1] #se
            ci <- matrix(c(z$fit$rep$log_SR_a[1],z$fit$rep$log_SR_b[1]) - c(z$truth$log_SR_a[1],z$truth$log_SR_b[1]),2,2) + 
              qnorm(0.975) * cbind(-out[,2],out[,2])
            out[,3] <- 0 >= ci[,1] & 0 <= ci[,2]
          }
        } else{
          out <- matrix(NA,1,3)
          if(length(z)) if(length(z$fit)) {
            out[1] <- exp(z$fit$rep$M_a[1])/exp(z$truth$M_a[1])
          }
          if(!is.null(z$fit$sdrep)){
            out[2] <- z$fit$sdrep$SE_par$M_a[1] #se
            ci <- z$fit$sdrep$Estimate_par$M_a[1] - z$truth$M_a[1] + qnorm(0.975) * c(-out[2],out[2])
            out[3] <- 0 >= ci[1] & 0 <= ci[2]
          }
        }
        return(out)
      })
      return(relres)
    })
    return(res)
  })
  return(all_out)
}

median_ci_fn <- function(x, alpha = 0.05){#, n.yrs = 1, n.sim = 100) {
  x <- x[which(!is.na(x))]
  n <- length(x)
  bnds <- qbinom(c(alpha/2,1-alpha/2), n, 0.5)/n # 95% CI bounds for median
  r <- quantile(x, probs = c(bnds[1], 0.5, bnds[2]))
  names(r) <- c("lo", "middle", "hi")
  r
}

make_plot_df <- function(om_type = "naa", res = naa_relSR_results, is_SE = FALSE, M_or_SR = "SR") {
  df.ems <- readRDS(here::here("Project_0","inputs", "df.ems.RDS"))
  if(om_type == "naa") {
    em_ind <- 1:20
  }
  if(om_type == "M") {
    em_ind <- 5:24
  }
  if(om_type == "Sel") {
    em_ind <- c(5:20,25:28)
  }
  if(om_type == "q") {
    em_ind <- c(5:20,29:32)
  }
  df.ems <- df.ems[em_ind,]
  if(M_or_SR == "SR") em_ind <- which(df.ems$SR_model == 3)
  if(M_or_SR == "M") em_ind <- which(df.ems$M_est)
  df.ems <- df.ems[em_ind,]

  res <- melt(res)
  names(res) <- c("par", "column", "value", "em", "sim","om") #em = 1: M fixed, em = 2: M estimated
  res <- cbind(df.ems[res$em,], res)
  res <- res %>% mutate(column = recode(column,
      "1" = "relerror",
      "2" = "cv",
      "3" = "in_ci"
    ))
  df <- res %>% pivot_wider(names_from = column, values_from = value) %>% as.data.frame
  if(is_SE) df <- filter(df, !is.na(cv))

  df$relerror = df$relerror - 1

  df <- df %>% group_by(om, em, par) %>%
    reframe(stats = median_ci_fn(relerror)) %>% as.data.frame
  df$type <- c("lo", "middle", "hi")
  if(om_type == "naa") {
    df.oms <- readRDS(here::here("Project_0", "inputs", "df.oms.RDS"))
  } else {
    df.oms <- readRDS(here::here("Project_0", "inputs", paste0("df.", om_type, ".oms.RDS")))
  }
  df <- df %>% pivot_wider(names_from = type, values_from = stats) %>% as.data.frame
  df <- cbind(df, df.oms[df$om,])
  df <- cbind(df, df.ems[df$em,])
  if(M_or_SR == "SR"){
    df <- df %>% mutate(par = recode(par,
      "1" = "italic(a)",
      "2" = "italic(b)"
    ))
    df <- df %>% mutate(M_config = if_else(M_est, "M estimated", "M known"))
  }
  if(M_or_SR == "M"){
    df <- df %>% mutate(par = recode(par, "1" = "italic(M)"))   
    df <- df %>% mutate(SR_model = recode(SR_model,
      "2" = "No SR",
      "3" = "SR"
    ))
  }
  df <- df %>%
    mutate(EM_process_error = recode(re_config,
      "rec" = "R",
      "rec+1" = "R+S",
      "M_re" = "R+M",
      "q_re" = "R+q",
      "sel_re" = "R+Sel"
    ))
  ind <- which(df$M_re_cor == "iid" | df$sel_re_cor == "iid" | df$q_re_cor == "iid")
  df$EM_process_error[ind] <- paste0(df$EM_process_error[ind], " (iid)")
  ind <- which(df$M_re_cor == "ar1_y" | df$sel_re_cor == "ar1_y" | df$q_re_cor == "ar1")
  df$EM_process_error[ind] <- paste0(df$EM_process_error[ind], " (AR1)")
  EM_process_error <- c("R","R+S","R+M (iid)","R+M (AR1)","R+Sel (iid)","R+Sel (AR1)","R+q (iid)","R+q (AR1)")
  df$EM_process_error <- factor(df$EM_process_error, levels = EM_process_error)

  df$correct_EM_PE <- "No"
  if(om_type == "naa") {
    df$NAA_sig[which(is.na(df$NAA_sig))] <- 0
    df$correct_EM_PE[df$NAA_sig == 0 & df$EM_process_error == "R"] <- "Yes"
    df$correct_EM_PE[df$NAA_sig >  0 & df$EM_process_error == "R+S"] <- "Yes"
    df <- df %>% mutate(NAA_sig = recode(NAA_sig,
        "0" = "sigma['2+'] == 0",
        "0.25" = "sigma['2+'] == 0.25",
        "0.5" = "sigma['2+'] == 0.5")) %>%
      mutate(R_sig = recode(R_sig,
        "0.5" = "sigma[italic(R)] == 0.5",
        "1.5" = "sigma[italic(R)] == 1.5"))
  }

  if(om_type == "M") {
    df$correct_EM_PE[df$M_cor == 0 & df$EM_process_error == "R+M (iid)"] <- "Yes"
    df$correct_EM_PE[df$M_cor >  0 & df$EM_process_error == "R+M (AR1)"] <- "Yes"
    df <- df %>% mutate(M_sig = recode(M_sig,
        "0.1" = "sigma[italic(M)] == 0.1",
        "0.5" = "sigma[italic(M)] == 0.5")) %>% 
      mutate(M_cor = recode(M_cor,
        "0" = "rho[italic(M)] == 0",
        "0.9" = "rho[italic(M)] == 0.9"))
  }
  if(om_type == "Sel") {
    df$correct_EM_PE[df$Sel_cor == 0 & df$EM_process_error == "R+Sel (iid)"] <- "Yes"
    df$correct_EM_PE[df$Sel_cor >  0 & df$EM_process_error == "R+Sel (AR1)"] <- "Yes"
    df <- df %>% mutate(Sel_sig = recode(Sel_sig,
        "0.1" = "sigma[Sel] == 0.1",
        "0.5" = "sigma[Sel] == 0.5")) %>% 
      mutate(Sel_cor = recode(Sel_cor,
        "0" = "rho[Sel] == 0",
        "0.9" = "rho[Sel] == 0.9"))
  }
  if(om_type == "q") {
    df$correct_EM_PE[df$q_cor == 0 & df$EM_process_error == "R+q (iid)"] <- "Yes"
    df$correct_EM_PE[df$q_cor >  0 & df$EM_process_error == "R+q (AR1)"] <- "Yes"
    df <- df %>% mutate(q_sig = recode(q_sig,
        "0.1" = "sigma[italic(q)] == 0.1",
        "0.5" = "sigma[italic(q)] == 0.5")) %>%
      mutate(q_cor = recode(q_cor,
        "0" = "rho[italic(q)] == 0",
        "0.9" = "rho[italic(q)] == 0.9"))  
  }
  df <- df %>%
    mutate(obs_error = recode(obs_error,
      "L" = "Low observation error",
      "H" = "High observation error"))
  df <- df %>% mutate(Fhist = recode(Fhist,
        "H-MSY"  = "2.5*italic(F)[MSY] %->% italic(F)[MSY]",
        "MSY" = "italic(F)[MSY]"))
  df <- df %>% as.data.frame
  facs <- c("om", "Fhist","NAA_sig", "R_sig", "M_sig", "M_cor", "Sel_sig", "Sel_cor", "q_sig", "q_cor", "obs_error", "correct_EM_PE")
  fac_names <- names(df)[names(df) %in% facs]  
  df[fac_names] <- lapply(df[fac_names], as.factor)
  df$correct_EM_PE[df$correct_EM_PE== "No"] <- NA
  return(df)
}
# temp <- make_plot_df(om_type = "M", res = M_rel_M_results, M_or_SR = "M")


df.ems <- readRDS(here::here("Project_0","inputs", "df.ems.RDS"))

make_results <- FALSE
if(make_results) {
  naa_rel_SR_results <- get_SR_M_estimates_fn("naa")
  saveRDS(naa_rel_SR_results, file = here("Project_0","results", "naa_relSR_results_all_EM_PE.RDS"))
  M_rel_SR_results <- get_SR_M_estimates_fn("M")
  saveRDS(M_rel_SR_results, file = here("Project_0","results", "M_relSR_results_all_EM_PE.RDS"))
  Sel_rel_SR_results <- get_SR_M_estimates_fn("Sel")
  saveRDS(Sel_rel_SR_results, file = here("Project_0","results", "Sel_relSR_results_all_EM_PE.RDS"))
  q_rel_SR_results <- get_SR_M_estimates_fn("q")
  saveRDS(q_rel_SR_results, file = here("Project_0","results", "q_relSR_results_all_EM_PE.RDS"))
}


naa_relSR_results <- readRDS(file = here("Project_0","results", "naa_relSR_results_all_EM_PE.RDS"))
M_relSR_results <- readRDS(file = here("Project_0","results", "M_relSR_results_all_EM_PE.RDS"))
Sel_relSR_results <- readRDS(file = here("Project_0","results", "Sel_relSR_results_all_EM_PE.RDS"))
q_relSR_results <- readRDS(file = here("Project_0","results", "q_relSR_results_all_EM_PE.RDS"))

# EM_process_error <- c("R","R+S","R+M (iid)","R+M (AR1)","R+Sel (iid)","R+Sel (AR1)","R+q (iid)","R+q (AR1)")
# R_S_df$EM_process_error <- factor(R_S_df$EM_process_error, levels = EM_process_error)


R_S_df <- make_plot_df(om_type = "naa", res = naa_relSR_results)
R_M_df <- make_plot_df(om_type = "M", res = M_relSR_results)
R_Sel_df <- make_plot_df(om_type = "Sel", res = Sel_relSR_results)
R_q_df <- make_plot_df(om_type = "q", res = q_relSR_results)


theme_set(theme_bw())
theme_update(strip.text = element_text(size = rel(2)), strip.placement = "outside", strip.background = element_rect(), #fill = "transparent"), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)))

naa_plt <- ggplot(R_S_df, aes(x = M_config, y = middle))  + 
    facet_nested(obs_error + Fhist ~ par + R_sig + NAA_sig, 
      labeller = labeller(Fhist = label_parsed, R_sig = label_parsed, NAA_sig = label_parsed, par = label_parsed)) +
    coord_cartesian(ylim = c(-1, 1)) + ylab("Median Relative error") + xlab("EM M assumption") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

naa_plt

M_plt <- ggplot(R_M_df, aes(x = M_config, y = middle))  + 
    facet_nested(obs_error + Fhist ~ par + M_sig + M_cor, 
      labeller = labeller(Fhist = label_parsed, M_sig = label_parsed, M_cor = label_parsed, par = label_parsed)) +
    coord_cartesian(ylim = c(-1, 1)) + ylab("Median Relative error") + xlab("EM M assumption") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
M_plt

Sel_plt <- ggplot(R_Sel_df, aes(x = M_config, y = middle))  + 
    facet_nested(obs_error + Fhist ~ par + Sel_sig + Sel_cor, 
      labeller = labeller(Fhist = label_parsed, Sel_sig = label_parsed, Sel_cor = label_parsed, par = label_parsed)) +
    coord_cartesian(ylim = c(-1, 1)) + ylab("Median Relative error") + xlab("EM M assumption") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
Sel_plt

q_plt <- ggplot(R_q_df, aes(x = M_config, y = middle))  + 
    facet_nested(obs_error + Fhist ~ par + q_sig + q_cor, 
      labeller = labeller(Fhist = label_parsed, q_sig = label_parsed, q_cor = label_parsed, par = label_parsed)) +
    coord_cartesian(ylim = c(-1, 1)) + ylab("Median Relative error") + xlab("EM M assumption") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
q_plt


Sel_plt_alt <- Sel_plt + 
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21, show.legend = TRUE) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) 

theme_set(theme_bw())
theme_update(strip.placement = "outside", strip.background = element_rect(), title = element_text(size = rel(1.5)))
cairo_pdf(here("Project_0","manuscript", "sr_bias_plots.pdf"), width = 30*2/3, height = 20*2/3)
design <- c(area(1,1,1,1), area(2,1,2,1), area(1,2,1,2), area(2,2,2,2))
(naa_plt +  xlab("") + labs(title = "A: R, R+S OMs") + theme(axis.title.x = element_blank(), strip.text.y=element_blank(), axis.text.x=element_blank(), legend.position="none")) +# , strip.text=element_text(margin=margin(b = 5)), plot.margin = margin(b = 1, t = 0))) + 
(M_plt + labs(title = "B: R+M OMs") + theme(legend.position="none", strip.text.y=element_blank())) + #axis.title.x = element_blank(), strip.text=element_text(margin=margin(b = 5)), axis.text.x=element_blank(), plot.margin = margin(b = 1, t = 0))) + 
(Sel_plt_alt + xlab("") + labs(title = "C: R+Sel OMs") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())) +# , strip.text=element_text(margin=margin(b = 5)), plot.margin = margin(b = 1, t = 0))) + 
(q_plt + labs(title = "D: R+q OMs") + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position = "none")) + 
  plot_layout(design = design, axis_titles = "collect")
dev.off()

#Natural Mortality

make_results <- FALSE
if(make_results) {
  naa_rel_M_results <- get_SR_M_estimates_fn("naa", M_or_SR = "M")
  M_rel_M_results <- get_SR_M_estimates_fn("M", M_or_SR = "M")
  Sel_rel_M_results <- get_SR_M_estimates_fn("Sel", M_or_SR = "M")
  q_rel_M_results <- get_SR_M_estimates_fn("q", M_or_SR = "M")
  saveRDS(naa_rel_M_results, file = here("Project_0","results", "naa_rel_M_results_all_EM_PE.RDS"))
  saveRDS(M_rel_M_results, file = here("Project_0","results", "M_rel_M_results_all_EM_PE.RDS"))
  saveRDS(Sel_rel_M_results, file = here("Project_0","results", "Sel_rel_M_results_all_EM_PE.RDS"))
  saveRDS(q_rel_M_results, file = here("Project_0","results", "q_rel_M_results_all_EM_PE.RDS"))
}

naa_rel_M_results <- readRDS(file = here("Project_0","results", "naa_rel_M_results_all_EM_PE.RDS"))
M_rel_M_results <- readRDS(file = here("Project_0","results", "M_rel_M_results_all_EM_PE.RDS"))
Sel_rel_M_results <- readRDS(file = here("Project_0","results", "Sel_rel_M_results_all_EM_PE.RDS"))
q_rel_M_results <- readRDS(file = here("Project_0","results", "q_rel_M_results_all_EM_PE.RDS"))

R_S_df <- make_plot_df(om_type = "naa", res = naa_rel_M_results, M_or_SR = "M")
R_M_df <- make_plot_df(om_type = "M", res = M_rel_M_results, M_or_SR = "M")
R_Sel_df <- make_plot_df(om_type = "Sel", res = Sel_rel_M_results, M_or_SR = "M")
R_q_df <- make_plot_df(om_type = "q", res = q_rel_M_results, M_or_SR = "M")
# R_S_df$EM_process_error <- factor(R_S_df$EM_process_error, levels = EM_process_error)

theme_set(theme_bw())
theme_update(strip.text = element_text(size = rel(2)), strip.placement = "outside", strip.background = element_rect(), #fill = "transparent"), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)))

naa_plt <- ggplot(R_S_df, aes(x = SR_model, y = middle))  + 
    facet_nested(obs_error + Fhist ~ R_sig + NAA_sig, 
      labeller = labeller(Fhist = label_parsed, R_sig = label_parsed, NAA_sig = label_parsed)) +
    coord_cartesian(ylim = c(-1, 1)) + ylab(expression(Median~relative~error~(italic(M)))) + xlab("EM SR assumption") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
naa_plt

M_plt <- ggplot(R_M_df, aes(x = SR_model, y = middle))  + 
    facet_nested(obs_error + Fhist ~ M_sig + M_cor, 
      labeller = labeller(Fhist = label_parsed, M_sig = label_parsed, M_cor = label_parsed)) +
    coord_cartesian(ylim = c(-1, 1)) + ylab(expression(Median~relative~error~(italic(M)))) + xlab("EM SR assumption") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
M_plt

Sel_plt <- ggplot(R_Sel_df, aes(x = SR_model, y = middle))  + 
    facet_nested(obs_error + Fhist ~ Sel_sig + Sel_cor, 
      labeller = labeller(Fhist = label_parsed, Sel_sig = label_parsed, Sel_cor = label_parsed)) +
    coord_cartesian(ylim = c(-1, 1)) + ylab(expression(Median~relative~error~(italic(M)))) + xlab("EM SR assumption") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
Sel_plt

q_plt <- ggplot(R_q_df, aes(x = SR_model, y = middle))  + 
    facet_nested(obs_error + Fhist ~ q_sig + q_cor, 
      labeller = labeller(Fhist = label_parsed, q_sig = label_parsed, q_cor = label_parsed)) +
    coord_cartesian(ylim = c(-1, 1)) + ylab(expression(Median~relative~error~(italic(M)))) + xlab("EM SR assumption") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
q_plt

Sel_plt_alt <- Sel_plt + 
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21, show.legend = TRUE) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) 

theme_set(theme_bw())
theme_update(strip.placement = "outside", strip.background = element_rect(), title = element_text(size = rel(1.5)))
cairo_pdf(here("Project_0","manuscript", "M_bias_plots.pdf"), width = 30*2/3, height = 20*2/3)
design <- c(area(1,1,1,1), area(2,1,2,1), area(1,2,1,2), area(2,2,2,2))
(naa_plt +  xlab("") + labs(title = "A: R, R+S OMs") + theme(axis.title.x = element_blank(), strip.text.y=element_blank(), axis.text.x=element_blank(), legend.position="none")) +
(M_plt + labs(title = "B: R+M OMs") + theme(legend.position="none", strip.text.y=element_blank())) + 
(Sel_plt_alt + xlab("") + labs(title = "C: R+Sel OMs") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())) +
(q_plt + labs(title = "D: R+q OMs") + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position = "none")) + 
  plot_layout(design = design, axis_titles = "collect")
dev.off()


########################################
#Get annual SSB, F relative error
########################################


make_annual_relerror_df <- function(om_type = "naa", res = all_naa_om_relssb, is_SE = FALSE) {
  df.ems <- readRDS(here::here("Project_0","inputs", "df.ems.RDS"))
  if(om_type == "naa") {
    em_ind <- 1:20
  }
  if(om_type == "M") {
    em_ind <- 5:24
  }
  if(om_type == "Sel") {
    em_ind <- c(5:20,25:28)
  }
  if(om_type == "q") {
    em_ind <- c(5:20,29:32)
  }
  df.ems <- df.ems[em_ind,]

  res <- melt(res)
  names(res) <- c("year", "column", "value", "em", "sim","om")
  res <- res %>% mutate(column = recode(column,
      "1" = "relerror",
      "2" = "cv",
      "3" = "in_ci",
    ))
  df <- res %>% tidyr::pivot_wider(names_from = column, values_from = value) %>% as.data.frame
  if(is_SE) df <- filter(df, !is.na(cv))
  
  df$relerror = df$relerror - 1

  df <- df %>% group_by(om, em, year) %>%
    reframe(stats = median_ci_fn(relerror)) %>% as.data.frame
  df$type <- c("lo", "middle", "hi")
  if(om_type == "naa") {
    df.oms <- readRDS(here::here("Project_0", "inputs", "df.oms.RDS"))
  } else {
    df.oms <- readRDS(here::here("Project_0", "inputs", paste0("df.", om_type, ".oms.RDS")))
  }
  df <- df %>% pivot_wider(names_from = type, values_from = stats) %>% as.data.frame
  df <- cbind(df, df.oms[df$om,])
  df <- cbind(df, df.ems[df$em,])
  df <- df %>%
    mutate(EM_process_error = recode(re_config,
      "rec" = "R",
      "rec+1" = "R+S",
      "M_re" = "R+M",
      "q_re" = "R+q",
      "sel_re" = "R+Sel"
    ))
  ind <- which(df$M_re_cor == "iid" | df$sel_re_cor == "iid" | df$q_re_cor == "iid")
  df$EM_process_error[ind] <- paste0(df$EM_process_error[ind], " (iid)")
  ind <- which(df$M_re_cor == "ar1_y" | df$sel_re_cor == "ar1_y" | df$q_re_cor == "ar1")
  df$EM_process_error[ind] <- paste0(df$EM_process_error[ind], " (AR1)")
  EM_process_error <- c("R","R+S","R+M (iid)","R+M (AR1)","R+Sel (iid)","R+Sel (AR1)","R+q (iid)","R+q (AR1)")
  df$EM_process_error <- factor(df$EM_process_error, levels = EM_process_error)
  
  df$correct_EM_PE <- "No"
  if(om_type == "naa") {
    df$NAA_sig[which(is.na(df$NAA_sig))] <- 0
    df$correct_EM_PE[df$NAA_sig == 0 & df$EM_process_error == "R"] <- "Yes"
    df$correct_EM_PE[df$NAA_sig >  0 & df$EM_process_error == "R+S"] <- "Yes"
    df <- df %>% mutate(NAA_sig = recode(NAA_sig,
        "0" = "sigma['2+'] == 0",
        "0.25" = "sigma['2+'] == 0.25",
        "0.5" = "sigma['2+'] == 0.5")) %>%
      mutate(R_sig = recode(R_sig,
        "0.5" = "sigma[italic(R)] == 0.5",
        "1.5" = "sigma[italic(R)] == 1.5"))
  }

  if(om_type == "M") {
    df$correct_EM_PE[df$M_cor == 0 & df$EM_process_error == "R+M (iid)"] <- "Yes"
    df$correct_EM_PE[df$M_cor >  0 & df$EM_process_error == "R+M (AR1)"] <- "Yes"
    df <- df %>% mutate(M_sig = recode(M_sig,
        "0.1" = "sigma[italic(M)] == 0.1",
        "0.5" = "sigma[italic(M)] == 0.5")) %>% 
      mutate(M_cor = recode(M_cor,
        "0" = "rho[italic(M)] == 0",
        "0.9" = "rho[italic(M)] == 0.9"))
  }
  if(om_type == "Sel") {
    df$correct_EM_PE[df$Sel_cor == 0 & df$EM_process_error == "R+Sel (iid)"] <- "Yes"
    df$correct_EM_PE[df$Sel_cor >  0 & df$EM_process_error == "R+Sel (AR1)"] <- "Yes"
    df <- df %>% mutate(Sel_sig = recode(Sel_sig,
        "0.1" = "sigma[Sel] == 0.1",
        "0.5" = "sigma[Sel] == 0.5")) %>% 
      mutate(Sel_cor = recode(Sel_cor,
        "0" = "rho[Sel] == 0",
        "0.9" = "rho[Sel] == 0.9"))
  }
  if(om_type == "q") {
    df$correct_EM_PE[df$q_cor == 0 & df$EM_process_error == "R+q (iid)"] <- "Yes"
    df$correct_EM_PE[df$q_cor >  0 & df$EM_process_error == "R+q (AR1)"] <- "Yes"
    df <- df %>% mutate(q_sig = recode(q_sig,
        "0.1" = "sigma[italic(q)] == 0.1",
        "0.5" = "sigma[italic(q)] == 0.5")) %>%
      mutate(q_cor = recode(q_cor,
        "0" = "rho[italic(q)] == 0",
        "0.9" = "rho[italic(q)] == 0.9"))  
  }
  df <- df %>%
    mutate(obs_error = recode(obs_error,
      "L" = "Low observation error",
      "H" = "High observation error"))
  df <- df %>% mutate(Fhist = recode(Fhist,
        "H-MSY"  = "2.5*italic(F)[MSY] %->% italic(F)[MSY]",
        "MSY" = "italic(F)[MSY]"))
  df <- df %>% as.data.frame
  facs <- c("om", "Fhist","NAA_sig", "R_sig", "M_sig", "M_cor", "Sel_sig", "Sel_cor", "q_sig", "q_cor", "obs_error", "correct_EM_PE")
  fac_names <- names(df)[names(df) %in% facs]
  df[fac_names] <- lapply(df[fac_names], as.factor)
  df$correct_EM_PE[df$correct_EM_PE== "No"] <- NA
  est_config <- c("M known, No SR", "M estimated, No SR", "M known, SR", "M estimated, SR")
  df$est_config <- est_config[1]
  df$est_config[df$M_est & df$SR_model == 2] <- est_config[2]
  df$est_config[!df$M_est & df$SR_model == 3] <- est_config[3]
  df$est_config[df$M_est & df$SR_model == 3] <- est_config[4]
  df$est_config <- factor(df$est_config, levels = est_config)
  return(df)
}



########################################
#terminal SSB
########################################

all_naa_om_relssb <- readRDS(file = here("Project_0","results", "all_naa_relssb_results.RDS"))

ylims <- c(-0.4,0.4)
plot.df <- make_annual_relerror_df()
plot.df$outside <- plot.df$hi < min(ylims) | plot.df$lo > max(ylims)

temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_S_df <- temp

all_M_relssb <- readRDS(file = here("Project_0","results", "all_M_relssb_results.RDS"))
plot.df <- make_annual_relerror_df("M", all_M_relssb)
plot.df$outside <- plot.df$hi < min(ylims) | plot.df$lo > max(ylims)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_M_df <- temp

all_Sel_relssb <- readRDS(file = here("Project_0","results", "all_Sel_relssb_results.RDS"))
plot.df <- make_annual_relerror_df("Sel", all_Sel_relssb)
plot.df$outside <- plot.df$hi < min(ylims) | plot.df$lo > max(ylims)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_Sel_df <- temp

all_q_relssb <- readRDS(file = here("Project_0","results", "all_q_relssb_results.RDS"))
plot.df <- make_annual_relerror_df("q", all_q_relssb)
plot.df$outside <- plot.df$hi < min(ylims) | plot.df$lo > max(ylims)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_q_df <- temp

theme_set(theme_bw())
theme_update(strip.text = element_text(size = rel(2)), strip.placement = "outside", strip.background = element_rect(), #fill = "transparent"), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)))

R_S_plt <- ggplot(filter(R_S_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ R_sig + NAA_sig, labeller = labeller(Fhist = label_parsed, R_sig = label_parsed, NAA_sig = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    coord_cartesian(ylim = ylims) + ylab("Median Relative Error (SSB)") + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_S_plt


R_M_plt <- ggplot(filter(R_M_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ M_sig + M_cor, labeller = labeller(Fhist = label_parsed, M_sig = label_parsed, M_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    coord_cartesian(ylim = ylims) + ylab("Median Relative Error (SSB)") + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_M_plt



R_Sel_plt <- ggplot(filter(R_Sel_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ Sel_sig + Sel_cor, labeller = labeller(Fhist = label_parsed, Sel_sig = label_parsed, Sel_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    coord_cartesian(ylim = ylims) + ylab("Median Relative Error (SSB)") + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_Sel_plt



R_q_plt <- ggplot(filter(R_q_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ q_sig + q_cor, labeller = labeller(Fhist = label_parsed, q_sig = label_parsed, q_cor = label_parsed)) +
    coord_cartesian(ylim = ylims) + ylab("Median Relative Error (SSB)") + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_q_plt


R_Sel_plt_alt <- R_Sel_plt + 
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21, show.legend = TRUE) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) 


theme_set(theme_bw())
theme_update(strip.placement = "outside", strip.background = element_rect(), title = element_text(size = rel(1.5)))
cairo_pdf(here("Project_0","manuscript", "term_SSB_bias_plots.pdf"), width = 30*2/3, height = 20*2/3)
design <- c(area(1,1,1,1), area(2,1,2,1), area(1,2,1,2), area(2,2,2,2))
(R_S_plt +  xlab("") + labs(title = "A: R, R+S OMs") + theme(axis.title.x = element_blank(), strip.text.y=element_blank(), axis.text.x=element_blank(), legend.position="none")) +
(R_M_plt + labs(title = "B: R+M OMs") + theme(legend.position="none", strip.text.y=element_blank())) + 
(R_Sel_plt_alt + xlab("") + labs(title = "C: R+Sel OMs") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())) +
(R_q_plt + labs(title = "D: R+q OMs") + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position = "none")) + 
  plot_layout(design = design, axis_titles = "collect")
dev.off()

########################################
#terminal F
########################################

all_naa_relF <- readRDS(file = here("Project_0","results", "all_naa_relF_results.RDS"))
plot.df <- make_annual_relerror_df("naa", all_naa_relF)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_S_df <- temp


all_M_relF <- readRDS(file = here("Project_0","results", "all_M_relF_results.RDS"))
plot.df <- make_annual_relerror_df("M", all_M_relF)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_M_df <- temp

all_Sel_relF <- readRDS(file = here("Project_0","results", "all_Sel_relF_results.RDS"))
plot.df <- make_annual_relerror_df("Sel", all_Sel_relF)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_Sel_df <- temp

all_q_relF <- readRDS(file = here("Project_0","results", "all_q_relF_results.RDS"))
plot.df <- make_annual_relerror_df("q", all_q_relF)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_q_df <- temp


R_S_plt <- ggplot(filter(R_S_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ R_sig + NAA_sig, labeller = labeller(Fhist = label_parsed, R_sig = label_parsed, NAA_sig = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    coord_cartesian(ylim = ylims) + ylab(bquote(Median~Relative~Error~(italic(F)))) + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_S_plt


R_M_plt <- ggplot(filter(R_M_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ M_sig + M_cor, labeller = labeller(Fhist = label_parsed, M_sig = label_parsed, M_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    coord_cartesian(ylim = ylims) + ylab(bquote(Median~Relative~Error~(italic(F)))) + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_M_plt



R_Sel_plt <- ggplot(filter(R_Sel_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ Sel_sig + Sel_cor, labeller = labeller(Fhist = label_parsed, Sel_sig = label_parsed, Sel_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    coord_cartesian(ylim = ylims) + ylab(bquote(Median~Relative~Error~(italic(F)))) + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_Sel_plt



R_q_plt <- ggplot(filter(R_q_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ q_sig + q_cor, labeller = labeller(Fhist = label_parsed, q_sig = label_parsed, q_cor = label_parsed)) +
    coord_cartesian(ylim = ylims) + ylab(bquote(Median~Relative~Error~(italic(F)))) + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_q_plt


R_Sel_plt_alt <- R_Sel_plt + 
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21, show.legend = TRUE) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) 


theme_set(theme_bw())
theme_update(strip.placement = "outside", strip.background = element_rect(), title = element_text(size = rel(1.5)))
cairo_pdf(here("Project_0","manuscript", "term_F_bias_plots.pdf"), width = 30*2/3, height = 20*2/3)
design <- c(area(1,1,1,1), area(2,1,2,1), area(1,2,1,2), area(2,2,2,2))
(R_S_plt +  xlab("") + labs(title = "A: R, R+S OMs") + theme(axis.title.x = element_blank(), strip.text.y=element_blank(), axis.text.x=element_blank(), legend.position="none")) +
(R_M_plt + labs(title = "B: R+M OMs") + theme(legend.position="none", strip.text.y=element_blank())) + 
(R_Sel_plt_alt + xlab("") + labs(title = "C: R+Sel OMs") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())) +
(R_q_plt + labs(title = "D: R+q OMs") + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position = "none")) + 
  plot_layout(design = design, axis_titles = "collect")
dev.off()


all_naa_relR <- readRDS(file = here("Project_0","results", "all_naa_relR_results.RDS"))
plot.df <- make_annual_relerror_df("naa", all_naa_relR)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_S_df <- temp


all_M_relR <- readRDS(file = here("Project_0","results", "all_M_relR_results.RDS"))
plot.df <- make_annual_relerror_df("M", all_M_relR)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_M_df <- temp

all_Sel_relR <- readRDS(file = here("Project_0","results", "all_Sel_relR_results.RDS"))
plot.df <- make_annual_relerror_df("Sel", all_Sel_relR)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_Sel_df <- temp

all_q_relR <- readRDS(file = here("Project_0","results", "all_q_relR_results.RDS"))
plot.df <- make_annual_relerror_df("q", all_q_relR)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_q_df <- temp


R_S_plt <- ggplot(filter(R_S_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ R_sig + NAA_sig, labeller = labeller(Fhist = label_parsed, R_sig = label_parsed, NAA_sig = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    coord_cartesian(ylim = ylims) + ylab(bquote(Median~Relative~Error~(italic(R)))) + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_S_plt


R_M_plt <- ggplot(filter(R_M_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ M_sig + M_cor, labeller = labeller(Fhist = label_parsed, M_sig = label_parsed, M_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    coord_cartesian(ylim = ylims) + ylab(bquote(Median~Relative~Error~(italic(R)))) + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_M_plt



R_Sel_plt <- ggplot(filter(R_Sel_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ Sel_sig + Sel_cor, labeller = labeller(Fhist = label_parsed, Sel_sig = label_parsed, Sel_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    coord_cartesian(ylim = ylims) + ylab(bquote(Median~Relative~Error~(italic(R)))) + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_Sel_plt



R_q_plt <- ggplot(filter(R_q_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ q_sig + q_cor, labeller = labeller(Fhist = label_parsed, q_sig = label_parsed, q_cor = label_parsed)) +
    coord_cartesian(ylim = ylims) + ylab(bquote(Median~Relative~Error~(italic(R)))) + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_q_plt


R_Sel_plt_alt <- R_Sel_plt + 
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21, show.legend = TRUE) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) 


theme_set(theme_bw())
theme_update(strip.placement = "outside", strip.background = element_rect(), title = element_text(size = rel(1.5)))
cairo_pdf(here("Project_0","manuscript", "term_R_bias_plots.pdf"), width = 30*2/3, height = 20*2/3)
design <- c(area(1,1,1,1), area(2,1,2,1), area(1,2,1,2), area(2,2,2,2))
(R_S_plt +  xlab("") + labs(title = "A: R, R+S OMs") + theme(axis.title.x = element_blank(), strip.text.y=element_blank(), axis.text.x=element_blank(), legend.position="none")) +
(R_M_plt + labs(title = "B: R+M OMs") + theme(legend.position="none", strip.text.y=element_blank())) + 
(R_Sel_plt_alt + xlab("") + labs(title = "C: R+Sel OMs") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())) +
(R_q_plt + labs(title = "D: R+q OMs") + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position = "none")) + 
  plot_layout(design = design, axis_titles = "collect")
dev.off()

