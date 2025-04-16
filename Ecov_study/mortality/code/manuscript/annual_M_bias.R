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
M_bias <- readRDS(file.path(here::here(),"Ecov_study","mortality", "results", "M_results.RDS"))
library(reshape2)
#res <- melt(ssb_bias)


#all EM PE assumptions
plot_df_fn <- function(df.ems, df.oms, Ecov_est = FALSE, M_est = FALSE) {
  custom_boxplot_stat <- function(x){#, n.yrs = 1, n.sim = 100) {
    x <- x[which(!is.na(x))]
    n <- length(x)
    bnds95 <- qbinom(c(0.025,0.975), n, 0.5)/n # 95% CI bounds for median (sims x years)
    bnds80 <- qbinom(c(0.1,0.9), n, 0.5)/n # 80% CI bounds for median (sims x years)
    r <- quantile(x, probs = c(bnds95[1], bnds80[1], 0.5, bnds80[2], bnds95[2]))
    names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
    r
  }
  for(h in 1:3) {
    re_mod <- c("rec", "rec+1", "rec+M")[h]
    #EM:  M fixed, mean_M estimated
    em_ind <- which(df.ems$Ecov_est== Ecov_est & df.ems$M_est == M_est)
    print("em_ind")
    print(em_ind)
    om_ind <- which(df.oms$NAA_M_re == re_mod) #om and em match
    res <- lapply(om_ind, function(x) {
      M_res <- lapply(1:40, function(y) {
        out <- matrix(NA,length(em_ind), 8)
        for(j in em_ind){
          r_e <- M_bias[[x]][,j,y,2]/M_bias[[x]][,j,y,1]-1
          log_M_e <- log(M_bias[[x]][,j,y,2]) - log(M_bias[[x]][,j,y,1])
          rmse <- sqrt(mean(log_M_e^2, na.rm = TRUE))
          out[which(em_ind==j),] <- c(median(r_e, na.rm = TRUE), sd(r_e, na.rm=T)/sqrt(sum(!is.na(r_e))), custom_boxplot_stat(r_e), rmse)
        }
        colnames(out) <- c("bias_est", "bias_se", "ymin", "lower", "middle", "upper", "ymax", "rmse")
        return(out)
      })
      return(M_res)
    })
    res <- reshape2::melt(res)
    colnames(res) <- c("em_config","type","value","year", "om")
    res$om <- om_ind[res$om]

    res <- cbind.data.frame(df.oms[res$om,], res)
    print(head(res))
    if(h == 1) {
      all_res <- res
    } else {
      all_res <- rbind.data.frame(all_res, res)
    }
  }
  all_res$type <- c("bias_est", "bias_se", "ymin", "lower", "middle", "upper", "ymax", "rmse")[all_res$type]
  all_res$em_config <- c("rec", "rec+1", "rec+M")[all_res$em_config]
  print(head(all_res))

  all_res<- all_res %>% tidyr::pivot_wider(names_from = type, values_from = value) %>% as.data.frame
  print(head(all_res))
  facs <- c("Ecov_obs_sig", "Fhist", "Ecov_re_sig","Ecov_re_cor", "obs_error", "NAA_M_re", "em_config")
  all_res[facs] <- lapply(all_res[facs], factor)
  all_res_mod <- all_res
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
    mutate(EM_process_error = recode(em_config,
      "rec" = "R",
      "rec+1" = "R+S",
      "rec+M" = "R+M"
    ))
  all_res_mod <- all_res_mod %>% mutate(Fhist = recode(Fhist,
      "H-MSY"  = "2.5*italic(F)[MSY] %->% italic(F)[MSY]",
      "MSY" = "italic(F)[MSY]"))
  return(all_res_mod)
}

all_res <- rbind(
  plot_df_fn(df.ems, df.oms, Ecov_est = FALSE, M_est = FALSE),
  plot_df_fn(df.ems, df.oms, Ecov_est = TRUE, M_est = FALSE),
  plot_df_fn(df.ems, df.oms, Ecov_est = FALSE, M_est = TRUE),
  plot_df_fn(df.ems, df.oms, Ecov_est = TRUE, M_est = TRUE))

n <- NROW(all_res)/4
all_res$M_est = rep(c(FALSE, TRUE), each = 2*n)
all_res$Ecov_est = rep(c(FALSE, TRUE,FALSE, TRUE), each = n)

temp <- subset(all_res, year %in% c(1,21,40))
# temp$M_est <- factor(temp$Ecov_effect)
# temp$Ecov_est <- factor(temp$Ecov_effect)
temp$year <- factor(temp$year)
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp <- subset(temp, obs_error == "Low observation error" & Fhist == "2.5*italic(F)[MSY] %->% italic(F)[MSY]" & year == "End")

theme_set(theme_bw())
theme_update(strip.text = element_text(size = rel(1.5)), strip.placement = "outside", strip.background = element_rect(), #fill = "transparent"), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)), legend.title.align=0.5, panel.grid.minor.x = element_blank())

temp1 <- subset(temp, M_est == TRUE & Ecov_est == TRUE)
plt <- ggplot(temp1, aes(x = Ecov_effect, y = bias_est, colour = EM_process_error)) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) + 
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    scale_x_continuous(breaks = c(0,0.25,0.5)) + 
    facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor  ~ OM_process_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
    coord_cartesian(ylim = c(-0.45, 0.45)) + 
    labs(colour = "EM process error") +
    ylab(bquote(Median~Relative~Error~of~hat(M))) + xlab(expression("True "*beta[Ecov])) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .01, linewidth = 1, position = position_dodge(0.1))
plt

cairo_pdf(here("Ecov_study","mortality","manuscript", "terminal_year_M_bias.pdf"), width = 30*2/3, height = 20*2/3)
print(plt)
dev.off()


temp1 <- subset(temp, M_est == TRUE & Ecov_est == FALSE)
plt <- ggplot(temp1, aes(x = Ecov_effect, y = bias_est, colour = EM_process_error)) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) + 
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    scale_x_continuous(breaks = c(0,0.25,0.5)) + 
    facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor  ~ OM_process_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
    coord_cartesian(ylim = c(-0.45, 0.45)) + 
    labs(colour = "EM process error") +
    ylab(bquote(Median~Relative~Error~of~hat(M))) + xlab(expression("True "*beta[Ecov])) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .01, linewidth = 1, position = position_dodge(0.1))
plt

cairo_pdf(here("Ecov_study","mortality","manuscript", "terminal_year_M_bias_ecov_beta_0.pdf"), width = 30*2/3, height = 20*2/3)
print(plt)
dev.off()


temp1 <- subset(temp, M_est == TRUE & Ecov_est == TRUE)
plt <- ggplot(temp1, aes(x = Ecov_effect, y = rmse, colour = EM_process_error)) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) + 
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    scale_x_continuous(breaks = c(0,0.25,0.5)) + 
    facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor  ~ OM_process_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
     coord_cartesian(ylim = c(0, 1)) + 
    labs(colour = "EM process error") +
    ylab(bquote(RMSE~log(M))) + xlab(expression("True "*beta[Ecov]))# + geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .01, linewidth = 1, position = position_dodge(0.1))
plt


temp1 <- subset(temp, M_est == TRUE & Ecov_est == FALSE)
temp1 <- subset(temp, M_est == FALSE & Ecov_est == FALSE)
plt <- ggplot(temp1, aes(x = Ecov_effect, y = rmse, colour = EM_process_error)) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) + 
    #geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    scale_x_continuous(breaks = c(0,0.25,0.5)) + 
    facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor  ~ OM_process_error, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, Ecov_obs_sig = label_parsed, Fhist = label_parsed)) +
     coord_cartesian(ylim = c(0, 1)) + 
    labs(colour = "EM process error") +
    ylab(bquote(RMSE~log(M))) + xlab(expression("True "*beta[Ecov]))# + geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .01, linewidth = 1, position = position_dodge(0.1))
plt


subset(df.oms, NAA_M_re == "rec+M" & Ecov_obs_sig == 0.1 & Ecov_re_sig == 0.5 & Ecov_re_cor == 0.5 & Ecov_effect == 0.5 & Fhist == "H-MSY" & obs_error == "L")
#em 11: M estimated, Ecov_beta == 0, EMPE = R+M
ind <- which(M_bias[[69]][,1,40,1]>0.2) #just realized terminal M> true median
median(M_bias[[69]][ind,11,40,2]/M_bias[[69]][ind,11,40,1] - 1, na.rm = TRUE)

ind <- which(M_bias[[69]][,1,40,1]<0.2) #just realized terminal M> true median
median(M_bias[[69]][ind,11,40,2]/M_bias[[69]][ind,11,40,1] - 1, na.rm = TRUE)


plot(log(M_bias[[69]][,9,40,1]), log(M_bias[[69]][,9,40,2]) - log(M_bias[[69]][,9,40,1]), ylim = c(-2,2))

om.ind <- which(df.oms$NAA_M_re == "rec+M" & df.oms$Ecov_effect == 0.5 & df.oms$Fhist == "H-MSY" & df.oms$obs_error == "L")
# df.oms.RM <- df.oms[om.ind,]
em.ind <- which(df.ems$M_est)
df.ems.RM <- subset(df.ems, M_est == TRUE)
bdf <- cbind.data.frame(logM = log(M_bias[[69]][,9,40,1]), bias = log(M_bias[[69]][,9,40,2]) - log(M_bias[[69]][,9,40,1]))

res <- readRDS(here::here("Ecov_study","mortality", "results", "aggregated_results", paste0("om_", 69, "_em_", 11, ".RDS")))
temp <- sapply(res, \(x) {
  out <- rep(NA,69)
  if(!is.null(x$fit$opt$par)) out[] <- x$fit$opt$par
  return(out)
})
res_median <- apply(temp,1,median, na.rm = TRUE)
names(res_median) <- colnames(temp) <- names(res[[1]]$fit$opt$par)
res_median <- res_median[-1] #mean recruitment but BH used for OMs
res_true <- c(res[[1]]$truth$logit_q, res[[1]]$truth$log_F1,res[[1]]$truth$F_devs,res[[1]]$truth$log_N1_pars, res[[1]]$truth$log_NAA_sigma,
  res[[1]]$truth$logit_selpars[,11:12], res[[1]]$truth$catch_paa_pars[1],res[[1]]$truth$index_paa_pars[1:2],
  res[[1]]$truth$M_a[1],res[[1]]$truth$M_repars[c(1,3)],res[[1]]$truth$Ecov_process_pars)
round(res_median - res_true,2)

res <- readRDS(here::here("Ecov_study","mortality", "results", "aggregated_results", paste0("om_", 69, "_em_", 5, ".RDS"))) #Ecov beta estimated
temp <- sapply(res, \(x) {
  out <- rep(NA,70)
  if(!is.null(x$fit$opt$par)) out[] <- x$fit$opt$par
  return(out)
})
res_median <- apply(temp,1,median, na.rm = TRUE)
names(res_median) <- names(res[[1]]$fit$opt$par)
res_median <- res_median[-1] #mean recruitment but BH used for OMs
res_true <- c(res[[1]]$truth$logit_q, res[[1]]$truth$log_F1,res[[1]]$truth$F_devs,res[[1]]$truth$log_N1_pars, res[[1]]$truth$log_NAA_sigma,
  res[[1]]$truth$logit_selpars[,11:12], res[[1]]$truth$catch_paa_pars[1],res[[1]]$truth$index_paa_pars[1:2],
  res[[1]]$truth$M_a[1],res[[1]]$truth$M_repars[c(1,3)],res[[1]]$truth$Ecov_beta[2,1,1,1], res[[1]]$truth$Ecov_process_pars)
round(res_median - res_true,2)


res <- readRDS(here::here("Ecov_study","mortality", "results", "aggregated_results", paste0("om_", 69, "_em_", 11, ".RDS")))

log(0.2) + res[[1]]$truth$M_re[,1] + 0.5*res[[1]]$truth$Ecov_x
log(res[[1]]$truth$MAA[,1])
temp <- sapply(res, \(x) {
  out <- rep(NA,40)
  if(!is.null(x$truth$MAA)) out[] <- log(x$truth$MAA[,1])
    return(out)
})
cor(t(temp), use = "complete.obs")[1:5,1:5]
cor(t(temp), use = "complete.obs")[cbind(2:40, 1:39)]
temp <- sapply(res, \(x) {
  out <- rep(NA,40)
  if(!is.null(x$fit$rep$Ecov_x)) out[] <- x$fit$rep$Ecov_x
    return(out)
})
cor(t(temp), use = "complete.obs")[1:5,1:5]
cor(t(temp), use = "complete.obs")[cbind(2:40, 1:39)]
cor(t(temp), use = "complete.obs")[1:5,1:5]
par(mfrow = c(3,3))
for(i in 1:9){
  plot(res[[i]]$truth$MAA[,1], ylim  = c(0,1))
  points(res[[i]]$fit$rep$MAA[,1], col = 'red')
}
temp <- sapply(res, \(x) {
  out <- NA
  if(!is.null(x$fit$rep$MAA)) out <- sqrt(mean((log(x$truth$MAA[,1]) - log(x$fit$rep$MAA[,1]))^2, na.rm = TRUE))
  return(out)
})
min(temp, na.rm = TRUE)

  plot(res[[23]]$truth$MAA[,1], ylim  = c(0,1))
  points(res[[23]]$fit$rep$MAA[,1], col = 'red')

temp <- sapply(res, \(x) {
  out <- rep(NA,2)
  if(!is.null(x$fit$opt$par)) {
    sd <- exp(x$truth$M_repars[1])
    rho <- -1 + 2/(1 + exp(-x$truth$M_repars[3]))
    out[1] <- sd/sqrt(1 - rho^2)
    sd <- exp(x$fit$opt$par[ind[1]])
    rho <- -1 + 2/(1 + exp(-x$fit$opt$par[ind[2]]))
    out[2] <- sd/sqrt(1 - rho^2)
  }
  return(out)
})
res_median <- apply(temp,1,median, na.rm = TRUE)


names(res_median) <- names(res[[1]]$fit$opt$par)
res_median <- res_median[-1] #mean recruitment but BH used for OMs
res_true <- c(res[[1]]$truth$logit_q, res[[1]]$truth$log_F1,res[[1]]$truth$F_devs,res[[1]]$truth$log_N1_pars, res[[1]]$truth$log_NAA_sigma,
  res[[1]]$truth$logit_selpars[,11:12], res[[1]]$truth$catch_paa_pars[1],res[[1]]$truth$index_paa_pars[1:2],
  res[[1]]$truth$M_a[1],res[[1]]$truth$M_repars[c(1,3)],res[[1]]$truth$Ecov_process_pars)
round(res_median - res_true,2)

temp <- sapply(res, \(x) {
  out <- rep(NA,40)
  if(!is.null(x$fit$opt)) out[] <- log(x$truth$MAA[,1])
  return(out)
})
annual_log_M_sd <- apply(temp,1,sd, na.rm = TRUE)
temp <- sapply(res, \(x) {
  out <- rep(NA,40)
  if(!is.null(x$fit$sdrep$Estimate_par$M_a)) if(!is.na(x$fit$sdrep$SE_par$M_a[1])) out[] <- log(x$fit$rep$MAA[,1])
  return(out)
})
annual_log_hat_M_sd <- apply(temp,1,sd, na.rm = TRUE)


ind <- which(names(res_median) == "M_repars")
exp(res_median[ind[1]])/sqrt(1 - (-1 + 2/(1 + exp(-res_median[ind[2]])))^2)

em.inputs <- readRDS(here("Ecov_study","mortality", "inputs", "em_inputs.RDS"))

bdfs <- lapply(om.ind, \(w) {
  bdfs <- lapply(em.ind, \(x) {
    res <- readRDS(here::here("Ecov_study","mortality", "results", "aggregated_results", paste0("om_", w, "_em_", x, ".RDS")))
    res <- lapply(res, \(z) {
      out <- list()
      out <- rep(NA, 11)
      names(out) <- c("Ecov_beta", "Ecov_beta_est", "Ecov_beta_se", "log_M", "log_M_est", "log_M_se", "term_M", "term_M_est", "log_M_sig","log_M_sig_est","log_M_sig_se")
      if(!is.null(z$truth$Ecov_beta[2,1,1,1])) out["Ecov_beta"] <- z$truth$Ecov_beta[2,1,1,1]
      if(!is.null(z$fit$sdrep$Estimate_par$Ecov_beta[2,1,1,1])) out["Ecov_beta_est"] <- z$fit$sdrep$Estimate_par$Ecov_beta[2,1,1,1]
      if(!is.null(z$fit$sdrep$SE_par$Ecov_beta[2,1,1,1])) out["Ecov_beta_se"] <- z$fit$sdrep$SE_par$Ecov_beta[2,1,1,1]
      if(!is.null(z$truth$M_a[1])) out["log_M"] <- z$truth$M_a[1]
      if(!is.null(z$fit$sdrep$Estimate_par$M_a[1])) out["log_M_est"] <-z$fit$sdrep$Estimate_par$M_a[1]
      if(!is.null(z$fit$sdrep$SE_par$M_a[1])) out["log_M_se"] <- z$fit$sdrep$SE_par$M_a[1]
      if(!is.null(z$truth$MAA[40,1])) out["term_M"] <- z$truth$MAA[40,1]
      if(!is.null(z$fit$rep$MAA[40,1])) out["term_M_est"] <-z$fit$rep$MAA[40,1]
      if(!is.null(z$truth$M_repars[1])) out["log_M_sig"] <- z$truth$M_repars[1]
      if(!is.null(z$fit$sdrep$Estimate_par$M_repars[1])) out["log_M_sig_est"] <- z$fit$sdrep$Estimate_par$M_repars[1]
      if(!is.null(z$fit$sdrep$SE_par$M_repars[1])) out["log_M_sig_se"] <- z$fit$sdrep$SE_par$M_repars[1]
      out <- cbind.data.frame(df.oms[w,], df.ems[x,], rbind(out))
      return(out)
    })
    res <- do.call(rbind, res)
    print(dim(res))
    return(res)
  })
  bdf <- do.call(rbind, bdfs)
})
bdf <- do.call(rbind, bdfs)
x <- subset(bdf, Ecov_est == FALSE)

sum((is.na(bdf$Ecov_beta_se ) & bdf$Ecov_est) | is.na(bdf$log_M_se))

ind <- which((is.na(bdf$Ecov_beta_se ) & bdf$Ecov_est) | is.na(bdf$log_M_se))
bdf[ind, c("Ecov_beta_est", "log_M_est", "term_M_est")] <- NA
bdf$Ecov_beta_bias <- bdf$Ecov_beta_est - bdf$Ecov_beta
bdf$log_M_bias <- bdf$log_M_est - bdf$log_M
bdf$log_term_M_bias <- log(bdf$term_M_est) - log(bdf$term_M)

ggplot(bdf, aes(x = log(term_M), y = log_term_M_bias)) + facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor~re_config+Ecov_est) +
  geom_point() + coord_cartesian(ylim = c(-2.5,2.5)) + 
  stat_smooth(method = "lm", col = "red")

x <- subset(bdf, !Ecov_est & Ecov_obs_sig == 0.1 & Ecov_re_sig == 0.5 & Ecov_re_cor == 0.5 & re_config == "rec+M")
median(x$log_M_sig_est, na.rm = TRUE)
median(x$log_M_est, na.rm = TRUE)
median(x$term_M_est, na.rm = TRUE)
plot(x$Ecov_beta_est, x$log_M_est)
summary(lm(log_M_est~Ecov_beta_est, data = x))


plot(log(M_bias[[69]][,11,40,1]), log(M_bias[[69]][,11,40,2]) - log(M_bias[[69]][,11,40,1]), ylim = c(-2,2))
bdf <- cbind.data.frame(x = log(M_bias[[69]][,11,40,1]), y = log(M_bias[[69]][,11,40,2]) - log(M_bias[[69]][,11,40,1]))
summary(lm(y~x, data = bdf))
