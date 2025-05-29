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
F_bias <- readRDS(file.path(here::here(),"Ecov_study","mortality", "results", "F_results.RDS"))
conv_res <- readRDS(here("Ecov_study","mortality", "results", "convergence_results.RDS"))
library(reshape2)
#res <- melt(F_bias)

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

custom_boxplot_stat <- function(x){#, n.yrs = 1, n.sim = 100) {
  x <- x[which(!is.na(x))]
  n <- length(x)
  bnds95 <- qbinom(c(0.025,0.975), n, 0.5)/n # 95% CI bounds for median (sims x years)
  bnds80 <- qbinom(c(0.1,0.9), n, 0.5)/n # 80% CI bounds for median (sims x years)
  r <- quantile(x, probs = c(bnds95[1], bnds80[1], 0.5, bnds80[2], bnds95[2]))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
plot_df_fn <- function(df.ems, df.oms, Ecov_est = FALSE, M_est = FALSE, conv_type = 3) {
  print(which(df.oms$NAA_M_re == "rec+1"))
  for(i in 1:3) {
    re_mod <- c("rec", "rec+1", "rec+M")[i]
    #EM:  M fixed, mean_M estimated
    em_ind <- which(df.ems$Ecov_est== Ecov_est & df.ems$M_est == M_est)
    print("em_ind")
    print(em_ind)
    om_ind <- which(df.oms$NAA_M_re == re_mod) #om and em match
    res <- lapply(om_ind, function(x) {
      F_res <- lapply(1:40, function(y) {
        out <- matrix(NA,length(em_ind), 8)
        for(j in em_ind){
          conv_ind <- 1:100 #all of them
          if(!is.null(conv_type)) conv_ind <- conv_fn(x,j,conv_res,Type = conv_type) #subset consistent with convergence results
          r_e <- F_bias[[x]][conv_ind,j,y,2]/F_bias[[x]][conv_ind,j,y,1]-1
          log_M_e <- log(F_bias[[x]][conv_ind,j,y,2]) - log(F_bias[[x]][conv_ind,j,y,1])
          rmse <- sqrt(mean((F_bias[[x]][conv_ind,j,y,2]-F_bias[[x]][conv_ind,j,y,1])^2, na.rm = TRUE))
          out[which(em_ind==j),] <- c(median(r_e, na.rm = TRUE), sd(r_e, na.rm=T)/sqrt(sum(!is.na(r_e))), custom_boxplot_stat(r_e), rmse)
        }
        colnames(out) <- c("bias_est", "bias_se", "ymin", "lower", "middle", "upper", "ymax", "rmse")
        return(out)
      })
    })
    res <- reshape2::melt(res)
    colnames(res) <- c("em_config","type","value","year", "om")
    res$om <- om_ind[res$om]

    res <- cbind.data.frame(df.oms[res$om,], res)
    if(i == 1) {
      all_res <- res
    } else {
      all_res <- rbind.data.frame(all_res, res)
    }
  }
  all_res$type <- c("bias_est", "bias_se", "ymin", "lower", "middle", "upper", "ymax", "rmse")[all_res$type]
  all_res$em_config <- c("rec", "rec+1", "rec+M")[all_res$em_config]
  all_res<- all_res %>% tidyr::pivot_wider(names_from = type, values_from = value) %>% as.data.frame

  facs <- c("Ecov_obs_sig", "Fhist", "Ecov_re_sig","Ecov_re_cor", "obs_error", "NAA_M_re", "em_config")
  
  all_res[facs] <- lapply(all_res[facs], factor)
  
  df <- all_res
  df <- all_res %>%
    mutate(Ecov_obs_sig = recode(Ecov_obs_sig,
      "0.1" = "sigma[italic(e)] == 0.1",
      "0.5" = "sigma[italic(e)] == 0.5"
    ))
  df <- df %>%
    mutate(Ecov_re_sig = recode(Ecov_re_sig,
      "0.1" = "sigma[italic(E)] == 0.1",
      "0.5" = "sigma[italic(E)] == 0.5"
    ))
  df <- df %>%
    mutate(Ecov_re_cor = recode(Ecov_re_cor,
      "0" = "rho[italic(E)] == 0",
      "0.5" = "rho[italic(E)] == 0.5"
    ))
  df <- df %>%
    mutate(beta_Ecov = recode(as.character(Ecov_est),
      "TRUE" = 'beta[italic(E)]*" estimated"',
      "FALSE" = 'beta[italic(E)]==0'
    ))
  df <- df %>%
    mutate(obs_error = recode(obs_error,
      "L" = "Low observation error",
      "H" = "High observation error"
    ))
  df <- df %>%
    mutate(oe = recode(obs_error,
        "Low observation error" = "Low OE",
        "High observation error" = "High OE"
    ))
  df$oe  <- factor(df$oe, levels = c("Low OE", "High OE"))
  df <- df %>%
    mutate(OM_process_error = recode(NAA_M_re,
      "rec" = "R OMs",
      "rec+1" = "R+S OMs",
      "rec+M" = "R+M OMs"
    ))
  df <- df %>%
    mutate(EM_process_error = recode(em_config,
      "rec" = "R",
      "rec+1" = "R+S",
      "rec+M" = "R+M"
    ))
  df <- df %>% mutate(Fhist = recode(Fhist,
      "H-MSY"  = "2.5*italic(F)[MSY] %->% italic(F)[MSY]",
      "MSY" = "italic(F)[MSY]"))
  df <- df %>%
    mutate(M = recode(as.character(M_est),
      "TRUE" = "Estimated", #'"Median "*italic(M)*" estimated"',
      "FALSE" = "Known" #'"Median "*italic(M)*" known"'
    ))
  return(df)
}

all_res <- rbind(
  plot_df_fn(df.ems, df.oms, Ecov_est = FALSE, M_est = FALSE, conv_type = 3),
  plot_df_fn(df.ems, df.oms, Ecov_est = TRUE, M_est = FALSE, conv_type = 3),
  plot_df_fn(df.ems, df.oms, Ecov_est = FALSE, M_est = TRUE, conv_type = 3),
  plot_df_fn(df.ems, df.oms, Ecov_est = TRUE, M_est = TRUE, conv_type = 3))

n <- NROW(all_res)/4
all_res$M_est = rep(c(FALSE, TRUE), each = 2*n)
all_res$Ecov_est = rep(c(FALSE, TRUE,FALSE, TRUE), each = n)

all_res_less <- subset(all_res, year %in% c(1,21,40))
all_res_less$year <- factor(all_res_less$year)
all_res_less <- all_res_less %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))

theme_set(theme_bw())
theme_update(strip.text = element_text(size = rel(1.5)), strip.placement = "outside", strip.background = element_rect(), #fill = "transparent"), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)), legend.title.align=0.5, panel.grid.minor.x = element_blank())

temp <- subset(all_res_less, obs_error == "Low observation error" & Fhist == "2.5*italic(F)[MSY] %->% italic(F)[MSY]" & year == "End")
plt <- ggplot(temp, aes(x = Ecov_effect, y = bias_est, colour = EM_process_error, shape = M)) + 
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) + 
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + 
    geom_point(position = position_dodge(0.1), size = 4) + 
    scale_x_continuous(breaks = c(0,0.25,0.5)) + 
    facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ OM_process_error + beta_Ecov, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, Ecov_obs_sig = label_parsed, beta_Ecov = label_parsed)) +
    coord_cartesian(ylim = c(-0.25, 0.25)) + 
    labs(colour = "EM process error", fill = "EM process error", shape = expression("Median "*italic(M)*" assumption")) +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", fill = "none") +
    ylab(bquote(MRE(hat(italic(F)))~"in"~terminal~year)) + 
    xlab(expression("True "*beta[italic(E)])) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = .01, linewidth = 1, position = position_dodge(0.1))
plt

cairo_pdf(here("Ecov_study","mortality","manuscript", "terminal_year_F_bias_main.pdf"), width = 30*2/3, height = 20*2/3)
print(plt)
dev.off()

temp_rmse <- temp
temp_rmse$rmse[which(temp_rmse$rmse==0)] <- NA
# temp_rmse <- filter(temp_rmse, rmse<10)
plt <- ggplot(temp_rmse, aes(x = Ecov_effect, y = rmse, colour = EM_process_error, shape = M)) + 
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + 
    geom_point(position = position_dodge(0.1), size = 4) + 
    scale_x_continuous(breaks = c(0,0.25,0.5)) + 
    facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ OM_process_error + beta_Ecov, 
      labeller = labeller(Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed, Ecov_obs_sig = label_parsed, beta_Ecov = label_parsed)) +
    labs(colour = "EM process error", fill = "EM process error", shape = expression("Median "*italic(M)*" assumption")) +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", fill = "none") +
    # coord_cartesian(ylim = c(1e3,1e5), clip = 'on') + 
    # scale_y_log10(breaks = c(1e3,1e4,1e5), labels = parse(text = c("10^3","10^4","10^5"))) + #ylim(c(0,10)) +
    ylab(bquote(RMSE(hat(italic(F)))~"in"~terminal~year)) + 
    xlab(expression("True "*beta[italic(E)]))
plt

cairo_pdf(here("Ecov_study","mortality","manuscript", "terminal_year_F_rmse_main.pdf"), width = 30*2/3, height = 20*2/3)
print(plt)
dev.off()


OMs <- levels(all_res_less$OM_process_error)
OMs_lab <- c("Rom","RSom","RMom")
for(i in 1:length(OMs)){
  temp <- subset(all_res_less, OM_process_error == OMs[i] & year == "End")
  plt <- ggplot(temp, aes(x = Ecov_effect, y = bias_est, colour = EM_process_error, shape = M)) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) + 
    geom_hline(aes(yintercept=0), linewidth = 2, linetype = "dashed", colour = "grey") +
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    scale_x_continuous(breaks = c(0,0.25,0.5)) + 
    facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ Fhist + oe + beta_Ecov,
      labeller = labeller(Ecov_obs_sig = label_parsed, Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed,  Fhist = label_parsed, beta_Ecov = label_parsed)) +
    coord_cartesian(ylim = c(-0.25, 0.25)) + ylab(bquote(MRE(hat(italic(F)))~"in"~terminal~year)) + xlab(expression("True "*beta[italic(E)])) +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", fill = "none") +
    labs(colour = "EM process error", fill = "EM process error", shape = expression("Median "*italic(M)*" assumption")) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = .05, position = position_dodge(0.1))
  plt
  cairo_pdf(here("Ecov_study","mortality","manuscript", paste0("terminal_year_F_bias_",OMs_lab[i],".pdf")), width = 30*2/3, height = 20*2/3)
  print(plt)
  dev.off()
  
  temp_rmse <- temp
  temp_rmse$rmse[which(temp_rmse$rmse==0)] <- NA
  plt <- ggplot(temp_rmse, aes(x = Ecov_effect, y = rmse, colour = EM_process_error, shape = M)) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) + 
      geom_line(position = position_dodge(0.1), linewidth = 1) + 
      geom_point(position = position_dodge(0.1), size = 4) + 
      scale_x_continuous(breaks = c(0,0.25,0.5)) + 
      facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ Fhist + oe + beta_Ecov,
        labeller = labeller(Ecov_obs_sig = label_parsed, Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed,  Fhist = label_parsed, beta_Ecov = label_parsed)) +
      ylab(bquote(RMSE(hat(italic(F)))~"in"~terminal~year)) + xlab(expression("True "*beta[italic(E)])) + 
      guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", fill = "none") +
      # coord_cartesian(ylim = c(1e3,1e5), clip = 'on') + 
      # scale_y_log10(breaks = c(1e3,1e4,1e5), labels = parse(text = c("10^3","10^4","10^5"))) +
      labs(colour = "EM process error", fill = "EM process error", shape = expression("Median "*italic(M)*" assumption"))
  plt
  cairo_pdf(here("Ecov_study","mortality","manuscript", paste0("terminal_year_F_rmse_",OMs_lab[i],".pdf")), width = 30*2/3, height = 20*2/3)
  print(plt)
  dev.off()
}
