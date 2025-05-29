library(here)
library(ggplot2)
library(dplyr)
library(ggh4x)
library(tidyr)
library(patchwork)
library(scales)
library(ggpattern)

df.ems = readRDS(here::here("Ecov_study","mortality","inputs", "df.ems.RDS"))
df.oms = readRDS(here::here("Ecov_study","mortality","inputs", "df.oms.RDS"))
aic_res <- readRDS(here::here("Ecov_study","mortality", "results", "aic_results.RDS"))

aic_fn <- function(all, est_ind, om_ind = NULL){
  if(!is.null(om_ind)) all <- all[om_ind]
  out <- sapply(all, function(res){
    tmp <- sapply(res[est_ind], function(x) return(x))
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
    om_ind <- which(df.oms$NAA_M_re == re_mod) #om and em match
    Mfixed_rec <- aic_fn(aic_res, em_ind, om_ind)
    res <- cbind(df.oms[rep(om_ind, each = length(em_ind)),], df.ems[rep(em_ind, length(om_ind)),], n = c(Mfixed_rec))
    if(i == 1) {
      df <- res
    } else {
      df <- rbind(df, res)
    }
  }
  facs <- c("Ecov_obs_sig", "Fhist", "Ecov_re_sig","Ecov_re_cor","obs_error", "NAA_M_re", "re_config", "Ecov_est")
  df[facs] <- lapply(df[facs], factor)
  df <- df %>%
    mutate(Ecov_obs_sig = recode(Ecov_obs_sig,
      "0.1" = "sigma[italic(e)] == 0.1",
      "0.5" = "sigma[italic(e)] == 0.5"
    ))
  df <- df %>%
    mutate(obs_error = recode(obs_error,
      "L" = "Low observation error",
      "H" = "High observation error"
    ))
  df <- df %>%
    mutate(OM_process_error = recode(NAA_M_re,
      "rec" = "R OMs",
      "rec+1" = "R+S OMs",
      "rec+M" = "R+M OMs"
    ))
  df <- df %>%
    mutate(EM_process_error = recode(re_config,
      "rec" = "R",
      "rec+1" = "R+S",
      "rec+M" = "R+M"
    ))
  df <- df %>% mutate(Fhist = recode(Fhist,
      "H-MSY"  = "2.5*italic(F)[MSY] %->% italic(F)[MSY]",
      "MSY" = "italic(F)[MSY]"))
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
    mutate(Ecov_assumption = recode(as.character(Ecov_est),
      "TRUE" = "beta[italic(E)]*' Estimated'",
      "FALSE" = "beta[italic(E)] == 0"
    ))
  df <- df %>%
    mutate(M_assumption = recode(as.character(M_est),
      "TRUE" = '"Median "*italic(M)*" estimated"',
      "FALSE" = '"Median "*italic(M)*" known"'
    ))
    # mutate(M_assumption = recode(as.character(M_est),
    #   "TRUE" = "Median M estimated",
    #   "FALSE" = "Median M known"
    # ))
  df <- df %>%
    mutate(oe = recode(obs_error,
        "Low observation error" = "Low OE",
        "High observation error" = "High OE"
    ))
  df$oe  <- factor(df$oe, levels = c("Low OE", "High OE"))

  #1: correct effect assumption and RE will be 1, 
  #2: correct RE, wrong effect assumption
  #3: wrong RE, correct effect assumption
  #4: wrong RE, wrong effect assumption
  df$correct <- 0 
  # df$correct_PE <- df$correct_Ecov_effect <- as.character(NA)
  df$correct_PE <- df$correct_Ecov_effect <- "No"
  df$correct_type <- "No"
  df$correct[df$re_config == df$NAA_M_re & df$Ecov_effect == 0 & df$Ecov_est == FALSE] <- 1
  df$correct[df$re_config == df$NAA_M_re & df$Ecov_effect > 0 & df$Ecov_est == TRUE] <- 1
  
  df$correct[df$re_config == df$NAA_M_re & df$Ecov_effect == 0 & df$Ecov_est == TRUE] <- 2
  df$correct[df$re_config == df$NAA_M_re & df$Ecov_effect > 0 & df$Ecov_est == FALSE] <- 2
  
  df$correct[df$re_config != df$NAA_M_re & df$Ecov_effect == 0 & df$Ecov_est == FALSE] <- 3
  df$correct[df$re_config != df$NAA_M_re & df$Ecov_effect > 0 & df$Ecov_est == TRUE] <- 3
  
  df$correct[df$re_config != df$NAA_M_re & df$Ecov_effect == 0 & df$Ecov_est == TRUE] <- 4
  df$correct[df$re_config != df$NAA_M_re & df$Ecov_effect > 0 & df$Ecov_est == FALSE] <- 4
  
  df$correct_PE[which(df$correct %in% c(1,2))] <- "Yes"
  df$correct_Ecov_effect[which(df$correct %in% c(1,3))] <- "Yes"

  df$correct_type[df$correct %in% c(1,2)] <- "PE"
  df$correct_type[df$correct %in% c(1,3)] <- "Ecov_effect"
  df$correct_type <- factor(df$correct_type)
  df$Ecov_effect <- factor(df$Ecov_effect)
 
  df <- df %>%
    mutate(correct = recode(as.character(correct),
      "1" = "Correct Effect, Correct PE",
      "2" = "Wrong Effect, Correct PE",
      "3" = "Correct Effect, Wrong PE",
      "4" = "Wrong Effect, Wrong PE"
    ))
  df$correct <- factor(df$correct)
  return(df)
}

all_res_mod <- rbind(make_df_fn(M_est = TRUE), make_df_fn(M_est = FALSE))
theme_set(theme_bw())
theme_update(strip.text = element_text(size = rel(2)), strip.placement = "outside", strip.background = element_rect(), #fill = "transparent"), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)), legend.title.align=0.5)

temp <- subset(all_res_mod, obs_error == "Low observation error" & Fhist == "2.5*italic(F)[MSY] %->% italic(F)[MSY]" & M_est)
temp <- subset(all_res_mod, obs_error == "Low observation error" &  M_est)
OMs <- levels(all_res_mod$OM_process_error)
OMs_lab <- c("Rom","RSom","RMom")
for(i in 1:length(OMs)){
  temp <- subset(all_res_mod, OM_process_error == OMs[i])
  plt <- ggplot(temp, aes(x = Ecov_effect, y = n, fill = correct)) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) + 
      geom_col(position = "fill") + labs(fill = "EM assumption") + 
      facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ M_assumption + Fhist + oe,
        labeller = labeller(Ecov_obs_sig = label_parsed, Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed,  Fhist = label_parsed, M_assumption = label_parsed)) +
      coord_cartesian(ylim = c(0, 1)) + ylab("Proportion with lowest AIC") + xlab(expression("True "*beta[italic(E)]))
  plt
  cairo_pdf(here("Ecov_study","mortality","manuscript", paste0("aic_",OMs_lab[i],".pdf")), width = 30*2/3, height = 20*2/3)
  print(plt)
  dev.off()
}

temp <- subset(all_res_mod, obs_error == "Low observation error" & Fhist == "2.5*italic(F)[MSY] %->% italic(F)[MSY]" & M_est)
  plt <- ggplot(temp, aes(x = Ecov_effect, y = n, fill = correct)) + scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) + 
      geom_col(position = "fill") + labs(fill = "EM assumption") + 
      facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ OM_process_error,
        labeller = labeller(Ecov_obs_sig = label_parsed, Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed)) +
      coord_cartesian(ylim = c(0, 1)) + ylab("Proportion with lowest AIC") + xlab(expression("True "*beta[italic(E)]))
  plt
cairo_pdf(here("Ecov_study","mortality","manuscript", "aic_main.pdf"), width = 30*2/3, height = 20*2/3)
print(plt)
dev.off()
