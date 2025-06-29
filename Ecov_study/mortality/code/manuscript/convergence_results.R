library(here)
library(ggplot2)
library(ggh4x)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

#convergence information:
#1: 0/1 model finished optimizing
#2: nlminb convergence flag
#3: number of NaNs in SEs for parameters, 0 = good invertible hessian
#4: max gradient value < 1e-6
#5: maximum non-NaN SE estimate < 10


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
      df <- res
    } else {
      df <- rbind(df, res)
    }
  }
  for(i in 1:5){
    df[[paste0("Type", i, "_p_pass")]] <- df[[paste0("Type", i, "_n_pass")]]/df$n_sim
    df[[paste0("Type", i, "_ci_lo")]] <- sapply(1:NROW(df), function(x) {
      binom.test(df[[paste0("Type", i, "_n_pass")]][x], df$n_sim[x], p = df[[paste0("Type", i, "_p_pass")]][x])$conf.int[1]
    })
    df[[paste0("Type", i, "_ci_hi")]]  <- sapply(1:NROW(df), function(x) {
      binom.test(df[[paste0("Type", i, "_n_pass")]][x], df$n_sim[x], p = df[[paste0("Type", i, "_p_pass")]][x])$conf.int[2]
    })
  }
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
    mutate(oe = recode(obs_error,
        "Low observation error" = "Low OE",
        "High observation error" = "High OE"
    ))
  df$oe  <- factor(df$oe, levels = c("Low OE", "High OE"))

  df <- df %>%
    mutate(M_assumption = recode(as.character(M_est),
      "TRUE" = "Estimated",
      "FALSE" = "Known"
    ))
  # df$M_assumption <- factor(df$M_assumption)
  df <- df %>%
    mutate(Ecov_assumption = recode(as.character(Ecov_est),
      "TRUE" = "beta[italic(E)]*' Estimated'",
      "FALSE" = "beta[italic(E)] == 0"
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
  df$correct_EM_PE <- "Yes"
  df$correct_EM_PE[df$OM_process_error != df$EM_process_error] <- NA
  df <- df %>% pivot_longer(
    cols = Type1_p_pass:Type5_ci_hi,
    names_to = c("Type", "locate"),
    names_pattern = "Type(.)_(.*)" #magic
  )
  df <- df %>% pivot_wider(names_from = locate, values_from = value) %>% as.data.frame
  df$Type <- as.numeric(df$Type)

  facs <- c("Ecov_obs_sig", "Fhist", "Ecov_re_sig","Ecov_re_cor", "obs_error", "NAA_M_re", "re_config", "OM_process_error", "EM_process_error", "correct_EM_PE")
  df[facs] <- lapply(df[facs], factor)
  df$OM_process_error <- factor(df$OM_process_error, levels= c("R OMs","R+S OMs", "R+M OMs"))
  df$EM_process_error <- factor(df$EM_process_error, levels= c("R","R+S", "R+M"))
  return(df)
}


df.ems = readRDS(here("Ecov_study","mortality","inputs", "df.ems.RDS"))
df.oms = readRDS(here("Ecov_study","mortality","inputs", "df.oms.RDS"))
conv_res <- readRDS(here("Ecov_study","mortality", "results", "convergence_results.RDS"))

df <- rbind(
  conv_res_plotting_fn(conv_res, M_est = FALSE, Ecov_est = TRUE),
  conv_res_plotting_fn(conv_res, M_est = TRUE, Ecov_est = TRUE),
  conv_res_plotting_fn(conv_res, M_est = FALSE, Ecov_est = FALSE),
  conv_res_plotting_fn(conv_res, M_est = TRUE, Ecov_est = FALSE))

theme_set(theme_bw())
theme_update(strip.text = element_text(size = rel(1.5)), strip.placement = "outside", strip.background = element_rect(), #fill = "transparent"), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)), legend.title.align=0.5)

levels(df$M_assumption)

# c(expression("Median "*italic(M)*" estimated"), expression("Median "*italic(M)*" known"))

# shape.labs <- list(expression("Median "*italic(M)*" estimated"), expression("Median "*italic(M)*" known"))
# library(latex2exp) #https://github.com/stefano-meschiari/latex2exp
# scale_color_viridis_d(begin = 0.2, end = 0.8, option = "turbo", labels=lapply(levels(df$model), TeX)) + 


temp <- subset(df, Type == 3 & obs_error == "Low observation error" & Fhist == "2.5*italic(F)[MSY] %->% italic(F)[MSY]")
plt <- ggplot(temp, aes(x = Ecov_effect, y = p_pass, colour = EM_process_error, fill = EM_process_error, shape = M_assumption)) + 
    facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ OM_process_error + Ecov_assumption, 
      labeller = labeller(Ecov_obs_sig = label_parsed, Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed,Ecov_assumption = label_parsed)) +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    coord_cartesian(ylim = c(0, 1)) + ylab("Probability of convergence") + xlab(expression("True "*beta[italic(E)])) +
    scale_x_continuous(breaks = c(0,0.25, 0.5)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + 
    geom_point(position = position_dodge(0.1), size = 4) + 
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi, colour = EM_process_error), width = 0, position = position_dodge(0.1), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error", shape = expression("Median "*italic(M)*" assumption")) +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", fill = "none")
plt


cairo_pdf(here("Ecov_study","mortality","manuscript", "convergence_main.pdf"), width = 30*2/3, height = 20*2/3)
print(plt)
dev.off()


OMs <- levels(df$OM_process_error)
OMs_lab <- c("Rom","RSom","RMom")
for(i in 1:length(OMs)){
  temp <- subset(df, Type == 3 & OM_process_error == OMs[i])
  plt <- ggplot(temp, aes(x = Ecov_effect, y = p_pass, colour = EM_process_error, fill = EM_process_error, shape = M_assumption)) + 
      facet_nested(Ecov_obs_sig+Ecov_re_sig+Ecov_re_cor ~ Fhist + oe + Ecov_assumption,
        labeller = labeller(Ecov_obs_sig = label_parsed, Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed,  Fhist = label_parsed, Ecov_assumption = label_parsed)) +
      scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
      scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
      coord_cartesian(ylim = c(0, 1)) + ylab("Probability of convergence") + xlab(expression("True "*beta[italic(E)])) +
      scale_x_continuous(breaks = c(0,0.25, 0.5)) + 
      geom_line(position = position_dodge(0.1), linewidth = 1) + 
      geom_point(position = position_dodge(0.1), size = 4) + 
      geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi, colour = EM_process_error), width = 0, position = position_dodge(0.1), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error", shape = expression("Median "*italic(M)*" assumption")) +
      guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", fill = "none")
  cairo_pdf(here("Ecov_study","mortality","manuscript", paste0("convergence_",OMs_lab[i],".pdf")), width = 30*2/3, height = 20*2/3)
  print(plt)
  dev.off()
}
