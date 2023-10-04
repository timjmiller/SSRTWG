tree_df_fn <- function(df.ems, df.oms, M_est = FALSE) {
  for(i in 1:3) {
    re_mod <- c("rec", "rec+1", "rec+M")[i]
    #EM:  M fixed, Ecov_beta estimated, OM and EM RE assumption match
    df.ems. <- df.ems[df.ems$Ecov_est,]
    em_ind <- which(df.ems.$M_est == M_est & df.ems.$re_config == re_mod)
    om_ind <- which(df.oms$NAA_M_re == re_mod) #om and em match
    res <- lapply(all_beta_bias[om_ind], function(x) {
      x <- x[[em_ind]]
      x <- x[which(!is.na(x[,2])),1] #remove fits without standard error estimates
      #x <- x[which(x[,2] < 100),1] #remove fits with bad (really big) standard error estimates
#      print(x)
      # out <- c(mean(x,na.rm=T), sd(x, na.rm=T)/sqrt(sum(!is.na(x))))
      # out <- c(out, out[1] + qt(0.025, sum(!is.na(x))) * out[2])
      # out <- c(out, out[1] + qt(0.975, sum(!is.na(x))) * out[2])
      # out <- c(median(x, na.rm = TRUE), sd(x, na.rm=T)/sqrt(sum(!is.na(x))), custom_boxplot_stat(x))#quantile(x, probs = c(0.025,0.975)))
      # return(out)
      return(x)
    })
    res <- reshape2::melt(res)
    res <- cbind.data.frame(df.oms[om_ind[res[,2]],], error = res[,1])
    if(i == 1) {
      all_res <- res
    } else {
      all_res <- rbind.data.frame(all_res, res)
    }
  }
  facs <- c("Ecov_obs_sig", "Fhist", "Ecov_re_sig","Ecov_re_cor", "obs_error", "NAA_M_re")
  all_res[facs] <- lapply(all_res[facs], factor)
  all_res_mod <- all_res %>%
    mutate(Ecov_obs_sig = recode(Ecov_obs_sig,
      "0.1" = "sigma[ecov] == 0.1",
      "0.5" = "sigma[ecov] == 0.5"
    ))
  # all_res_mod <- all_res_mod %>%
  #   mutate(Ecov_re_sig = recode(Ecov_re_sig,
  #     "0.1" = "sigma(Ecov) = 0.1",
  #     "0.5" = "sigma(Ecov) = 0.5"
  #   ))
  # all_res_mod <- all_res_mod %>%
  #   mutate(Ecov_re_cor = recode(Ecov_re_cor,
  #     "0" = "rho(Ecov) = 0",
  #     "0.5" = "rho(Ecov) = 0.5"
  #   ))
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
  all_res_mod <- all_res_mod %>% mutate(Fhist = recode(Fhist,
      "H-MSY" = "F history: High->FMSY",
      "MSY" = "F history: FMSY"))
  return(all_res_mod)
}

all_res <- tree_df_fn(df.ems, df.oms, M_est = FALSE)

fit <- rpart::rpart(error ~ NAA_M_re + Ecov_obs_sig + Ecov_re_sig + Ecov_re_cor + Ecov_effect + Fhist + obs_error, data=all_res, control=rpart.control(cp=0.01))
rpart.plot(fit)

fit <- rpart::rpart(error ~ NAA_M_re + Ecov_obs_sig + Ecov_re_sig + Ecov_re_cor + Ecov_effect + Fhist + obs_error, data=all_res, control=rpart.control(cp=0.001))
rpart.plot(fit)

fit <- randomForest(error ~ NAA_M_re + Ecov_obs_sig + Ecov_re_sig + Ecov_re_cor + Ecov_effect + Fhist + obs_error, data=all_res)
plot(fit)
getTree(fit, 1, labelVar = TRUE)

setwd(here::here("Ecov_study","mortality","paper"))
rmarkdown::render("paper.Rmd", output_file = "paper.pdf")
