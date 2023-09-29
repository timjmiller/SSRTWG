library(here)
library(ggplot2)
df.ems = readRDS(here("Ecov_study","mortality","inputs", "df.ems.RDS"))
df.oms = readRDS(here("Ecov_study","mortality","inputs", "df.oms.RDS"))

#x <- readRDS(here::here("Ecov_study","mortality", "results", paste0("om", 1), paste0("sim", 3, "_em", 4, ".RDS")))
#x <- readRDS(here::here("Ecov_study","mortality", "results", "aggregated_results", paste0("om_", 1, "_em_", 4, ".RDS")))
conv_res = lapply(1:NROW(df.oms), function(y){
#conv_res = lapply(1:1, function(y){
  lapply(1:NROW(df.ems), function(x) {
  #lapply(4:4, function(x) {
    res <- readRDS(here::here("Ecov_study","mortality", "results", "aggregated_results", paste0("om_", y, "_em_", x, ".RDS")))
    print(paste(y, x))
    conv_res = t(sapply(res, function(z){
      out <- rep(NA,5) #convergence type  1, and 2
      if(!is.character(z)) {
#        print(names(z))
        if(!is.null(z$fit$opt)) {
          out[1] <- 1
          out[2] <- z$fit$opt$conv
        }
        if(!is.null(z$fit$sdrep)) out[3] <- as.integer(sum(sapply(z$fit$sdrep$SE_par, function(g) any(is.nan(g)))))
        if(!is.null(z$fit$final_gradient)) out[4] <- max(abs(z$fit$final_gradient))
        if(!is.null(z$fit$sdrep)) out[5] <- max(sapply(z$fit$sdrep$SE_par, function(g) max(g,na.rm=TRUE)))
      }
      return(out)
    }))
    return(conv_res)
  })
})
saveRDS(conv_res, file = file.path(here(),"Ecov_study","mortality", "results", "convergence_results.RDS"))

aic_res = lapply(1:NROW(df.oms), function(y){
  lapply(1:NROW(df.ems), function(x) {
  #lapply(4:4, function(x) {
    res <- readRDS(here::here("Ecov_study","mortality", "results", "aggregated_results", paste0("om_", y, "_em_", x, ".RDS")))
    print(paste(y, x))
    aic = sapply(res, function(z){
      out <- NA
      if(!is.character(z)) if(!is.null(z$fit$opt)) out = 2*(z$fit$opt$obj + length(z$fit$opt$par))
      return(out)
    })
    return(aic)
  })
})
saveRDS(aic_res, file = file.path(here(),"Ecov_study","mortality", "results", "aic_results.RDS"))

#fit = try(readRDS(here::here("Ecov_study","mortality", "results", paste0("om", 1), paste0("sim", 1, "_em", 1, ".RDS"))), silent = TRUE)
x <- readRDS(here::here("Ecov_study","mortality", "results", "aggregated_results", paste0("om_", 53, "_em_", 4, ".RDS")))
x <- readRDS(here::here("Ecov_study","mortality", "results", "aggregated_results", paste0("om_", 53, "_em_", 10, ".RDS")))

#res <- readRDS(here::here("Ecov_study","mortality", "results", "aggregated_results", paste0("om_", 1, "_em_", 1, ".RDS")))

#bias of ecov beta
ecov_beta_bias_res = lapply(1:NROW(df.oms), function(y){
  lapply(which(df.ems$Ecov_est), function(x) { #only the estimating models that estimate Ecov_beta
  #lapply(4:4, function(x) {
    res <- readRDS(here::here("Ecov_study","mortality", "results", "aggregated_results", paste0("om_", y, "_em_", x, ".RDS")))
    print(paste(y, x))
  #res = sapply(1:20, function(x) {
    bias <- t(sapply(res, function(z) { 
      #fit = try(readRDS(here::here("Ecov_study","mortality", "results", paste0("om", y), paste0("sim", x, "_em", z, ".RDS"))), silent = TRUE)
      out <- rep(NA,3) #bias, se, true_in_ci (0/1)
      if(!is.character(z)) if(!is.null(z$fit$opt)) {
        out[1] = z$fit$opt$par["Ecov_beta"] - df.oms$Ecov_effect[y]
        if(!is.null(z$fit$sdrep)){
          ind <- which(!is.na(z$fit$sdrep$SE_par$Ecov_beta))
          if(length(ind)){
            out[2] <- z$fit$sdrep$SE_par$Ecov_beta[ind[1]] #se
            ci = z$fit$opt$par["Ecov_beta"] + c(-1,1)*qnorm(0.975) * out[2]
            out[3] <- as.integer(df.oms$Ecov_effect[y] >= ci[1] & df.oms$Ecov_effect[y] <= ci[2])
          }
        }
      }
      return(out)
    }))
    return(bias)
  })
})
saveRDS(ecov_beta_bias_res, file = here::here("Ecov_study","mortality", "results", "ecov_beta_bias_results.RDS"))

#bias of beta_M (mean M)
mean_M_bias_res = lapply(1:NROW(df.oms), function(y){
  lapply(which(df.ems$Ecov_est), function(x) { #only the estimating models that estimate Ecov_beta
  #lapply(4:4, function(x) {
    res <- readRDS(here::here("Ecov_study","mortality", "results", "aggregated_results", paste0("om_", y, "_em_", x, ".RDS")))
    print(paste(y, x))
  #res = sapply(1:20, function(x) {
    bias <- t(sapply(res, function(z) { 
      #fit = try(readRDS(here::here("Ecov_study","mortality", "results", paste0("om", y), paste0("sim", x, "_em", z, ".RDS"))), silent = TRUE)
      out <- rep(NA,3) #bias, se, true_in_ci (0/1)
      if(!is.character(z)) if(!is.null(z$fit$opt)) {
        out[1] = z$fit$opt$par["M_a"] - z$truth$M_a[1]
        if(!is.null(z$fit$sdrep)){
          ind <- which(!is.na(z$fit$sdrep$SE_par$M_a))[1]
          if(length(ind)){
            out[2] <- z$fit$sdrep$SE_par$M_a[ind] #se
            ci = z$fit$opt$par["M_a"] + c(-1,1)*qnorm(0.975) * out[2]
            out[3] <- as.integer(z$truth$M_a[1] >= ci[1] & z$truth$M_a[1] <= ci[2])
          }
        }
      }
      return(out)
    }))
    return(bias)
  })
})
saveRDS(mean_M_bias_res, file = here::here("Ecov_study","mortality", "results", "mean_M_bias_results.RDS"))

#temp = try(readRDS(here::here("Ecov_study","mortality", "results", paste0("om", 1), paste0("sim", 1, "_em", 1, ".RDS"))), silent = TRUE)

M_res = lapply(1:NROW(df.oms), function(i){
  res_i <- array(NA, dim = c(100, NROW(df.ems), 40, 2)) #nsims, nems, nyears, 2=(true, est)
  for(j in 1:NROW(df.ems)){
    print(paste(i, j))
    res <- readRDS(here::here("Ecov_study","mortality", "results", "aggregated_results", paste0("om_", i, "_em_", j, ".RDS")))
    for(k in 1:dim(res_i)[1]) if(!is.character(res[[k]])) if(!is.null(res[[k]]$fit$opt)) {
      res_i[k,j,,1] <-res[[k]]$truth$MAA[,1]
      res_i[k,j,,2] <- res[[k]]$fit$rep$MAA[,1]
    }
  }
  return(res_i)
})
saveRDS(M_res, file = here::here("Ecov_study","mortality", "results", "M_results.RDS"))


ssb_res = lapply(1:NROW(df.oms), function(i){
  res_i <- array(NA, dim = c(100, NROW(df.ems), 40, 2)) #nsims, nems, nyears, 2=(true, est)
  for(j in 1:NROW(df.ems)){
    print(paste(i, j))
    res <- readRDS(here::here("Ecov_study","mortality", "results", "aggregated_results", paste0("om_", i, "_em_", j, ".RDS")))
    for(k in 1:dim(res_i)[1]) if(!is.character(res[[k]])) if(!is.null(res[[k]]$fit$opt)) {
      res_i[k,j,,1] <-res[[k]]$truth$SSB
      res_i[k,j,,2] <- res[[k]]$fit$rep$SSB
    }
  }
  return(res_i)
})
saveRDS(ssb_res, file = here::here("Ecov_study","mortality", "results", "ssb_results.RDS"))

F_res = lapply(1:NROW(df.oms), function(i){
  res_i <- array(NA, dim = c(100, NROW(df.ems), 40, 2)) #nsims, nems, nyears, 2=(true, est)
  for(j in 1:NROW(df.ems)){
    print(paste(i, j))
    res <- readRDS(here::here("Ecov_study","mortality", "results", "aggregated_results", paste0("om_", i, "_em_", j, ".RDS")))
    for(k in 1:dim(res_i)[1]) if(!is.character(res[[k]])) if(!is.null(res[[k]]$fit$opt)) {
      res_i[k,j,,1] <-res[[k]]$truth$F
      res_i[k,j,,2] <- res[[k]]$fit$rep$F
      # print(res_i[k,j,,])
      # stop()
    }
  }
  return(res_i)
})
saveRDS(F_res, file = here::here("Ecov_study","mortality", "results", "F_results.RDS"))

################################################################
#depletion (SSB terminal/SSB in year 1)
depletion_error <- data.frame(om = integer(), em = integer(), sim = integer(), depletion = numeric())
for(h in 1:length(ssb_res)) for(i in 1:100)for(j in 1:NROW(df.ems))  {
  print(paste(h,i,j))
  true <- ssb_res[[h]][i,j,40,1]/ssb_res[[h]][i,j,1,1] #true depletion
  est <- ssb_res[[h]][i,j,40,2]/ssb_res[[h]][i,j,1,2] #estimated depletion
  depletion_error = rbind.data.frame(depletion_error, c(h, j, i, (est-true)/true))
}
colnames(depletion_error) <- c("om", "em", "sim", "depletion")
saveRDS(depletion_error, file = here::here("Ecov_study","mortality", "results", "depletion_error_results.RDS"))

#########################################################
#bias in terminal SSB

ssb_term_bias_res = lapply(1:length(ssb_res), function(i){
  res_i <- array(NA, dim = c(100, NROW(df.ems))) #nsims, nems, nyears, 2=(true, est)
  for(j in 1:NROW(df.ems)){
    print(paste(i, j))
    res <- readRDS(here::here("Ecov_study","mortality", "results", "aggregated_results", paste0("om_", i, "_em_", j, ".RDS")))
    for(k in 1:dim(res_i)[1]) if(!is.character(res[[k]])) if(!is.null(res[[k]]$fit$opt)) {
      true <- res[[k]]$truth$SSB[40]
      est <- res[[k]]$fit$rep$SSB[40]
      res_i[k,j] <- (est-true)/true
    }
  }
  return(res_i)
})

saveRDS(ssb_term_bias_res, file = here::here("Ecov_study","mortality", "results", "ssb_term_bias_res.RDS"))

res <- readRDS(here::here("Ecov_study","mortality", "results", "aggregated_results", paste0("om_", 1, "_em_", 1, ".RDS")))
ind <- which(sapply(res, length) != 0)
sapply(ind, function(x) max(abs(res[[x]]$fit$final_gradient)))
which(sapply(res[[50]]$fit$sdrep$SE_par, function(x) any(is.nan(x))))
res[[50]]$fit$opt


res <- readRDS(here::here("Ecov_study","mortality", "results", "aggregated_results", paste0("om_", 184, "_em_", 6, ".RDS")))
ind <- which(sapply(res, length) != 0)
sapply(ind, function(x) max(res[[x]]$fit$final_gradient))

sapply(res[[4]]$fit$sdrep$SE_par, function(x) any(is.nan(x)))
sapply(res[[13]]$fit$sdrep$SE_par, function(x) any(is.nan(x)))
which(sapply(res[[4]]$fit$sdrep$SE_par, function(x) any(is.nan(x))))
res[[13]]$fit$sdrep$SE_par$Ecov_process_pars
res[[15]]$fit$sdrep$SE_par$Ecov_process_pars
