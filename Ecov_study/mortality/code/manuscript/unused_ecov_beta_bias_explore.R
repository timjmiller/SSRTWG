# tail(sort(dfs[["mean_M"]]$error),20)
# 
# head(sort(dfs[["mean_M"]]$error),20)
# tail(sort(dfs[["mean_M"]]$error),20)
# head(dfs[["ecov_beta"]])
# head(sort(dfs[["ecov_beta"]]$error),20)
# tail(sort(dfs[["ecov_beta"]]$error),20)
# hist(dfs[["mean_M"]]$error)
# hist(dfs[["ecov_beta"]]$error)
# temp <- subset(dfs[["ecov_beta"]], conv == "True")
# temp <- subset(temp, error > -5 & error < 5)
# hist(temp[[1]]$error)
# hist(dfs[["ecov_beta"]]$error)
# 
# head(sort(temp$error),20)
# tail(sort(temp$error),20)
# temp <- list(ecov_beta = temp)
# facs <- factors[!factors %in% c("EM_beta_ecov", "conv")]
# fits <- get_bias_reg_fits(factors = facs, dfs = temp, type = "ecov_beta")
# prd_table <- get_bias_PRD_tables(glm_fits=fits, factors = factors[-1])
# prd_table*100
# round(prd_table/max(prd_table),2)
# round(PRD.tables[["ecov_beta"]][["converged"]]/max(PRD.tables[["ecov_beta"]][["converged"]]),2)
# ind <- which(temp[["ecov_beta"]]$error < 172.8)
# temp[["ecov_beta"]]$error
# x <- sapply(fits[["R"]], AIC)
# x-min(x)


##################################################################
#quantile (median) regression
library(quantreg)
qr_fits <- list()
for(i in c("mean_M", "ecov_beta")){
  if(i == "mean_M") facs <- factors[factors != "EM_M"]
  if(i == "ecov_beta") facs <- factors[factors != "EM_beta_ecov"]
  qr_fits[[i]] <- list()
  qr_fits[[i]][["complete"]] <- get_bias_reg_fits(factors = facs, dfs = dfs, type = i, method = "qr")
  qr_fits[[i]][["converged"]] <- get_bias_reg_fits(factors = facs, dfs = dfs, type = i, method = "qr", is_SE = TRUE)
}

#R1 is analogous to R^2 which is reduction in deviance for gaussian regression
R1.tables <- list()
for(i in c("mean_M", "ecov_beta")){
  R1.tables[[i]] <- list()
  for(j in c("complete", "converged")){
    R1.tables[[i]][[j]] <- get_bias_PRD_tables(glm_fits=qr_fits[[i]][[j]], factors = factors[-1], method = "qr")
  }
}
R1.tables[["ecov_beta"]]


############################################################
#AIC tables
AIC.tables <- list()
for(i in c("mean_M", "ecov_beta")){
  AIC.tables[[i]] <- list()
  for(j in c("complete", "converged", "se", "ci_coverage")){
    AIC.tables[[i]][[j]] <- get_bias_PRD_tables(glm_fits=glm_fits[[i]][[j]], factors = factors[-1], aic = TRUE)
    for(k in 1:NCOL(AIC.tables[[i]][[j]])) AIC.tables[[i]][[j]][,k] <- AIC.tables[[i]][[j]][,k]- min(AIC.tables[[i]][[j]][,k])
  }
}

  
#########################################################
#what is the best deviance we could hope for with different means for all interactions of factors
factors <- factors[which(factors != "conv")]
factors <- factors[factors != "EM_beta_ecov"]
df <- dfs[["ecov_beta"]]
df <- subset(df, conv == "True")
for(i in c("R","R+S", "R+M")){
  temp <- subset(df, OM_PE == i)
  glm_fits[["ecov_beta"]][["converged"]][[i]][["full"]] <- glm(as.formula(paste("error ~ ", paste(factors[-1],collapse = "*"))), family = "gaussian",
    data = temp)
}


#Percent possible reduction in deviance (reduction for a given model divided by reduction for "full" model with interactions of all factors)

x <- PRD.tables[["ecov_beta"]][["converged"]]
for(i in c("R","R+S", "R+M")){
  x[,i] <- x[,i]/(1 - glm_fits[["ecov_beta"]][["converged"]][[i]][["full"]]$deviance/glm_fits[["ecov_beta"]][["converged"]][[i]][["full"]]$null.deviance)
}
sapply(c("R","R+S", "R+M"), \(i) (1 - glm_fits[["ecov_beta"]][["converged"]][[i]][["full"]]$deviance/glm_fits[["ecov_beta"]][["converged"]][[i]][["full"]]$null.deviance))
