for(SR_par in c("a","b","M", "SSB")) {
  print(SR_par)
  par_type <- ifelse(SR_par %in% c("a","b"), "SR", SR_par)
  print(par_type)
  glm_fits[[SR_par]] <- dev.tables[[SR_par]] <- PRD.tables[[SR_par]] <- full.trees[[SR_par]] <- list()
  for(OM_type in names(factors[[par_type]])){
    print(OM_type)
    facs <- factors[[par_type]][[OM_type]]
    dfs <- obs_dfs[[par_type]]

    if(OM_type == "R") temp <- subset(dfs[["naa"]], OM_NAA_sigma == 0 )
    if(OM_type == "R+S") temp <- subset(dfs[["naa"]], OM_NAA_sigma != 0)
    if(OM_type == "R+M") temp <- dfs[["M"]]
    if(OM_type == "R+Sel") temp <- dfs[["Sel"]]
    if(OM_type == "R+q") temp <- dfs[["q"]]
    if(SR_par %in% c("a","b","M")) temp <- subset(temp, par == paste0("italic(",SR_par,")"))
    print(dim(temp))
    if(!is.na(cv_limit)) temp <- subset(temp, cv < cv_limit) #delta-method based cv, not log-normal
    print(dim(temp))

    temp$relerror_trans <- log(temp$relerror + 1)
    temp$relerror_trans[which(is.infinite(temp$relerror_trans))] <- NA
    temp$relerror_trans2 <- log(abs(temp$relerror_trans)) # log of absolute errors on log scale (higher values are differences further from 0)
    glm_fits[[SR_par]][[OM_type]] <- list()
    dev.tables[[SR_par]][[OM_type]] <- list()
    for(i in facs){
      glm_fits[[SR_par]][[OM_type]][[i]] <- glm(as.formula(paste("relerror_trans", "~", i)), family = gaussian, data = temp)
    }
    sapply(glm_fits[[SR_par]][[OM_type]][facs[-1]], \(x) anova(x, test = "LRT")[[5]][2])
    glm_fits[[SR_par]][[OM_type]][["all"]] <- glm(as.formula(paste("relerror_trans", "~", paste(facs,collapse = "+"))), family = gaussian, data = temp)
    glm_fits[[SR_par]][[OM_type]][["all2"]] <- glm(as.formula(paste("relerror_trans", "~ (", paste(facs[-1],collapse = "+"), ")^2")), family = gaussian, data = temp)
    glm_fits[[SR_par]][[OM_type]][["all3"]] <- glm(as.formula(paste("relerror_trans", "~ (", paste(facs[-1],collapse = "+"), ")^3")), family = gaussian, data = temp)

    #percent reduction in deviance
    dev.tables[[SR_par]][[OM_type]] <- sapply(glm_fits[[SR_par]][[OM_type]][facs], \(x) 1 - x$deviance/glm_fits[[SR_par]][[OM_type]][[1]]$null.deviance)

    #Regression trees
    form <- as.formula(paste("relerror_trans ~", paste(facs, collapse = "+")))# (EM_PE + OM_R_SD + EM_M + EM_SR + OM_Obs._Error + OM_F_History)
    full.trees[[SR_par]][[OM_type]] <- rpart(form, data=temp, method = "anova", control=rpart.control(cp=0, xval = 100), model = TRUE)#, roundint = FALSE)
    # full.trees[[SR_par]][[OM_type]] <- add_to_frame(full.trees[[SR_par]][[OM_type]], temp)
    print("OM_type done")
  }
}
dev.tables[["a"]]
plot.prune(full.trees[["a"]][["R"]], cp = 0.005, type = "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0))
