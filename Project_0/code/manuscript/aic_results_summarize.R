library(here)
library(dplyr)
library(tidyr)
library(mgcv)
library(gratia)
library(Hmisc)
library(VGAM)
#modified rpart.plot package to use plotmath in split labels...
pkgload::load_all(file.path(here(),"../../rpart.plot"))
# pkgload::load_all("/home/tmiller/FSAM_research/aug_backup/work/rpart.plot")
# pkgload::load_all("c:/work/rpart.plot")

df.oms = list(naa = readRDS(here("Project_0", "inputs", "df.oms.RDS")))
df.oms$M = readRDS(here("Project_0", "inputs", "df.M.oms.RDS"))
df.oms$Sel = readRDS(here("Project_0", "inputs", "df.Sel.oms.RDS"))
df.oms$q = readRDS(here("Project_0", "inputs", "df.q.oms.RDS"))
om.em.rows <- list(naa = 1:20, M = 5:24, Sel = c(5:20, 25:28), q = c(5:20, 29:32))
df.ems = readRDS(here("Project_0", "inputs", "df.ems.RDS"))

types <- c("naa", "M", "Sel", "q")
types.plt <- c("R, R+S","R+M", "R+Sel", "R+q")
est_config <- c("M known, No SR", "M estimated, No SR", "M known, SR", "M estimated, SR")
EM_process_error <- c("R","R+S","R+M (iid)","R+M (AR1)","R+Sel (iid)","R+Sel (AR1)","R+q (iid)","R+q (AR1)")

om_aic_res <- list()
for(i in types) om_aic_res[[i]] <- readRDS(file = here("Project_0","results", paste0("all_", i, "_aic_results.RDS")))

aic_fn <- function(all, em_rows, om_type, df.ems, om_ind = NULL){

  #df.ems. is the subset of df.ems for the current om_type
  #est_ind indexes ROWs of df.ems.
  #20 EMs for each OM type, but the EMs are not all the same.
  #df.ems has all of them
  df.ems. <- df.ems[em_rows[[om_type]],]
  #est_ind <- which(df.ems.$SR_model == ifelse(SR_est,3,2) & df.ems.$M_est == M_est)
  # if(is.null(est_ind)) est_ind <- 1:NROW(df.ems.)

  est_ind <- list(
    SR_no_M_no = which(df.ems.$SR_model == 2 & !df.ems.$M_est),
    SR_yes_M_no = which(df.ems.$SR_model == 3 & !df.ems.$M_est),
    SR_no_M_yes = which(df.ems.$SR_model == 2 & df.ems.$M_est),
    SR_yes_M_yes = which(df.ems.$SR_model == 3 & df.ems.$M_est))
  print(est_ind)

  all <- all[[om_type]]
  if(is.null(om_ind)) om_ind <- 1:length(all)
  all <- all[om_ind]
  res_by_om <- lapply(1:length(all), \(om_i){
    # print(paste("om", om_i, "begin"))
    res_by_sim <- lapply(1:NCOL(all[[om_i]]), \(sim_j) {
      # if(om_i == 7) print(paste("sim", sim_j, "begin"))
      res_by_sim_est_ind <- lapply(1:length(est_ind), \(k){
      # if(om_i == 7 & sim_j == 12) print(paste("k", k, "begin"))
        res <- all[[om_i]][est_ind[[k]],sim_j]
      if(om_i == 7 & sim_j == 12) print(res)
        if(any(!is.na(res))) {
          if(length(which(res == min(res,na.rm=T)))>1) {
            print(res)
            stop()
          }
          best <- as.integer(res == min(res,na.rm=T))
          best_ind <- est_ind[[k]][which(best == 1)]
          out <- cbind.data.frame(EM_assumption = names(est_ind)[k], best = df.ems.$re_config[best_ind], df.ems.[best_ind,])
        } else {
          best <- rep(NA, length(res))
          best_ind <- NA
          na_df <- do.call(cbind, lapply(names(df.ems.), \(cols) cbind.data.frame(as.character(NA))))
          colnames(na_df) <- names(df.ems.)
          out <- cbind.data.frame(EM_assumption = names(est_ind)[k], best = NA, na_df)
        }
        # if(om_i == 7 & sim_j == 12) print(out)
        # if(om_i == 7 & sim_j == 12) print(paste("k", k, "end"))
        return(out)
      })
      out <- cbind.data.frame(sim = sim_j, do.call(rbind,res_by_sim_est_ind))
      # if(om_i == 7) print(paste("sim", sim_j, "begin"))
      return(out)
    })
    out <- do.call(rbind, res_by_sim)
    # print(paste("om", om_i, "end"))
    return(cbind(om = om_i, out))
  })
  out <- do.call(rbind, res_by_om)
  return(out)
}
  # aic_df[[i]] <- aic_fn(om_aic_res, em_rows = om.em.rows, om_type = i, df.ems = df.ems)

# table(aic_df[[i]]$best)
# unique(aic_df[[i]]$M_re_cor)

aic_df <- list()
for(i in names(om.em.rows)) {
  aic_df[[i]] <- aic_fn(om_aic_res, em_rows = om.em.rows, om_type = i, df.ems = df.ems)
  aic_df[[i]] <- cbind(df.oms[[i]][aic_df[[i]][,"om"],], aic_df[[i]])
  aic_df[[i]]$Model <- paste0(i,"_",aic_df[[i]]$Model)
  aic_df[[i]]$OM <- paste0(i,"_",aic_df[[i]]$om)
  aic_df[[i]]$OM_sim  <- paste0(i, "_", (aic_df[[i]][["om"]]-1)*100 + aic_df[[i]][["sim"]])
  aic_df[[i]]$OM_type <- i
  print(i)
  print(head(aic_df[[i]]))
  if(i == "naa"){
    aic_df[[i]]$NAA_sig[is.na(aic_df[[i]]$NAA_sig)] <- 0
    other_oms <- names(om.em.rows)[names(om.em.rows) != i]
    for(j in other_oms){
      aic_df[[i]][[paste0(j, "_sig")]] <- NA
      aic_df[[i]][[paste0(j, "_cor")]] <- NA
    }
  } else{
    aic_df[[i]]$NAA_sig <- NA
    aic_df[[i]][[paste0(i, "_sig")]] <- as.factor(aic_df[[i]][[paste0(i, "_sig")]])
    aic_df[[i]][[paste0(i, "_cor")]] <- as.factor(aic_df[[i]][[paste0(i, "_cor")]])
    other_oms <- names(om.em.rows)[!names(om.em.rows) %in% c(i,"naa")]
    for(j in other_oms){
      aic_df[[i]][[paste0(j, "_sig")]] <- NA
      aic_df[[i]][[paste0(j, "_cor")]] <- NA
    }
  }
  aic_df[[i]] <- aic_df[[i]] %>%
    mutate(best_PE = recode(best,
      "rec" = "R",
      "rec+1" = "R+S",
      "M_re" = "R+M",
      "q_re" = "R+q",
      "sel_re" = "R+Sel"
    ))
  aic_df[[i]]$best_PE_cor <- aic_df[[i]]$best_PE
  for(re in c("M","sel","q")){
    ind <- which(aic_df[[i]]$re_config == paste0(re,"_re") & aic_df[[i]][[paste0(re,"_re_cor")]] == "iid")
    aic_df[[i]]$best_PE_cor[ind] <- paste0(aic_df[[i]]$best_PE_cor[ind], "(iid)")
    ind <- which(aic_df[[i]]$re_config == paste0(re,"_re") & sapply(aic_df[[i]][[paste0(re,"_re_cor")]], \(x) length(grep("ar1", x, fixed = TRUE)) == 1))
    aic_df[[i]]$best_PE_cor[ind] <- paste0(aic_df[[i]]$best_PE_cor[ind], "(AR1)")
  }
  # ind <- which(aic_df[[i]]$M_re_cor == "iid" | aic_df[[i]]$sel_re_cor == "iid" | aic_df[[i]]$q_re_cor == "iid")
  # aic_df[[i]]$EM_PE[ind] <- paste0(aic_df[[i]]$EM_PE[ind], " (iid)")
  # ind <- which(aic_df[[i]]$M_re_cor == "ar1_y" | aic_df[[i]]$sel_re_cor == "ar1_y" | aic_df[[i]]$q_re_cor == "ar1")
  # aic_df[[i]]$EM_PE[ind] <- paste0(aic_df[[i]]$EM_PE[ind], " (AR1)")
  

  
  aic_df[[i]] <- aic_df[[i]] %>%
    mutate(OM_Obs._Error = recode(obs_error,
      "L" = "Low",
      "H" = "High"
    ))
  
  aic_df[[i]] <- aic_df[[i]] %>%
    mutate(OM_F_History = recode(Fhist,
      "H-MSY"  = "2.5*italic(F)[MSY] %->% italic(F)[MSY]",
      "MSY" = "italic(F)[MSY]"))
  
  aic_df[[i]] <- aic_df[[i]] %>%
  mutate(EM_SR = recode(as.character(SR_model),
                        "2" = "None",
                        "3" = "Estimated"
  ))
  
  aic_df[[i]] <- aic_df[[i]] %>%
    mutate(EM_M = recode(as.character(M_est),
                       "TRUE" = "Estimated",
                       "FALSE" = "Known"
  ))
  
  aic_df[[i]] <- aic_df[[i]][sort(names(aic_df[[i]]))]

}
aic_df <- do.call(rbind.data.frame, aic_df)
aic_df$EM_assumption <- as.factor(aic_df$EM_assumption)
aic_df$OM_F_History <- as.factor(aic_df$OM_F_History)
aic_df$EM_M <- as.factor(aic_df$EM_M)
aic_df$EM_SR <- as.factor(aic_df$EM_SR)
aic_df$OM_R_SD <- as.factor(aic_df$R_sig)
aic_df$OM_NAA_SD <- as.factor(aic_df$NAA_sig)
aic_df$sim <- as.factor(aic_df$sim)
aic_df$OM_sim <- as.factor(aic_df$OM_sim)
aic_df$OM_type <- as.factor(aic_df$OM_type)
aic_df$OM_Obs._Error <- as.factor(aic_df$OM_Obs._Error)
aic_df$OM_q_sig <- as.factor(aic_df$q_sig)
aic_df$OM_q_cor <- as.factor(aic_df$q_cor)
aic_df$OM_M_sig <- as.factor(aic_df$M_sig)
aic_df$OM_M_cor <- as.factor(aic_df$M_cor)
aic_df$OM_Sel_sig <- as.factor(aic_df$Sel_sig)
aic_df$OM_Sel_cor <- as.factor(aic_df$Sel_cor)
aic_df$best_PE <- as.factor(aic_df$best_PE)
aic_df$best_PE_cor <- as.factor(aic_df$best_PE_cor)

vglm_fits <- list()
dev.tables <- list(PRD = list())
factors <- list()
# null <- gam(conv ~ 1, family = binomial, method = "REML", data = temp)
factors[["R"]] <- c("1","EM_M", "EM_SR", "OM_Obs._Error", "OM_F_History", "OM_R_SD")
factors[["R+S"]] <- c("1", "EM_M", "EM_SR", "OM_Obs._Error", "OM_F_History", "OM_R_SD","OM_NAA_SD")
factors[["R+M"]] <- c("1", "EM_M", "EM_SR", "OM_Obs._Error", "OM_F_History","OM_M_sig", "OM_M_cor")
factors[["R+Sel"]] <- c("1", "EM_M", "EM_SR", "OM_Obs._Error", "OM_F_History","OM_Sel_sig", "OM_Sel_cor")
factors[["R+q"]] <- c("1", "EM_M", "EM_SR", "OM_Obs._Error", "OM_F_History","OM_q_sig", "OM_q_cor")

for(type in names(factors)){
  if(type == "R") temp <- subset(aic_df, OM_type == "naa" & OM_NAA_SD == 0)
  if(type == "R+S") temp <- subset(aic_df, OM_type == "naa" & OM_NAA_SD != 0)
  if(type == "R+M") temp <- subset(aic_df, OM_type == "M")
  if(type == "R+Sel") temp <- subset(aic_df, OM_type == "Sel")
  if(type == "R+q") temp <- subset(aic_df, OM_type == "q")
  temp$best_PE <- as.factor(temp$best_PE)
  temp$best_PE_cor <- as.factor(temp$best_PE_cor)
  vglm_fits[[type]] <- list()
  dev.tables[[type]] <- list()
  for(i in factors[[type]]){
    vglm_fits[[type]][[i]] <- vglm(as.formula(paste("best_PE", "~", i)), family = multinomial, data = temp)
  }
  sapply(vglm_fits[[type]][factors[[type]][-1]], \(x) lrtest(vglm_fits[[type]][[factors[[type]][1]]], x)@Body[[5]][2])
  vglm_fits[[type]][["all"]] <- vglm(as.formula(paste("best_PE", "~", paste(factors[[type]],collapse = "+"))), family = multinomial, data = temp)
  vglm_fits[[type]][["all2"]] <- vglm(as.formula(paste("best_PE", "~ (", paste(factors[[type]][-1],collapse = "+"), ")^2")), family = multinomial, data = temp)
  vglm_fits[[type]][["all3"]] <- vglm(as.formula(paste("best_PE", "~ (", paste(factors[[type]][-1],collapse = "+"), ")^3")), family = multinomial, data = temp)
  #precent reduction in deviance
  dev.tables[["PRD"]][[type]] <- sapply(vglm_fits[[type]], \(x) 1 - x@criterion$deviance/vglm_fits[[type]][[1]]@criterion$deviance)
}

All.facs <- c("EM_M", "EM_SR", "OM_Obs._Error", "OM_F_History","OM_R_SD","OM_NAA_SD","OM_M_sig", "OM_M_cor","OM_Sel_sig", "OM_Sel_cor","OM_q_sig", "OM_q_cor")
PRD.table <- matrix(NA, length(All.facs)+3,5)
colnames(PRD.table) <- c("R","R+S","R+M","R+Sel","R+q")
rnames <- gsub("_", " ", All.facs, fixed = TRUE)
rnames <- gsub("assumption", "Assumption", rnames, fixed = TRUE)
rnames <- gsub("PE", "Process Error", rnames, fixed = TRUE)
rnames <- gsub("F ", "$F$ ", rnames, fixed = TRUE)
rnames <- gsub("EM M", "EM $M$ Assumption", rnames, fixed = TRUE)
rnames <- gsub("EM SR", "EM SR Assumption", rnames, fixed = TRUE)
rnames <- gsub("NAA SD", "$\\sigma_{2+}$ ", rnames, fixed = TRUE)
rnames <- gsub("R SD", "$\\sigma_R$", rnames, fixed = TRUE)
rnames <- gsub("q sig", "$\\sigma_q$", rnames, fixed = TRUE)
rnames <- gsub("M sig", "$\\sigma_M$", rnames, fixed = TRUE)
rnames <- gsub("Sel sig", "$\\sigma_{Sel}$", rnames, fixed = TRUE)
rnames <- gsub("M cor", "$\\rho_R$", rnames, fixed = TRUE)
rnames <- gsub("Sel cor", "$\\rho_{Sel}$", rnames, fixed = TRUE)
rnames <- gsub("q cor", "$\\rho_q$", rnames, fixed = TRUE)
rownames(PRD.table) <- c(rnames, "All factors", "+ All Two Way", "+ All Three Way")
for(i in c("R","R+S","R+M","R+Sel","R+q")) PRD.table[,i] <- c(dev.tables[["PRD"]][[i]][match(All.facs, names(dev.tables[["PRD"]][[i]]))], dev.tables[["PRD"]][[i]][paste0("all",c("",2:3))])
x <- PRD.table
x[] <- format(round(100*x,2), nsmall = 2)
dim(x)
x[which(is.na(as.numeric(x)))] <- "--"
x[which(as.numeric(x) == 0)] <- "< 0.01"
x <- latex(x, file = here("Project_0","manuscript","AIC_PE_PRD_table.tex"), 
  table.env = FALSE, col.just = rep("r", dim(x)[2]), rowlabel = "Factor", rowlabel.just = "l")#, rowname = NULL)


#########################################
#Regression/classification trees
# install.packages("rpart.plot")
library(rpart)
#modified rpart.plot package to use plotmath in split labels...
#library(rpart.plot)
# pkgload::load_all("c:/work/rpart.plot")
split.fun <- function(type = "R") {
  # replace commas with spaces (needed for strwrap)
  if(!type %in% c("R","R+S")){
    fn <-function(x, labs, digits, varlen, faclen){
      labs <- gsub("_", " ", labs, fixed = TRUE)
      labs <- gsub(" = ", "==", labs, fixed = TRUE)
      labs <- gsub(",", "*','*", labs, fixed = TRUE)
      others <- c("R+M","R+Sel","R+q")[c("R+M","R+Sel","R+q") != type]
      for(i in others) labs <- gsub(paste(i,"(iid)"), i, labs, fixed = TRUE)
      labs <- gsub(" (iid)","(iid)", labs, fixed = TRUE)
      labs <- gsub(" (AR1)","(AR1)", labs, fixed = TRUE)
      labs <- gsub("q sig", "sigma[q]", labs, fixed = TRUE)
      labs <- gsub(" M sig", " sigma[M]", labs, fixed = TRUE)
      labs <- gsub("Sel sig", "sigma[Sel]", labs, fixed = TRUE)
      labs <- gsub("M cor", "rho[M]", labs, fixed = TRUE)
      labs <- gsub("Sel cor", "rho[Sel]", labs, fixed = TRUE)
      labs <- gsub("q cor", "rho[q]", labs, fixed = TRUE)
      labs <- gsub(" M", "*phantom(0)*italic(M)", labs, fixed = TRUE) #italic M
      labs <- gsub("PE", "Process*phantom(0)*Error", labs, fixed = TRUE)
      labs <- gsub("F ", "italic(F)*phantom(0)*", labs, fixed = TRUE)
      labs <- gsub("H->MSY", "2.5*italic(F)[MSY]%->%phantom(0)*italic(F)[MSY]", labs, fixed = TRUE)
      labs <- gsub("MSY", "italic(F)[MSY]", labs, fixed = TRUE)
      labs <- gsub("EM SR", "EM*phantom(0)*SR*phantom(0)*Assumption", labs, fixed = TRUE)
      labs <- gsub("NAA SD", "sigma['2+']", labs, fixed = TRUE) #sigma_2+
      labs <- gsub("R SD", "sigma['R']", labs, fixed = TRUE) #sigma_R
      labs <- gsub(" ", "*phantom(0)*", labs, fixed = TRUE)
      labs
    }
  } else {
    fn <-function(x, labs, digits, varlen, faclen){
      labs <- gsub("_", " ", labs, fixed = TRUE)
      labs <- gsub(" (iid)", "", labs, fixed = TRUE)
      labs <- gsub(" = ", "==", labs, fixed = TRUE)
      labs <- gsub(",", "*','*", labs, fixed = TRUE)
      labs <- gsub("q sig", "sigma[q]", labs, fixed = TRUE)
      labs <- gsub(" M sig", "sigma[M]", labs, fixed = TRUE)
      labs <- gsub("Sel sig", "sigma[Sel]", labs, fixed = TRUE)
      labs <- gsub("M cor", "rho[M]", labs, fixed = TRUE)
      labs <- gsub("Sel cor", "rho[Sel]", labs, fixed = TRUE)
      labs <- gsub("q cor", "rho[q]", labs, fixed = TRUE)
      labs <- gsub(" M", "*phantom(0)*italic(M)", labs, fixed = TRUE) #italic M
      labs <- gsub("PE", "Process*phantom(0)*Error", labs, fixed = TRUE)
      labs <- gsub("F ", "italic(F)*phantom(0)*", labs, fixed = TRUE)
      labs <- gsub("H->MSY", "2.5*italic(F)[MSY]%->%phantom(0)*italic(F)[MSY]", labs, fixed = TRUE)
      labs <- gsub("MSY", "italic(F)[MSY]", labs, fixed = TRUE)
      labs <- gsub("EM SR", "EM*phantom(0)*SR*phantom(0)*Assumption", labs, fixed = TRUE)
      labs <- gsub("NAA SD", "sigma['2+']", labs, fixed = TRUE) #sigma_2+
      labs <- gsub("R SD", "sigma['R']", labs, fixed = TRUE) #sigma_R
      labs <- gsub(" ", "*phantom(0)*", labs, fixed = TRUE)
      labs
    }
  }
  return(fn)
}
# plot.prune(full.trees[["R+q"]],cp = 0, type = "R+q", tweak = 1.2, mar = c(0,0,5,0))

node.fun <- function(x, labs, digits, varlen)
{
  labs <- sapply(1:NROW(x$frame$yval2), \(y) {
    levs <- attr(x, "ylevels")
    nlevs <- length(levs)
    vals <- x$frame$yval2[y,1+nlevs + 1:nlevs]
    ind <- which(vals>0.1)
    lab <- paste0(levs[ind], ": ", format(round(vals[ind],2), nsmall = 2), collapse = ", ")
    if(length(ind)<length(vals)) lab <- paste0(lab, ", Others < 0.1")
    lab
  }) 
  paste0(labs, "\n", x$frame$n)
}

#modified from part package to allow more flexibility in pruning in plots
prune.rpart <- function(tree, cp, factor = "complexity", ...)
{
    ff <- tree$frame
    id <- as.integer(row.names(ff))
    toss <- id[ff[[factor]] <= cp & ff$var != "<leaf>"] #not a leaf
    if (length(toss) == 0L) return(tree)   # all the tree is retained
    newx <- rpart:::snip.rpart(tree, toss)
    newx
}

plot.prune<- function(mod, extra = 7, cp, type, factor = "complexity", ...) {
  rpart.plot(prune.rpart(mod, cp = cp, factor = factor),
    yesno=FALSE,
    type=4, 
    clip.right.labs=TRUE, 
    xcompact = FALSE,
    ycompact = FALSE,
    extra = extra, 
    node.fun = node.fun, 
    split.fun = split.fun(type), 
    box.palette = "RdGn", 
    #box.palette = "Greens", 
    branch = 0.2, 
    fallen.leaves = FALSE, ...)
}

full.trees <- list()
for(type in names(factors)){

  if(type == "R") {
    temp_df <- subset(aic_df, OM_type == "naa" & OM_NAA_SD == 0)
    for(resp in c("best_PE","best_PE_cor")){
      levs <- levels(factor(temp_df[[resp]]))
      temp_df[[resp]] <- factor(temp_df[[resp]], levels = c(levs[levs != type],type))
    }
  }
  if(type == "R+S") {
    temp_df <- subset(aic_df, OM_type == "naa" & OM_NAA_SD != 0)
    for(resp in c("best_PE","best_PE_cor")){
      levs <- levels(factor(temp_df[[resp]]))
      temp_df[[resp]] <- factor(temp_df[[resp]], levels = c(levs[levs != type],type))
    }
  }
  if(type %in% c("R+M","R+q", "R+Sel")){
    if(type == "R+M") temp_df <- subset(aic_df, OM_type == "M")
    if(type == "R+Sel") temp_df <- subset(aic_df, OM_type == "Sel")
    if(type == "R+q") temp_df <- subset(aic_df, OM_type == "q")
    resp <- "best_PE"
    levs <- levels(factor(temp_df[[resp]]))
    temp_df[[resp]] <- factor(temp_df[[resp]], levels = c(levs[levs != type],type))
  }
  form <- as.formula(paste("best_PE ~", paste(factors[[type]], collapse = "+")))# (EM_assumption + OM_R_SD + OM_Obs._Error + OM_F_History)
  full.trees[[type]] <- rpart(form, data=temp_df, method = "class", control=rpart.control(cp=0, xval = 100))
  full.trees[[type]]$correct_level <- type
}
pkgload::load_all(file.path(here(),"../../rpart.plot"))
# pkgload::load_all("/home/tmiller/FSAM_research/aug_backup/work/rpart.plot")

cairo_pdf(here("Project_0","manuscript", paste0("AIC_PE_classification_plots.pdf")), width = 30*2/3, height = 15*2/3)
# x <- matrix(c(1,2,3,4), 2, 2, byrow = TRUE)
anova.palette.sd <- 0.25
x <- matrix(c(1,2,2,3,3,0,4,4,5,5), 2, 5, byrow = TRUE)
layout.x <- layout(x)
par(oma = c(0,2,0,0))
plot.prune(full.trees[["R"]],cp = 0, type = "R", tweak = 1.4, mar = c(0,0,5,0))
mtext("R OMs", side = 3, line = 0, cex = 2)
plot.prune(full.trees[["R+S"]],cp = 0.1, type = "R+S", tweak = 1.4, mar = c(0,0,5,5), legend.x = NA)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
plot.prune(full.trees[["R+M"]],cp = 0.04, type = "R+M", tweak = 1.4, mar = c(0,0,5,5), legend.x = NA)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
plot.prune(full.trees[["R+Sel"]],cp = 0.1, type = "R+Sel", tweak = 1.4, mar = c(0,0,5,5), legend.x = NA)
mtext("R+Sel OMs", side = 3, line = 0, cex = 2)
plot.prune(full.trees[["R+q"]],cp = 800, type = "R+q", factor = "dev", tweak = 1.4, mar = c(0,0,5,5), legend.x = NA)
mtext("R+q OMs", side = 3, line = 0, cex = 2)
dev.off()
