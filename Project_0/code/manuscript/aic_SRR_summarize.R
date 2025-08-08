library(here)
library(dplyr)
library(tidyr)
library(Hmisc)

aic_fn <- function(x, all_aic, df.oms, df.ems, rec_mod, M_est) {
  if(!is.null(df.oms$NAA_sig)){
    if(is.na(df.oms$NAA_sig[x])) re_config <- "rec"
    else re_config <- "rec+1"  
    em_ind <- which(df.ems$SR_model == rec_mod & df.ems$M_est == M_est & df.ems$re_config == re_config)
  }
  if(!is.null(df.oms$M_sig)){
    if(df.oms$M_cor[x] == 0) re_config <- "iid"
    else re_config <- "ar1_y"  
    em_ind <- which(df.ems$SR_model == rec_mod & df.ems$M_est == M_est & df.ems$M_re_cor == re_config)
  }
  if(!is.null(df.oms$Sel_sig)){
    if(df.oms$Sel_cor[x] == 0) re_config <- "iid"
    else re_config <- "ar1_y"  
    em_ind <- which(df.ems$SR_model == rec_mod & df.ems$M_est == M_est & df.ems$sel_re_cor == re_config)
  }
  if(!is.null(df.oms$q_sig)){
    if(df.oms$q_cor[x] == 0) re_config <- "iid"
    else re_config <- "ar1"  
    em_ind <- which(df.ems$SR_model == rec_mod & df.ems$M_est == M_est & df.ems$q_re_cor == re_config)
  }
  return(all_aic[[x]][em_ind,]) #single EM/row
}

df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))
# em_ind <- c(5:20,25:28)
# use.df.ems <- df.ems[em_ind,]
types <- c("naa", "M", "Sel", "q")
types.plt <- c("R, R+S","R+M", "R+Sel", "R+q")
est_conds <- cbind(c(FALSE,TRUE))
em_inds <- cbind(1:20, 5:24, c(5:20,25:28), c(5:20, 29:32))
pred_dfs <- list()
obs_dfs <- list()
for(i in 1:4) {
  all_aic <- readRDS(file = here("Project_0","results", paste0("all_", types[i], "_aic_results.RDS")))
  df.oms = readRDS(here("Project_0","inputs", paste0("df.", ifelse(i==1, "", paste0(types[i],".")),"oms.RDS")))
  use.df.ems <- df.ems[em_inds[,i],]
  df.oms$OM <- 1:NROW(df.oms)
  sd_log_SSB <- readRDS(file = here("Project_0","results", paste0("all_", types[i], "_sd_log_SSB.RDS")))
  sd_log_SSB <- t(matrix(unname(unlist(sd_log_SSB)), nrow = 100))
#  colnames(sd_log_SSB) <- paste0("sim", 1:100)
  df <- cbind(df.oms, sd_log_SSB)
  pred_df <- cbind(df.oms, t(sapply(1:NROW(df.oms), \(x) seq(min(sd_log_SSB[x,]), max(sd_log_SSB[x,]),length.out = 100))))
  df <- df %>% pivot_longer(cols = paste0(1:100), names_to = "sim", values_to = "sd_log_SSB")
  pred_df <- pred_df %>% pivot_longer(cols = paste0(1:100), names_to = NULL, values_to = "sd_log_SSB")
  df$aic_BH_ME <- c(sapply(1:length(all_aic), aic_fn, all_aic = all_aic, df.oms = df.oms, df.ems = use.df.ems, 3, TRUE))
  df$aic_R_ME <- c(sapply(1:length(all_aic), aic_fn, all_aic = all_aic, df.oms = df.oms, df.ems = use.df.ems, 2, TRUE))
  df$BH_best_ME <- as.integer(df$aic_BH_ME < df$aic_R_ME)
  df$aic_BH_MF <- c(sapply(1:length(all_aic), aic_fn, all_aic = all_aic, df.oms = df.oms, df.ems = use.df.ems, 3, FALSE))
  df$aic_R_MF <- c(sapply(1:length(all_aic), aic_fn, all_aic = all_aic, df.oms = df.oms, df.ems = use.df.ems, 2, FALSE))
  df$BH_best_MF <- as.integer(df$aic_BH_MF < df$aic_R_MF)
  df <- df %>% pivot_longer(cols = c("BH_best_MF","BH_best_ME"), names_to = "EM_M", values_to = "BH_best", names_prefix = "BH_best_") %>% as.data.frame
  if(i==1) {
    df$NAA_sig[which(is.na(df$NAA_sig))] <- 0
    df$OM_NAA_sigma <- df$NAA_sig
    df$OM_R_sigma <- df$R_sig
  }
  if(i ==2){
    df$OM_M_sigma <- df$M_sig
    df$OM_M_rho <- df$M_cor
  }
  if(i == 3){
    df$OM_Sel_sigma <- df$Sel_sig
    df$OM_Sel_rho <- df$Sel_cor
  }
  if(i == 4){
    df$OM_q_sigma <- df$q_sig
    df$OM_q_rho <- df$q_cor
  }
  df$log_sd_log_SSB <- log(df$sd_log_SSB)
  df <- df %>% mutate(OM_F_History = recode(Fhist,
        "H-MSY"  = "H->MSY",
        "MSY" = "MSY"))
  df <- df %>% mutate(OM_Obs._Error = recode(obs_error,
      "L" = "Low",
      "H" = "High"))
  df <- df %>% mutate(EM_M = recode(EM_M,
      "ME" = "Estimated",
      "MF" = "Known"))
  df <- df %>% as.data.frame
  facs <- c("OM", "OM_F_History","OM_NAA_sigma", "OM_R_sigma", "OM_M_sigma", "OM_M_rho", "OM_Sel_sigma", "OM_Sel_rho", "OM_q_sigma", "OM_q_rho", "OM_Obs._Error","EM_M")
  fac_names <- names(df)[names(df) %in% facs]
  df[fac_names] <- lapply(df[fac_names], as.factor)
  obs_dfs[[types[i]]] <- df
}

glm_fits <- list()
dev.tables <- list(PRD = list())
factors <- list()
# null <- gam(conv ~ 1, family = binomial, method = "REML", data = temp)
factors[["R"]] <- c("1", "OM_R_sigma","EM_M", "OM_Obs._Error", "OM_F_History", "log_sd_log_SSB")
factors[["R+S"]] <- c("1", "OM_R_sigma","EM_M", "OM_Obs._Error", "OM_F_History","OM_NAA_sigma", "log_sd_log_SSB")
factors[["R+M"]] <- c("1", "EM_M","OM_Obs._Error", "OM_F_History","OM_M_sigma", "OM_M_rho", "log_sd_log_SSB")
factors[["R+Sel"]] <- c("1", "EM_M","OM_Obs._Error", "OM_F_History","OM_Sel_sigma", "OM_Sel_rho", "log_sd_log_SSB")
factors[["R+q"]] <- c("1", "EM_M","OM_Obs._Error", "OM_F_History","OM_q_sigma", "OM_q_rho", "log_sd_log_SSB")

for(type in names(factors)){
  if(type == "R") temp <- subset(obs_dfs[["naa"]], OM_NAA_sigma == 0)
  if(type == "R+S") temp <- subset(obs_dfs[["naa"]], OM_NAA_sigma != 0)
  if(type == "R+M") temp <- obs_dfs[["M"]]
  if(type == "R+Sel") temp <- obs_dfs[["Sel"]]
  if(type == "R+q") temp <- obs_dfs[["q"]]
  glm_fits[[type]] <- list()
  dev.tables[[type]] <- list()
  for(i in factors[[type]]){
    glm_fits[[type]][[i]] <- glm(as.formula(paste("BH_best", "~", i)), family = binomial, data = temp)
  }
  sapply(glm_fits[[type]][factors[[type]][-1]], \(x) anova(x, test = "LRT")[[5]][2])
  glm_fits[[type]][["all"]] <- glm(as.formula(paste("BH_best", "~", paste(factors[[type]],collapse = "+"))), family = binomial, data = temp)
  glm_fits[[type]][["all2"]] <- glm(as.formula(paste("BH_best", "~ (", paste(factors[[type]][-1],collapse = "+"), ")^2")), family = binomial, data = temp)
  glm_fits[[type]][["all3"]] <- glm(as.formula(paste("BH_best", "~ (", paste(factors[[type]][-1],collapse = "+"), ")^3")), family = binomial, data = temp)

  #precent reduction in deviance
  dev.tables[["PRD"]][[type]] <- sapply(glm_fits[[type]][factors[[type]]], \(x) 1 - x$deviance/glm_fits[[type]][[1]]$null.deviance)
}

interactions.dev.table <- sapply(names(factors), \(x) 1 - c(glm_fits[[x]][["all"]]$deviance,glm_fits[[x]][["all2"]]$deviance,glm_fits[[x]][["all3"]]$deviance)/glm_fits[[x]][[1]]$null.deviance )
x <- as.data.frame(round(100*interactions.dev.table,2))
x <- cbind(Model = c("No Interactions", "+ All Two Way", "+ All Three Way"), x)

All.facs <- c("EM_M", "OM_Obs._Error", "OM_F_History","OM_R_sigma","OM_NAA_sigma","OM_M_sigma", "OM_M_rho","OM_Sel_sigma", "OM_Sel_rho","OM_q_sigma", "OM_q_rho", "log_sd_log_SSB")
PRD.table <- matrix(NA, length(All.facs),5)
colnames(PRD.table) <- c("R","R+S","R+M","R+Sel","R+q")
rnames <- gsub("_", " ", All.facs, fixed = TRUE)
rnames <- gsub("F ", "$F$ ", rnames, fixed = TRUE)
rnames <- gsub("EM M", "EM $M$ Assumption", rnames, fixed = TRUE)
rnames <- gsub("log sd log SSB", "$\\log\\left(SD_\\text{SSB}]\\right)$", rnames, fixed = TRUE)
rnames <- gsub("NAA sigma", "$\\sigma_{2+}$ ", rnames, fixed = TRUE)
rnames <- gsub("R sigma", "$\\sigma_R$", rnames, fixed = TRUE)
rnames <- gsub("q sigma", "$\\sigma_q$", rnames, fixed = TRUE)
rnames <- gsub("M sigma", "$\\sigma_M$", rnames, fixed = TRUE)
rnames <- gsub("Sel sigma", "$\\sigma_{Sel}$", rnames, fixed = TRUE)
rnames <- gsub("M rho", "$\\rho_R$", rnames, fixed = TRUE)
rnames <- gsub("Sel rho", "$\\rho_{Sel}$", rnames, fixed = TRUE)
rnames <- gsub("q rho", "$\\rho_q$", rnames, fixed = TRUE)

rownames(PRD.table) <- rnames
for(i in c("R","R+S","R+M","R+Sel","R+q")) PRD.table[,i] <- dev.tables[["PRD"]][[i]][match(All.facs, names(dev.tables[["PRD"]][[i]]))]
y <- interactions.dev.table
rownames(y) <- c("All factors", "+ All Two Way", "+ All Three Way")
PRD.table <- rbind(PRD.table,y)

x[] <- format(round(100*PRD.table,2), nsmall = 2)
dim(x)
x[which(is.na(as.numeric(x)))] <- "--"
x[which(as.numeric(x) == 0)] <- "< 0.01"
x <- latex(x, file = here("Project_0","manuscript","AIC_SRR_PRD_table.tex"), 
  table.env = FALSE, col.just = rep("r", dim(x)[2]), rowlabel = "Factor", rowlabel.just = "l")#, rowname = NULL)


#########################################
#Regression/classification trees
# install.packages("rpart.plot")
library(rpart)
#modified rpart.plot package to use plotmath in split labels...
#library(rpart.plot)
pkgload::load_all("c:/work/rpart.plot")
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
      labs <- gsub(" M", "*phantom(0)*italic(M)", labs, fixed = TRUE) #italic M
      labs <- gsub("F ", "italic(F)*phantom(0)*", labs, fixed = TRUE)
      labs <- gsub("H->MSY", "2.5*italic(F)[MSY]%->%phantom(0)*italic(F)[MSY]", labs, fixed = TRUE)
      labs <- gsub("MSY", "italic(F)[MSY]", labs, fixed = TRUE)
      labs <- gsub("R sigma", "sigma[R]", labs, fixed = TRUE) #sigma_2+
      labs <- gsub("NAA sigma", "sigma['2+']", labs, fixed = TRUE) #sigma_2+
      labs <- gsub("log sd log SSB", "log(SD[SSB])", labs, fixed = TRUE)
      labs <- gsub(" ", "*phantom(0)*", labs, fixed = TRUE)
      labs <- gsub("*<", "<", labs, fixed = TRUE)
      labs <- gsub("*>", ">", labs, fixed = TRUE)
      labs <- gsub("=*", "=", labs, fixed = TRUE)
      labs <- gsub(">*", ">", labs, fixed = TRUE)
      labs <- gsub("<*", "<", labs, fixed = TRUE)
      labs <- gsub("^\\*", "", labs)
      #can't do this with math expressions
      # for(i in 1:length(labs)) {
      #   # split labs[i] into multiple lines
      #   # labs[i] <- paste(strwrap(labs[i], width = 15), collapse = "\n")
      # }
      labs
    }
  } else {
    fn <-function(x, labs, digits, varlen, faclen){
      labs <- gsub("_", " ", labs, fixed = TRUE)
      labs <- gsub(" (iid)", "", labs, fixed = TRUE)
      labs <- gsub(" = ", "==", labs, fixed = TRUE)
      labs <- gsub(",", "*','*", labs, fixed = TRUE)
      labs <- gsub(" M", "*phantom(0)*italic(M)", labs, fixed = TRUE) #italic M
      labs <- gsub("F ", "italic(F)*phantom(0)*", labs, fixed = TRUE)
      labs <- gsub("H->MSY", "2.5*italic(F)[MSY]%->%phantom(0)*italic(F)[MSY]", labs, fixed = TRUE)
      labs <- gsub("MSY", "italic(F)[MSY]", labs, fixed = TRUE)
      labs <- gsub("R sigma", "sigma[R]", labs, fixed = TRUE) #sigma_2+
      labs <- gsub("NAA sigma", "sigma['2+']", labs, fixed = TRUE) #sigma_2+
      labs <- gsub("log sd log SSB", "log(SD[SSB])", labs, fixed = TRUE)
      labs <- gsub(" ", "*phantom(0)*", labs, fixed = TRUE)
      labs <- gsub("*<", "<", labs, fixed = TRUE)
      labs <- gsub("*>", ">", labs, fixed = TRUE)
      labs <- gsub("=*", "=", labs, fixed = TRUE)
      labs <- gsub(">*", ">", labs, fixed = TRUE)
      labs <- gsub("<*", "<", labs, fixed = TRUE)
      labs <- gsub("^\\*", "", labs)
      #can't do this with math expressions
      # for(i in 1:length(labs)) {
      #   split labs[i] into multiple lines
      #   labs[i] <- paste(strwrap(labs[i], width = 15), collapse = "\n")
      # }
      labs
    }
  }
  return(fn)
}

node.fun <- function(x, labs, digits, varlen)
{
  paste0("Prop. Correct = ", format(round(x$frame$yval2[,5],3), nsmall = 3), "\nn = ", x$frame$n)
}

#modified from part package to allow more flexibility in pruning in plots
prune.rpart <- function(tree, cp, factor = "complexity", ...)
{
    ff <- tree$frame
    id <- as.integer(row.names(ff))
    toss <- id[ff[[factor]] <= cp & ff$var != "<leaf>"] #not a leaf
    if (length(toss) == 0L) return(tree)   # all the tree is retained
    newx <- rpart:::snip.rpart(tree, toss)
    ## Now cut down the CP table
    # temp <- pmax(tree$cptable[, 1L], cp)
    # keep <- match(unique(temp), temp)
    # newx$cptable <- tree$cptable[keep, , drop = FALSE]
    # newx$cptable[length(keep), 1L] <- cp
    # # Reset the variable importance
    # newx$variable.importance <- rpart:::importance(newx)
    newx
}

plot.prune<- function(mod,cp, type, factor = "complexity", ...) {
  rpart.plot(prune.rpart(mod, cp = cp, factor = factor),
    yesno=FALSE,
    type=4, 
    clip.right.labs=TRUE, 
    xcompact = FALSE,
    ycompact = FALSE,
    extra = 7, 
    node.fun = node.fun, 
    split.fun = split.fun(type), 
    box.palette = "RdGn", 
    branch = 0.2, 
    fallen.leaves = FALSE, ...)
}

full.trees <- list()
for(type in names(factors)){

  if(type == "R") temp <- subset(obs_dfs[["naa"]], OM_NAA_sigma == 0)
  if(type == "R+S") temp <- subset(obs_dfs[["naa"]], OM_NAA_sigma != 0)
  if(type == "R+M") temp <- obs_dfs[["M"]]
  if(type == "R+Sel") temp <- obs_dfs[["Sel"]]
  if(type == "R+q") temp <- obs_dfs[["q"]]
  temp$BH_best_fac <- as.factor(temp$BH_best)
  
  form <- as.formula(paste("BH_best_fac ~", paste(factors[[type]], collapse = "+")))# (EM_PE + OM_R_SD + EM_M + EM_SR + OM_Obs._Error + OM_F_History)
  full.trees[[type]] <- rpart(form, data=temp, method = "class", control=rpart.control(cp=0, xval = 100), roundint = FALSE)
}
plot.prune(full.trees[["R"]], cp = 0.1, type = "R", roundint = FALSE)

cairo_pdf(here("Project_0","manuscript", paste0("AIC_SRR_classification_plots.pdf")), width = 30*2/3, height = 20*2/3)
x <- matrix(c(1,1,2,2,3,3,0,4,4,5,5,0), 2, 6, byrow = TRUE)
layout.x <- layout(x) 
par(mar = c(0,0,5,0), oma = c(0,0,0,0))
plot.prune(full.trees[["R"]],0.1, "R", tweak = 1.2, mar = c(0,0,5,0), roundint = FALSE)
mtext("R OMs", side = 3, line = 0, cex = 2)
plot.prune(full.trees[["R+S"]],0.03, "R+S", tweak = 1.2, mar = c(0,0,5,0), roundint = FALSE)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
x <- prune(full.trees[["R+M"]],0.005)
plot.prune(x,cp = 100, factor = "dev", type = "R+M", tweak= 1.2, mar = c(0,0,5,0), roundint= FALSE)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
plot.prune(full.trees[["R+Sel"]],0.005, "R+Sel", tweak = 1.2, mar = c(0,0,5,0), roundint= FALSE)
mtext("R+Sel OMs", side = 3, line = 0, cex = 2)
x <- prune(full.trees[["R+q"]],0.01)
plot.prune(x,100, factor = "dev", type = "R+q", tweak = 1.2, mar = c(0,0,5,0), roundint= FALSE)
mtext("R+q OMs", side = 3, line = 0, cex = 2)
dev.off()
