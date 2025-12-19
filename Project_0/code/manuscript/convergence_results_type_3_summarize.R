library(here)
library(dplyr)
library(tidyr)
library(mgcv)
library(gratia)
library(Hmisc)
#modified rpart.plot package to use plotmath in split labels...
pkgload::load_all(file.path(here(),"../../rpart.plot"))


df.oms = list(naa = readRDS(here("Project_0", "inputs", "df.oms.RDS")))
df.oms$M = readRDS(here("Project_0", "inputs", "df.M.oms.RDS"))
df.oms$Sel = readRDS(here("Project_0", "inputs", "df.Sel.oms.RDS"))
df.oms$q = readRDS(here("Project_0", "inputs", "df.q.oms.RDS"))
om.em.rows <- list(naa = 1:20, M = 5:24, Sel = c(5:20, 25:28), q = c(5:20, 29:32))
df.ems = readRDS(here("Project_0", "inputs", "df.ems.RDS"))

#convergence information:
#1: 0/1 model finished optimizing
#2: nlminb convergence flag
#3: max gradient value
#4: number of NaNs in SEs for parameters, 0 = good invertible hessian
#5: maximum non-NaN SE estimate

om_conv_res <- list(naa = readRDS(here("Project_0", "results", "naa_om_convergence_results.RDS")))
# naa_om_conv_res <- readRDS(here("Project_0", "results", "naa_om_convergence_results.RDS"))
om_conv_res$M <- readRDS(here("Project_0", "results", "M_om_convergence_results.RDS"))
om_conv_res$Sel <- readRDS(here("Project_0", "results", "Sel_om_convergence_results.RDS"))
om_conv_res$q <- readRDS(here("Project_0", "results", "q_om_convergence_results.RDS"))
est_config <- c("M known, No SR", "M estimated, No SR", "M known, SR", "M estimated, SR")
EM_process_error <- c("R","R+S","R+M (iid)","R+M (AR1)","R+Sel (iid)","R+Sel (AR1)","R+q (iid)","R+q (AR1)")

conv_fn <- function(all, em_rows, om_type, df.ems, om_ind = NULL, est_ind = NULL, type = 3){
  
  #df.ems. is the subset of df.ems for the current om_type
  #est_ind indexes ROWs of df.ems.
  #20 EMs for each OM type, but the EMs are not all the same.
  #df.ems has all of them
  df.ems. <- df.ems[em_rows[[om_type]],]
  if(is.null(est_ind)) est_ind <- 1:NROW(df.ems.)
  
  all <- all[[om_type]]
  if(is.null(om_ind)) om_ind <- 1:length(all)
  all <- all[om_ind]
  
  #all has results for a given om type 
  res_by_om <- lapply(1:length(all), \(om_i){
    #length(res) = nsims for a given simulation
    res_by_sim <- lapply(1:length(all[[om_i]]), \(sim_j) {
      #w has results for a given sim dim(w) = n_ems x 5 (columns have convergence info)
      res_by_sim_em <- lapply(est_ind, \(k){
        simres <- all[[om_i]][[sim_j]][k,]
        if(type == 1) conv <- !is.na(simres[1]) #catastrophic failure
        if(type == 2) conv <- !is.na(simres[2]) & simres[2] == 0  #nlminb convergence flag
        if(type == 3) conv <- !is.na(simres[4]) & simres[4] < 1e-6 #gradient
        if(type == 4) conv <- !is.na(simres[3]) & simres[3] == 0 #hessian
        if(type == 5) conv <- !is.na(simres[5]) & simres[3] == 0 & simres[5] < 100 #all SEs less than 100
        return(cbind(est_ind = k, conv = as.integer(conv)))
      })
      out <- do.call(rbind, res_by_sim_em)
      return(cbind(sim = sim_j, out))
    })
    out <- do.call(rbind, res_by_sim)
    return(cbind(om = om_i, out))
  })
  out <- do.call(rbind, res_by_om)
  return(out)
}

conv_df <- list()
for(i in names(om.em.rows)) {
  conv_df[[i]] <- conv_fn(om_conv_res, em_rows = om.em.rows, om_type = i, df.ems)
  conv_df[[i]] <- cbind(df.oms[[i]][conv_df[[i]][,"om"],], df.ems[om.em.rows[[i]][conv_df[[i]][,"est_ind"]],], conv_df[[i]])
  conv_df[[i]]$Model <- paste0(i,"_",conv_df[[i]]$Model)
  conv_df[[i]]$OM <- paste0(i,"_",conv_df[[i]]$om)
  conv_df[[i]]$OM_sim  <- paste0(i, "_", (conv_df[[i]][["om"]]-1)*100 + conv_df[[i]][["sim"]])
  conv_df[[i]]$OM_type <- i
  print(i)
  print(head(conv_df[[i]]))
  if(i == "naa"){
    conv_df[[i]]$NAA_sig[is.na(conv_df[[i]]$NAA_sig)] <- 0
    other_oms <- names(om.em.rows)[names(om.em.rows) != i]
    for(j in other_oms){
      conv_df[[i]][[paste0(j, "_sig")]] <- NA
      conv_df[[i]][[paste0(j, "_cor")]] <- NA
    }
  } else{
    conv_df[[i]]$NAA_sig <- NA
    conv_df[[i]][[paste0(i, "_sig")]] <- as.factor(conv_df[[i]][[paste0(i, "_sig")]])
    conv_df[[i]][[paste0(i, "_cor")]] <- as.factor(conv_df[[i]][[paste0(i, "_cor")]])
    other_oms <- names(om.em.rows)[!names(om.em.rows) %in% c(i,"naa")]
    for(j in other_oms){
      conv_df[[i]][[paste0(j, "_sig")]] <- NA
      conv_df[[i]][[paste0(j, "_cor")]] <- NA
    }
  }
  conv_df[[i]] <- conv_df[[i]] %>%
    mutate(EM_PE = recode(re_config,
                          "rec" = "R",
                          "rec+1" = "R+S",
                          "M_re" = "R+M",
                          "q_re" = "R+q",
                          "sel_re" = "R+Sel"
    ))
  ind <- which(conv_df[[i]]$M_re_cor == "iid" | conv_df[[i]]$sel_re_cor == "iid" | conv_df[[i]]$q_re_cor == "iid")
  conv_df[[i]]$EM_PE[ind] <- paste0(conv_df[[i]]$EM_PE[ind], " (iid)")
  ind <- which(conv_df[[i]]$M_re_cor == "ar1_y" | conv_df[[i]]$sel_re_cor == "ar1_y" | conv_df[[i]]$q_re_cor == "ar1")
  conv_df[[i]]$EM_PE[ind] <- paste0(conv_df[[i]]$EM_PE[ind], " (AR1)")
  conv_df[[i]] <- conv_df[[i]] %>%
    mutate(Obs._Error = recode(obs_error,
                               "L" = "Low",
                               "H" = "High"
    ))
  conv_df[[i]] <- conv_df[[i]] %>%
    mutate(EM_M = recode(as.character(M_est),
                         "TRUE" = "Estimated",
                         "FALSE" = "Known"
    ))
  conv_df[[i]] <- conv_df[[i]] %>%
    mutate(EM_SR = recode(as.character(SR_model),
                          "2" = "None",
                          "3" = "Estimated"
    ))
  conv_df[[i]] <- conv_df[[i]] %>%
    mutate(F_History = recode(as.character(Fhist),
      "H-MSY"  = "2.5*italic(F)[MSY] %->% italic(F)[MSY]",
      "MSY" = "italic(F)[MSY]"
    ))
  
  conv_df[[i]] <- conv_df[[i]][sort(names(conv_df[[i]]))]
  
}
conv_df <- do.call(rbind.data.frame, conv_df)
conv_df$OM_F_History <- as.factor(conv_df$F_History)
conv_df$OM_R_SD <- as.factor(conv_df$R_sig)
conv_df$OM_NAA_SD <- as.factor(conv_df$NAA_sig)
conv_df$EM_SR <- as.factor(conv_df$EM_SR)
conv_df$EM_M <- as.factor(conv_df$EM_M)
conv_df$EM_PE <- as.factor(conv_df$EM_PE)
conv_df$M_re_cor <- as.factor(conv_df$M_re_cor)
conv_df$q_re_cor <- as.factor(conv_df$q_re_cor)
conv_df$sel_re_cor <- as.factor(conv_df$sel_re_cor)
conv_df$sim <- as.factor(conv_df$sim)
conv_df$OM_sim <- as.factor(conv_df$OM_sim)
conv_df$OM_type <- as.factor(conv_df$OM_type)
conv_df$OM_Obs._Error <- as.factor(conv_df$Obs._Error)
conv_df$OM_q_sig <- as.factor(conv_df$q_sig)
conv_df$OM_q_cor <- as.factor(conv_df$q_cor)
conv_df$OM_M_sig <- as.factor(conv_df$M_sig)
conv_df$OM_M_cor <- as.factor(conv_df$M_cor)
conv_df$OM_Sel_sig <- as.factor(conv_df$Sel_sig)
conv_df$OM_Sel_cor <- as.factor(conv_df$Sel_cor)
conv_df$conv_txt <- as.factor(c("Success","Fail")[match(conv_df$conv,1:0)])

glm_fits <- list()
dev.tables <- list(PRD = list(), PPRD = list())
factors <- list()
# null <- gam(conv ~ 1, family = binomial, method = "REML", data = temp)
factors[["R"]] <- c("1", "EM_PE","OM_R_SD","EM_M","EM_SR", "OM_Obs._Error", "OM_F_History")
factors[["R+S"]] <- c("1", "EM_PE","OM_R_SD","EM_M","EM_SR", "OM_Obs._Error", "OM_F_History","OM_NAA_SD")
factors[["R+M"]] <- c("1", "EM_PE","EM_M","EM_SR", "OM_Obs._Error", "OM_F_History","OM_M_sig", "OM_M_cor")
factors[["R+Sel"]] <- c("1", "EM_PE","EM_M","EM_SR", "OM_Obs._Error", "OM_F_History","OM_Sel_sig", "OM_Sel_cor")
factors[["R+q"]] <- c("1", "EM_PE","EM_M","EM_SR", "OM_Obs._Error", "OM_F_History","OM_q_sig", "OM_q_cor")

for(type in names(factors)){
  if(type == "R") temp <- subset(conv_df, OM_type == "naa" & OM_NAA_SD == 0)
  if(type == "R+S") temp <- subset(conv_df, OM_type == "naa" & OM_NAA_SD != 0)
  if(type == "R+M") temp <- subset(conv_df, OM_type == "M")
  if(type == "R+Sel") temp <- subset(conv_df, OM_type == "Sel")
  if(type == "R+q") temp <- subset(conv_df, OM_type == "q")
  glm_fits[[type]] <- list()
  dev.tables[[type]] <- list()
  for(i in factors[[type]]){
    glm_fits[[type]][[i]] <- glm(as.formula(paste("conv", "~", i)), family = binomial, data = temp)
  }
  sapply(glm_fits[[type]][factors[[type]][-1]], \(x) anova(x, test = "LRT")[[5]][2])
  glm_fits[[type]][["all"]] <- glm(as.formula(paste("conv", "~", paste(factors[[type]],collapse = "+"))), family = binomial, data = temp)
  glm_fits[[type]][["all2"]] <- glm(as.formula(paste("conv", "~ (", paste(factors[[type]][-1],collapse = "+"), ")^2")), family = binomial, data = temp)
  glm_fits[[type]][["all3"]] <- glm(as.formula(paste("conv", "~ (", paste(factors[[type]][-1],collapse = "+"), ")^3")), family = binomial, data = temp)
  
  #precent reduction in deviance
  dev.tables[["PRD"]][[type]] <- sapply(glm_fits[[type]][factors[[type]]], \(x) 1 - x$deviance/glm_fits[[type]][[1]]$null.deviance)
  #precent "possible" reduction in deviance: measure reduction in deviance relative to that with all factors in the model (no interactions)
  dev.tables[["PPRD"]][[type]] <- (glm_fits[[type]][[1]]$null.deviance - sapply(glm_fits[[type]][factors[[type]]], \(x) x$deviance))/(glm_fits[[type]][[1]]$null.deviance - glm_fits[[type]][["all"]]$deviance)
}

interactions.dev.table <- sapply(names(factors), \(x) 1 - c(glm_fits[[x]][["all"]]$deviance,glm_fits[[x]][["all2"]]$deviance,glm_fits[[x]][["all3"]]$deviance)/glm_fits[[x]][[1]]$null.deviance )
x <- as.data.frame(round(100*interactions.dev.table,2))
x <- cbind(Model = c("No Interactions", "+ All Two Way", "+ All Three Way"), x)

# x <- latex(x, file = here("Project_0","manuscript","convergence_glm_interactions_PRD_table.tex"), 
#   table.env = FALSE, col.just = c("l",rep("r", dim(x)[2]-1)), rowname = NULL)


All.facs <- c("EM_PE","EM_M","EM_SR", "OM_Obs._Error", "OM_F_History","OM_R_SD","OM_NAA_SD","OM_M_sig", "OM_M_cor","OM_Sel_sig", "OM_Sel_cor","OM_q_sig", "OM_q_cor")
PRD.table <- matrix(NA, length(All.facs),5)
colnames(PRD.table) <- c("R","R+S","R+M","R+Sel","R+q")
rnames <- gsub("_", " ", All.facs, fixed = TRUE)
rnames <- gsub("PE", "Process Error", rnames, fixed = TRUE)
rnames <- gsub("F ", "$F$ ", rnames, fixed = TRUE)
rnames <- gsub("EM M", "EM $M$ Assumption", rnames, fixed = TRUE)
rnames <- gsub("EM SR", "EM SR Assumption", rnames, fixed = TRUE)
rnames <- gsub("NAA SD", "$\\sigma_{2+}$ ", rnames, fixed = TRUE)
rnames <- gsub("R SD", "$\\sigma_R$", rnames, fixed = TRUE)
rnames <- gsub("q sig", "$\\sigma_q$", rnames, fixed = TRUE)
rnames <- gsub("M sig", "$\\sigma_M$", rnames, fixed = TRUE)
rnames <- gsub("Sel sig", "$\\sigma_{\\text{Sel}}$", rnames, fixed = TRUE)
rnames <- gsub("M cor", "$\\rho_M$", rnames, fixed = TRUE)
rnames <- gsub("Sel cor", "$\\rho_{\\text{Sel}}$", rnames, fixed = TRUE)
rnames <- gsub("q cor", "$\\rho_q$", rnames, fixed = TRUE)
rownames(PRD.table) <- rnames
for(i in c("R","R+S","R+M","R+Sel","R+q")) PRD.table[,i] <- dev.tables[["PRD"]][[i]][match(All.facs, names(dev.tables[["PRD"]][[i]]))]
y <- interactions.dev.table
rownames(y) <- c("All factors", "+ All Two Way", "+ All Three Way")
x <- rbind(PRD.table,y)

x[] <- format(round(100*x,2), nsmall = 2)
dim(x)
x[which(is.na(as.numeric(x)))] <- "--"
x[which(as.numeric(x) == 0)] <- "< 0.01"
x <- latex(x, file = here("Project_0","manuscript","convergence_gradient_PRD_table.tex"), 
           table.env = FALSE, col.just = rep("r", dim(x)[2]), rowlabel = "Factor", rowlabel.just = "l")#, rowname = NULL)


#########################################
#Regression/classification trees
# install.packages("rpart.plot")
library(rpart)
#modified rpart.plot package to use plotmath in split labels...
#library(rpart.plot)
pkgload::load_all(file.path(here(),"../../rpart.plot"))
#pkgload::load_all("c:/work/rpart.plot")
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
      labs <- gsub("PE", "Process*phantom(0)*Error", labs, fixed = TRUE)
      labs <- gsub("F ", "italic(F)*phantom(0)*", labs, fixed = TRUE)
      # labs <- gsub("H->MSY", "2.5*italic(F)[MSY]%->%phantom(0)*italic(F)[MSY]", labs, fixed = TRUE)
      # labs <- gsub("MSY", "italic(F)[MSY]", labs, fixed = TRUE)
      labs <- gsub(" %->% ", "*phantom(0)%->%phantom(0)*", labs, fixed = TRUE)
      labs <- gsub("OM ", "OM*phantom(0)*", labs, fixed = TRUE)
      labs <- gsub("Obs. Error", "Obs.*phantom(0)*Error", labs, fixed = TRUE)
      labs <- gsub("EM Process", "EM*phantom(0)*Process", labs, fixed = TRUE)
      labs <- gsub("EM SR", "EM*phantom(0)*SR*phantom(0)*Assumption", labs, fixed = TRUE)
      labs <- gsub("NAA SD", "sigma['2+']", labs, fixed = TRUE) #sigma_2+
      labs <- gsub(" ", "*phantom(0)*", labs, fixed = TRUE)
      #can't do this with math expressions
      # for(i in 1:length(labs)) {
      # 	# split labs[i] into multiple lines
      # 	# labs[i] <- paste(strwrap(labs[i], width = 15), collapse = "\n")
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
      labs <- gsub("PE", "Process*phantom(0)*Error", labs, fixed = TRUE)
      labs <- gsub("F ", "italic(F)*phantom(0)*", labs, fixed = TRUE)
      # labs <- gsub("H->MSY", "2.5*italic(F)[MSY]%->%phantom(0)*italic(F)[MSY]", labs, fixed = TRUE)
      # labs <- gsub("MSY", "italic(F)[MSY]", labs, fixed = TRUE)
      labs <- gsub(" %->% ", "*phantom(0)%->%phantom(0)*", labs, fixed = TRUE)
      labs <- gsub("OM ", "OM*phantom(0)*", labs, fixed = TRUE)
      labs <- gsub("Obs. Error", "Obs.*phantom(0)*Error", labs, fixed = TRUE)
      labs <- gsub("EM Process", "EM*phantom(0)*Process", labs, fixed = TRUE)
      labs <- gsub("EM SR", "EM*phantom(0)*SR*phantom(0)*Assumption", labs, fixed = TRUE)
      labs <- gsub("NAA SD", "sigma['2+']", labs, fixed = TRUE) #sigma_2+
      labs <- gsub(" ", "*phantom(0)*", labs, fixed = TRUE)
      #can't do this with math expressions
      # for(i in 1:length(labs)) {
      # 	split labs[i] into multiple lines
      # 	labs[i] <- paste(strwrap(labs[i], width = 15), collapse = "\n")
      # }
      labs
    }
  }
  return(fn)
}

node.fun <- function(x, labs, digits, varlen)
{
  # paste0("Conv. Rate = ", format(round(x$frame$yval2[,5],3), nsmall = 3), "\nn = ", x$frame$n)
  n <- x$frame$n
  n.print<- character()
  n.print[which(nchar(n)<5)]<- format(n[which(nchar(n)<5)])
  n.print[which(nchar(n)>4)] <- format(n, big.mark = " ")[which(nchar(n)>4)]
  paste0(format(round(x$frame$yval2[,5]*100,1), nsmall = 1),"%\n", n.print)
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
  
  if(type == "R") temp_df <- subset(conv_df, OM_type == "naa" & OM_NAA_SD == 0)
  if(type == "R+S") temp_df <- subset(conv_df, OM_type == "naa" & OM_NAA_SD != 0)
  if(type == "R+M") temp_df <- subset(conv_df, OM_type == "M")
  if(type == "R+Sel") temp_df <- subset(conv_df, OM_type == "Sel")
  if(type == "R+q") temp_df <- subset(conv_df, OM_type == "q")
  form <- as.formula(paste("conv_txt ~", paste(factors[[type]], collapse = "+")))# (EM_PE + OM_R_SD + EM_M + EM_SR + OM_Obs._Error + OM_F_History)
  full.trees[[type]] <- rpart(form, data=temp_df, method = "class", control=rpart.control(cp=0, xval = 100))
}
anova.palette.sd <- 0.25

cairo_pdf(here("Project_0","manuscript", paste0("convergence_gradient_classification_plots.pdf")), width = 30*2/3, height = 20*2/3)
x <- matrix(c(1,1,2,2,3,3,0,4,4,5,5,0), 2, 6, byrow = TRUE)
layout.x <- layout(x) 
par(oma = c(0,0,0,0))
plot.prune(full.trees[["R"]],0.01, "R", tweak = 1.6, mar = c(0,0,5,0))
mtext("R OMs", side = 3, line = 0, cex = 2)
plot.prune(full.trees[["R+S"]],0.015, "R+S", tweak = 1.6, mar = c(0,0,5,0))
mtext("R+S OMs", side = 3, line = 0, cex = 2)
plot.prune(full.trees[["R+M"]],0.02, "R+M", tweak= 1.6, mar = c(0,0,5,0))
mtext("R+M OMs", side = 3, line = 0, cex = 2)
plot.prune(full.trees[["R+Sel"]],0.02, "R+Sel", tweak = 1.6, mar = c(0,0,5,0))
mtext("R+Sel OMs", side = 3, line = 0, cex = 2)
plot.prune(full.trees[["R+q"]],0.03, "R+q", tweak = 1.6, mar = c(0,0,5,0))
mtext("R+q OMs", side = 3, line = 0, cex = 2)
dev.off()
