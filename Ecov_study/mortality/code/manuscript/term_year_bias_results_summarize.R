library(here)
library(dplyr)
library(tidyr)
library(Hmisc)
library(rpart)
#modified rpart.plot package to use plotmath in split labels...
pkgload::load_all(file.path(here(),"../../rpart.plot"))
df.ems = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.ems.RDS"))
df.oms = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.oms.RDS"))
# ssb_bias <- readRDS(file.path(here::here(),"Ecov_study","mortality", "results", "ssb_results.RDS"))
conv_res <- readRDS(here("Ecov_study","mortality", "results", "convergence_results.RDS"))
#library(reshape2)
#res <- melt(ssb_bias)

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

#all EM PE assumptions
plot_df_fn <- function(df.ems, df.oms, Ecov_est = FALSE, M_est = FALSE, conv_type = 3, years = 40,bias_res, conv_res, FE = NULL) {
  all_res <- list()
  for(i in 1:3) {
    re_mod <- c("rec", "rec+1", "rec+M")[i] #OM re config
    print("re_mod")
    print(re_mod)
    em_ind <- which(df.ems$Ecov_est== Ecov_est & df.ems$M_est == M_est) #all three PE configs for the ems
    if(!is.null(FE)) {
      if(FE == "mean_M") em_ind <- which(df.ems$M_est & df.ems$Ecov_est == Ecov_est) #all three PE configs for the ems
      if(FE == "ecov_beta") em_ind <- which(df.ems$Ecov_est & df.ems$M_est == M_est) #all three PE configs for the ems
    }
    print("em_ind")
    print(em_ind)
    om_ind <- which(df.oms$NAA_M_re == re_mod) 
    res <- lapply(om_ind, function(om) {
      if(is.null(FE)) {
        if(is.null(years)) years <- 1:40
        yr_res <- lapply(years, function(y) {
        #out <- matrix(NA,length(em_ind), 8)
        out <- lapply(em_ind, function(em){
          conv_ind <- 1:100 #all of them
          if(!is.null(conv_type)) conv_ind <- conv_fn(om,em,conv_res,Type = conv_type) #subset consistent with convergence results
          r_e <- bias_res[[om]][,em,y,2]/bias_res[[om]][,em,y,1]-1
          conv <- rep(0,100)
          conv[conv_ind] <- 1
          if(!length(conv_ind)) {
            print(paste0("i: ", i," om: ", om, " em: ", em))
            r_e <- NA
          }
          return(cbind.data.frame(year = y, om=om,em=em,conv = conv, re=r_e))
        })
        return(do.call(rbind,out))
      })
      out <- do.call(rbind,yr_res)
      } else{
        print(FE)
        if(FE == "mean_M") em_ind_ <- match(em_ind,which(df.ems[["M_est"]]))
        if(FE == "ecov_beta") em_ind_ <- match(em_ind,which(df.ems[["Ecov_est"]]))
        print(em_ind_)
        out <- lapply(1:length(em_ind_), function(j){
          conv_ind <- 1:100 #all of them
          if(!is.null(conv_type)) conv_ind <- conv_fn(om,em_ind[j],conv_res,Type = conv_type) #subset consistent with convergence results
          error <- bias_res[[om]][[em_ind_[j]]][,1] #em indexes among only 6 ems fitted (either with ecov_beta or with median M estimated
          true_sd <- sd(bias_res[[om]][[em_ind_[j]]][,1], na.rm=T)
          sd_error <- bias_res[[om]][[em_ind_[j]]][,2] - true_sd
          sd_rel_error <- sd_error/true_sd
          conv <- rep(0,100)
          conv[conv_ind] <- 1
          if(!length(conv_ind)) {
            print(paste0("i: ", i," om: ", om, " em: ", em_ind[j]))
          }
          return(cbind.data.frame(om=om,em=em_ind[j],conv = conv, error=error, sd_error = sd_error, sd_rel_error = sd_rel_error))
        })
        out <- do.call(rbind,out)
        return(out)
      }
      return(out)
    })
    res <- do.call(rbind,res)
    res <- cbind(df.ems[res$em,, drop = FALSE], res)
    print(dim(res))
    all_res[[i]] <- cbind(df.oms[res$om,, drop = FALSE], res)
  }
  df <- do.call(rbind, all_res)
  print(dim(df))
  print(head(df))

  df <- df %>%
    mutate(EM_beta_ecov = recode(as.character(Ecov_est),
      "TRUE" = "Estimated",
      "FALSE" = "0"
    ))
  df <- df %>%
    mutate(obs_error = recode(obs_error,
      "L" = "Low",
      "H" = "High"
    ))
  df <- df %>%
    mutate(OM_PE = recode(NAA_M_re,
      "rec" = "R",
      "rec+1" = "R+S",
      "rec+M" = "R+M"
    ))
  df <- df %>%
    mutate(EM_PE = recode(re_config,
      "rec" = "R",
      "rec+1" = "R+S",
      "rec+M" = "R+M"
    ))
  df <- df %>% mutate(OM_F_history = recode(Fhist,
      "H-MSY"  = "2.5*italic(F)[MSY]*phantom(0)%->%phantom(0)*italic(F)[MSY]",
      "MSY" = "italic(F)[MSY]"))
  df <- df %>%
    mutate(EM_M = recode(as.character(M_est),
      "TRUE" = "Estimated",
      "FALSE" = "Known"
    ))
  df <- df %>%
    mutate(conv = recode(conv,
      "1" = "True",
      "0" = "False"
    ))
  return(df)
}


get_bias_reg_fits <- function(factors, dfs, type = "SSB"){
  glm_fits <- dev.tables <- PRD.tables <- list()
  df <- dfs[[type]]
  for(OM_type in levels(df$OM_PE)){
    print(OM_type)
    temp <- subset(df, OM_PE == OM_type)
    print(dim(temp))
    glm_fits[[OM_type]] <- list()
    dev.tables[[OM_type]] <- list()
    if(type %in% c("SSB","F","M")) response <- "relerror_trans"
    if(type %in% c("mean_M", "ecov_beta")) response <- "error"
    for(i in factors){
      print(i)
      glm_fits[[OM_type]][[i]] <- glm(as.formula(paste(response, "~", i)), family = gaussian, data = temp)
    }
    glm_fits[[OM_type]][["all"]] <- glm(as.formula(paste(response, "~", paste(factors,collapse = "+"))), family = gaussian, data = temp)
    glm_fits[[OM_type]][["all2"]] <- glm(as.formula(paste(response, "~ (", paste(factors[-1],collapse = "+"), ")^2")), family = gaussian,
      data = temp)
    glm_fits[[OM_type]][["all3"]] <- glm(as.formula(paste(response, "~ (", paste(factors[-1],collapse = "+"), ")^3")), family = gaussian,
      data = temp)
    #percent reduction in deviance
  }
  return(glm_fits)
}

get_bias_PRD_tables <- function(glm_fits, factors){
  dev.tables <- PRD.table <- list()
  factors <- factors[factors %in% names(glm_fits[[1]])]
  for(OM_type in names(glm_fits)){
    dev.tables[[OM_type]] <- sapply(glm_fits[[OM_type]][c(factors,paste0("all",c("",2:3)))], 
      \(x) 1 - x$deviance/glm_fits[[OM_type]][[1]]$null.deviance)
  }
  print(dev.tables[[1]])
  PRD.table <- do.call(cbind,dev.tables)
  rnames <- gsub("_", " ", factors, fixed = TRUE)
  rnames <- gsub("PE", "Process Error", rnames, fixed = TRUE)
  rnames <- gsub("conv", "Convergence", rnames, fixed = TRUE)
  rnames <- gsub("F ", "$F$ ", rnames, fixed = TRUE)
  rnames <- gsub("EM M", "EM $M$ assumption", rnames, fixed = TRUE)
  rnames <- gsub("obs error", "OM Obs. Error", rnames, fixed = TRUE)
  rnames <- gsub("Ecov obs sig", "$OM \\sigma_e$", rnames, fixed = TRUE)
  rnames <- gsub("Ecov re sig", "$OM \\sigma_E$", rnames, fixed = TRUE)
  rnames <- gsub("Ecov re cor", "$OM \\rho_{E}$", rnames, fixed = TRUE)
  rnames <- gsub("Ecov effect", "OM $\\beta_{E}$", rnames, fixed = TRUE)
  rnames <- gsub("EM beta ecov", "$EM \\beta_{E}$ assumption", rnames, fixed = TRUE)
  rnames <- c(rnames,"All factors", "+ All Two Way", "+ All Three Way")
  rownames(PRD.table) <- rnames
  print(PRD.table)
  return(PRD.table)
}

factors <- c("1", "conv", "OM_F_history", "obs_error", "Ecov_obs_sig", "Ecov_re_sig","Ecov_re_cor", "Ecov_effect", "EM_PE","EM_beta_ecov", "EM_M")


dfs <- list()
for(i in c("SSB", "F", "M")){
  bias_res <- readRDS(file.path(here::here(),"Ecov_study","mortality", "results", paste0(i,"_results.RDS")))
  df <- rbind(
    plot_df_fn(df.ems, df.oms, Ecov_est = FALSE, M_est = FALSE, conv_type = 3, bias_res = bias_res, conv_res = conv_res),
    plot_df_fn(df.ems, df.oms, Ecov_est = TRUE, M_est = FALSE, conv_type = 3, bias_res = bias_res, conv_res = conv_res),
    plot_df_fn(df.ems, df.oms, Ecov_est = FALSE, M_est = TRUE, conv_type = 3, bias_res = bias_res, conv_res = conv_res),
    plot_df_fn(df.ems, df.oms, Ecov_est = TRUE, M_est = TRUE, conv_type = 3, bias_res = bias_res, conv_res = conv_res))
  df$relerror_trans <- log(df$re + 1)
  df$relerror_trans[which(is.infinite(df$relerror_trans))] <- NA
  
  df[factors[-1]] <- lapply(df[factors[-1]], as.factor)
  df$OM_PE <- factor(df$OM_PE, levels= c("R","R+S", "R+M"))
  df$EM_PE <- factor(df$EM_PE, levels= c("R","R+S", "R+M"))
  dfs[[i]] <- df  
}

for(i in c("mean_M", "ecov_beta")){
  bias_res <- readRDS(file.path(here::here(),"Ecov_study","mortality", "results", paste0(i,"_bias_results.RDS")))
  if(i == "mean_M") {
    df <- rbind(
      plot_df_fn(df.ems, df.oms, Ecov_est = FALSE, M_est = TRUE, conv_type = 3, bias_res = bias_res, conv_res = conv_res, FE = i),
      plot_df_fn(df.ems, df.oms, Ecov_est = TRUE, M_est = TRUE, conv_type = 3, bias_res = bias_res, conv_res = conv_res, FE = i))
  } else {
    df <- rbind(
      plot_df_fn(df.ems, df.oms, Ecov_est = TRUE, M_est = TRUE, conv_type = 3, bias_res = bias_res, conv_res = conv_res, FE = i),
      plot_df_fn(df.ems, df.oms, Ecov_est = TRUE, M_est = FALSE, conv_type = 3, bias_res = bias_res, conv_res = conv_res, FE = i))
  }
  df$sd_rel_error_trans <- log(df$sd_rel_error + 1)
  df$sd_rel_error_trans[which(is.infinite(df$sd_rel_error_trans))] <- NA
  
  df[factors[-1]] <- lapply(df[factors[-1]], as.factor)
  df$OM_PE <- factor(df$OM_PE, levels= c("R","R+S", "R+M"))
  df$EM_PE <- factor(df$EM_PE, levels= c("R","R+S", "R+M"))
  dfs[[i]] <- df  
}

for(i in c("SSB", "F", "M")[-1]){
  glm_fits <- get_bias_reg_fits(factors = factors, dfs = dfs, type = i)
  PRD.table <- get_bias_PRD_tables(glm_fits=glm_fits, factors = factors[-1])
  x <- PRD.table
  x[] <- format(round(100*x,2), nsmall = 2)
  dim(x)
  x[which(is.na(as.numeric(x)))] <- "--"
  x[which(as.numeric(x) == 0)] <- "< 0.01"
  x <- latex(x, file = here("Ecov_study","mortality","manuscript",paste0("bias_",i,"_PRD_table.tex")), 
    table.env = FALSE, col.just = rep("r", dim(x)[2]), rowlabel = "Factor", rowlabel.just = "l")#, rowname = NULL)
}

for(i in c("mean_M", "ecov_beta")){
  if(i == "mean_M") facs <- factors[factors != "EM_M"]
  if(i == "ecov_beta") facs <- factors[factors != "EM_beta_ecov"]
  glm_fits <- get_bias_reg_fits(factors = facs, dfs = dfs, type = i)
  PRD.table <- get_bias_PRD_tables(glm_fits=glm_fits, factors = factors[-1])
  x <- PRD.table
  x[] <- format(round(100*x,2), nsmall = 2)
  dim(x)
  x[which(is.na(as.numeric(x)))] <- "--"
  x[which(as.numeric(x) == 0)] <- "< 0.01"
  x <- latex(x, file = here("Ecov_study","mortality","manuscript",paste0("bias_",i,"_PRD_table.tex")), 
    table.env = FALSE, col.just = rep("r", dim(x)[2]), rowlabel = "Factor", rowlabel.just = "l")#, rowname = NULL)
}

#Regression trees
#########################################
#Regression/classification trees plotting functions
split.fun <- function(type = "R") {
  # replace commas with spaces (needed for strwrap)
  fn <-function(x, labs, digits, varlen, faclen){
    #labs <- gsub("_", " ", labs, fixed = TRUE)
    labs <- gsub(" = ", "==", labs, fixed = TRUE)
    labs <- gsub(" == ", "*phantom(0)==phantom(0)*", labs, fixed = TRUE)
    labs <- gsub(",", "*','*", labs, fixed = TRUE)
    labs <- gsub("EM_M", "EM*phantom(0)*italic(M)", labs, fixed = TRUE) #italic M
    labs <- gsub("EM_beta_ecov", "EM*phantom(0)*beta[italic(E)]", labs, fixed = TRUE)
    labs <- gsub("Ecov_effect", "OM*phantom(0)*beta[italic(E)]", labs, fixed = TRUE)
    labs <- gsub("EM_PE", "EM*phantom(0)*Process*phantom(0)*Error", labs, fixed = TRUE)
    labs <- gsub("OM_F_history", "OM*phantom(0)*italic(F)*phantom(0)*History", labs, fixed = TRUE)
    labs <- gsub("Ecov_re_sig", "sigma[italic(E)]", labs, fixed = TRUE)
    labs <- gsub("Ecov_re_cor", "rho[italic(E)]", labs, fixed = TRUE)
    labs <- gsub("Ecov_obs_sig", "sigma[italic(e)]", labs, fixed = TRUE)
    labs <- gsub("obs_error", "OM*phantom(0)*Pop.*phantom(0)*OE", labs, fixed = TRUE)
    labs <- gsub("conv==", "EM*phantom(0)*Converged==", labs, fixed = TRUE)
    print(labs)
    labs
  }
  return(fn)
}

prune <- function(tree, toss){
  
  newx <- tree
  ff <- newx$frame
  id <- as.integer(row.names(ff))
  toss_not_leaves <- any(id %in% toss & ff$var != "<leaf>")
  while(toss_not_leaves) {
    kids <- list(toss)
    while(length(kids[[length(kids)]])){
      possible <- c(2*kids[[length(kids)]], 2*kids[[length(kids)]] + 1)
      kids[[length(kids)+1]] <- id[which(id %in% possible & ff$var != "<leaf>")]
    }
    newx <- rpart::snip.rpart(tree, kids[[length(kids)-1]])
    ff <- newx$frame
    id <- as.integer(row.names(ff))
    toss_not_leaves <- any(id %in% toss & ff$var != "<leaf>")
  }
  return(newx)
}


#modified from part package to allow more flexibility in pruning in plots
prune.rpart <- function(tree, cp, factor = "complexity", ...)
{
  ff <- tree$frame
  id <- as.integer(row.names(ff))
  toss <- id[ff[[factor]] <= cp & ff$var != "<leaf>"] #not a leaf
  if (length(toss) == 0L) return(tree)   # all the tree is retained
  newx <- rpart::snip.rpart(tree, toss)
  newx
}

plot.prune<- function(mod,cp, type, factor = "complexity", extra = 7, ...) {
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
             branch = 0.2, 
             fallen.leaves = FALSE, ...)
}

add_to_frame <- function(obj, data){
  
  origx <- newx <- obj
  origx$frame$mean_RE <- NA
  origx$frame$mean_abs_RE <- NA
  origx$frame$median_RE <- NA
  origx$frame$median_abs_RE <- NA
  origx$frame$n_test <- NA
  origx$frame$yval_test <- NA
  origx$frame$median_yval <- NA
  all_node_names <- node_names <- as.numeric(rownames(newx$frame))
  toss <- max(node_names %/% 2) #node furthest out from trunk (2 leaves are 2*toss and 2*toss +1)

  #do calculations for all leaves in full model
  leaf_rows_done <- leaf_rows <- sort(which(newx$frame$var == "<leaf>"))
  leaf_names_done <- leaf_names <- all_node_names[leaf_rows_done]
  
  nloop <- 0
  data <- subset(data, !is.na(relerror_trans))
  k = 0
  while(toss > -1){
    for(i in leaf_rows){
      origx_index <- which(all_node_names == node_names[i])
      origx$frame$n_test[origx_index] <- NROW(data[which(newx$where == i),])
      origx$frame$yval_test[origx_index] <- mean(data$relerror_trans[which(newx$where == i)])
      origx$frame$median_yval[origx_index] <- median(data$relerror_trans[which(newx$where == i)])
      x <- data$relerror[which(newx$where == i)]
      origx$frame$mean_RE[origx_index] <- mean(data$relerror[which(newx$where == i)])
      origx$frame$median_RE[origx_index] <- median(data$relerror[which(newx$where == i)])
      origx$frame$mean_abs_RE[origx_index] <- mean(abs(data$relerror[which(newx$where == i)]))
      origx$frame$median_abs_RE[origx_index] <- median(abs(data$relerror[which(newx$where == i)]))
    }
    newx <- suppressWarnings(rpart::snip.rpart(newx, toss))
    node_names <- as.numeric(rownames(newx$frame))
    leaf_names <- toss
    leaf_rows <- which(node_names == toss)
    if(max(node_names %/% 2) == toss & toss == 0) break
    toss <- max(node_names %/% 2)
    k <- k + 1
  }
  return(origx)
}

node.fun <- function(x, labs, digits, varlen){
  out <- rep("", NROW(x$frame))
  if(!is.null(x$frame$median_yval)) out <- paste0(format(round(x$frame$median_yval,3), nsmall = 3))
  paste0(out, "\n", x$frame$n)
}

full.trees <- list()
for(i in c("SSB", "F", "M")) {
  print(i)
  full.trees[[i]] <- list()
  for(OM_type in levels(dfs[[i]]$OM_PE)){
    print(OM_type)
    temp <- subset(dfs[[i]], OM_PE == OM_type)
    form <- as.formula(paste("relerror_trans ~", paste(factors, collapse = "+")))
    full.trees[[i]][[OM_type]] <- rpart(form, data=temp, method = "anova", control=rpart.control(cp=0, xval = 100), model = TRUE)
    print("OM_type done")
  }
}

for(i in c("SSB", "F", "M")) {
  print(i)
  for(OM_type in levels(dfs[[i]]$OM_PE)){
    print(OM_type)
    temp <- subset(dfs[[i]], OM_PE == OM_type)
    full.trees[[i]][[OM_type]] <- add_to_frame(full.trees[[i]][[OM_type]], temp)
  }
}

anova.palette.sd <- 0.15

cairo_pdf(here("Ecov_study","mortality","manuscript", paste0("term_SSB_bias_regtree_plots.pdf")), width = 30*2/3, height = 15*2/3)
x <- matrix(1:3, 1, 3, byrow = TRUE)
layout.x <- layout(x) 
par(oma = c(0,0,0,0))
#plot.prune(full.trees[["SSB"]][["R"]], cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["SSB"]][["R"]],c(15)), cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
mtext("R OMs", side = 3, line = 0, cex = 2)
#plot.prune(full.trees[["SSB"]][["R+S"]], cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["SSB"]][["R+S"]],c(14,30,31)), cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
# plot.prune(full.trees[["SSB"]][["R+M"]], cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["SSB"]][["R+M"]],c(62,63)), cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
dev.off()

cairo_pdf(here("Ecov_study","mortality","manuscript", paste0("term_F_bias_regtree_plots.pdf")), width = 30*2/3, height = 15*2/3)
x <- matrix(1:3, 1, 3, byrow = TRUE)
layout.x <- layout(x) 
par(oma = c(0,0,0,0))
#plot.prune(full.trees[["F"]][["R"]], cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["F"]][["R"]],c(8)), cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
mtext("R OMs", side = 3, line = 0, cex = 2)
#plot.prune(full.trees[["F"]][["R+S"]], cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["F"]][["R+S"]],c(9,16)), cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
#plot.prune(full.trees[["F"]][["R+M"]], cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["F"]][["R+M"]],c(32,33)), cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
dev.off()

cairo_pdf(here("Ecov_study","mortality","manuscript", paste0("term_M_bias_regtree_plots.pdf")), width = 30*2/3, height = 15*2/3)
x <- matrix(1:3, 1, 3, byrow = TRUE)
layout.x <- layout(x) 
par(oma = c(0,0,0,0))
#plot.prune(full.trees[["M"]][["R"]], cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["M"]][["R"]],c(8,18)), cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
mtext("R OMs", side = 3, line = 0, cex = 2)
#plot.prune(full.trees[["M"]][["R+S"]], cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["M"]][["R+S"]],c(8,9)), cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
#plot.prune(full.trees[["M"]][["R+M"]], cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["M"]][["R+M"]],c(8,36)), cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
dev.off()
