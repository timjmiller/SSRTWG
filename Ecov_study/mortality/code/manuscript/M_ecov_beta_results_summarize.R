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
plot_df_fn <- function(df.ems, df.oms, Ecov_est = FALSE, M_est = FALSE, conv_type = 3, bias_res, conv_res, FE = "mean_M") {
  all_res <- list()
  for(i in 1:3) {
    re_mod <- c("rec", "rec+1", "rec+M")[i] #OM re config
    print("re_mod")
    print(re_mod)
    if(FE == "mean_M") em_ind <- which(df.ems$M_est & df.ems$Ecov_est == Ecov_est) #all three PE configs for the ems
    if(FE == "ecov_beta") em_ind <- which(df.ems$Ecov_est & df.ems$M_est == M_est) #all three PE configs for the ems
    print("em_ind")
    print(em_ind)
    om_ind <- which(df.oms$NAA_M_re == re_mod) 
    res <- lapply(om_ind, function(om) {
      print(FE)
      if(FE == "mean_M") em_ind_ <- match(em_ind,which(df.ems[["M_est"]]))
      if(FE == "ecov_beta") em_ind_ <- match(em_ind,which(df.ems[["Ecov_est"]]))
      print(em_ind_)
      out <- lapply(1:length(em_ind_), function(j){
        conv_ind <- 1:100 #all of them
        if(!is.null(conv_type)) conv_ind <- conv_fn(om,em_ind[j],conv_res,Type = conv_type) #subset consistent with convergence results
        error <- bias_res[[om]][[em_ind_[j]]][,1] #em indexes among only 6 ems fitted (either with ecov_beta or with median M estimated
        sehat <- bias_res[[om]][[em_ind_[j]]][,2]
        true_se <- sd(bias_res[[om]][[em_ind_[j]]][,1], na.rm=T)
        se_error <- sehat - true_se
        se_rel_error <- se_error/true_se
        ci_coverage_ind <- bias_res[[om]][[em_ind_[j]]][,3]
        conv <- rep(0,100)
        conv[conv_ind] <- 1
        # if(any(conv == 0 & !is.na(se_error))) {
        #   print(paste0("i: ", i," om: ", om, " em: ", em_ind[j]))
        #   print(which(conv == 0 & !is.na(se_error)))
        #   stop()
        # }
        if(!length(conv_ind)) {
          print(paste0("i: ", i," om: ", om, " em: ", em_ind[j]))
        }
        return(cbind.data.frame(om=om,em=em_ind[j],conv = conv, error=error, se = sehat, se_error = se_error, se_rel_error = se_rel_error, ci_coverage = ci_coverage_ind))
      })
      out <- do.call(rbind,out)
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
  # df <- df %>%
  #   mutate(ci_coverage = recode(ci_coverage,
  #     "1" = "Yes",
  #     "0" = "No"
  #   ))
  return(df)
}


get_bias_reg_fits <- function(factors, dfs, type = "mean_M", response = "error", is_SE = FALSE, family = "gaussian"){
  glm_fits <- dev.tables <- PRD.tables <- list()
  df <- dfs[[type]]
  if(is_SE) {
    df <- subset(df, conv == "True")
    factors <- factors[which(factors != "conv")]
  }
  for(OM_type in levels(df$OM_PE)){
    print(OM_type)
    temp <- subset(df, OM_PE == OM_type)
    print(dim(temp))
    glm_fits[[OM_type]] <- list()
    dev.tables[[OM_type]] <- list()
    for(i in factors){
      print(i)
      glm_fits[[OM_type]][[i]] <- glm(as.formula(paste(response, "~", i)), family = family, data = temp)
    }
    glm_fits[[OM_type]][["all"]] <- glm(as.formula(paste(response, "~", paste(factors,collapse = "+"))), family = family, data = temp)
    glm_fits[[OM_type]][["all2"]] <- glm(as.formula(paste(response, "~ (", paste(factors[-1],collapse = "+"), ")^2")), family = family,
      data = temp)
    glm_fits[[OM_type]][["all3"]] <- glm(as.formula(paste(response, "~ (", paste(factors[-1],collapse = "+"), ")^3")), family = family,
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
  rnames <- gsub("conv", "EM Convergence", rnames, fixed = TRUE)
  rnames <- gsub("F ", "$F$ ", rnames, fixed = TRUE)
  rnames <- gsub("EM M", "EM $M$ assumption", rnames, fixed = TRUE)
  rnames <- gsub("obs error", "OM Obs. Error", rnames, fixed = TRUE)
  rnames <- gsub("Ecov obs sig", "OM $\\sigma_e$", rnames, fixed = TRUE)
  rnames <- gsub("Ecov re sig", "OM $\\sigma_E$", rnames, fixed = TRUE)
  rnames <- gsub("Ecov re cor", "OM $\\rho_{E}$", rnames, fixed = TRUE)
  rnames <- gsub("Ecov effect", "OM $\\beta_{E}$", rnames, fixed = TRUE)
  rnames <- gsub("EM beta ecov", "EM $\\beta_{E}$ assumption", rnames, fixed = TRUE)
  rnames <- c(rnames,"All factors", "+ All Two Way", "+ All Three Way")
  rownames(PRD.table) <- rnames
  print(PRD.table)
  return(PRD.table)
}

factors <- c("1", "OM_F_history", "obs_error", "Ecov_obs_sig", "Ecov_re_sig","Ecov_re_cor", "Ecov_effect", "conv", "EM_PE","EM_beta_ecov", "EM_M")


dfs <- list()
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
  df$se_rel_error_trans <- log(df$se_rel_error + 1)
  df$se_rel_error_trans[which(is.infinite(df$se_rel_error_trans))] <- NA
  
  df[factors[-1]] <- lapply(df[factors[-1]], as.factor)
  df$OM_PE <- factor(df$OM_PE, levels= c("R","R+S", "R+M"))
  df$EM_PE <- factor(df$EM_PE, levels= c("R","R+S", "R+M"))
  dfs[[i]] <- df  
}

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

glm_fits <- list()
for(i in c("mean_M", "ecov_beta")){
  if(i == "mean_M") facs <- factors[factors != "EM_M"]
  if(i == "ecov_beta") facs <- factors[factors != "EM_beta_ecov"]
  glm_fits[[i]] <- list()
  glm_fits[[i]][["complete"]] <- get_bias_reg_fits(factors = facs, dfs = dfs, type = i)
  glm_fits[[i]][["converged"]] <- get_bias_reg_fits(factors = facs, dfs = dfs, type = i, is_SE = TRUE)
  glm_fits[[i]][["se"]] <- get_bias_reg_fits(factors = facs, dfs = dfs, type = i, response = "se_rel_error_trans", is_SE = TRUE)
  glm_fits[[i]][["ci_coverage"]] <- get_bias_reg_fits(factors = facs, dfs = dfs, type = i, response = "ci_coverage", is_SE = TRUE, family = "binomial")
}

PRD.tables <- list()
for(i in c("mean_M", "ecov_beta")){
  PRD.tables[[i]] <- list()
  for(j in c("complete", "converged", "se", "ci_coverage")){
    PRD.tables[[i]][[j]] <- get_bias_PRD_tables(glm_fits=glm_fits[[i]][[j]], factors = factors[-1])
  }
}
PRD.tables[["ecov_beta"]]

for(i in c("mean_M", "ecov_beta")){
  for(j in c("complete", "converged", "se", "ci_coverage")){
    x <- PRD.tables[[i]][[j]]
    x[] <- format(round(100*x,2), nsmall = 2)
    dim(x)
    x[which(is.na(as.numeric(x)))] <- "--"
    x[which(as.numeric(x) == 0)] <- "< 0.01"
    j_name <- j
    if(j %in% c("complete","converged")) j_name <- paste0("bias_",j)
    x <- latex(x, file = here("Ecov_study","mortality","manuscript",paste0(i,"_", j_name, "_PRD_table.tex")), 
      table.env = FALSE, col.just = rep("r", dim(x)[2]), rowlabel = "Factor", rowlabel.just = "l")#, rowname = NULL)
  }
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
  fn <- node.fun
  print(names(mod))
  if(mod$method == "class") fn <- ci.coverage.node.fun
  rpart.plot(prune.rpart(mod, cp = cp, factor = factor),
             yesno=FALSE,
             type=4, 
             clip.right.labs=TRUE, 
             xcompact = FALSE,
             ycompact = FALSE,
             extra = extra, 
             node.fun = fn, 
             split.fun = split.fun(type), 
             box.palette = "RdGn", 
             branch = 0.2, 
             fallen.leaves = FALSE, ...)
}

add_to_frame <- function(obj, data, response = "error"){
  
  origx <- newx <- obj
  origx$frame$n_test <- NA
  origx$frame$yval_test <- NA
  origx$frame$median_yval <- NA
  origx$frame$mean_abs_yval <- NA
  origx$frame$median_abs_yval <- NA
  all_node_names <- node_names <- as.numeric(rownames(newx$frame))
  toss <- max(node_names %/% 2) #node furthest out from trunk (2 leaves are 2*toss and 2*toss +1)

  #do calculations for all leaves in full model
  leaf_rows_done <- leaf_rows <- sort(which(newx$frame$var == "<leaf>"))
  leaf_names_done <- leaf_names <- all_node_names[leaf_rows_done]
  
  nloop <- 0
  data <- subset(data, !is.na(data[[response]]))
  k = 0
  while(toss > -1){
    for(i in leaf_rows){
      origx_index <- which(all_node_names == node_names[i])
      origx$frame$n_test[origx_index] <- NROW(data[which(newx$where == i),])
      origx$frame$yval_test[origx_index] <- mean(data[[response]][which(newx$where == i)])
      origx$frame$median_yval[origx_index] <- median(data[[response]][which(newx$where == i)])
      #origx$frame$mean_RE[origx_index] <- mean(data[[response]][which(newx$where == i)])
      #origx$frame$median_RE[origx_index] <- median(data[[response]][which(newx$where == i)])
      origx$frame$mean_abs_yval[origx_index] <- mean(abs(data[[response]][which(newx$where == i)]))
      origx$frame$median_abs_yval[origx_index] <- median(abs(data[[response]][which(newx$where == i)]))
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
  n <- x$frame$n
  n.print<- character()
  n.print[which(nchar(n)<5)]<- format(n[which(nchar(n)<5)])
  n.print[which(nchar(n)>4)] <- format(n, big.mark = " ")[which(nchar(n)>4)]
  paste0(out, "\n", n.print)
}

ci.coverage.node.fun <- function(x, labs, digits, varlen) {
  # paste0("Conv. Rate = ", format(round(x$frame$yval2[,5],3), nsmall = 3), "\nn = ", x$frame$n)
  
  n <- x$frame$n
  n.print<- character()
  n.print[which(nchar(n)<5)]<- format(n[which(nchar(n)<5)])
  n.print[which(nchar(n)>4)] <- format(n, big.mark = " ")[which(nchar(n)>4)]
  paste0(format(round(x$frame$yval2[,5]*100,1), nsmall = 1),"%\n", n.print)
}


full.trees <- list()
for(i in c("mean_M", "ecov_beta")){
  full.trees[[i]] <- list()
  for(OM_type in levels(dfs[[i]]$OM_PE)){
    if(i == "mean_M") facs <- factors[factors != "EM_M"]
    if(i == "ecov_beta") facs <- factors[factors != "EM_beta_ecov"]
    print(OM_type)
    full.trees[[i]][[OM_type]] <- list()
    temp <- subset(dfs[[i]], OM_PE == OM_type)
    temp$ci_cov_txt <- as.factor(c("Yes","No")[match(temp$ci_coverage,1:0)])
    response <- "error"
    form <- as.formula(paste(response,"~", paste(facs, collapse = "+")))
    full.trees[[i]][[OM_type]][["complete"]] <- rpart(form, data=temp, method = "anova", control=rpart.control(cp=0, xval = 100), model = TRUE)
    facs <- facs[which(facs != "conv")]
    form <- as.formula(paste(response,"~", paste(facs, collapse = "+")))
    temp <- subset(temp, conv == "True")
    full.trees[[i]][[OM_type]][["converged"]] <- rpart(form, data=temp, method = "anova", control=rpart.control(cp=0, xval = 100), model = TRUE)
    response <- "se_rel_error_trans"
    form <- as.formula(paste(response,"~", paste(facs, collapse = "+")))
    full.trees[[i]][[OM_type]][["se"]] <- rpart(form, data=temp, method = "anova", control=rpart.control(cp=0, xval = 100), model = TRUE)
    response <- "ci_cov_txt"
    form <- as.formula(paste(response,"~", paste(facs, collapse = "+")))
    full.trees[[i]][[OM_type]][["ci_coverage"]] <- rpart(form, data=temp, method = "class", control=rpart.control(cp=0, xval = 100), model = TRUE)
    print("OM_type done")
  }
}

# for(i in c("mean_M", "ecov_beta")){
#   if(i == "mean_M") facs <- factors[factors != "EM_M"]
#   if(i == "ecov_beta") facs <- factors[factors != "EM_beta_ecov"]
#   for(OM_type in levels(dfs[[i]]$OM_PE)){
#     temp <- subset(dfs[[i]], OM_PE == OM_type)
#     temp$ci_cov_txt <- as.factor(c("Yes","No")[match(temp$ci_coverage,1:0)])
#     facs <- facs[which(facs != "conv")]
#     response <- "ci_cov_txt"
#     form <- as.formula(paste(response,"~", paste(facs, collapse = "+")))
#     full.trees[[i]][[OM_type]][["ci_coverage"]] <- rpart(form, data=temp, method = "class", control=rpart.control(cp=0, xval = 100), model = TRUE)
#   }
# }


for(i in c("mean_M", "ecov_beta")){
  print(i)
  for(OM_type in levels(dfs[[i]]$OM_PE)){
    print(OM_type)
    temp <- subset(dfs[[i]], OM_PE == OM_type)
    full.trees[[i]][[OM_type]][["complete"]] <- add_to_frame(full.trees[[i]][[OM_type]][["complete"]], temp, response = "error")
    temp <- subset(temp, conv == "True") #median error for each node
    full.trees[[i]][[OM_type]][["converged"]] <- add_to_frame(full.trees[[i]][[OM_type]][["converged"]], temp, response = "error") #median error for each node
    full.trees[[i]][[OM_type]][["se"]] <- add_to_frame(full.trees[[i]][[OM_type]][["se"]], temp, response = "se_rel_error") #median RE for each node
  }
}

anova.palette.sd <- 0.15

#median error
cairo_pdf(here("Ecov_study","mortality","manuscript", paste0("median_M_bias_regtree_plots.pdf")), width = 30*2/3, height = 15*2/3)
x <- matrix(1:3, 1, 3, byrow = TRUE)
layout.x <- layout(x) 
par(oma = c(0,0,0,0))
#par(mfrow = c(1,2))
# plot.prune(full.trees[["mean_M"]][["R"]][["complete"]], cp = 1e-9, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
# plot.prune(prune(full.trees[["mean_M"]][["R"]][["complete"]],c(64,66,34,36,37,80,48,49,25,13,7)), cp = 1e-9, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["mean_M"]][["R"]][["converged"]], c(3,20,32,33,17,36,37,19)), cp = 0.0001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.6, split.yshift = -3)
mtext("R OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["mean_M"]][["R+S"]][["converged"]],c(4,20,42,22,6,14)), cp = 0.0001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.6, split.yshift = -3)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["mean_M"]][["R+M"]][["converged"]],c(8,36,38,10,11,3)), cp = 0.0001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.6, split.yshift = -3)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
dev.off()

#median error
cairo_pdf(here("Ecov_study","mortality","manuscript", paste0("ecov_beta_bias_regtree_plots.pdf")), width = 30*2/3, height = 15*2/3)
x <- matrix(1:3, 1, 3, byrow = TRUE)
layout.x <- layout(x) 
par(oma = c(0,0,0,0))
#par(mfrow = c(1,2))
#plot.prune(full.trees[["ecov_beta"]][["R"]][["converged"]], cp = 0.0001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["ecov_beta"]][["R"]][["converged"]], c(4,3)), cp = 0.0001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.6, split.yshift = -3)
mtext("R OMs", side = 3, line = 0, cex = 2)
#plot.prune(full.trees[["ecov_beta"]][["R+S"]][["converged"]], cp = 0.0001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["ecov_beta"]][["R+S"]][["converged"]], c(16,34,35,24)), cp = 0.0001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.6, split.yshift = -3)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
#plot.prune(full.trees[["ecov_beta"]][["R+M"]][["converged"]], cp = 0.0001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["ecov_beta"]][["R+M"]][["converged"]],c(64,65,66,34)), cp = 0.0001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.6, split.yshift = -1)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
dev.off()

#median SE error
cairo_pdf(here("Ecov_study","mortality","manuscript", paste0("median_M_SE_bias_regtree_plots.pdf")), width = 30*2/3, height = 15*2/3)
#cairo_pdf(here("Ecov_study","mortality","manuscript", paste0("temp.pdf")), width = 30*2/3, height = 15*2/3)
x <- matrix(c(1,2,3), 1, 3, byrow = TRUE)
layout.x <- layout(x) 
par(oma = c(0,1,0,0))
#plot.prune(full.trees[["mean_M"]][["R"]][["se"]], cp = 0, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
#plot.prune(prune(full.trees[["mean_M"]][["R"]][["se"]], c(8,5,12,15,28)), cp = 0, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.6)
#plot.prune(prune(full.trees[["mean_M"]][["R"]][["se"]], c(8,5,12,15,28,59,26,54,18,38,78,79,55,116,117)), cp = 0, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.6)
plot.prune(prune(full.trees[["mean_M"]][["R"]][["se"]], c(2,8,5,12,15,28,59,26,54,38,78,79,55,116,117)), cp = 0, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.4)
#plot.prune(prune(full.trees[["mean_M"]][["R"]][["se"]], c(8,18,38,78,12,15)), cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.6)
#plot.prune(prune(full.trees[["mean_M"]][["R"]][["se"]], c(8,9,12,18,38,78,5,25,15)), cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.6)
mtext("R OMs", side = 3, line = 0, cex = 2)
#plot.prune(full.trees[["mean_M"]][["R+S"]][["se"]], cp = 0, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["mean_M"]][["R+S"]][["se"]],c(4,10,22,3,94,190,191,184,370,371,186,374,92)), cp = 0, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.4, split.yshift = 0.5)
#plot.prune(prune(full.trees[["mean_M"]][["R+S"]][["se"]],c(4,10,6,7)), cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.6, split.yshift = -3)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
#plot.prune(full.trees[["mean_M"]][["R+M"]][["se"]], cp = 0, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["mean_M"]][["R+M"]][["se"]], c(32,33,34,70,71,18,38,39,5,3)), cp = 0, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.4, split.yshift = 0.5)
#plot.prune(prune(full.trees[["mean_M"]][["R+M"]][["se"]], c(4:7)), cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.6, split.yshift = -3)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
dev.off()

#median SE error
cairo_pdf(here("Ecov_study","mortality","manuscript", paste0("ecov_beta_SE_bias_regtree_plots.pdf")), width = 30*2/3, height = 15*2/3)
x <- matrix(c(1,1,2,3,3), 1, 5, byrow = TRUE)
layout.x <- layout(x) 
par(oma = c(0,0,0,0))
#plot.prune(full.trees[["ecov_beta"]][["R"]][["se"]], cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["ecov_beta"]][["R"]][["se"]], c(8,18,19,5,12,13,14,15)), cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.6, split.yshift = -3)
mtext("R OMs", side = 3, line = 0, cex = 2)
#plot.prune(full.trees[["ecov_beta"]][["R+S"]][["se"]], cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["ecov_beta"]][["R+S"]][["se"]],c(2,6,14,30)), cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.6, split.yshift = -3)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
#plot.prune(full.trees[["ecov_beta"]][["R+M"]][["se"]], cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["ecov_beta"]][["R+M"]][["se"]], c(2,4,10,22,46,47,12,26,14,30)), cp = 0.001, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.6, split.yshift = -3)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
dev.off()

anova.palette.sd <- 0.25

anova.palette.sd <- 0.1
cairo_pdf(here("Ecov_study","mortality","manuscript", paste0("median_M_CI_coverage_regtree_plots.pdf")), width = 30*2/3, height = 15*2/3)
x <- matrix(1:3, 1, 3, byrow = TRUE)
layout.x <- layout(x) 
par(oma = c(0,0,0,0))
#plot.prune(full.trees[["mean_M"]][["R"]][["ci_coverage"]], cp = 0, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["mean_M"]][["R"]][["ci_coverage"]],c(4)), cp = 0, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.6, split.yshift = -3)
mtext("R OMs", side = 3, line = 0, cex = 2)
#plot.prune(full.trees[["mean_M"]][["R+S"]][["ci_coverage"]], cp = 0, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["mean_M"]][["R+S"]][["ci_coverage"]],c(4,10,12)), cp = 0, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.6, split.yshift = -3)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
#plot.prune(full.trees[["mean_M"]][["R+M"]][["ci_coverage"]], cp = 0, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["mean_M"]][["R+M"]][["ci_coverage"]],c(4)), cp = 0, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.6, split.yshift = -3)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
dev.off()


anova.palette.sd <- 0.1
cairo_pdf(here("Ecov_study","mortality","manuscript", paste0("ecov_beta_CI_coverage_regtree_plots.pdf")), width = 30*2/3, height = 15*2/3)
x <- matrix(1:3, 1, 3, byrow = TRUE)
layout.x <- layout(x) 
par(oma = c(0,0,0,0))
#plot.prune(full.trees[["ecov_beta"]][["R"]][["ci_coverage"]], cp = 0, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["ecov_beta"]][["R"]][["ci_coverage"]],c(8,18,38)), cp = 0, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.6, split.yshift = -3)
mtext("R OMs", side = 3, line = 0, cex = 2)
#plot.prune(full.trees[["ecov_beta"]][["R+S"]][["ci_coverage"]], cp = 0, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["ecov_beta"]][["R+S"]][["ci_coverage"]],c(16,34,18,38)), cp = 0, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.6, split.yshift = -3)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
#plot.prune(full.trees[["ecov_beta"]][["R+M"]][["ci_coverage"]], cp = 0, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1)
plot.prune(prune(full.trees[["ecov_beta"]][["R+M"]][["ci_coverage"]],c(4,12,28)), cp = 0, "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.6, split.yshift = -3)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
dev.off()

# library(dplyr)
# median_summary <- dfs[["ecov_beta"]] %>% 
#   filter(OM_PE == "R" & conv == "True") %>%
#   group_by(Ecov_obs_sig,Ecov_effect) %>%
#   summarise(Median_Value = median(se_error, na.rm = TRUE)) %>% # na.rm = TRUE handles missing values
#   as.data.frame
# median_summary  
# median_summary <- dfs[["ecov_beta"]] %>% 
#   filter(conv == "True" & Ecov_est) %>%
#   group_by(OM_PE,OM_F_history, obs_error, Ecov_effect, Ecov_obs_sig, Ecov_re_sig, Ecov_re_cor, EM_PE, EM_M) %>%
#   summarise(Median_Value = median(se_error, na.rm = TRUE)) %>% # na.rm = TRUE handles missing values
#   as.data.frame
# median_summary  
# temp <- subset(median_summary, Median_Value>0)
# table(temp$OM_PE)
# subset(temp, OM_PE == "R+M")
# 
# median_summary <- dfs[["ecov_beta"]] %>% 
#   filter(conv == "True" & OM_PE == "R+M" & Ecov_est) %>%
#   group_by(Ecov_effect, obs_error, Ecov_re_sig) %>%
#   summarise(Median_Value = median(se_error, na.rm = TRUE)) %>% # na.rm = TRUE handles missing values
#   as.data.frame
# median_summary  
# 
# median_summary <- dfs[["ecov_beta"]] %>% 
#   filter(conv == "True" & OM_PE == "R" & Ecov_est) %>%
#   group_by(Ecov_effect, Ecov_obs_sig) %>%
#   summarise(median = median(se_error, na.rm = TRUE)) %>% # na.rm = TRUE handles missing values
#   as.data.frame
# median_summary  
# 
# CI_summary <- dfs[["ecov_beta"]] %>% 
#   filter(conv == "True" & OM_PE == "R" & Ecov_est) %>%
#   group_by(Ecov_effect, Ecov_obs_sig) %>%
#   summarise(cov_rate = mean(ci_coverage, na.rm = TRUE),cov_rate2 = mean((error < (qnorm(0.975)*se)) & (error>(qnorm(0.025)*se)), na.rm = TRUE)) %>% # na.rm = TRUE handles missing values
#   as.data.frame
# CI_summary  
