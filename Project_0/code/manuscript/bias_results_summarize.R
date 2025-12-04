library(here)
library(dplyr)
library(tidyr)
library(Hmisc)
library(reshape2)
library(rpart)
library(viridis)
#modified rpart.plot package to use plotmath in split labels...
pkgload::load_all(file.path(here(),"../../rpart.plot"))
# pkgload::load_all("/home/tmiller/FSAM_research/aug_backup/work/rpart.plot")
# pkgload::load_all("c:/work/rpart.plot")


make_plot_df <- function(om_type = "naa", res = naa_relSR_results, is_SE = FALSE, M_or_SR = "SR", year = NULL) {
  df.ems <- readRDS(here::here("Project_0","inputs", "df.ems.RDS"))
  if(om_type == "naa") {
    em_ind <- 1:20
  }
  if(om_type == "M") {
    em_ind <- 5:24
  }
  if(om_type == "Sel") {
    em_ind <- c(5:20,25:28)
  }
  if(om_type == "q") {
    em_ind <- c(5:20,29:32)
  }
  df.ems <- df.ems[em_ind,]
  if(!is.null(M_or_SR)) {
    if(M_or_SR == "SR") em_ind <- which(df.ems$SR_model == 3)
    if(M_or_SR == "M") em_ind <- which(df.ems$M_est)
    df.ems <- df.ems[em_ind,]
  }

  res <- melt(res)
  names(res) <- c("par", "column", "value", "em", "sim","om") #em = 1: M fixed, em = 2: M estimated
  if(!is.null(year)) res <- filter(res, par == year)

  if(!(all(unique(res$em) %in% 1:NROW(df.ems)) & length(unique(res$em)) == NROW(df.ems))) stop("df.ems does not seem to be difined correctly for the number of ems in res.")
  res <- cbind(df.ems[res$em,], res)
  res <- res %>% mutate(column = recode(column,
      "1" = "relerror",
      "2" = "cv",
      "3" = "in_ci"
    ))


  df <- res %>% pivot_wider(names_from = column, values_from = value) %>% as.data.frame
  if(is_SE) df <- filter(df, !is.na(cv))

  df$relerror = df$relerror - 1

  if(om_type == "naa") {
    df.oms <- readRDS(here::here("Project_0", "inputs", "df.oms.RDS"))
  } else {
    df.oms <- readRDS(here::here("Project_0", "inputs", paste0("df.", om_type, ".oms.RDS")))
  }
  # df <- df %>% pivot_wider(names_from = type, values_from = stats) %>% as.data.frame
  print(dim(df))
  print(head(df))
  df <- cbind(df, df.oms[df$om,])
  if(!is.null(M_or_SR)) if(M_or_SR == "SR"){
    sd_log_SSB <- readRDS(file = here("Project_0","results", paste0("all_", om_type, "_sd_log_SSB.RDS")))
    sd_log_SSB <- t(matrix(unname(unlist(sd_log_SSB)), nrow = 100))
    print(dim(sd_log_SSB))
    df$sd_log_SSB <- NA
    for(i in 1:NROW(df.oms)) for(j in 1:100){
      df$sd_log_SSB[df$om == i & df$sim == j] <- sd_log_SSB[i,j]
    }
    print(dim(df))
    print(head(df))
    df$log_sd_log_SSB <- log(df$sd_log_SSB)
  }
    
#  df <- cbind(df, df.oms[df$om,])
  print(head(df))
  if(!is.null(M_or_SR)) {
    if(M_or_SR == "SR"){
      df <- df %>% mutate(par = recode(par,
        "1" = "italic(a)",
        "2" = "italic(b)"
      ))
      if(!is.null(year)) stop("Doing SR parameters, should not specify year in call to make_plot_df")
      
      #df <- df %>% mutate(M_config = if_else(M_est, "M estimated", "M known"))
    }
    if(M_or_SR == "M"){
      df <- df %>% mutate(par = recode(par, "1" = "italic(M)"))   
      df <- df %>% mutate(SR_model = recode(SR_model,
        "2" = "No SR",
        "3" = "SR"
      ))
    }
  }
  if(!is.null(year)){
    df <- df %>% mutate(SR_model = recode(SR_model,
      "2" = "No SR",
      "3" = "SR"
    ))
  }
  df <- df %>%
    mutate(EM_process_error = recode(re_config,
      "rec" = "R",
      "rec+1" = "R+S",
      "M_re" = "R+M",
      "q_re" = "R+q",
      "sel_re" = "R+Sel"
    ))
  ind <- which(df$M_re_cor == "iid" | df$sel_re_cor == "iid" | df$q_re_cor == "iid")
  df$EM_process_error[ind] <- paste0(df$EM_process_error[ind], "(iid)")
  ind <- which(df$M_re_cor == "ar1_y" | df$sel_re_cor == "ar1_y" | df$q_re_cor == "ar1")
  df$EM_process_error[ind] <- paste0(df$EM_process_error[ind], "(AR1)")
  EM_process_error <- c("R","R+S","R+M(iid)","R+M(AR1)","R+Sel(iid)","R+Sel(AR1)","R+q(iid)","R+q(AR1)")
  df$EM_process_error <- factor(df$EM_process_error, levels = EM_process_error)

  df$correct_EM_PE <- "No"
  if(om_type == "naa") {
    df$NAA_sig[which(is.na(df$NAA_sig))] <- 0
    df$correct_EM_PE[df$NAA_sig == 0 & df$EM_process_error == "R"] <- "Yes"
    df$correct_EM_PE[df$NAA_sig >  0 & df$EM_process_error == "R+S"] <- "Yes"
    df$OM_NAA_sigma <- df$NAA_sig
    df$OM_R_sigma <- df$R_sig
  }

  if(om_type == "M") {
    df$correct_EM_PE[df$M_cor == 0 & df$EM_process_error == "R+M(iid)"] <- "Yes"
    df$correct_EM_PE[df$M_cor >  0 & df$EM_process_error == "R+M(AR1)"] <- "Yes"
    df$OM_M_sigma <- df$M_sig
    df$OM_M_rho <- df$M_cor
  }
  if(om_type == "Sel") {
    df$correct_EM_PE[df$Sel_cor == 0 & df$EM_process_error == "R+Sel(iid)"] <- "Yes"
    df$correct_EM_PE[df$Sel_cor >  0 & df$EM_process_error == "R+Sel(AR1)"] <- "Yes"
    df$OM_Sel_sigma <- df$Sel_sig
    df$OM_Sel_rho <- df$Sel_cor
  }
  if(om_type == "q") {
    df$correct_EM_PE[df$q_cor == 0 & df$EM_process_error == "R+q(iid)"] <- "Yes"
    df$correct_EM_PE[df$q_cor >  0 & df$EM_process_error == "R+q(AR1)"] <- "Yes"
    df$OM_q_sigma <- df$q_sig
    df$OM_q_rho <- df$q_cor
  }
  df <- df %>% mutate(EM_M = recode(as.character(M_est),
      "TRUE" = "Estimated",
      "FALSE" = "Known"))
  df <- df %>%
    mutate(OM_Obs._Error = recode(obs_error,
      "L" = "Low",
      "H" = "High"))
  df <- df %>% mutate(OM_F_History = recode(Fhist,
        "H-MSY"  = "2.5*italic(F)[MSY] %->% italic(F)[MSY]",
        "MSY" = "italic(F)[MSY]"))
  df <- df %>% as.data.frame
  facs <- c("OM", "OM_F_History","OM_NAA_sigma", "OM_R_sigma", "OM_M_sigma", "OM_M_rho", "OM_Sel_sigma", "OM_Sel_rho", "OM_q_sigma", "OM_q_rho", "OM_Obs._Error","EM_M", "correct_EM_PE")
  fac_names <- names(df)[names(df) %in% facs]  
  df[fac_names] <- lapply(df[fac_names], as.factor)
  # df$correct_EM_PE[df$correct_EM_PE== "No"] <- NA
  # fac_names <- names(df)[names(df) %in% facs]  
  return(df)
}
# obs_dfs$SR$naa <- make_plot_df(om_type = "naa", res = naa_relSR_results)
# obs_dfs$SSB$naa <- make_plot_df(om_type = "naa", res = all_naa_relssb, M_or_SR = NULL, year = 40)
# temp <- make_plot_df(om_type = "M", res = M_rel_M_results, M_or_SR = "M")

#########################################
#Regression/classification trees plotting functions
split.fun <- function(type = "R") {
  # replace commas with spaces (needed for strwrap)
  if(!type %in% c("R","R+S")){
    fn <-function(x, labs, digits, varlen, faclen){
      labs <- gsub("_", " ", labs, fixed = TRUE)
      labs <- gsub(" = ", "==", labs, fixed = TRUE)
      labs <- gsub(",", "*','*", labs, fixed = TRUE)
      labs <- gsub("+", "*'+'*", labs, fixed = TRUE)
      others <- c("M","Sel","q")[c("R+M","R+Sel","R+q") != type]
      for(i in others) labs <- gsub(paste0(i,"(iid)"), i, labs, fixed = TRUE)
      labs <- gsub("EM SR", "EM*phantom(0)*SR*phantom(0)*Assumption", labs, fixed = TRUE)
      labs <- gsub("Sel sigma", "sigma['Sel']", labs, fixed = TRUE) #sigma_Sel
      labs <- gsub("q sigma", "sigma['q']", labs, fixed = TRUE) #sigma_q
      labs <- gsub(" M sigma", " sigma['M']", labs, fixed = TRUE) #sigma_M
      labs <- gsub("Sel rho", "rho['Sel']", labs, fixed = TRUE) #sigma_Sel
      labs <- gsub("q rho", "rho['q']", labs, fixed = TRUE) #sigma_q
      labs <- gsub(" M rho", " rho['M']", labs, fixed = TRUE) #sigma_M
      labs <- gsub("R sigma", "sigma[R]", labs, fixed = TRUE) #sigma_R
      labs <- gsub("NAA sigma", "sigma['2+']", labs, fixed = TRUE) #sigma_2+
      labs <- gsub("log sd log SSB", "log(SD[SSB])", labs, fixed = TRUE)
      labs <- gsub(" M", "*phantom(0)*italic(M)", labs, fixed = TRUE) #italic M
      labs <- gsub("process error", "Process*phantom(0)*Error", labs, fixed = TRUE)
      labs <- gsub("F ", "italic(F)*phantom(0)*", labs, fixed = TRUE)
      labs <- gsub(" ", "*phantom(0)*", labs, fixed = TRUE)
      labs <- gsub("*%", "%", labs, fixed = TRUE)
      labs <- gsub("%*", "%", labs, fixed = TRUE)
      labs <- gsub("*<", "<", labs, fixed = TRUE)
      labs <- gsub("*>", ">", labs, fixed = TRUE)
      labs <- gsub("=*", "=", labs, fixed = TRUE)
      labs <- gsub(">*", ">", labs, fixed = TRUE)
      labs <- gsub("<*", "<", labs, fixed = TRUE)
      labs <- gsub("^\\*", "", labs)
      labs
    }
  } else {
    fn <-function(x, labs, digits, varlen, faclen){
      labs <- gsub("_", " ", labs, fixed = TRUE)
      labs <- gsub("(iid)", "", labs, fixed = TRUE)
      labs <- gsub(" = ", "==", labs, fixed = TRUE)
      labs <- gsub(",", "*','*", labs, fixed = TRUE)
      labs <- gsub("+", "*'+'*", labs, fixed = TRUE)
      labs <- gsub(" M", "*phantom(0)*italic(M)", labs, fixed = TRUE) #italic M
      labs <- gsub("log sd log SSB", "log(SD[SSB])", labs, fixed = TRUE)
      labs <- gsub("process error", "Process*phantom(0)*Error", labs, fixed = TRUE)
      labs <- gsub("F ", "italic(F)*phantom(0)*", labs, fixed = TRUE)
      labs <- gsub("EM SR", "EM*phantom(0)*SR*phantom(0)*Assumption", labs, fixed = TRUE)
      labs <- gsub("R sigma", "sigma[R]", labs, fixed = TRUE) #sigma_2+
      labs <- gsub("NAA sigma", "sigma['2+']", labs, fixed = TRUE) #sigma_2+
      labs <- gsub(" ", "*phantom(0)*", labs, fixed = TRUE)
      labs <- gsub("*%", "%", labs, fixed = TRUE)
      labs <- gsub("%*", "%", labs, fixed = TRUE)
      labs <- gsub("*<", "<", labs, fixed = TRUE)
      labs <- gsub("*>", ">", labs, fixed = TRUE)
      labs <- gsub("=*", "=", labs, fixed = TRUE)
      labs <- gsub(">*", ">", labs, fixed = TRUE)
      labs <- gsub("<*", "<", labs, fixed = TRUE)
      labs <- gsub("^\\*", "", labs)
      labs
    }
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
    ## Now cut down the CP table
    # temp <- pmax(tree$cptable[, 1L], cp)
    # keep <- match(unique(temp), temp)
    # newx$cptable <- tree$cptable[keep, , drop = FALSE]
    # newx$cptable[length(keep), 1L] <- cp
    # # Reset the variable importance
    # newx$variable.importance <- rpart:::importance(newx)
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
  # print(which((node_names %/% 2) == toss))
  # print(toss)
  #do calculations for all leaves in full model
  leaf_rows_done <- leaf_rows <- sort(which(newx$frame$var == "<leaf>"))
  leaf_names_done <- leaf_names <- all_node_names[leaf_rows_done]

  nloop <- 0
  data <- subset(data, !is.na(relerror_trans))
  k = 0
  while(toss > -1){
    # print(paste0("n leaves: ", length(leaf_rows)))
    for(i in leaf_rows){
      #origx_index <- which(all_node_names == node_names[i])
      # if(k == 1) print(i)
      # if(k == 1) print(node_names[i])
      origx_index <- which(all_node_names == node_names[i])
      # if(k == 1) print(origx_index)
      origx$frame$n_test[origx_index] <- NROW(data[which(newx$where == i),])
      origx$frame$yval_test[origx_index] <- mean(data$relerror_trans[which(newx$where == i)])
      origx$frame$median_yval[origx_index] <- median(data$relerror_trans[which(newx$where == i)])
      x <- data$relerror[which(newx$where == i)]
      origx$frame$mean_RE[origx_index] <- mean(data$relerror[which(newx$where == i)])
      origx$frame$median_RE[origx_index] <- median(data$relerror[which(newx$where == i)])
      origx$frame$mean_abs_RE[origx_index] <- mean(abs(data$relerror[which(newx$where == i)]))
      origx$frame$median_abs_RE[origx_index] <- median(abs(data$relerror[which(newx$where == i)]))
      # if(k == 1) print(origx$frame[origx_index,])
      #print(paste0("end i: ", i))
    }
    # print(paste0("toss: ", toss))
    newx <- suppressWarnings(rpart::snip.rpart(newx, toss))
    # print("newx done")
    node_names <- as.numeric(rownames(newx$frame))
    #print(tail(sort(node_names)))
    # print(newx$frame[which(node_names == toss),])
    leaf_names <- toss
    # print(leaf_names)
    leaf_rows <- which(node_names == toss)
    # print(leaf_rows)
    # print(paste0("n node_names: ", length(node_names)))
    # leaf_rows <- which(newx$frame$var == "<leaf>")
    # leaf_names <- node_names[leaf_rows]
    # leaf_rows <- leaf_rows[which(!leaf_names %in% leaf_names_done)]
    # leaf_names <- node_names[leaf_rows]
    if(max(node_names %/% 2) == toss & toss == 0) break
    toss <- max(node_names %/% 2)
    # print(paste0("toss: ", toss))
    k <- k + 1
    #if(k == 2) stop()
    # stop()
  }
  return(origx)
}
# full.trees[[SR_par]][[OM_type]] <- add_to_frame(full.trees[[SR_par]][[OM_type]], temp)
# OM_type = "R+Sel"
# SR_par = "a"
# temp <- get_small_data(SR_par, OM_type, obs_dfs, factors, cv_limit)
# head(full.trees[[SR_par]][[OM_type]]$frame)

node.fun <- function(x, labs, digits, varlen){
  # out <- paste0("Mean log(|RE|) = ", format(round(x$frame$yval,3), nsmall = 3))
  # out <- paste0("Mean |RE| = ", format(round(x$frame$mean_abs_RE,3), nsmall = 3))
  out <- rep("", NROW(x$frame))
  # if(!is.null(x$frame$median_abs_RE)) out <- paste0("Median |RE| = ", format(round(x$frame$median_abs_RE,3), nsmall = 3))
  # if(!is.null(x$frame$median_RE)) out <- paste0("Median RE = ",format(round(x$frame$median_RE,3), nsmall = 3))
  # out <- paste0(out, "\nMean response = ", format(round(x$frame$yval,3), nsmall = 3))
  if(!is.null(x$frame$median_yval)) out <- paste0(format(round(x$frame$median_yval,3), nsmall = 3))
  paste0(out, "\n", x$frame$n)
}

get_small_data <- function(SR_par, OM_type, obs_dfs, factors, cv_limit){
  print("get_small_data")
  print(paste0("SR_par:", SR_par))
  par_type <- ifelse(SR_par %in% c("a","b"), "SR", SR_par)
  print(paste0("par_type: ", par_type))
  print(paste0("OM_type: ", OM_type))
  facs <- factors[[par_type]][[OM_type]]
  dfs <- obs_dfs[[par_type]]
  
  if(OM_type == "R") temp <- subset(dfs[["naa"]], OM_NAA_sigma == 0 )
  if(OM_type == "R+S") temp <- subset(dfs[["naa"]], OM_NAA_sigma != 0)
  if(OM_type == "R+M") temp <- dfs[["M"]]
  if(OM_type == "R+Sel") temp <- dfs[["Sel"]]
  if(OM_type == "R+q") temp <- dfs[["q"]]
  if(SR_par %in% c("a","b","M")) temp <- subset(temp, par == paste0("italic(",SR_par,")"))
  if(!is.na(cv_limit)) temp <- subset(temp, cv < cv_limit) #delta-method based cv, not log-normal
  temp$relerror_trans <- log(temp$relerror + 1)
  temp$relerror_trans[which(is.infinite(temp$relerror_trans))] <- NA
  return(temp)
}

########################################
#Stock-recruit pars

naa_relSR_results <- readRDS(file = here("Project_0","results", "naa_relSR_results_all_EM_PE.RDS"))
M_relSR_results <- readRDS(file = here("Project_0","results", "M_relSR_results_all_EM_PE.RDS"))
Sel_relSR_results <- readRDS(file = here("Project_0","results", "Sel_relSR_results_all_EM_PE.RDS"))
q_relSR_results <- readRDS(file = here("Project_0","results", "q_relSR_results_all_EM_PE.RDS"))

# EM_process_error <- c("R","R+S","R+M (iid)","R+M (AR1)","R+Sel (iid)","R+Sel (AR1)","R+q (iid)","R+q (AR1)")
# R_S_df$EM_process_error <- factor(R_S_df$EM_process_error, levels = EM_process_error)

obs_dfs <- list(SR = list(), M = list())
obs_dfs$SR$naa <- make_plot_df(om_type = "naa", res = naa_relSR_results)
obs_dfs$SR$M <- make_plot_df(om_type = "M", res = M_relSR_results)
obs_dfs$SR$Sel <- make_plot_df(om_type = "Sel", res = Sel_relSR_results)
obs_dfs$SR$q <- make_plot_df(om_type = "q", res = q_relSR_results)

# temp <- make_plot_df(om_type = "M", res = M_relSR_results)
# dim(temp)
# dim(obs_dfs$SR$M)
# temp <- make_plot_df(om_type = "Sel", res = Sel_relSR_results)
# dim(temp)
# dim(obs_dfs$SR$Sel)
# temp <- make_plot_df(om_type = "q", res = q_relSR_results)
# dim(temp)
# dim(obs_dfs$SR$q)

########################################
#Natural Mortality

naa_rel_M_results <- readRDS(file = here("Project_0","results", "naa_rel_M_results_all_EM_PE.RDS"))
M_rel_M_results <- readRDS(file = here("Project_0","results", "M_rel_M_results_all_EM_PE.RDS"))
Sel_rel_M_results <- readRDS(file = here("Project_0","results", "Sel_rel_M_results_all_EM_PE.RDS"))
q_rel_M_results <- readRDS(file = here("Project_0","results", "q_rel_M_results_all_EM_PE.RDS"))

obs_dfs$M$naa <- make_plot_df(om_type = "naa", res = naa_rel_M_results, M_or_SR = "M")
obs_dfs$M$M <- make_plot_df(om_type = "M", res = M_rel_M_results, M_or_SR = "M")
obs_dfs$M$Sel <- make_plot_df(om_type = "Sel", res = Sel_rel_M_results, M_or_SR = "M")
obs_dfs$M$q <- make_plot_df(om_type = "q", res = q_rel_M_results, M_or_SR = "M")

# temp <- make_plot_df(om_type = "naa", res = naa_rel_M_results, M_or_SR = "M")
# dim(temp)
# dim(obs_dfs$M$naa)
# temp <- make_plot_df(om_type = "M", res = M_relS_M_results, M_or_SR = "M")
# dim(temp)
# dim(obs_dfs$M$M)
# temp <- make_plot_df(om_type = "Sel", res = Sel_rel_M_results, M_or_SR = "M")
# dim(temp)
# dim(obs_dfs$M$Sel)
# temp <- make_plot_df(om_type = "q", res = q_rel_M_results, M_or_SR = "M")
# dim(temp)
# dim(obs_dfs$M$q)

########################################
#make dfs for SSB,F

all_naa_relssb <- readRDS(file = here("Project_0","results", "all_naa_relssb_results.RDS"))
all_M_relssb <- readRDS(file = here("Project_0","results", "all_M_relssb_results.RDS"))
all_Sel_relssb <- readRDS(file = here("Project_0","results", "all_Sel_relssb_results.RDS"))
all_q_relssb <- readRDS(file = here("Project_0","results", "all_q_relssb_results.RDS"))

obs_dfs$SSB <- list()
obs_dfs$SSB$naa <- make_plot_df(om_type = "naa", res = all_naa_relssb, M_or_SR = NULL, year = 40)
obs_dfs$SSB$M <- make_plot_df(om_type = "M", res = all_M_relssb, M_or_SR = NULL, year = 40)
obs_dfs$SSB$Sel <- make_plot_df(om_type = "Sel", res = all_Sel_relssb, M_or_SR = NULL, year = 40)
obs_dfs$SSB$q <- make_plot_df(om_type = "q", res = all_q_relssb, M_or_SR = NULL, year = 40)

all_naa_relF <- readRDS(file = here("Project_0","results", "all_naa_relF_results.RDS"))
all_M_relF <- readRDS(file = here("Project_0","results", "all_M_relF_results.RDS"))
all_Sel_relF <- readRDS(file = here("Project_0","results", "all_Sel_relF_results.RDS"))
all_q_relF <- readRDS(file = here("Project_0","results", "all_q_relF_results.RDS"))

obs_dfs$F <- list()
obs_dfs$F$naa <- make_plot_df(om_type = "naa", res = all_naa_relF, M_or_SR = NULL, year = 40)
obs_dfs$F$M <- make_plot_df(om_type = "M", res = all_M_relF, M_or_SR = NULL, year = 40)
obs_dfs$F$Sel <- make_plot_df(om_type = "Sel", res = all_Sel_relF, M_or_SR = NULL, year = 40)
obs_dfs$F$q <- make_plot_df(om_type = "q", res = all_q_relF, M_or_SR = NULL, year = 40)

all_naa_relR <- readRDS(file = here("Project_0","results", "all_naa_relR_results.RDS"))
all_M_relR <- readRDS(file = here("Project_0","results", "all_M_relR_results.RDS"))
all_Sel_relR <- readRDS(file = here("Project_0","results", "all_Sel_relR_results.RDS"))
all_q_relR <- readRDS(file = here("Project_0","results", "all_q_relR_results.RDS"))

obs_dfs$R <- list()
obs_dfs$R$naa <- make_plot_df(om_type = "naa", res = all_naa_relR, M_or_SR = NULL, year = 40)
obs_dfs$R$M <- make_plot_df(om_type = "M", res = all_M_relR, M_or_SR = NULL, year = 40)
obs_dfs$R$Sel <- make_plot_df(om_type = "Sel", res = all_Sel_relR, M_or_SR = NULL, year = 40)
obs_dfs$R$q <- make_plot_df(om_type = "q", res = all_q_relR, M_or_SR = NULL, year = 40)

########################################
#define factors for fits
factors <- list(SR = list(), M = list())
# null <- gam(conv ~ 1, family = binomial, method = "REML", data = temp)
factors[["SR"]][["R"]] <- c("1", "EM_process_error","OM_R_sigma","EM_M", "OM_Obs._Error", "OM_F_History")
factors[["SR"]][["R+S"]] <- c("1", "EM_process_error","OM_R_sigma","EM_M", "OM_Obs._Error", "OM_F_History","OM_NAA_sigma")
factors[["SR"]][["R+M"]] <- c("1", "EM_process_error","EM_M","OM_Obs._Error", "OM_F_History","OM_M_sigma", "OM_M_rho")
factors[["SR"]][["R+Sel"]] <- c("1", "EM_process_error","EM_M","OM_Obs._Error", "OM_F_History","OM_Sel_sigma", "OM_Sel_rho")
factors[["SR"]][["R+q"]] <- c("1", "EM_process_error","EM_M","OM_Obs._Error", "OM_F_History","OM_q_sigma", "OM_q_rho")

factors[["M"]][["R"]] <- c("1", "EM_process_error","OM_R_sigma","SR_model", "OM_Obs._Error", "OM_F_History")
factors[["M"]][["R+S"]] <- c("1", "EM_process_error","OM_R_sigma","SR_model", "OM_Obs._Error", "OM_F_History","OM_NAA_sigma")
factors[["M"]][["R+M"]] <- c("1", "EM_process_error","SR_model","OM_Obs._Error", "OM_F_History","OM_M_sigma", "OM_M_rho")
factors[["M"]][["R+Sel"]] <- c("1", "EM_process_error","SR_model","OM_Obs._Error", "OM_F_History","OM_Sel_sigma", "OM_Sel_rho")
factors[["M"]][["R+q"]] <- c("1", "EM_process_error","SR_model","OM_Obs._Error", "OM_F_History","OM_q_sigma", "OM_q_rho")

factors[["SSB"]] <- list()
factors[["SSB"]][["R"]] <- c("1", "EM_process_error","EM_M","OM_R_sigma","SR_model", "OM_Obs._Error", "OM_F_History")
factors[["SSB"]][["R+S"]] <- c("1", "EM_process_error","EM_M","OM_R_sigma","SR_model", "OM_Obs._Error", "OM_F_History","OM_NAA_sigma")
factors[["SSB"]][["R+M"]] <- c("1", "EM_process_error","EM_M","SR_model","OM_Obs._Error", "OM_F_History","OM_M_sigma", "OM_M_rho")
factors[["SSB"]][["R+Sel"]] <- c("1", "EM_process_error","EM_M","SR_model","OM_Obs._Error", "OM_F_History","OM_Sel_sigma", "OM_Sel_rho")
factors[["SSB"]][["R+q"]] <- c("1", "EM_process_error","EM_M","SR_model","OM_Obs._Error", "OM_F_History","OM_q_sigma", "OM_q_rho")

factors[["F"]] <- factors[["R"]] <- factors[["SSB"]]
########################################

# glm_fits <- dev.tables <- PRD.tables <- list()
get_bias_reg_fits <- function(pars = c("a","b","M", "SSB", "F", "R"), factors, obs_dfs, cv_limit=NA){
#cv_limit <- NA
#for(SR_par in c("a","b","M", "SSB", "F", "R")) {
  for(SR_par in pars) {
    print(SR_par)
    par_type <- ifelse(SR_par %in% c("a","b"), "SR", SR_par)
    print(par_type)
    glm_fits[[SR_par]] <- dev.tables[[SR_par]] <- PRD.tables[[SR_par]] <- list()
    for(OM_type in names(factors[[par_type]])){
      print(OM_type)
      facs <- factors[[par_type]][[OM_type]]
      temp <- get_small_data(SR_par, OM_type, obs_dfs, factors, cv_limit)
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
    }
  }
  return(glm_fits)
}

get_bias_PRD_tables <- function(pars = c("a","b","M", "SSB", "F", "R"), glm_fits, factors){
  dev.tables <- PRD.tables <- list()
  for(SR_par in pars) {
    print(SR_par)
    par_type <- ifelse(SR_par %in% c("a","b"), "SR", SR_par)
    print(par_type)
    dev.tables[[SR_par]] <- PRD.tables[[SR_par]] <- list()
    for(OM_type in names(glm_fits[[SR_par]])){
      facs <- factors[[par_type]][[OM_type]]
      dev.tables[[SR_par]][[OM_type]] <- sapply(glm_fits[[SR_par]][[OM_type]][facs], \(x) 1 - x$deviance/glm_fits[[SR_par]][[OM_type]][[1]]$null.deviance)
    }
    interactions.dev.table <- sapply(names(factors[[par_type]]), \(x) 1 - c(glm_fits[[SR_par]][[x]][["all"]]$deviance, glm_fits[[SR_par]][[x]][["all2"]]$deviance,
      glm_fits[[SR_par]][[x]][["all3"]]$deviance)/glm_fits[[SR_par]][[x]][[1]]$null.deviance)
    x <- as.data.frame(round(100*interactions.dev.table,2))
    x <- cbind(Model = c("No Interactions", "+ All Two Way", "+ All Three Way"), x)
  
    All.facs <- c("EM_process_error","OM_Obs._Error", "OM_F_History","OM_R_sigma","OM_NAA_sigma","OM_M_sigma", "OM_M_rho","OM_Sel_sigma", "OM_Sel_rho", 
      "OM_q_sigma", "OM_q_rho")
    if(SR_par %in% c("a","b")) All.facs <- c("EM_M", All.facs)
    if(SR_par %in% c("M")) All.facs <- c("SR_model", All.facs)
    if(SR_par %in% c("SSB", "F", "R")) All.facs <- c("EM_M", "SR_model",All.facs)
    PRD.table <- matrix(NA, length(All.facs),5)
    colnames(PRD.table) <- c("R","R+S","R+M","R+Sel","R+q")
    rnames <- gsub("_", " ", All.facs, fixed = TRUE)
    rnames <- gsub("process error", "Process Error", rnames, fixed = TRUE)
    rnames <- gsub("F ", "$F$ ", rnames, fixed = TRUE)
    rnames <- gsub("EM M", "EM $M$ Assumption", rnames, fixed = TRUE)
    rnames <- gsub("log sd log SSB", "$\\log\\left(SD_\\text{SSB}]\\right)$", rnames, fixed = TRUE)
    rnames <- gsub("NAA sigma", "$\\sigma_{2+}$ ", rnames, fixed = TRUE)
    rnames <- gsub("R sigma", "$\\sigma_R$", rnames, fixed = TRUE)
    rnames <- gsub("q sigma", "$\\sigma_q$", rnames, fixed = TRUE)
    rnames <- gsub("M sigma", "$\\sigma_M$", rnames, fixed = TRUE)
    rnames <- gsub("Sel sigma", "$\\sigma_{\\text{Sel}}$", rnames, fixed = TRUE)
    rnames <- gsub("M rho", "$\\rho_M$", rnames, fixed = TRUE)
    rnames <- gsub("Sel rho", "$\\rho_{\\text{Sel}}$", rnames, fixed = TRUE)
    rnames <- gsub("q rho", "$\\rho_q$", rnames, fixed = TRUE)
    rnames <- gsub("SR model", "EM SR assumption", rnames, fixed = TRUE)
  
    rownames(PRD.table) <- rnames
    print(dim(PRD.table))
    for(i in c("R","R+S","R+M","R+Sel","R+q")) {
      print(names(dev.tables[[SR_par]][[i]]))
      print(All.facs)
      print(length(dev.tables[[SR_par]][[i]][match(All.facs, names(dev.tables[[SR_par]][[i]]))]))
      PRD.table[,i] <- dev.tables[[SR_par]][[i]][match(All.facs, names(dev.tables[[SR_par]][[i]]))]
    }
    y <- interactions.dev.table
    rownames(y) <- c("All factors", "+ All Two Way", "+ All Three Way")
    PRD.tables[[SR_par]] <- rbind(PRD.table,y)
  }
  return(PRD.tables)
}

glm_fits <- get_bias_reg_fits(factors = factors, obs_dfs = obs_dfs)
PRD.tables <- get_bias_PRD_tables(glm_fits=glm_fits, factors = factors)

saveRDS(glm_fits, here::here("Project_0","results", "glm_fits_bias.RDS"))


# # examine which combination of factors are important for bias of SRR pars for R and R+S OMs : what is PRD if we leave one out?
# x <- factors[["SR"]][["R"]]
# tfits <- lapply(x[-1], \(y) glm(as.formula(paste("relerror_trans", "~", paste(x[x!=y],collapse = "+"))), family = gaussian, data = get_small_data("a", "R", obs_dfs, factors, NA)))
# names(tfits) <- paste0("not_",x[2:6])
# sapply(tfits, \(x) 1 - x$deviance/glm_fits[["a"]][["R"]][[1]]$null.deviance)
# 
# x <- factors[["SR"]][["R+S"]]
# tfits <- lapply(x[2:8], \(y) glm(as.formula(paste("relerror_trans", "~", paste(x[x!=y],collapse = "+"))), family = gaussian, data = get_small_data("a", "R+S", obs_dfs, factors, NA)))
# names(tfits) <- paste0("not_",x[2:7])
# tfits[["R_sigma+SD_SSB"]] <- glm(as.formula(paste("relerror_trans", "~", paste(x[c(3,7:8)],collapse = "+"))), family = gaussian, data = get_small_data("a", "R+S", obs_dfs, factors, NA))
# sapply(tfits, \(x) 1 - x$deviance/glm_fits[["a"]][["R+S"]][[1]]$null.deviance)
####################################################

x <- cbind(PRD.tables[["a"]],PRD.tables[["b"]])
x[] <- format(round(100*x,2), nsmall = 2)
dim(x)
x[which(is.na(as.numeric(x)))] <- "--"
x[which(as.numeric(x) == 0)] <- "< 0.01"
x <- latex(x, file = here("Project_0","manuscript","bias_SR_pars_PRD_table.tex"), 
  cgroup = c("Beverton-Holt $a$","Beverton-Holt $b$"), n.cgroup = c(5,5),
  table.env = FALSE, col.just = rep("r", dim(x)[2]), rowlabel = "Factor", rowlabel.just = "l")#, rowname = NULL)

x <- PRD.tables[["M"]]
x[] <- format(round(100*x,2), nsmall = 2)
dim(x)
x[which(is.na(as.numeric(x)))] <- "--"
x[which(as.numeric(x) == 0)] <- "< 0.01"
x <- latex(x, file = here("Project_0","manuscript","bias_median_M_PRD_table.tex"), 
  table.env = FALSE, col.just = rep("r", dim(x)[2]), rowlabel = "Factor", rowlabel.just = "l")#, rowname = NULL)

for(i in c("SSB", "F", "R")){
  x <- PRD.tables[[i]]
  x[] <- format(round(100*x,2), nsmall = 2)
  dim(x)
  x[which(is.na(as.numeric(x)))] <- "--"
  x[which(as.numeric(x) == 0)] <- "< 0.01"
  x <- latex(x, file = here("Project_0","manuscript",paste0("bias_",i,"_PRD_table.tex")), 
             table.env = FALSE, col.just = rep("r", dim(x)[2]), rowlabel = "Factor", rowlabel.just = "l")#, rowname = NULL)
}
#Regression trees

full.trees <- list()#a = list(), b = list(), M = list())
for(SR_par in c("a","b","M", "SSB", "F", "R")) {
  print(SR_par)
  par_type <- ifelse(SR_par %in% c("a","b"), "SR", SR_par)
  print(par_type)
  full.trees[[SR_par]] <- list()
  for(OM_type in names(factors[[par_type]])){
    print(OM_type)
    temp <- get_small_data(SR_par, OM_type, obs_dfs, factors, cv_limit)
    facs <- factors[[par_type]][[OM_type]]
    form <- as.formula(paste("relerror_trans ~", paste(facs, collapse = "+")))# (EM_PE + OM_R_SD + EM_M + EM_SR + OM_Obs._Error + OM_F_History)
    full.trees[[SR_par]][[OM_type]] <- rpart(form, data=temp, method = "anova", control=rpart.control(cp=0, xval = 100), model = TRUE)#, roundint = FALSE)
    print("OM_type done")
  }
}

for(SR_par in c("a","b","M", "SSB", "F", "R")) {
  print(SR_par)
  par_type <- ifelse(SR_par %in% c("a","b"), "SR", SR_par)
  print(par_type)
  for(OM_type in names(factors[[par_type]])){
    print(OM_type)
    temp <- get_small_data(SR_par, OM_type, obs_dfs, factors, cv_limit)
    full.trees[[SR_par]][[OM_type]] <- add_to_frame(full.trees[[SR_par]][[OM_type]], temp)
  }
}

saveRDS(full.trees, here::here("Project_0","results", "reg_trees_bias.RDS"))

#this is needed inside my modified handle.anova.palette function
anova.palette.sd <- 0.15
# full.trees <- readRDS(here::here("Project_0","results", "reg_trees_bias.RDS"))

cairo_pdf(here("Project_0","manuscript", paste0("SR_a_bias_regtree_plots.pdf")), width = 30*2/3, height = 15*2/3)
x <- matrix(c(1,1,2,2,3,3,0,4,4,5,5,0), 2, 6, byrow = TRUE)
# x <- matrix(c(1,1,2,3,1,1,4,5), 2, 4, byrow = TRUE)
layout.x <- layout(x) 
par(oma = c(0,0,0,0), bg = NA)
#par(mfrow = c(1,2))
plot.prune(prune(full.trees[["a"]][["R"]],3), cp = 0.005, type = "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.4, split.yshift = 0)
mtext("R OMs", side = 3, line = 0, cex = 2)
plot.prune(full.trees[["a"]][["R+S"]], cp = 0.01, type = "R+S", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.4, split.yshift = 0)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["a"]][["R+M"]], c(13,7)), cp = 0.001, type = "R+M", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 3, split.yshift = 1)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["a"]][["R+Sel"]],c(4,5,3)), cp = 0.0001, type = "R+Sel", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.4, split.yshift = 0)
mtext("R+Sel OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["a"]][["R+q"]], c(4,5,3)), cp = 0.0001, type = "R+q", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.4, split.yshift = 0)
mtext("R+q OMs", side = 3, line = 0, cex = 2)
dev.off()


# pkgload::load_all("c:/work/rpart.plot")
# show_col(viridis::viridis_pal(option = "turbo", begin = 0.2, end = 0.8)(4))
round(PRD.tables$b,2)*100
cairo_pdf(here("Project_0","manuscript", paste0("SR_b_bias_regtree_plots.pdf")), width = 30*2/3, height = 15*2/3)
x <- matrix(c(1,3,2,4,5,5), 2, 3)
layout.x <- layout(x) 
par(oma = c(0,2,0,0))
plot.prune(prune(full.trees[["b"]][["R"]], c(2,12,13,7)), cp = 0.0001, type = "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2, split.yshift = 0)
mtext("R OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["b"]][["R+S"]],c(8,9,5,3)), cp = 0.0001, type = "R+S", roundint = FALSE, extra = 1, mar = c(0,1,5,0), tweak = 1.4, split.yshift = 0)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["b"]][["R+M"]],c(2,12,7)), cp = 0.001, type = "R+M", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.4, split.yshift = 0)
mtext("R+M OMs", side = 3, line = 1, cex = 2)
plot.prune(full.trees[["b"]][["R+Sel"]], cp = 0.02,type = "R+Sel", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.4, split.yshift = -3)
mtext("R+Sel OMs", side = 3, line = 1, cex = 2)
plot.prune(prune(full.trees[["b"]][["R+q"]], c(3,5,9)), cp = 0.0001,type = "R+q", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.4, split.yshift = -3)
mtext("R+q OMs", side = 3, line = 1, cex = 2)
dev.off()


# pkgload::load_all("c:/work/rpart.plot")
# show_col(viridis::viridis_pal(option = "turbo", begin = 0.2, end = 0.8)(4))
round(PRD.tables$M,2)*100
cairo_pdf(here("Project_0","manuscript", paste0("med_M_bias_regtree_plots.pdf")), width = 30*2/3, height = 15*2/3)
x <- matrix(c(rep(1,4),rep(2,4),rep(3,4),4,5), 2, 7)
layout.x <- layout(x) 
par(oma = c(0,0,0,1))
plot.prune(prune(full.trees[["M"]][["R"]], c(3,5,16:17,18:19)), cp = 100, factor = "n", type = "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.1, split.yshift = -6)
mtext("R OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["M"]][["R+S"]],c(2,12,14:15,27,28:29,52,53)), cp = 0, type = "R+S", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.1, split.yshift = -6)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["M"]][["R+M"]], c(8,18,38:39,5,3)), cp = 0, type = "R+M", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.1, split.yshift = -6)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
# par(mfrow = c(1,1), oma = c(0,0,0,0))
plot.prune(prune(full.trees[["M"]][["R+Sel"]],2), cp = 0.01,type = "R+Sel", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2, split.yshift = -6)
mtext("R+Sel OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["M"]][["R+q"]],4), cp = 0.03, type = "R+q", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2, split.yshift = -3)
mtext("R+q OMs", side = 3, line = 0, cex = 2)
dev.off()

PRD.tables$SSB
cairo_pdf(here("Project_0","manuscript", paste0("SSB_bias_regtree_plots.pdf")), width = 30*2/3, height = 15*2/3)
x <- matrix(c(1,3,rep(2,4),rep(4,4), rep(5,4)), 2, 7)
layout.x <- layout(x) 
par(mar = c(0,0,5,0), oma = c(0,4,0,1))
plot.prune(prune(full.trees[["SSB"]][["R"]],c(5,15)), cp = 0.01, type = "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["SSB"]][["R+S"]],c(4,10,15)), cp = 0.01, type = "R+S", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.4)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["SSB"]][["R+M"]],c(4,7,9,10,15)), cp = 0.01, type = "R+M", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.4)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["SSB"]][["R+Sel"]],c(6,15,29)), cp = 0.01,type = "R+Sel", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.42)
mtext("R+Sel OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["SSB"]][["R+q"]],c(19,31,60)), cp = 0.01, type = "R+q", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R+q OMs", side = 3, line = 0, cex = 2)
dev.off()

cairo_pdf(here("Project_0","manuscript", paste0("F_bias_regtree_plots.pdf")), width = 30*2/3, height = 15*2/3)
x <- matrix(c(1,3,rep(2,4),rep(4,4), rep(5,4)), 2, 7)
layout.x <- layout(x) 
par(mar = c(0,0,5,0), oma = c(0,4,0,1))
plot.prune(prune(full.trees[["F"]][["R"]],c(5,8,15)), cp = 0.01, type = "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["F"]][["R+S"]],c(4,10,15)), cp = 0.01, type = "R+S", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.4)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["F"]][["R+M"]],c(4,7,9,10,15)), cp = 0.01, type = "R+M", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.4)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["F"]][["R+Sel"]],c(6,8,15,18,29)), cp = 0.01,type = "R+Sel", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.42)
mtext("R+Sel OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["F"]][["R+q"]],c(16,19,31,34,35,60)), cp = 0.01, type = "R+q", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R+q OMs", side = 3, line = 0, cex = 2)
dev.off()


cairo_pdf(here("Project_0","manuscript", paste0("R_bias_regtree_plots.pdf")), width = 30*2/3, height = 15*2/3)
x <- matrix(c(1,1,2,2,3,3,4,4,5,5), 2, 5)
layout.x <- layout(x) 
par(mar = c(0,0,5,0), oma = c(0,4,0,1))
plot.prune(prune(full.trees[["R"]][["R"]],c(5,15)), cp = 0.01, type = "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["R"]][["R+S"]],c(4,10,15)), cp = 0.01, type = "R+S", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.4)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["R"]][["R+M"]], c(15,24)), cp = 0.001, type = "R+M", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.4)
#plot.prune(prune(full.trees[["R"]][["R+M"]],c(4,7,9,10,15)), cp = 0.01, type = "R+M", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.4)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["R"]][["R+Sel"]],c(6,15,29)), cp = 0.01,type = "R+Sel", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.42)
mtext("R+Sel OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["R"]][["R+q"]],c(19,31,60)), cp = 0.01, type = "R+q", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R+q OMs", side = 3, line = 0, cex = 2)
dev.off()
