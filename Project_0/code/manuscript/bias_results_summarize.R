library(here)
library(dplyr)
library(tidyr)
library(Hmisc)
library(reshape2)
library(rpart)
library(viridis)
#modified rpart.plot package to use plotmath in split labels...
pkgload::load_all("c:/work/rpart.plot")


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
  print("df.ems")
  print(df.ems)

  res <- melt(res)
  names(res) <- c("par", "column", "value", "em", "sim","om") #em = 1: M fixed, em = 2: M estimated
  if(!is.null(year)) res <- filter(res, par == year)

  print("unique(res$em)")
  print(unique(res$em))
  print(all(unique(res$em) %in% 1:NROW(df.ems)))
  print(length(unique(res$em)))
  print(length(unique(res$em)))
  print(NROW(df.ems))
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
  df <- cbind(df, df.oms[df$om,])

  if(!is.null(M_or_SR)) {
    if(M_or_SR == "SR"){
      df <- df %>% mutate(par = recode(par,
        "1" = "italic(a)",
        "2" = "italic(b)"
      ))
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
  facs <- c("OM", "OM_F_History","OM_NAA_sigma", "OM_R_sigma", "OM_M_sigma", "OM_M_rho", "OM_Sel_sigma", "OM_Sel_rho", "OM_q_sigma", "OM_q_rho", "OM_Obs._Error","EM_M", "make_plot_df", "correct_EM_PE")
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
      labs <- gsub("NAA SD", "sigma['2+']", labs, fixed = TRUE) #sigma_2+
      labs <- gsub(" M", "*phantom(0)*italic(M)", labs, fixed = TRUE) #italic M
      labs <- gsub("process error", "Process*phantom(0)*Error", labs, fixed = TRUE)
      labs <- gsub("F ", "italic(F)*phantom(0)*", labs, fixed = TRUE)
      labs <- gsub(" ", "*phantom(0)*", labs, fixed = TRUE)
      labs <- gsub("*%", "%", labs, fixed = TRUE)
      labs <- gsub("%*", "%", labs, fixed = TRUE)
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
      labs <- gsub("process error", "Process*phantom(0)*Error", labs, fixed = TRUE)
      labs <- gsub("F ", "italic(F)*phantom(0)*", labs, fixed = TRUE)
      labs <- gsub("EM SR", "EM*phantom(0)*SR*phantom(0)*Assumption", labs, fixed = TRUE)
      labs <- gsub("R sigma", "sigma[R]", labs, fixed = TRUE) #sigma_2+
      labs <- gsub("NAA SD", "sigma['2+']", labs, fixed = TRUE) #sigma_2+
      labs <- gsub(" ", "*phantom(0)*", labs, fixed = TRUE)
      labs <- gsub("*%", "%", labs, fixed = TRUE)
      labs <- gsub("%*", "%", labs, fixed = TRUE)
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
  all_node_names <- node_names <- as.numeric(rownames(newx$frame))
  toss <- max(node_names %/% 2)
  nloop <- 0
  data <- subset(data, !is.na(relerror_trans2))
  while(toss > -1){
    leaf_rows <- sort(which(newx$frame$var == "<leaf>"))
    for(i in leaf_rows){
      origx_index <- which(all_node_names == node_names[i])
      origx$frame$n_test[origx_index] <- NROW(data[which(newx$where == i),])
      origx$frame$yval_test[origx_index] <- mean(data$relerror_trans2[which(newx$where == i)])
      x <- data$relerror[which(newx$where == i)]
      origx$frame$mean_RE[origx_index] <- mean(data$relerror[which(newx$where == i)])
      origx$frame$median_RE[origx_index] <- median(data$relerror[which(newx$where == i)])
      origx$frame$mean_abs_RE[origx_index] <- mean(abs(data$relerror[which(newx$where == i)]))
      origx$frame$median_abs_RE[origx_index] <- median(abs(data$relerror[which(newx$where == i)]))
    }
    if(max(node_names %/% 2) == toss & toss == 0) break
    toss <- max(node_names %/% 2)
    newx <- suppressWarnings(rpart::snip.rpart(newx, toss))
    node_names <- as.numeric(rownames(newx$frame))
  }
  return(origx)
}
  
node.fun <- function(x, labs, digits, varlen)
{
  # out <- paste0("Mean log(|RE|) = ", format(round(x$frame$yval,3), nsmall = 3))
  # out <- paste0("Mean |RE| = ", format(round(x$frame$mean_abs_RE,3), nsmall = 3))
  out <- ""
  if(!is.null(x$frame$median_abs_RE)) out <- paste0(out, "Median |RE| = ", format(round(x$frame$median_abs_RE,3), nsmall = 3))
  if(!is.null(x$frame$median_RE)) out <- paste0(out, "\nMedian RE = ",format(round(x$frame$median_RE,3), nsmall = 3))
  if(out == "") out <- paste0("yval = ", format(round(x$frame$yval,3), nsmall = 3))
  paste0(out, "\nn = ", x$frame$n)
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

temp <- make_plot_df(om_type = "M", res = M_relSR_results)
dim(temp)
dim(obs_dfs$SR$M)
temp <- make_plot_df(om_type = "Sel", res = Sel_relSR_results)
dim(temp)
dim(obs_dfs$SR$Sel)
temp <- make_plot_df(om_type = "q", res = q_relSR_results)
dim(temp)
dim(obs_dfs$SR$q)

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

temp <- make_plot_df(om_type = "naa", res = naa_rel_M_results, M_or_SR = "M")
dim(temp)
dim(obs_dfs$M$naa)
temp <- make_plot_df(om_type = "M", res = M_relS_M_results, M_or_SR = "M")
dim(temp)
dim(obs_dfs$M$M)
temp <- make_plot_df(om_type = "Sel", res = Sel_rel_M_results, M_or_SR = "M")
dim(temp)
dim(obs_dfs$M$Sel)
temp <- make_plot_df(om_type = "q", res = q_rel_M_results, M_or_SR = "M")
dim(temp)
dim(obs_dfs$M$q)

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


########################################

glm_fits <- list(a = list(), b = list(), M = list())
dev.tables <- list(a = list(), b = list(), M = list())
PRD.tables <- list(a = list(), b = list(), M = list())
full.trees <- list(a = list(), b = list(), M = list())

cv_limit <- NA
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
      glm_fits[[SR_par]][[OM_type]][[i]] <- glm(as.formula(paste("relerror_trans2", "~", i)), family = gaussian, data = temp)
    }
    sapply(glm_fits[[SR_par]][[OM_type]][facs[-1]], \(x) anova(x, test = "LRT")[[5]][2])
    glm_fits[[SR_par]][[OM_type]][["all"]] <- glm(as.formula(paste("relerror_trans2", "~", paste(facs,collapse = "+"))), family = gaussian, data = temp)
    glm_fits[[SR_par]][[OM_type]][["all2"]] <- glm(as.formula(paste("relerror_trans2", "~ (", paste(facs[-1],collapse = "+"), ")^2")), family = gaussian, data = temp)
    glm_fits[[SR_par]][[OM_type]][["all3"]] <- glm(as.formula(paste("relerror_trans2", "~ (", paste(facs[-1],collapse = "+"), ")^3")), family = gaussian, data = temp)

    #percent reduction in deviance
    dev.tables[[SR_par]][[OM_type]] <- sapply(glm_fits[[SR_par]][[OM_type]][facs], \(x) 1 - x$deviance/glm_fits[[SR_par]][[OM_type]][[1]]$null.deviance)

    #Regression trees
    form <- as.formula(paste("relerror_trans2 ~", paste(facs, collapse = "+")))# (EM_PE + OM_R_SD + EM_M + EM_SR + OM_Obs._Error + OM_F_History)
    full.trees[[SR_par]][[OM_type]] <- rpart(form, data=temp, method = "anova", control=rpart.control(cp=0, xval = 100), model = TRUE)#, roundint = FALSE)
    full.trees[[SR_par]][[OM_type]] <- add_to_frame(full.trees[[SR_par]][[OM_type]], temp)
    print("OM_type done")
  }

  interactions.dev.table <- sapply(names(factors[[par_type]]), \(x) 1 - c(glm_fits[[SR_par]][[x]][["all"]]$deviance,glm_fits[[SR_par]][[x]][["all2"]]$deviance,glm_fits[[SR_par]][[x]][["all3"]]$deviance)/glm_fits[[SR_par]][[x]][[1]]$null.deviance)
  x <- as.data.frame(round(100*interactions.dev.table,2))
  x <- cbind(Model = c("No Interactions", "+ All Two Way", "+ All Three Way"), x)

  All.facs <- c("EM_process_error","OM_Obs._Error", "OM_F_History","OM_R_sigma","OM_NAA_sigma","OM_M_sigma", "OM_M_rho","OM_Sel_sigma", "OM_Sel_rho","OM_q_sigma", "OM_q_rho")
  if(SR_par %in% c("a","b")) All.facs <- c("EM_M", All.facs)
  if(SR_par %in% c("M")) All.facs <- c("SR_model", All.facs)
  if(SR_par %in% c("SSB")) All.facs <- c("EM_M", "SR_model",All.facs)
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
  rnames <- gsub("Sel sigma", "$\\sigma_{Sel}$", rnames, fixed = TRUE)
  rnames <- gsub("M rho", "$\\rho_R$", rnames, fixed = TRUE)
  rnames <- gsub("Sel rho", "$\\rho_{Sel}$", rnames, fixed = TRUE)
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

saveRDS(glm_fits, here::here("Project_0","results", "glm_fits_bias.RDS"))
saveRDS(full.trees, here::here("Project_0","results", "reg_trees_bias.RDS"))

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

x <- PRD.tables[["SSB"]]
x[] <- format(round(100*x,2), nsmall = 2)
dim(x)
x[which(is.na(as.numeric(x)))] <- "--"
x[which(as.numeric(x) == 0)] <- "< 0.01"
x <- latex(x, file = here("Project_0","manuscript","bias_SSB_PRD_table.tex"), 
  table.env = FALSE, col.just = rep("r", dim(x)[2]), rowlabel = "Factor", rowlabel.just = "l")#, rowname = NULL)

#this is needed inside my modified handle.anova.palette function
anova.palette.sd <- 0.3

cairo_pdf(here("Project_0","manuscript", paste0("SR_a_bias_regtree_plots.pdf")), width = 30*2/3, height = 20*2/3)
x <- matrix(c(1,1,2,2,3,3,0,4,4,5,5,0), 2, 6, byrow = TRUE)
layout.x <- layout(x) 
par(oma = c(0,0,0,0))
plot.prune(full.trees[["a"]][["R"]], cp = 0.005, type = "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0))
mtext("R OMs", side = 3, line = 0, cex = 2)
plot.prune(full.trees[["a"]][["R+S"]], cp = 0.005, type = "R+S", roundint = FALSE, extra = 1, mar = c(0,0,5,0))
mtext("R+S OMs", side = 3, line = 0, cex = 2)
plot.prune(full.trees[["a"]][["R+M"]], cp = 0.015, type = "R+M", roundint = FALSE, extra = 1, mar = c(0,0,5,0))
mtext("R+M OMs", side = 3, line = 0, cex = 2)
plot.prune(full.trees[["a"]][["R+Sel"]], cp = 0.015, type = "R+Sel", roundint = FALSE, extra = 1, mar = c(0,0,5,0))
mtext("R+Sel OMs", side = 3, line = 0, cex = 2)
plot.prune(full.trees[["a"]][["R+q"]], cp = 0.03, type = "R+q", roundint = FALSE, extra = 1, mar = c(0,0,5,0))
mtext("R+q OMs", side = 3, line = 0, cex = 2)
dev.off()

# pkgload::load_all("c:/work/rpart.plot")
# show_col(viridis::viridis_pal(option = "turbo", begin = 0.2, end = 0.8)(4))
PRD.tables$b
cairo_pdf(here("Project_0","manuscript", paste0("SR_b_bias_regtree_plots.pdf")), width = 30*2/3, height = 20*2/3)
x <- matrix(c(1,1,2,2,3,3,0,4,4,5,5,0), 2, 6, byrow = TRUE)
layout.x <- layout(x) 
par(oma = c(0,0,0,0))
plot.prune(full.trees[["b"]][["R"]], cp = 0.01, type = "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0))
mtext("R OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["b"]][["R+S"]],5), cp = 0.01, type = "R+S", roundint = FALSE, extra = 1, mar = c(0,0,5,0))
mtext("R+S OMs", side = 3, line = 0, cex = 2)
plot.prune(full.trees[["b"]][["R+M"]], cp = 0.01, type = "R+M", roundint = FALSE, extra = 1, mar = c(0,0,5,0))
mtext("R+M OMs", side = 3, line = 0, cex = 2)
plot.prune(full.trees[["b"]][["R+Sel"]], cp = 0.01,type = "R+Sel", roundint = FALSE, extra = 1, mar = c(0,0,5,0))
mtext("R+Sel OMs", side = 3, line = 0, cex = 2)
plot.prune(full.trees[["b"]][["R+q"]], cp = 0.01,type = "R+q", roundint = FALSE, extra = 1, mar = c(0,0,5,0))
mtext("R+q OMs", side = 3, line = 0, cex = 2)
dev.off()


# pkgload::load_all("c:/work/rpart.plot")
# show_col(viridis::viridis_pal(option = "turbo", begin = 0.2, end = 0.8)(4))
PRD.tables$M
cairo_pdf(here("Project_0","manuscript", paste0("med_M_bias_regtree_plots.pdf")), width = 30*2/3, height = 20*2/3)
x <- matrix(c(1,1,2,2,3,3,0,4,4,5,5,0), 2, 6, byrow = TRUE)
layout.x <- layout(x) 
par(mar = c(0,0,5,0), oma = c(0,0,0,0))
plot.prune(full.trees[["M"]][["R"]], cp = 0.01, type = "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0))
mtext("R OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["M"]][["R+S"]],5), cp = 0.01, type = "R+S", roundint = FALSE, extra = 1, mar = c(0,0,5,0))
mtext("R+S OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["M"]][["R+M"]],c(4:6)), cp = 0.01, type = "R+M", roundint = FALSE, extra = 1, mar = c(0,0,5,0))
mtext("R+M OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["M"]][["R+Sel"]], c(4,6)), cp = 0.01,type = "R+Sel", roundint = FALSE, extra = 1, mar = c(0,0,5,0))
mtext("R+Sel OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["M"]][["R+q"]], c(9,13)), cp = 0.01, type = "R+q", roundint = FALSE, extra = 1, mar = c(0,0,5,0))
mtext("R+q OMs", side = 3, line = 0, cex = 2)
dev.off()

PRD.tables$SSB
cairo_pdf(here("Project_0","manuscript", paste0("SSB_bias_regtree_plots.pdf")), width = 30*2/3, height = 20*2/3)
x <- matrix(c(1,1,2,2,3,3,4,4,5,5,5,5), 2, 6, byrow = TRUE)
layout.x <- layout(x) 
par(mar = c(0,0,5,0), oma = c(0,0,0,0))
plot.prune(prune(full.trees[["SSB"]][["R"]],5), cp = 0.01, type = "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["SSB"]][["R+S"]],c(4,10)), cp = 0.01, type = "R+S", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.4)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["SSB"]][["R+M"]],c(4,9,10)), cp = 0.01, type = "R+M", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.4)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["SSB"]][["R+Sel"]],6), cp = 0.01,type = "R+Sel", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.4)
mtext("R+Sel OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["SSB"]][["R+q"]],19), cp = 0.01, type = "R+q", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R+q OMs", side = 3, line = 0, cex = 2)
dev.off()

########################################
#Get annual SSB, F relative error
########################################


make_annual_relerror_df <- function(om_type = "naa", res = all_naa_om_relssb, is_SE = FALSE) {
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

  res <- melt(res)
  names(res) <- c("year", "column", "value", "em", "sim","om")
  res <- res %>% mutate(column = recode(column,
      "1" = "relerror",
      "2" = "cv",
      "3" = "in_ci",
    ))
  df <- res %>% tidyr::pivot_wider(names_from = column, values_from = value) %>% as.data.frame
  if(is_SE) df <- filter(df, !is.na(cv))
  
  df$relerror = df$relerror - 1

  df <- df %>% group_by(om, em, year) %>%
    reframe(stats = median_ci_fn(relerror)) %>% as.data.frame
  df$type <- c("lo", "middle", "hi")
  if(om_type == "naa") {
    df.oms <- readRDS(here::here("Project_0", "inputs", "df.oms.RDS"))
  } else {
    df.oms <- readRDS(here::here("Project_0", "inputs", paste0("df.", om_type, ".oms.RDS")))
  }
  df <- df %>% pivot_wider(names_from = type, values_from = stats) %>% as.data.frame
  df <- cbind(df, df.oms[df$om,])
  df <- cbind(df, df.ems[df$em,])
  df <- df %>%
    mutate(EM_process_error = recode(re_config,
      "rec" = "R",
      "rec+1" = "R+S",
      "M_re" = "R+M",
      "q_re" = "R+q",
      "sel_re" = "R+Sel"
    ))
  ind <- which(df$M_re_cor == "iid" | df$sel_re_cor == "iid" | df$q_re_cor == "iid")
  df$EM_process_error[ind] <- paste0(df$EM_process_error[ind], " (iid)")
  ind <- which(df$M_re_cor == "ar1_y" | df$sel_re_cor == "ar1_y" | df$q_re_cor == "ar1")
  df$EM_process_error[ind] <- paste0(df$EM_process_error[ind], " (AR1)")
  EM_process_error <- c("R","R+S","R+M (iid)","R+M (AR1)","R+Sel (iid)","R+Sel (AR1)","R+q (iid)","R+q (AR1)")
  df$EM_process_error <- factor(df$EM_process_error, levels = EM_process_error)
  
  df$correct_EM_PE <- "No"
  if(om_type == "naa") {
    df$NAA_sig[which(is.na(df$NAA_sig))] <- 0
    df$correct_EM_PE[df$NAA_sig == 0 & df$EM_process_error == "R"] <- "Yes"
    df$correct_EM_PE[df$NAA_sig >  0 & df$EM_process_error == "R+S"] <- "Yes"
    df <- df %>% mutate(NAA_sig = recode(NAA_sig,
        "0" = "sigma['2+'] == 0",
        "0.25" = "sigma['2+'] == 0.25",
        "0.5" = "sigma['2+'] == 0.5")) %>%
      mutate(R_sig = recode(R_sig,
        "0.5" = "sigma[italic(R)] == 0.5",
        "1.5" = "sigma[italic(R)] == 1.5"))
  }

  if(om_type == "M") {
    df$correct_EM_PE[df$M_cor == 0 & df$EM_process_error == "R+M (iid)"] <- "Yes"
    df$correct_EM_PE[df$M_cor >  0 & df$EM_process_error == "R+M (AR1)"] <- "Yes"
    df <- df %>% mutate(M_sig = recode(M_sig,
        "0.1" = "sigma[italic(M)] == 0.1",
        "0.5" = "sigma[italic(M)] == 0.5")) %>% 
      mutate(M_cor = recode(M_cor,
        "0" = "rho[italic(M)] == 0",
        "0.9" = "rho[italic(M)] == 0.9"))
  }
  if(om_type == "Sel") {
    df$correct_EM_PE[df$Sel_cor == 0 & df$EM_process_error == "R+Sel (iid)"] <- "Yes"
    df$correct_EM_PE[df$Sel_cor >  0 & df$EM_process_error == "R+Sel (AR1)"] <- "Yes"
    df <- df %>% mutate(Sel_sig = recode(Sel_sig,
        "0.1" = "sigma[Sel] == 0.1",
        "0.5" = "sigma[Sel] == 0.5")) %>% 
      mutate(Sel_cor = recode(Sel_cor,
        "0" = "rho[Sel] == 0",
        "0.9" = "rho[Sel] == 0.9"))
  }
  if(om_type == "q") {
    df$correct_EM_PE[df$q_cor == 0 & df$EM_process_error == "R+q (iid)"] <- "Yes"
    df$correct_EM_PE[df$q_cor >  0 & df$EM_process_error == "R+q (AR1)"] <- "Yes"
    df <- df %>% mutate(q_sig = recode(q_sig,
        "0.1" = "sigma[italic(q)] == 0.1",
        "0.5" = "sigma[italic(q)] == 0.5")) %>%
      mutate(q_cor = recode(q_cor,
        "0" = "rho[italic(q)] == 0",
        "0.9" = "rho[italic(q)] == 0.9"))  
  }
  df <- df %>%
    mutate(obs_error = recode(obs_error,
      "L" = "Low observation error",
      "H" = "High observation error"))
  df <- df %>% mutate(Fhist = recode(Fhist,
        "H-MSY"  = "2.5*italic(F)[MSY] %->% italic(F)[MSY]",
        "MSY" = "italic(F)[MSY]"))
  df <- df %>% as.data.frame
  facs <- c("om", "Fhist","NAA_sig", "R_sig", "M_sig", "M_cor", "Sel_sig", "Sel_cor", "q_sig", "q_cor", "obs_error", "correct_EM_PE")
  fac_names <- names(df)[names(df) %in% facs]
  df[fac_names] <- lapply(df[fac_names], as.factor)
  df$correct_EM_PE[df$correct_EM_PE== "No"] <- NA
  est_config <- c("M known, No SR", "M estimated, No SR", "M known, SR", "M estimated, SR")
  df$est_config <- est_config[1]
  df$est_config[df$M_est & df$SR_model == 2] <- est_config[2]
  df$est_config[!df$M_est & df$SR_model == 3] <- est_config[3]
  df$est_config[df$M_est & df$SR_model == 3] <- est_config[4]
  df$est_config <- factor(df$est_config, levels = est_config)
  return(df)
}



########################################
#terminal SSB
########################################

all_naa_om_relssb <- readRDS(file = here("Project_0","results", "all_naa_relssb_results.RDS"))

ylims <- c(-0.4,0.4)
plot.df <- make_annual_relerror_df()
plot.df$outside <- plot.df$hi < min(ylims) | plot.df$lo > max(ylims)

temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_S_df <- temp

all_M_relssb <- readRDS(file = here("Project_0","results", "all_M_relssb_results.RDS"))
plot.df <- make_annual_relerror_df("M", all_M_relssb)
plot.df$outside <- plot.df$hi < min(ylims) | plot.df$lo > max(ylims)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_M_df <- temp

all_Sel_relssb <- readRDS(file = here("Project_0","results", "all_Sel_relssb_results.RDS"))
plot.df <- make_annual_relerror_df("Sel", all_Sel_relssb)
plot.df$outside <- plot.df$hi < min(ylims) | plot.df$lo > max(ylims)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_Sel_df <- temp

all_q_relssb <- readRDS(file = here("Project_0","results", "all_q_relssb_results.RDS"))
plot.df <- make_annual_relerror_df("q", all_q_relssb)
plot.df$outside <- plot.df$hi < min(ylims) | plot.df$lo > max(ylims)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_q_df <- temp

theme_set(theme_bw())
theme_update(strip.text = element_text(size = rel(2)), strip.placement = "outside", strip.background = element_rect(), #fill = "transparent"), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)))

R_S_plt <- ggplot(filter(R_S_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ R_sig + NAA_sig, labeller = labeller(Fhist = label_parsed, R_sig = label_parsed, NAA_sig = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    coord_cartesian(ylim = ylims) + ylab("Median Relative Error (SSB)") + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_S_plt


R_M_plt <- ggplot(filter(R_M_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ M_sig + M_cor, labeller = labeller(Fhist = label_parsed, M_sig = label_parsed, M_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    coord_cartesian(ylim = ylims) + ylab("Median Relative Error (SSB)") + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_M_plt



R_Sel_plt <- ggplot(filter(R_Sel_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ Sel_sig + Sel_cor, labeller = labeller(Fhist = label_parsed, Sel_sig = label_parsed, Sel_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    coord_cartesian(ylim = ylims) + ylab("Median Relative Error (SSB)") + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_Sel_plt



R_q_plt <- ggplot(filter(R_q_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ q_sig + q_cor, labeller = labeller(Fhist = label_parsed, q_sig = label_parsed, q_cor = label_parsed)) +
    coord_cartesian(ylim = ylims) + ylab("Median Relative Error (SSB)") + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_q_plt


R_Sel_plt_alt <- R_Sel_plt + 
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21, show.legend = TRUE) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) 


theme_set(theme_bw())
theme_update(strip.placement = "outside", strip.background = element_rect(), title = element_text(size = rel(1.5)))
cairo_pdf(here("Project_0","manuscript", "term_SSB_bias_plots.pdf"), width = 30*2/3, height = 20*2/3)
design <- c(area(1,1,1,1), area(2,1,2,1), area(1,2,1,2), area(2,2,2,2))
(R_S_plt +  xlab("") + labs(title = "A: R, R+S OMs") + theme(axis.title.x = element_blank(), strip.text.y=element_blank(), axis.text.x=element_blank(), legend.position="none")) +
(R_M_plt + labs(title = "B: R+M OMs") + theme(legend.position="none", strip.text.y=element_blank())) + 
(R_Sel_plt_alt + xlab("") + labs(title = "C: R+Sel OMs") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())) +
(R_q_plt + labs(title = "D: R+q OMs") + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position = "none")) + 
  plot_layout(design = design, axis_titles = "collect")
dev.off()

########################################
#terminal F
########################################

all_naa_relF <- readRDS(file = here("Project_0","results", "all_naa_relF_results.RDS"))
plot.df <- make_annual_relerror_df("naa", all_naa_relF)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_S_df <- temp


all_M_relF <- readRDS(file = here("Project_0","results", "all_M_relF_results.RDS"))
plot.df <- make_annual_relerror_df("M", all_M_relF)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_M_df <- temp

all_Sel_relF <- readRDS(file = here("Project_0","results", "all_Sel_relF_results.RDS"))
plot.df <- make_annual_relerror_df("Sel", all_Sel_relF)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_Sel_df <- temp

all_q_relF <- readRDS(file = here("Project_0","results", "all_q_relF_results.RDS"))
plot.df <- make_annual_relerror_df("q", all_q_relF)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_q_df <- temp


R_S_plt <- ggplot(filter(R_S_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ R_sig + NAA_sig, labeller = labeller(Fhist = label_parsed, R_sig = label_parsed, NAA_sig = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    coord_cartesian(ylim = ylims) + ylab(bquote(Median~Relative~Error~(italic(F)))) + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_S_plt


R_M_plt <- ggplot(filter(R_M_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ M_sig + M_cor, labeller = labeller(Fhist = label_parsed, M_sig = label_parsed, M_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    coord_cartesian(ylim = ylims) + ylab(bquote(Median~Relative~Error~(italic(F)))) + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_M_plt



R_Sel_plt <- ggplot(filter(R_Sel_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ Sel_sig + Sel_cor, labeller = labeller(Fhist = label_parsed, Sel_sig = label_parsed, Sel_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    coord_cartesian(ylim = ylims) + ylab(bquote(Median~Relative~Error~(italic(F)))) + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_Sel_plt



R_q_plt <- ggplot(filter(R_q_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ q_sig + q_cor, labeller = labeller(Fhist = label_parsed, q_sig = label_parsed, q_cor = label_parsed)) +
    coord_cartesian(ylim = ylims) + ylab(bquote(Median~Relative~Error~(italic(F)))) + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_q_plt


R_Sel_plt_alt <- R_Sel_plt + 
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21, show.legend = TRUE) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) 


theme_set(theme_bw())
theme_update(strip.placement = "outside", strip.background = element_rect(), title = element_text(size = rel(1.5)))
cairo_pdf(here("Project_0","manuscript", "term_F_bias_plots.pdf"), width = 30*2/3, height = 20*2/3)
design <- c(area(1,1,1,1), area(2,1,2,1), area(1,2,1,2), area(2,2,2,2))
(R_S_plt +  xlab("") + labs(title = "A: R, R+S OMs") + theme(axis.title.x = element_blank(), strip.text.y=element_blank(), axis.text.x=element_blank(), legend.position="none")) +
(R_M_plt + labs(title = "B: R+M OMs") + theme(legend.position="none", strip.text.y=element_blank())) + 
(R_Sel_plt_alt + xlab("") + labs(title = "C: R+Sel OMs") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())) +
(R_q_plt + labs(title = "D: R+q OMs") + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position = "none")) + 
  plot_layout(design = design, axis_titles = "collect")
dev.off()


all_naa_relR <- readRDS(file = here("Project_0","results", "all_naa_relR_results.RDS"))
plot.df <- make_annual_relerror_df("naa", all_naa_relR)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_S_df <- temp


all_M_relR <- readRDS(file = here("Project_0","results", "all_M_relR_results.RDS"))
plot.df <- make_annual_relerror_df("M", all_M_relR)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_M_df <- temp

all_Sel_relR <- readRDS(file = here("Project_0","results", "all_Sel_relR_results.RDS"))
plot.df <- make_annual_relerror_df("Sel", all_Sel_relR)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_Sel_df <- temp

all_q_relR <- readRDS(file = here("Project_0","results", "all_q_relR_results.RDS"))
plot.df <- make_annual_relerror_df("q", all_q_relR)
temp <- filter(plot.df, year %in% c(1,21,40))
temp <- temp %>% mutate(year = recode(year,
    "1" = "Start",
    "21" = "Middle",
    "40" = "End"))
temp$year <- factor(temp$year, levels = c("Start", "Middle", "End"))
R_q_df <- temp


R_S_plt <- ggplot(filter(R_S_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ R_sig + NAA_sig, labeller = labeller(Fhist = label_parsed, R_sig = label_parsed, NAA_sig = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    coord_cartesian(ylim = ylims) + ylab(bquote(Median~Relative~Error~(italic(R)))) + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_S_plt


R_M_plt <- ggplot(filter(R_M_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ M_sig + M_cor, labeller = labeller(Fhist = label_parsed, M_sig = label_parsed, M_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    coord_cartesian(ylim = ylims) + ylab(bquote(Median~Relative~Error~(italic(R)))) + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_M_plt



R_Sel_plt <- ggplot(filter(R_Sel_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ Sel_sig + Sel_cor, labeller = labeller(Fhist = label_parsed, Sel_sig = label_parsed, Sel_cor = label_parsed)) +#, obs_error = label_wrap_gen(width=30))) + #labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    coord_cartesian(ylim = ylims) + ylab(bquote(Median~Relative~Error~(italic(R)))) + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_Sel_plt



R_q_plt <- ggplot(filter(R_q_df, year == "End"), aes(x = est_config, y = middle)) + 
    facet_nested(obs_error + Fhist ~ q_sig + q_cor, labeller = labeller(Fhist = label_parsed, q_sig = label_parsed, q_cor = label_parsed)) +
    coord_cartesian(ylim = ylims) + ylab(bquote(Median~Relative~Error~(italic(R)))) + xlab("EM M and Stock Recruit Assumptions") +
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) + 
    scale_size_manual(values = c('No'=0, 'Yes'=6)) + 
    geom_errorbar(aes(ymin = lo, ymax = hi, colour = EM_process_error), width = 0, position = position_dodge(0.7), linewidth = 1) +
    labs(colour = "EM process error", fill = "EM process error") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0), order = 1), size = "none", shape = "none", fill = "none") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
R_q_plt


R_Sel_plt_alt <- R_Sel_plt + 
    geom_point(position = position_dodge(0.7), size =  2, mapping = aes(colour = EM_process_error, fill = EM_process_error), shape = 21, show.legend = TRUE) + 
    geom_point(mapping = aes(fill = EM_process_error, size = correct_EM_PE), shape = 1, position = position_dodge(0.7), na.rm = TRUE) 


theme_set(theme_bw())
theme_update(strip.placement = "outside", strip.background = element_rect(), title = element_text(size = rel(1.5)))
cairo_pdf(here("Project_0","manuscript", "term_R_bias_plots.pdf"), width = 30*2/3, height = 20*2/3)
design <- c(area(1,1,1,1), area(2,1,2,1), area(1,2,1,2), area(2,2,2,2))
(R_S_plt +  xlab("") + labs(title = "A: R, R+S OMs") + theme(axis.title.x = element_blank(), strip.text.y=element_blank(), axis.text.x=element_blank(), legend.position="none")) +
(R_M_plt + labs(title = "B: R+M OMs") + theme(legend.position="none", strip.text.y=element_blank())) + 
(R_Sel_plt_alt + xlab("") + labs(title = "C: R+Sel OMs") + theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())) +
(R_q_plt + labs(title = "D: R+q OMs") + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position = "none")) + 
  plot_layout(design = design, axis_titles = "collect")
dev.off()

