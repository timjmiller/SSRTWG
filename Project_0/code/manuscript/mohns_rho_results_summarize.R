library(here)
library(dplyr)
library(tidyr)
library(Hmisc)
library(reshape2)
library(rpart)
library(viridis)
library(gratia)
#modified rpart.plot package to use plotmath in split labels...
pkgload::load_all(file.path(here(),"../../rpart.plot"))
# pkgload::load_all("/home/tmiller/FSAM_research/aug_backup/work/rpart.plot")


make_plot_df <- function(om_type = "naa", res = naa_relSR_results) {
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
  print("df.ems")
  print(df.ems)
  
  res <- melt(res)
  names(res) <- c("rho", "em", "sim","om")
  res$Type = c("SSB","F","R")
  print(res[1:20,])

  print("unique(res$em)")
  print(unique(res$em))
  print(all(unique(res$em) %in% 1:NROW(df.ems)))
  print(length(unique(res$em)))
  print(length(unique(res$em)))
  print(NROW(df.ems))
  if(!(all(unique(res$em) %in% 1:NROW(df.ems)) & length(unique(res$em)) == NROW(df.ems))) stop("df.ems does not seem to be difined correctly for the number of ems in res.")
  res <- cbind(df.ems[res$em,], res)

  df <- res

  if(om_type == "naa") {
    df.oms <- readRDS(here::here("Project_0", "inputs", "df.oms.RDS"))
  } else {
    df.oms <- readRDS(here::here("Project_0", "inputs", paste0("df.", om_type, ".oms.RDS")))
  }
  df <- cbind(df, df.oms[df$om,])
  
  df <- df %>%
    mutate(EM_process_error = recode(re_config,
                                     "rec" = "R",
                                     "rec+1" = "R+S",
                                     "M_re" = "R+M",
                                     "q_re" = "R+q",
                                     "sel_re" = "R+Sel"
  ))
  
  df <- df %>% mutate(SR_model = recode(SR_model,
                                        "2" = "None",
                                        "3" = "Estimated"
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


all_naa_mohns_rho <-  readRDS(file = here("Project_0","results", "all_naa_mohns_rho_results.RDS"))
all_Sel_mohns_rho <-  readRDS(file = here("Project_0","results", "all_Sel_mohns_rho_results.RDS"))
all_M_mohns_rho <-  readRDS(file = here("Project_0","results", "all_M_mohns_rho_results.RDS"))
all_q_mohns_rho <-  readRDS(file = here("Project_0","results", "all_q_mohns_rho_results.RDS"))

obs_dfs <- list()
obs_dfs$naa <- make_plot_df(om_type = "naa", res = all_naa_mohns_rho)
obs_dfs$M <- make_plot_df(om_type = "M", res = all_M_mohns_rho)
obs_dfs$Sel <- make_plot_df(om_type = "Sel", res = all_Sel_mohns_rho)
obs_dfs$q <- make_plot_df(om_type = "q", res = all_q_mohns_rho)


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
      labs <- gsub("NAA sigma", "sigma['2+']", labs, fixed = TRUE) #sigma_2+
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
  origx$frame$median_yval <- NA
  all_node_names <- node_names <- as.numeric(rownames(newx$frame))
  toss <- max(node_names %/% 2)
  nloop <- 0
  data <- subset(data, !is.na(relerror_trans))
  while(toss > -1){
    leaf_rows <- sort(which(newx$frame$var == "<leaf>"))
    for(i in leaf_rows){
      origx_index <- which(all_node_names == node_names[i])
      origx$frame$n_test[origx_index] <- NROW(data[which(newx$where == i),])
      origx$frame$yval_test[origx_index] <- mean(data$relerror_trans[which(newx$where == i)])
      origx$frame$median_yval[origx_index] <- median(data$relerror_trans[which(newx$where == i)])
      x <- data$rho[which(newx$where == i)]
      origx$frame$mean_RE[origx_index] <- mean(data$rho[which(newx$where == i)])
      origx$frame$median_RE[origx_index] <- median(data$rho[which(newx$where == i)])
      origx$frame$mean_abs_RE[origx_index] <- mean(abs(data$rho[which(newx$where == i)]))
      origx$frame$median_abs_RE[origx_index] <- median(abs(data$rho[which(newx$where == i)]))
    }
    if(max(node_names %/% 2) == toss & toss == 0) break
    toss <- max(node_names %/% 2)
    newx <- suppressWarnings(rpart::snip.rpart(newx, toss))
    node_names <- as.numeric(rownames(newx$frame))
  }
  return(origx)
}

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

get_small_data <- function(Type = "SSB", OM_type, obs_dfs, factors){
  print(OM_type)
  facs <- factors[[Type]][[OM_type]]
  dfs <- obs_dfs
  if(OM_type == "R") temp <- subset(dfs[["naa"]], OM_NAA_sigma == 0 )
  if(OM_type == "R+S") temp <- subset(dfs[["naa"]], OM_NAA_sigma != 0)
  if(OM_type == "R+M") temp <- dfs[["M"]]
  if(OM_type == "R+Sel") temp <- dfs[["Sel"]]
  if(OM_type == "R+q") temp <- dfs[["q"]]
  # temp <- subset(temp, Type == Type)
  print(head(temp))
  print(dim(temp))
  temp$relerror_trans <- log(temp$rho + 1)
  temp$relerror_trans[which(is.infinite(temp$relerror_trans))] <- NA
  return(temp)
}

########################################
#define factors for fits
factors <- list()
factors[["SSB"]] <- list()
factors[["SSB"]][["R"]] <- c("1", "EM_process_error","EM_M","OM_R_sigma","SR_model", "OM_Obs._Error", "OM_F_History")
factors[["SSB"]][["R+S"]] <- c("1", "EM_process_error","EM_M","OM_R_sigma","SR_model", "OM_Obs._Error", "OM_F_History","OM_NAA_sigma")
factors[["SSB"]][["R+M"]] <- c("1", "EM_process_error","EM_M","SR_model","OM_Obs._Error", "OM_F_History","OM_M_sigma", "OM_M_rho")
factors[["SSB"]][["R+Sel"]] <- c("1", "EM_process_error","EM_M","SR_model","OM_Obs._Error", "OM_F_History","OM_Sel_sigma", "OM_Sel_rho")
factors[["SSB"]][["R+q"]] <- c("1", "EM_process_error","EM_M","SR_model","OM_Obs._Error", "OM_F_History","OM_q_sigma", "OM_q_rho")

########################################

glm_fits <- dev.tables <- PRD.tables <- full.trees <-list(SSB = list(),F = list(),R = list())

for(OM_type in names(factors[["SSB"]])){
  print(OM_type)
  facs <- factors[["SSB"]][[OM_type]]
  temp <- get_small_data(Type = "SSB", OM_type, obs_dfs, factors)
  for(par_type in c("F","R","SSB")){
    temp_p <- subset(temp, Type == par_type)
    glm_fits[[par_type]][[OM_type]] <- list()
    dev.tables[[par_type]][[OM_type]] <- list()
    for(i in facs){
      glm_fits[[par_type]][[OM_type]][[i]] <- glm(as.formula(paste("relerror_trans", "~", i)), family = gaussian, data = temp_p)
    }
    sapply(glm_fits[[par_type]][[OM_type]][facs[-1]], \(x) anova(x, test = "LRT")[[5]][2])
    glm_fits[[par_type]][[OM_type]][["all"]] <- glm(as.formula(paste("relerror_trans", "~", paste(facs,collapse = "+"))), family = gaussian, data = temp_p)
    glm_fits[[par_type]][[OM_type]][["all2"]] <- glm(as.formula(paste("relerror_trans", "~ (", paste(facs[-1],collapse = "+"), ")^2")), family = gaussian, data = temp_p)
    glm_fits[[par_type]][[OM_type]][["all3"]] <- glm(as.formula(paste("relerror_trans", "~ (", paste(facs[-1],collapse = "+"), ")^3")), family = gaussian, data = temp_p)
    
    #percent reduction in deviance
    dev.tables[[par_type]][[OM_type]] <- sapply(glm_fits[[par_type]][[OM_type]][facs], \(x) 1 - x$deviance/glm_fits[[par_type]][[OM_type]][[1]]$null.deviance)
    
    #Regression trees
    form <- as.formula(paste("relerror_trans ~", paste(facs, collapse = "+")))# (EM_PE + OM_R_SD + EM_M + EM_SR + OM_Obs._Error + OM_F_History)
    full.trees[[par_type]][[OM_type]] <- rpart(form, data=temp_p, method = "anova", control=rpart.control(cp=0, xval = 100), model = TRUE)#, roundint = FALSE)
    print("par_type done")
  }
  print("OM_type done")
}
saveRDS(glm_fits, here::here("Project_0","results", "glm_fits_mohns_rho.RDS"))

for(par_type in c("F","R","SSB")){
  interactions.dev.table <- sapply(names(factors[["SSB"]]), \(x) 1 - c(glm_fits[[par_type]][[x]][["all"]]$deviance,glm_fits[[par_type]][[x]][["all2"]]$deviance,glm_fits[[par_type]][[x]][["all3"]]$deviance)/glm_fits[[par_type]][[x]][[1]]$null.deviance)
  x <- as.data.frame(round(100*interactions.dev.table,2))
  x <- cbind(Model = c("No Interactions", "+ All Two Way", "+ All Three Way"), x)
  
  All.facs <- c("EM_process_error","OM_Obs._Error", "OM_F_History","OM_R_sigma","OM_NAA_sigma","OM_M_sigma", "OM_M_rho","OM_Sel_sigma", "OM_Sel_rho","OM_q_sigma", "OM_q_rho")
  All.facs <- c("EM_M", "SR_model",All.facs)
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
    print(names(dev.tables[[par_type]][[i]]))
    print(All.facs)
    print(length(dev.tables[[par_type]][[i]][match(All.facs, names(dev.tables[[par_type]][[i]]))]))
    PRD.table[,i] <- dev.tables[[par_type]][[i]][match(All.facs, names(dev.tables[[par_type]][[i]]))]
  }
  y <- interactions.dev.table
  rownames(y) <- c("All factors", "+ All Two Way", "+ All Three Way")
  PRD.tables[[par_type]] <- rbind(PRD.table,y)
}

  
for(OM_type in names(factors[["SSB"]])){
  print(OM_type)
  temp <- get_small_data(Type = "SSB", OM_type, obs_dfs, factors)
  for(par_type in c("F","R","SSB")){
    temp_p <- subset(temp, Type == par_type)
    full.trees[[par_type]][[OM_type]] <- add_to_frame(full.trees[[par_type]][[OM_type]], temp_p)
  }
}
saveRDS(full.trees, here::here("Project_0","results", "reg_trees_mohns_rho.RDS"))


x <- PRD.tables[["SSB"]]
x[] <- format(round(100*x,2), nsmall = 2)
dim(x)
x[which(is.na(as.numeric(x)))] <- "--"
x[which(as.numeric(x) == 0)] <- "< 0.01"
x <- latex(x, file = here("Project_0","manuscript","mohns_rho_SSB_PRD_table.tex"), 
           table.env = FALSE, col.just = rep("r", dim(x)[2]), rowlabel = "Factor", rowlabel.just = "l")#, rowname = NULL)

x <- PRD.tables[["F"]]
x[] <- format(round(100*x,2), nsmall = 2)
dim(x)
x[which(is.na(as.numeric(x)))] <- "--"
x[which(as.numeric(x) == 0)] <- "< 0.01"
x <- latex(x, file = here("Project_0","manuscript","mohns_rho_F_PRD_table.tex"), 
           table.env = FALSE, col.just = rep("r", dim(x)[2]), rowlabel = "Factor", rowlabel.just = "l")#, rowname = NULL)

x <- PRD.tables[["R"]]
x[] <- format(round(100*x,2), nsmall = 2)
dim(x)
x[which(is.na(as.numeric(x)))] <- "--"
x[which(as.numeric(x) == 0)] <- "< 0.01"
x <- latex(x, file = here("Project_0","manuscript","mohns_rho_R_PRD_table.tex"), 
           table.env = FALSE, col.just = rep("r", dim(x)[2]), rowlabel = "Factor", rowlabel.just = "l")#, rowname = NULL)

#this is needed inside my modified handle.anova.palette function
anova.palette.sd <- 0.15
anova.palette.sd <- 0.05

PRD.tables$SSB
cairo_pdf(here("Project_0","manuscript", paste0("SSB_mohns_rho_regtree_plots.pdf")), width = 30*2/3, height = 15*2/3)
x <- matrix(c(1,5,1,5,2,2,2,2,3,3,3,3,4,4,4,4), 2, 8)
layout.x <- layout(x) 
par(oma = c(0,4,0,1))
plot.prune(prune(full.trees[["SSB"]][["R"]],c(6,14,15)), cp = 0.001, type = "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["SSB"]][["R+S"]],c(6)), cp = 0.001, type = "R+S", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["SSB"]][["R+M"]],c(3)), cp = 0.001, type = "R+M", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["SSB"]][["R+Sel"]],c(14,30,63)), cp = 0.001,type = "R+Sel", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R+Sel OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["SSB"]][["R+q"]],c(14,15)), cp = 0.001, type = "R+q", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R+q OMs", side = 3, line = 0, cex = 2)
dev.off()

cairo_pdf(here("Project_0","manuscript", paste0("F_mohns_rho_regtree_plots.pdf")), width = 30*2/3, height = 15*2/3)
x <- matrix(c(1,0,1,4,2,4,2,5,3,5,3,0), 2, 6)
layout.x <- layout(x) 
par(mar = c(0,0,5,0), oma = c(0,4,0,1))
plot.prune(prune(full.trees[["F"]][["R"]], c(4,12,14,15)), cp = 0.0001, type = "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["F"]][["R+S"]], c(15)), cp = 0.001, type = "R+S", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["F"]][["R+M"]],c(30,31)), cp = 0.0005, type = "R+M", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["F"]][["R+Sel"]],c(16)), cp = 0.001, type = "R+Sel", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R+Sel OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["F"]][["R+q"]],c(12,15)), cp = 0.0005, type = "R+q", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R+q OMs", side = 3, line = 0, cex = 2)
dev.off()

par(mfrow = c(1,2))
plot.prune(full.trees[["R"]][["R+S"]], cp = 0.0001, type = "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
 
cairo_pdf(here("Project_0","manuscript", paste0("R_mohns_rho_regtree_plots.pdf")), width = 30*2/3, height = 15*2/3)
x <- matrix(c(1,0,1,4,2,4,2,5,3,5,3,0), 2, 6)
layout.x <- layout(x) 
par(mar = c(0,0,5,0), oma = c(0,4,0,1))
plot.prune(prune(full.trees[["R"]][["R"]],c(14:15)), cp = 0.0001, type = "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
#plot.prune(prune(full.trees[["R"]][["R"]], c(5,14,15)), cp = 0.0001, type = "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["R"]][["R+S"]], c(11,29)), cp = 0.001, type = "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["R"]][["R+M"]],31), cp = 0.001, type = "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["R"]][["R+Sel"]],c(14,15)), cp = 0.001, type = "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R+Sel OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["R"]][["R+q"]],c(30,31)), cp = 0.001, type = "R", roundint = FALSE, extra = 1, mar = c(0,0,5,0), tweak = 1.2)
mtext("R+q OMs", side = 3, line = 0, cex = 2)
dev.off()
