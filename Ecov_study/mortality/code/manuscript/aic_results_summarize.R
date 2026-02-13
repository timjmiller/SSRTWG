library(here)
library(dplyr)
library(tidyr)
library(mgcv)
library(gratia)
library(Hmisc)
library(VGAM)

df.ems = readRDS(here::here("Ecov_study","mortality","inputs", "df.ems.RDS"))
df.oms = readRDS(here::here("Ecov_study","mortality","inputs", "df.oms.RDS"))
aic_res <- readRDS(here::here("Ecov_study","mortality", "results", "aic_results.RDS"))

aic_fn <- function(all, est_ind, om_ind = NULL){
  if(!is.null(om_ind)) all <- all[om_ind]
  out <- lapply(all, function(res){
    tmp <- sapply(res[est_ind], function(x) return(x))
    tmp <- apply(tmp,1, function(x) {
      if(any(!is.na(x))) {
        # return(x == min(x,na.rm=T))
        return(est_ind[which(x == min(x,na.rm=T))])
      } else {
      	# return(rep(NA, length(x)))
      	return(NA)
      }
    })
    # tmp[] <- as.integer(tmp)
    out <- cbind.data.frame(em_ind = tmp, sim = 1:length(tmp))
    return(out)
  })
  return(out)
}

#for a given M assumption this will produced 28800 rows: #OMs x 100 Simulations: One best EM (out of 6 with same M assumption) per simulation.
make_df_fn <- function(M_est = TRUE){
  #all EM PE assumptions
  df <- list()
  df <- lapply(1:3, \(i) {
  # for(i in 1:3) {
    re_mod <- c("rec", "rec+1", "rec+M")[i]
    #EM:  M estimated
    em_ind <- which(df.ems$M_est == M_est) #all EM PE assumptions
    om_ind <- which(df.oms$NAA_M_re == re_mod) #om and em match
    res <- aic_fn(aic_res, em_ind, om_ind)
    print(length(res))
    print(dim(res[[1]]))
    res <- lapply(1:length(om_ind), \(x) cbind.data.frame(om = om_ind[x], res[[x]]))
    print(dim(res[[1]]))
    stop()
    print(length(res))
    res <- do.call(rbind, res)
    print(dim(res))
    res <- cbind(df.oms[res$om,], df.ems[res$em_ind,], res)
    print(dim(res))
    # print(res[1:20,])
    #res <- cbind(df.oms[rep(om_ind, each = length(em_ind)),], df.ems[rep(em_ind, length(om_ind)),], n = c(Mfixed_rec))
    # if(i == 1) {
    #   df <- res
    # } else {
    #   df <- rbind(df, res)
    # }
  })
  df <- do.call(rbind, df)
  print(dim(df))
  # df <- df %>%
  #   mutate(Ecov_obs_sig = recode(Ecov_obs_sig,
  #     "0.1" = "sigma[italic(e)] == 0.1",
  #     "0.5" = "sigma[italic(e)] == 0.5"
  #   ))
  # df <- df %>%
  #   mutate(Ecov_re_sig = recode(Ecov_re_sig,
  #     "0.1" = "sigma[italic(E)] == 0.1",
  #     "0.5" = "sigma[italic(E)] == 0.5"
  #   ))
  # df <- df %>%
  #   mutate(Ecov_re_cor = recode(Ecov_re_cor,
  #     "0" = "rho[italic(E)] == 0",
  #     "0.5" = "rho[italic(E)] == 0.5"
  #   ))
  df <- df %>%
    mutate(Ecov_assumption = recode(as.character(Ecov_est),
      "TRUE" = "beta[italic(E)]*' Estimated'",
      "FALSE" = "beta[italic(E)] == 0"
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
    mutate(EM_beta_ecov = recode(as.character(Ecov_est),
      "TRUE" = "Estimated",
      "FALSE" = "0"
    ))
  df <- df %>%
    mutate(EM_M = recode(as.character(M_est),
      "TRUE" = "Estimated",
      "FALSE" = "Known"
    ))
  df$best <- paste0(df$EM_PE, ", ", ifelse(df$Ecov_est, "Yes","No"))
  df$PE_correct <- factor(ifelse(df$EM_PE == df$OM_PE, "Yes", "No"))
  df$Ecov_correct <- as.character(NA)
  df$Ecov_correct[which(!df$Ecov_est & df$Ecov_effect != 0)] <- "No"
  df$Ecov_correct[which(!df$Ecov_est & df$Ecov_effect == 0)] <- "Yes"
  df$Ecov_correct[which(df$Ecov_est & df$Ecov_effect != 0)] <- "Yes"
  df$Ecov_correct[which(df$Ecov_est & df$Ecov_effect == 0)] <- "No"
  facs <- c("Ecov_obs_sig", "Ecov_re_sig","Ecov_re_cor","OM_F_history", "obs_error", "OM_PE", "EM_PE", "best", "EM_beta_ecov", "Ecov_effect", "PE_correct", "Ecov_correct")#, "Ecov_est")
  df[facs] <- lapply(df[facs], factor)
  df$obs_error  <- factor(df$obs_error, levels = c("Low", "High"))
  return(df)
}
x <- make_df_fn(M_est = TRUE)

get_PRD_tables <- function(glm_fits, factors){
  dev.tables <- PRD.table <- list()
  factors <- factors[factors %in% names(glm_fits[[1]])]
  for(OM_type in names(glm_fits)){
    dev.tables[[OM_type]] <- sapply(glm_fits[[OM_type]][c(factors,paste0("all",c("",2:3)))], 
      \(x) {
        # print(class(x))
      	if(any(class(x) == "vglm")) return(1 - x@criterion$deviance/glm_fits[[OM_type]][[1]]@criterion$deviance)
      	else return(1 - x$deviance/glm_fits[[OM_type]][[1]]$null.deviance)
     	}
    )
  }
  print(dev.tables[[1]])
  PRD.table <- do.call(cbind,dev.tables)
  rnames <- gsub("_", " ", factors, fixed = TRUE)
  rnames <- gsub("PE", "Process Error", rnames, fixed = TRUE)
  rnames <- gsub("conv", "Convergence", rnames, fixed = TRUE)
  rnames <- gsub("F ", "$F$ ", rnames, fixed = TRUE)
  rnames <- gsub("EM M", "EM $M$ Assumption", rnames, fixed = TRUE)
  rnames <- gsub("obs error", "OM Obs. Error", rnames, fixed = TRUE)
  rnames <- gsub("Ecov obs sig", "OM $\\sigma_e$", rnames, fixed = TRUE)
  rnames <- gsub("Ecov re sig", "OM $\\sigma_E$", rnames, fixed = TRUE)
  rnames <- gsub("Ecov re cor", "OM $\\rho_{E}$", rnames, fixed = TRUE)
  rnames <- gsub("Ecov effect", "OM $\\beta_{E}$", rnames, fixed = TRUE)
  rnames <- gsub("Ecov correct", "EM $\\beta_{E}$ Assumption Correct", rnames, fixed = TRUE)
  rnames <- gsub("Process Error correct", "EM Process Error Assumption Correct", rnames, fixed = TRUE)
  rnames <- gsub("EM beta ecov", "EM $\\beta_{E}$ Assumption", rnames, fixed = TRUE)
  rnames <- c(rnames,"All factors", "+ All Two Way", "+ All Three Way")
  rownames(PRD.table) <- rnames
  print(PRD.table)
  return(PRD.table)
}

all_res_mod <- rbind(make_df_fn(M_est = TRUE), make_df_fn(M_est = FALSE))
all_res_mod$best[which(all_res_mod$best == "NA, NA")] <- NA
all_res_mod$best <- factor(all_res_mod$best, levels = c("R, No",    "R, Yes",   "R+M, No",  "R+M, Yes", "R+S, No",  "R+S, Yes"))
levels(all_res_mod$best)

factors <- c("1", "OM_F_history", "obs_error", "Ecov_obs_sig", "Ecov_re_sig","Ecov_re_cor", "Ecov_effect", "EM_M", "Ecov_correct")
vglm_fits <- list()
for(type in levels(all_res_mod$OM_PE)){
	vglm_fits[[type]] <- list()
	temp <- subset(all_res_mod, OM_PE == type)
	for(i in factors){
		vglm_fits[[type]][[i]] <- vglm(as.formula(paste("EM_PE", "~", i)), family = multinomial, data = temp)
	}
  vglm_fits[[type]][["all"]] <- vglm(as.formula(paste("EM_PE", "~", paste(factors,collapse = "+"))), family = multinomial, data = temp)
  vglm_fits[[type]][["all2"]] <- vglm(as.formula(paste("EM_PE", "~ (", paste(factors[-1],collapse = "+"), ")^2")), family = multinomial, data = temp)
  vglm_fits[[type]][["all3"]] <- vglm(as.formula(paste("EM_PE", "~ (", paste(factors[-1],collapse = "+"), ")^3")), family = multinomial, data = temp)
  # #precent reduction in deviance
  # dev.tables[[type]] <- sapply(vglm_fits[[type]], \(x) 1 - x@criterion$deviance/vglm_fits[[type]][[1]]@criterion$deviance)
}

factors_ecov_assumption <- c("1", "OM_F_history", "obs_error", "OM_PE", "Ecov_obs_sig", "Ecov_re_sig","Ecov_re_cor", "EM_M", "PE_correct")
glm_fits <- list()
for(type in levels(all_res_mod$Ecov_effect)){
  temp <- subset(all_res_mod, Ecov_effect == type)
  glm_fits[[type]] <- list()
  for(i in factors_ecov_assumption){
    glm_fits[[type]][[i]] <- glm(as.formula(paste("Ecov_correct", "~", i)), family = binomial, data = temp)
  }
  glm_fits[[type]][["all"]] <- glm(as.formula(paste("Ecov_correct", "~", paste(factors_ecov_assumption,collapse = "+"))), family = binomial, data = temp)
  glm_fits[[type]][["all2"]] <- glm(as.formula(paste("Ecov_correct", "~ (", paste(factors_ecov_assumption[-1],collapse = "+"), ")^2")), family = binomial, data = temp)
  glm_fits[[type]][["all3"]] <- glm(as.formula(paste("Ecov_correct", "~ (", paste(factors_ecov_assumption[-1],collapse = "+"), ")^3")), family = binomial, data = temp)
}

PRD.table <- get_PRD_tables(glm_fits=vglm_fits, factors = factors[-1])
PRD.table <- PRD.table[,c("R","R+S","R+M")]
x <- PRD.table
x[] <- format(round(100*x,2), nsmall = 2)
dim(x)
x[which(is.na(as.numeric(x)))] <- "--"
x[which(as.numeric(x) == 0)] <- "< 0.01"
x <- latex(x, file = here("Ecov_study","mortality","manuscript","aic_PE_PRD_table.tex"), 
  table.env = FALSE, col.just = rep("r", dim(x)[2]), rowlabel = "Factor", rowlabel.just = "l")#, rowname = NULL)

PRD.table <- get_PRD_tables(glm_fits=glm_fits, factors = factors_ecov_assumption[-1])
colnames(PRD.table) <- paste0("$\\beta_E = ", c(0,0.25,0.5), "$")
x <- PRD.table
x[] <- format(round(100*x,2), nsmall = 2)
dim(x)
x[which(is.na(as.numeric(x)))] <- "--"
x[which(as.numeric(x) == 0)] <- "< 0.01"
x <- latex(x, file = here("Ecov_study","mortality","manuscript","aic_Ecov_effect_PRD_table.tex"), 
  table.env = FALSE, col.just = rep("r", dim(x)[2]), rowlabel = "Factor", rowlabel.just = "l")#, rowname = NULL)

#Regression trees
#########################################
#Regression/classification trees plotting functions
split.fun <- function(type = "R") {
  # replace commas with spaces (needed for strwrap)
  fn <-function(x, labs, digits, varlen, faclen){
    print(labs)
    # stop()
    #labs <- gsub("_", " ", labs, fixed = TRUE)
    labs <- gsub(" = ", "==", labs, fixed = TRUE)
    labs <- gsub(" == ", "*phantom(0)==phantom(0)*", labs, fixed = TRUE)
    labs <- gsub(",", "*','*", labs, fixed = TRUE)
    labs <- gsub("EM_M", "EM*phantom(0)*italic(M)", labs, fixed = TRUE) #italic M
    labs <- gsub("EM_beta_ecov", "EM*phantom(0)*beta[italic(E)]", labs, fixed = TRUE)
    labs <- gsub("Ecov_effect", "OM*phantom(0)*beta[italic(E)]", labs, fixed = TRUE)
    labs <- gsub("EM_PE", "EM*phantom(0)*Process*phantom(0)*Error", labs, fixed = TRUE)
    labs <- gsub("OM_PE", "OM*phantom(0)*Process*phantom(0)*Error", labs, fixed = TRUE)
    labs <- gsub("PE_correct", "Correct*phantom(0)*EM*phantom(0)*Process*phantom(0)*Error", labs, fixed = TRUE)
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

node.fun <- function(x, labs, digits, varlen)
{
  levs <- attr(x, "ylevels")
  nlevs <- length(levs)
  n <- x$frame$n
  n.print<- character()
  n.print[which(nchar(n)<5)]<- format(n[which(nchar(n)<5)])
  n.print[which(nchar(n)>4)] <- format(n, big.mark = " ")[which(nchar(n)>4)]
  labs <- sapply(1:NROW(x$frame$yval2), \(y) {
    if(nlevs>2){
      vals <- 100*x$frame$yval2[y,1+1:nlevs]/sum(x$frame$yval2[y,1+1:nlevs])
#      vals <- x$frame$yval2[y,1+nlevs + 1:nlevs]
      ind <- which(vals>10)
      lab <- paste0(levs[ind], ": ", format(round(vals[ind],1), nsmall = 1), collapse = "%,\n")
      if(length(ind)<length(vals)) lab <- paste0(lab, "%,\nOthers < 10%")
    } else{
      lab <- paste0(format(round(x$frame$yval2[y,5]*100,1), nsmall = 1),"%")
    }
    return(lab)
  }) 
  paste0(labs, "\n", n.print)
}


pkgload::load_all(file.path(here(),"../../rpart.plot"))

full.trees <- list()
for(type in levels(all_res_mod$OM_PE)){
	temp_df <- subset(all_res_mod, OM_PE == type)
  form <- as.formula(paste("EM_PE ~", paste(factors, collapse = "+")))
  full.trees[[type]] <- rpart(form, data=temp_df, method = "class", control=rpart.control(cp=0, xval = 100), model= TRUE)
  full.trees[[type]]$correct_level <- type
}

#adjust prior probabilities for EM_PE so that we get some branches
table(temp_df$EM_PE)/sum(table(temp_df$EM_PE))
type = "R+S"
temp_df <- subset(all_res_mod, OM_PE == type & !is.na(EM_PE))
full.trees[[type]] <- rpart(form, data=temp_df, method = "class", control=rpart.control(cp=0, xval = 100), model= TRUE, parms = list(prior =  c(0.1,0.1,0.8)))
full.trees[[type]]$correct_level <- type

anova.palette.sd <- 0.25

cairo_pdf(here("Ecov_study","mortality","manuscript", "AIC_PE_classification_plots.pdf"), width = 30*2/3, height = 15*2/3)
x <- matrix(1:3, 1, 3, byrow = TRUE)
layout.x <- layout(x) 
par(oma = c(0,0,0,0))
#plot.prune(full.trees[["R"]],cp = 0, type = "R", tweak = 1.4, mar = c(0,0,5,0))
plot.prune(prune(full.trees[["R"]], 14),cp = 0, type = "R", tweak = 1.6, mar = c(0,1,5,0.5), legend.x = NA, split.yshift = -3)
mtext("R OMs", side = 3, line = 0, cex = 2)
#plot.prune(full.trees[["R+S"]],cp = 0, type = "R+S", tweak = 1.6, mar = c(0,0,5,5), legend.x = NA, split.yshift = -3, parms = list(prior =  c(0.1,0.1,0.8)))
plot.prune(prune(full.trees[["R+S"]], c(4,5,6,7)),cp = 0, type = "R+S", tweak = 1.6, mar = c(0,0.5,5,0.5), legend.x = NA, split.yshift = -3)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
#plot.prune(full.trees[["R+M"]],cp = 0.001, type = "R+M", tweak = 1.4, mar = c(0,0,5,5), legend.x = NA)
plot.prune(prune(full.trees[["R+M"]], 7),cp = 0.001, type = "R+M", tweak = 1.6, mar = c(0,0,5,1), legend.x = NA, split.yshift = -3)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
dev.off()


#full.trees <- list()
for(type in levels(all_res_mod$Ecov_effect)){
  temp_df <- subset(all_res_mod, Ecov_effect == type)
  form <- as.formula(paste("Ecov_correct ~", paste(factors_ecov_assumption, collapse = "+")))
  full.trees[[type]] <- rpart(form, data=temp_df, method = "class", control=rpart.control(cp=0, xval = 100), model= TRUE)
  full.trees[[type]]$correct_level <- type
}

cairo_pdf(here("Ecov_study","mortality","manuscript", "AIC_Ecov_effect_classification_plots.pdf"), width = 30*2/3, height = 15*2/3)
x <- matrix(1:3, 1, 3, byrow = TRUE)
layout.x <- layout(x) 
par(oma = c(0,2,0,2))
#plot.prune(full.trees[["0"]],cp = 0, type = "R", tweak = 1.6, mar = c(0,0,5,0), roundint = FALSE)
plot.prune(prune(full.trees[["0"]],c(9)),cp = 0, type = "R", tweak = 1.6, mar = c(0,0,5,0), roundint = FALSE, split.yshift = -3)
mtext(expression(beta[E]==0*" OMs"), side = 3, line = 0, cex = 2)
#plot.prune(full.trees[["0.25"]],cp = 0, type = "R", tweak = 1.6, mar = c(0,0,5,0), roundint = FALSE)
plot.prune(prune(full.trees[["0.25"]],c(2,6,14,62)),cp = 0, type = "R", tweak = 1.6, mar = c(0,0,5,0), roundint = FALSE, split.yshift = -3)
mtext(expression(beta[E]==0.25*" OMs"), side = 3, line = 0, cex = 2)
#plot.prune(full.trees[["0.5"]],cp = 0, type = "R", tweak = 1.6, mar = c(0,0,5,0), roundint= FALSE)
plot.prune(prune(full.trees[["0.5"]],c(2,12,13,28,58)),cp = 0, type = "R", tweak = 1.6, mar = c(0,0,5,0), roundint= FALSE, split.yshift = -3)
mtext(expression(beta[E]==0.5*" OMs"), side = 3, line = 0, cex = 2)
dev.off()

x <- subset(all_res_mod, OM_PE == "R+S")
table(x$EM_PE)
table(x$EM_PE)/sum(table(x$EM_PE))
x <- subset(all_res_mod, OM_PE == "R+S" & obs_error == "Low")
table(x$EM_PE)
table(x$EM_PE)/sum(table(x$EM_PE))
x <- subset(all_res_mod, OM_PE == "R+S" & obs_error == "High")
table(x$EM_PE)
table(x$EM_PE)/sum(table(x$EM_PE))

alt.trees <- list()
form <- as.formula(paste("EM_PE ~", paste(factors, collapse = "+")))
for(type in levels(all_res_mod$OM_PE)){
  temp_df <- subset(all_res_mod, OM_PE == type)
  alt.trees[[type]] <- rpart(form, data=temp_df, method = "class", control=rpart.control(cp=0, xval = 100), model= TRUE, parms = list(prior = rep(1,3)/3))
  alt.trees[[type]]$correct_level <- type
}

# y <- rpart(EM_PE ~ obs_error, data=subset(all_res_mod, OM_PE == "R+S"), method = "class", control=rpart.control(cp=0, xval = 100), model= TRUE, parms = list(split = "information"))
# y <- rpart(EM_PE ~ obs_error, data=subset(all_res_mod, OM_PE == "R+S"), method = "class", control=rpart.control(cp=0, xval = 0, minsplit = 1), model= TRUE, parms = list(prior = rep(1,3)/3, split = "information"))
# y <- rpart(EM_PE ~ obs_error, data=subset(all_res_mod, OM_PE == "R+S"), method = "class", control=rpart.control(cp=0, xval = 0, minsplit = 1), model= TRUE, parms = list(prior = rep(1,3)/3))
# y <- rpart(form, data=subset(all_res_mod, OM_PE == "R+S"), method = "class", control=rpart.control(cp=0, xval = 100), model= TRUE, parms = list(prior = rep(1,3)/3))
# 
# y <- rpart(EM_PE ~ obs_error, data=subset(all_res_mod, OM_PE == "R+S"), method = "class", control=rpart.control(cp=0, xval = 100), model= TRUE, parms = list(prior =  rep(1,3)/3))
# 
# pkgload::load_all(file.path(here(),"../../rpart.plot"))
# 
# type = "R+S"
# temp_df <- subset(all_res_mod, OM_PE == type)
# set.seed(123)
# temp_df <- temp_df[sample(1:NROW(temp_df), 10000, replace = FALSE),]
# alt.trees[[type]] <- rpart(form, data=temp_df, method = "class", control=rpart.control(cp=0, xval = 100), model= TRUE)
# alt.trees[[type]]$correct_level <- type
# 
# type = "R+S"
# temp_df <- subset(all_res_mod, OM_PE == type)
# temp_df <- rbind(temp_df,temp_df,temp_df)
# alt.trees[[type]] <- rpart(form, data=temp_df, method = "class", control=rpart.control(cp=0, xval = 100), model= TRUE)
# alt.trees[[type]]$correct_level <- type
# 
# type = "R+S"
# temp_df <- subset(all_res_mod, OM_PE == type)
# alt.trees[[type]] <- rpart(form, data=temp_df, weights = rep(0.1, NROW(temp_df)), method = "class", control=rpart.control(cp=0, xval = 100), model= TRUE)
# alt.trees[[type]]$correct_level <- type
# 
# type = "R+S"
# temp_df <- subset(all_res_mod, OM_PE == type)
# alt.trees[[type]] <- rpart(form, data=temp_df, weights = rep(0.1, NROW(temp_df)), method = "class", control=rpart.control(cp=0, xval = 100), model= TRUE, 
#   parms = list(loss = full.trees[[type]]$parms$loss*0.1))
# alt.trees[[type]]$correct_level <- type
# 
# library(rpartScore)
# type = "R+S"
# temp_df <- subset(all_res_mod, OM_PE == type)
# alt.trees[[type]] <- rpartScore(form, data=temp_df, split = "quad", control=rpart.control(cp=0, xval = 100), model= TRUE)
# alt.trees[[type]]$correct_level <- type
# 
# type = "R+S"
# temp_df <- subset(all_res_mod, OM_PE == type)
# set.seed(123)
# temp_df <- temp_df[sample(1:NROW(temp_df), 15000, replace = TRUE),]
# alt.trees[[type]] <- rpart(form, data=temp_df, method = "class", control=rpart.control(cp=0, xval = 100), model= TRUE)
# alt.trees[[type]]$correct_level <- type
# 
# library(party)
# type = "R+S"
# temp_df <- subset(all_res_mod, OM_PE == type & !is.na(EM_PE))
# temp_df$EM_M <- factor(temp_df$EM_M)
# alt.trees[[type]] <- ctree(form, data=temp_df)
# alt.trees[[type]]$correct_level <- type
# 
# 
# plot(alt.trees[[type]])
# 
# table(temp_df$EM_PE)/sum(table(temp_df$EM_PE))
# type = "R+S"
# temp_df <- subset(all_res_mod, OM_PE == type & !is.na(EM_PE))
# alt.trees[[type]] <- rpart(form, data=temp_df, method = "class", control=rpart.control(cp=0, xval = 100), model= TRUE, parms = list(prior =  c(0.1,0.1,0.8)))
# alt.trees[[type]]$correct_level <- type
# 
# 
# type = "R"
# temp_df <- subset(all_res_mod, OM_PE == type)
# temp_df <- rbind(temp_df,temp_df,temp_df)
# alt.trees[[type]] <- rpart(form, data=temp_df, method = "class", control=rpart.control(cp=0, xval = 100), model= TRUE)
# alt.trees[[type]]$correct_level <- type
# 
# par(mfrow = c(2,1))
# plot.prune(prune(full.trees[[type]], 7),cp = 0, type = "R", tweak = 1, mar = c(0,0,5,0), legend.x = NA, split.yshift = -3)
# plot.prune(prune(alt.trees[[type]], c(4,5,6,7)),cp = 0, type = "R", tweak = 1, mar = c(0,0,5,0), legend.x = NA, split.yshift = -3)
