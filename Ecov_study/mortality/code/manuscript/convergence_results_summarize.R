library(here)
library(dplyr)
library(tidyr)
library(mgcv)
library(gratia)
library(Hmisc)
#modified rpart.plot package to use plotmath in split labels...
pkgload::load_all(file.path(here(),"../../rpart.plot"))

#convergence information:
#1: 0/1 model finished optimizing
#2: nlminb convergence flag
#3: max gradient value
#4: number of NaNs in SEs for parameters, 0 = good invertible hessian
#5: maximum non-NaN SE estimate

conv_fn <- function(all, em_ind, om_ind = NULL, type = 4){
  if(is.null(om_ind)) om_ind <- 1:length(all)
  # if(!is.null(om_ind)) all <- all[om_ind]
  out <- lapply(om_ind, function(om){
    tmp <- lapply(em_ind, function(em) {
      simres <- all[[om]][[em]]
      if(type == 1) conv <- !is.na(simres[,1]) #catastrophic failure
      if(type == 2) conv <- !is.na(simres[,2]) & simres[,2] == 0  #nlminb convergence flag
      if(type == 3) conv <- !is.na(simres[,4]) & simres[,4] < 1e-6 #gradient
      if(type == 4) conv <- !is.na(simres[,3]) & simres[,3] == 0 #hessian
      if(type == 5) conv <- !is.na(simres[,5]) & simres[,3] == 0 & simres[,5] < 100 #all SEs less than 100
      return(cbind(em_ind = em, conv = as.integer(conv)))
    })
    tmp <- cbind(om_ind = om, do.call(rbind, tmp))
    return(tmp)
  })
  out <- do.call(rbind, out)
  return(out)
}

# re_mods <- c("rec", "rec+1", "rec+M")
# i = j = 1
# om_ind <- which(df.oms$NAA_M_re == re_mods[i])
# em_ind <- which(df.ems$Ecov_est == TRUE & df.ems$M_est == TRUE & df.ems$re_config == re_mods[j])
# conv <- conv_fn(conv_res, em_ind, om_ind)

conv_res_plotting_fn <- function(conv_res, M_est = TRUE, Ecov_est = TRUE){
  re_mods <- c("rec", "rec+1", "rec+M")
  df <- lapply(1:3, \(i) { 
    out <- lapply(1:3, \(j) {
      om_ind <- which(df.oms$NAA_M_re == re_mods[i])
      em_ind <- which(df.ems$Ecov_est == Ecov_est & df.ems$M_est == M_est & df.ems$re_config == re_mods[j])
      conv <- as.data.frame(conv_fn(conv_res, em_ind, om_ind))
      conv <- cbind(df.oms[conv$om_ind,], df.ems[conv$em_ind,], conv)
      return(conv)
    })
    out <- do.call(rbind,out)
    print(dim(out))
    return(out)
  })
  df <- do.call(rbind, df)
  print(dim(df))
  print(head(df))
  df <- df %>%
    mutate(obs_error = recode(obs_error,
                              "L" = "Low",
                              "H" = "High"
    ))

  df <- df %>%
    mutate(EM_M = recode(as.character(M_est),
                                 "TRUE" = "Estimated",
                                 "FALSE" = "Known"
    ))
  df <- df %>%
    mutate(EM_beta_ecov = recode(as.character(Ecov_est),
                                    "TRUE" = "Estimated",
                                    "FALSE" = "0"
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
  facs <- c("Ecov_obs_sig", "OM_F_history", "Ecov_re_sig","Ecov_re_cor", "obs_error", "OM_PE", "EM_PE","EM_beta_ecov", "EM_M", "Ecov_effect")
  df[facs] <- lapply(df[facs], factor)
  df$OM_PE <- factor(df$OM_PE, levels= c("R","R+S", "R+M"))
  df$EM_PE <- factor(df$EM_PE, levels= c("R","R+S", "R+M"))
  return(df)
}

df.ems = readRDS(here("Ecov_study","mortality","inputs", "df.ems.RDS"))
df.oms = readRDS(here("Ecov_study","mortality","inputs", "df.oms.RDS"))
conv_res <- readRDS(here("Ecov_study","mortality", "results", "convergence_results.RDS"))

df <- rbind(
  conv_res_plotting_fn(conv_res, M_est = FALSE, Ecov_est = TRUE),
  conv_res_plotting_fn(conv_res, M_est = TRUE, Ecov_est = TRUE),
  conv_res_plotting_fn(conv_res, M_est = FALSE, Ecov_est = FALSE),
  conv_res_plotting_fn(conv_res, M_est = TRUE, Ecov_est = FALSE))

glm_fits <- list()
dev.tables <- list(PRD = list(), PPRD = list())
# null <- gam(conv ~ 1, family = binomial, method = "REML", data = temp)
#factors <- c("1", "Ecov_obs_sig", "OM_F_history", "Ecov_re_sig","Ecov_re_cor", "obs_error", "EM_PE","EM_beta_ecov", "EM_M", "Ecov_effect")
factors <- c("1", "OM_F_history", "obs_error", "Ecov_obs_sig", "Ecov_re_sig","Ecov_re_cor", "Ecov_effect", "EM_PE","EM_beta_ecov", "EM_M")

for(type in levels(df$OM_PE)){
	temp <- subset(df, OM_PE == type)
	glm_fits[[type]] <- list()
	dev.tables[[type]] <- list()
	for(i in factors){
		glm_fits[[type]][[i]] <- glm(as.formula(paste("conv", "~", i)), family = binomial, data = temp)
	}
	sapply(glm_fits[[type]][factors[-1]], \(x) anova(x, test = "LRT")[[5]][2])
	glm_fits[[type]][["all"]] <- glm(as.formula(paste("conv", "~", paste(factors,collapse = "+"))), family = binomial, data = temp)
	glm_fits[[type]][["all2"]] <- glm(as.formula(paste("conv", "~ (", paste(factors[-1],collapse = "+"), ")^2")), family = binomial, data = temp)
	glm_fits[[type]][["all3"]] <- glm(as.formula(paste("conv", "~ (", paste(factors[-1],collapse = "+"), ")^3")), family = binomial, data = temp)

	#precent reduction in deviance
	dev.tables[["PRD"]][[type]] <- sapply(glm_fits[[type]][factors], \(x) 1 - x$deviance/glm_fits[[type]][[1]]$null.deviance)
	#precent "possible" reduction in deviance: measure reduction in deviance relative to that with all factors in the model (no interactions)
	dev.tables[["PPRD"]][[type]] <- (glm_fits[[type]][[1]]$null.deviance - sapply(glm_fits[[type]][factors], \(x) x$deviance))/(glm_fits[[type]][[1]]$null.deviance - glm_fits[[type]][["all"]]$deviance)
}

interactions.dev.table <- sapply(levels(df$OM_PE), \(x) 1 - c(glm_fits[[x]][["all"]]$deviance,glm_fits[[x]][["all2"]]$deviance,glm_fits[[x]][["all3"]]$deviance)/glm_fits[[x]][[1]]$null.deviance )
x <- as.data.frame(round(100*interactions.dev.table,2))
x <- cbind(Model = c("No Interactions", "+ All Two Way", "+ All Three Way"), x)

facs <- factors[-1]

PRD.table <- matrix(NA, length(facs),length(levels(df$OM_PE)))
colnames(PRD.table) <- levels(df$OM_PE)
rnames <- gsub("_", " ", facs, fixed = TRUE)
rnames <- gsub("PE", "Process Error", rnames, fixed = TRUE)
rnames <- gsub("F ", "$F$ ", rnames, fixed = TRUE)
rnames <- gsub("EM M", "EM $M$ Assumption", rnames, fixed = TRUE)
rnames <- gsub("obs error", "OM Obs. Error", rnames, fixed = TRUE)
rnames <- gsub("Ecov obs sig", "OM $\\sigma_e$", rnames, fixed = TRUE)
rnames <- gsub("Ecov re sig", "OM $\\sigma_E$", rnames, fixed = TRUE)
rnames <- gsub("Ecov re cor", "OM $\\rho_{E}$", rnames, fixed = TRUE)
rnames <- gsub("Ecov effect", "OM $\\beta_{E}$", rnames, fixed = TRUE)
rnames <- gsub("EM beta ecov", "EM $\\beta_{E}$ assumption", rnames, fixed = TRUE)
rownames(PRD.table) <- rnames
for(i in levels(df$OM_PE)) PRD.table[,i] <- dev.tables[["PRD"]][[i]][match(facs, names(dev.tables[["PRD"]][[i]]))]
y <- interactions.dev.table
rownames(y) <- c("All factors", "+ All Two Way", "+ All Three Way")
x <- rbind(PRD.table,y)

x[] <- format(round(100*x,2), nsmall = 2)
dim(x)
x[which(is.na(as.numeric(x)))] <- "--"
x[which(as.numeric(x) == 0)] <- "< 0.01"
x <- latex(x, file = here("Ecov_study","mortality","manuscript","convergence_PRD_table.tex"), 
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
		labs
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

df$conv_txt <- as.factor(c("Success","Fail")[match(df$conv,1:0)])


full.trees <- list()
for(type in levels(df$OM_PE)){

	temp_df <- subset(df, OM_PE == type)
  form <- as.formula(paste("conv_txt ~", paste(factors, collapse = "+")))
	full.trees[[type]] <- rpart(form, data=temp_df, method = "class", control=rpart.control(cp=0, xval = 100))
}
anova.palette.sd <- 0.25

cairo_pdf(here("Ecov_study","mortality","manuscript", "convergence_classification_plots.pdf"), width = 30*2/3, height = 20*2/3)
x <- matrix(1:3, 1, 3, byrow = TRUE)
layout.x <- layout(x) 
par(oma = c(0,0,0,0))
#plot.prune(full.trees[["R"]], cp = 0.001, "R", tweak = 1.2, mar = c(0,0,5,0), split.yshift = -3)
plot.prune(prune(full.trees[["R"]],c(12, 56)), cp = 0.001, "R", tweak = 1.6, mar = c(0,0,5,0), split.yshift = -3)
mtext("R OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["R+S"]],24),0.01, "R+S", tweak = 1.6, mar = c(0,0,5,0), split.yshift = -3)
mtext("R+S OMs", side = 3, line = 0, cex = 2)
plot.prune(prune(full.trees[["R+M"]],c(10,22,24)),cp = 0.001, "R+M", tweak= 1.6, mar = c(0,0,5,0), split.yshift = -3)
mtext("R+M OMs", side = 3, line = 0, cex = 2)
dev.off()

# "\U1D440"#italic M
# "\U1D5B3"#subscript S
# "\U2091" #subscript e
# "\U2097" #subscript l
# "\U1D70E\U1D5B3\U2091\U2097" #sigma_Sel
# "\U1D70C" #rho
# "\U1D70E\U2082\U208A" #sigma_2+
# "\U2192" #right arrow

