library(rpart)
library(rpart.plot)
library(wham)
library(randomForest)
##--Read data--###############

dir <- "~/dropbox/working/state_space_assessments/cluster/01_25_2022/"
#dir <- "~/dropbox/working/state_space_assessments/SSRTWG/Ecov_study/recruitment/results/"
em_fits <- readRDS(file.path(dir,'em_fits_GLB_recruitment_doparallel.RDS'))
df.mods <- readRDS(file.path(dir,'om_sim_inputs_GLB_recruitment_doparallel.RDS'))
n.mods <- nrow(df.mods)
nsim <- length(em_fits[[1]])

#check convergence
#model 13 giving NaNs for unknown reason
conv <- unlist(lapply(1:n.mods, function(y) unlist(lapply(1:nsim, function(x) em_fits[[1]][[1]]$opt$convergence))))
sum(conv)

fin <- 40 #number of years to extract last year

#names of input variables
nms <- c('dssb','d2ssb','dr','d2r','dF','d2F')#,'dssb_f','d2ssb_f','dr_f','d2r_f','dF_f','d2F_f')

temp <- list()
for(m in 1:n.mods){
    temp[[m]] <- lapply(1:nsim, function(x){
    ##SSB
    ssb_sim <- em_fits[[m]][[x]]$input$data$SSB
    ssb_est <- em_fits[[m]][[x]]$rep$SSB
    
    r_sim <- em_fits[[m]][[x]]$input$data$NAA[,1]
    r_est <- em_fits[[m]][[x]]$rep$NAA[,1]

    F_sim <- em_fits[[m]][[x]]$input$data$F
    F_est <- em_fits[[m]][[x]]$rep$F
    
    res <- list(mean(ssb_est - ssb_sim),
                sqrt(mean((ssb_est - ssb_sim)^2)),
                mean(r_est - r_sim),
                sqrt(mean((r_est - r_sim)^2)),
                mean(F_est - F_sim),
                sqrt(mean((F_est - F_sim)^2))
                #ssb_est[fin] - ssb_sim[fin],
                #sqrt((ssb_est[fin] - ssb_sim[fin])^2),
                #r_est[fin] - r_sim[fin],
                #sqrt((r_est[fin] - r_sim[fin])^2),
                #F_est[fin] - F_sim[fin],
                #sqrt((F_est[fin] - F_sim[fin])^2)
                ); names(res) <- nms
    return(res)
  })
}
results <- unlist(temp)

##--Combine into dataframe--#############
df <- df.mods[rep(1:nrow(df.mods),each=nsim),]
for(i in 1:length(nms)){
  df <- cbind(df,as.numeric(results[names(results)==nms[i]]))
}         
colnames(df)[10:15] <- nms


###############################################
## PLOT RAW DATA ##############################
###############################################
par(mfrow=c(2,3),mar=c(2,2,2,2),oma=c(2,2,2,2))
for(i in 1:length(nms)){
  x <- unlist(select(df,nms[i]) %>% filter(!is.na(.)))
  ii <- which(x < mean(x) + 3*sd(x) & x>mean(x)-3*sd(x))
  hist(x[ii],main=nms[i])
  abline(v=mean(x[ii]),lwd=2,col='red')
}


################################################
## MODEL BUILDING ##############################
################################################
##--Regression trees--###########
fit <- rpart(d2r ~ Ecov_sig + Ecov_phi + beta + obs_sig + NAA_sig, data=df)


##--Linear regression--#########
fit <- lm(dr ~ (Ecov_sig + Ecov_phi + beta + obs_sig + NAA_sig)^2, data=df)
summary(fit)

fit <- lm(d2r ~ (Ecov_sig + Ecov_phi + beta + obs_sig + NAA_sig)^2, data=df)
summary(fit)


##--Random forest--###########
fit <- randomForest(d2r ~ Ecov_sig + Ecov_phi + beta + obs_sig + NAA_sig, data=df, na.action=na.omit)
summary(fit)

cor(df$d2r[!is.na(df$d2r)],predict(fit))^2


#visualization
plot_min_depth_distribution(fit)

#ale plot
X <- cbind(Ecov_phi=df$Ecov_phi,
           Ecov_sig=df$Ecov_sig,
           beta=df$beta,
           obs_sig=df$obs_sig,
           NAA_sig=df$NAA_sig)

yhat <- function(X.model,newdata) as.numeric(predict(X.model,newdata))
ALEPlot(X,fit,pred.fun=yhat,J=2)

#partial dependence plot
par(mar=c(3,3,3,3))
PDPlot(X,fit,pred.fun=yhat,J=c(2,3),K=10)








