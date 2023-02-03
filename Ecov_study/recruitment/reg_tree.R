library(rpart)
library(rpart.plot)
library(wham)
library(randomForest)
##--Read data--###############

dir <- "~/dropbox/working/state_space_assessments/cluster/01_31_2023/"
#dir <- "~/dropbox/working/state_space_assessments/SSRTWG/Ecov_study/recruitment/results/"
df.mods <- readRDS(file.path(dir,'om_sim_inputs_GLB_recruitment_doparallel.RDS'))
n.mods <- nrow(df.mods)
nsim <- 100

#em_fits <- list()
#for(i in 1:n.mods){
#print(i)
#  em_fits[[i]] <- readRDS(paste0(dir,"em_fits_m_",i))  
#}

#em_fits <- readRDS(file.path(dir,'em_fits_GLB_recruitment_doparallel.RDS'))
#nsim <- length(em_fits[[1]])

#check convergence
#model 13 giving NaNs for unknown reason
#conv <- unlist(lapply(1:n.mods, function(y) unlist(lapply(1:nsim, function(x) em_fits[[1]][[1]]$opt$convergence))))
#sum(conv)

#fin <- 40 #number of years to extract last year

#names of input variables
nms <- c('dssb','d2ssb','rssb','dr','d2r','rr','dF','d2F','rF')#,'dssb_f','d2ssb_f','dr_f','d2r_f','dF_f','d2F_f')

temp <- list()
for(m in 1:n.mods){
#for(m in c(1:14,16:27)){
print(m)
    em_fits <- readRDS(paste0(dir,"em_fits_m_",m))
    
    temp[[m]] <- lapply(1:nsim, function(x){

    ssb_sim <- em_fits[[x]]$input$data$SSB
    ssb_est <- em_fits[[x]]$rep$SSB
    
    r_sim <- em_fits[[x]]$input$data$NAA[,1]
    r_est <- em_fits[[x]]$rep$NAA[,1]

    F_sim <- em_fits[[x]]$input$data$F
    F_est <- em_fits[[x]]$rep$F
    
    res <- list(mean(ssb_est - ssb_sim),
                sqrt(mean((ssb_est - ssb_sim)^2)),
                mean(ssb_est/ssb_sim),
                mean(r_est - r_sim),
                sqrt(mean((r_est - r_sim)^2)),
                mean(r_est/r_sim),
                mean(F_est - F_sim),
                sqrt(mean((F_est - F_sim)^2)),
                mean(F_est/F_sim)
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
#df <- df.mods[rep(c(1:14,16:27),each=nsim),]

for(i in 1:length(nms)){
  df <- cbind(df,as.numeric(results[names(results)==nms[i]]))
}         
colnames(df)[10:ncol(df)] <- nms

###############################################
## PLOT RAW DATA ##############################
###############################################
pdf('~/Dropbox/Working/state_space_assessments/SSRTWG/Ecov_study/recruitment/results/hists.pdf',
    height=5,width=7)
par(mfrow=c(3,3),mar=c(2,2,2,2),oma=c(2,2,2,2))
for(i in 1:length(nms)){
  x <- unlist(select(df,nms[i]) %>% filter(!is.na(.)))
  ii <- which(x < mean(x) + 3*sd(x) & x>mean(x)-3*sd(x))
  hist(x[ii],main='');mtext(side=1,line=2.5,nms[i],cex=0.8)
  abline(v=mean(x[ii]),lwd=2,col='red')
}
dev.off()

################################################
## MODEL BUILDING ##############################
################################################
##--Regression trees--###########
cps <- c(0.01,0.005,0.001)
pdf('~/Dropbox/Working/state_space_assessments/SSRTWG/Ecov_study/recruitment/results/trees2.pdf',
    height=5,width=7)
par(mfrow=c(3,3),oma=c(2,2,2,2))
for(j in 1:3){
for(i in 1:length(nms)){
  form <- paste0(nms[i],' ~ Ecov_sig + Ecov_phi + beta + obs_sig + NAA_sig')
  fit <- rpart(form, data=df, control=rpart.control(cp=cps[j]))
  rpart.plot(fit); mtext(side=1,line=5,nms[i],cex=0.6)
}
  mtext(outer=TRUE,paste0('cp = ',cps[j]*10))
}
dev.off()


##--Linear regression--#########
plotlm <- function(fit,i,nm){
  nn <- nrow(summary(fit)$coefficients)
  coefs <- summary(fit)$coefficients[2:nn,1]
  ses   <- summary(fit)$coefficients[2:nn,2]
  maxx <- max(abs(coefs+2*ses))
  plot(1:(nn-1),coefs,ylim=c(-maxx,maxx),bty='n',xaxt='n',pch=19)
  if(i%in%7:9){axis(side=1,at=1:(nn-1),labels=names(coefs),las=2)}else{
    axis(side=1,at=1:(nn-1),labels=NA)}
  segments(x0=1:(nn-1),x1=1:(nn-1),y0=coefs-2*ses,y1=coefs+2*ses)
  abline(h=0)
  mtext(nm)
}

df2 <- df
df2[,4:8] <- scale(df2[,4:8])

pdf('~/Dropbox/Working/state_space_assessments/SSRTWG/Ecov_study/recruitment/results/coefs_noint.pdf',
    height=5,width=7)
par(mfrow=c(3,3),oma=c(7,2,2,2),mar=c(2,2,2,2),cex.lab=0.7)
for(i in 1:length(nms)){
  #form <- paste0(nms[i],' ~ (Ecov_sig + Ecov_phi + beta + obs_sig + NAA_sig)^2')
  form <- paste0(nms[i],' ~ Ecov_sig + Ecov_phi + beta + obs_sig + NAA_sig')
  fit <- lm(form, data=df2)
  plotlm(fit,i,nm=nms[i])
}
dev.off()


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








