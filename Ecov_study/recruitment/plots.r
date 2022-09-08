library(here)

#################################################
##--LOAD OMs, EMs, FITS--########################
#################################################
dir <- file.path(here(),"Ecov_study","recruitment","results")
  
sim_input <- readRDS(file.path(dir, "om_sim_data_GLB_recruitment.RDS"))
em_input  <- readRDS(file.path(dir, "em_input_GLB_recruitment.RDS"))
em_fits   <- readRDS(file.path(dir, "em_fits_GLB_recruitment.RDS"))
df.mods   <- readRDS(file.path(dir, "om_sim_inputs_GLB_recruitment.RDS"))
n.mods    <- df.mods$nsim
nsim      <- n.mods[1]

yrs <- em_fits[[1]][[1]]$years

##--check which models did not converge--#########################
lapply(1:n.mods, function(y) unlist(lapply(1:nsim, function(x) check_convergence(em_fits[[y]][[x]],ret=TRUE)$convergence)))
lapply(1:n.mods, function(y) unlist(lapply(1:nsim, function(x) check_convergence(em_fits[[y]][[x]],ret=TRUE)$is_sdrep)))

#############################################################################
## ECOV ANALYSIS ############################################################
#############################################################################

##--EXAMPLE PLOTS OF TRUE VS SIMULATED ECOV REs--############################
j <- 10
i <- 1

plot(-999,xlim=range(yrs), ylim=c(-1,1))
for(i in 1:2){
  sim <- sim_input[[j]][[i]]$data$Ecov_re
  est <- em_fits[[j]][[i]]$parList$Ecov_re
  
  lines(yrs,  sim, lty=i)
  lines(yrs,est,col='red', lty=i)
}

##_-PLOT SIMULATED AND ESTIMATION ECOV REs--##############################
pdf('ecov_re_diff.pdf',width=9)
par(mfrow=c(2,2))
for(j in c(1,8,16,20)){
  plot(-999,xlim=range(yrs), ylim=c(-1,1))
  for(i in 1:4){
    sim <- sim_input[[j]][[i]]$data$Ecov_re
    est <- em_fits[[j]][[i]]$parList$Ecov_re
    
    lines(yrs,  sim)
    lines(yrs,est,col='red')
  }
}
dev.off()


#############################################################################
## BETA ANALYSIS ############################################################
#############################################################################

##--PLOT BETAS WITH TRUE VALUE ACROSS OBS_SIG, ECOV_SIG, BETA--#############################
betas    <- unique(df.mods$beta)
ecov_sigs <- unique(df.mods$Ecov_sig)
obs.sigs  <- unique(df.mods$obs_sig)
l <- which(row.names(em_fits[[1]][[1]]$sdrep)=="Ecov_beta")

design <- expand.grid(obs.sigs,ecov_sigs,betas)

pdf('beta_ecov_hists_obs_sig_ecov_sig_beta.pdf',height=10,width=10)
par(mfrow=c(4,3),mar=c(2,3,2,2),oma=c(2,2,2,2))
for(i in 1:nrow(design)){
    whic <- which(df.mods$Ecov_sig==design[i,2] & df.mods$obs_sig==design[i,1] & df.mods$beta==design[i,3])
    hist(unlist(lapply(whic, function(y) lapply(1:nsim, function(x) em_fits[[y]][[x]]$sdrep[l,1]))),
           main=paste0("obs_sig=",obs.sigs[j],"  Ecov_sig=",ecov_sigs[i]), xlab='', xlim=c(-10,10),breaks=10000)
    abline(v=betas[i],lwd=2,col='red')
}
mtext(outer=TRUE,'beta',side=1)
dev.off()

##--PLOT BETAS WITH TRUE VALUE ACROSS OBS_SIG, ECOV_SIG, COLLAPSE BETA--#############################
ecov_sigs <- unique(df.mods$Ecov_sig)
obs.sigs  <- unique(df.mods$obs_sig)
l <- which(row.names(em_fits[[1]][[1]]$sdrep)=="Ecov_beta")
design <- expand.grid(obs.sigs,ecov_sigs)

pdf('beta_ecov_hists_obs_sig_ecov_sig_across_beta.pdf',height=7,width=7)
par(mfrow=c(2,2),mar=c(2,3,2,2),oma=c(2,2,2,2))
for(i in 1:nrow(design)){
    whic <- which(df.mods$Ecov_sig==design[i,2] & df.mods$obs_sig==design[i,1])
    hist(unlist(lapply(whic, function(y) lapply(1:nsim, function(x) em_fits[[y]][[x]]$sdrep[l,1]))),
           main=paste0("obs_sig=",obs.sigs[j],"  Ecov_sig=",ecov_sigs[i]), xlab='', xlim=c(-20,20),breaks=10000)
}
mtext(outer=TRUE,'beta',side=1)
dev.off()

##--INDIVIDUAL BETA ESTIMATES--#######
j <- which(row.names(em_fits[[1]][[1]]$sdrep)=="Ecov_beta")

pdf('beta_ecov_hists_individual.pdf',height=9,width=6)
par(mfrow=c(5,2),mar=c(2,2,2,2),oma=c(2,2,2,2))
for(p in 1:n.mods){
  #p <- 23
  beta_true    <- df.mods$beta[p]
  obs_sig_true <- df.mods$obs_sig[p] 
  for(i in 1:nsim){
      plot(-999, xlim=c(-5,5), ylim=c(0,1),yaxt='n') 
        beta_hat <- em_fits[[p]][[i]]$sdrep[j,1] 
        beta_se  <- em_fits[[p]][[i]]$sdrep[j,2]
        abline(v=beta_hat)
        abline(v=beta_hat + 2*c(-1,1)*em_fits[[p]][[i]]$sdrep[j,2],lty=2)
        abline(v=beta_true,col='red')
        mtext(outer=TRUE,'beta',side=1)
        #mtext(outer=TRUE,paste0('beta=',beta_true,"  obs.sig=",obs_sig_true),side=3)
        mtext(outer=TRUE,paste0('ecov_sig=',df.mods[p,4],'  ecov_phi=',df.mods[p,5],'  beta=',df.mods[p,7],'  obs_sig=',df.mods[p,9]),side=3)
  }
}
dev.off()

###############################################################################
## RECRUITMENT ################################################################
###############################################################################
betas    <- unique(df.mods$beta)
ecov_sigs <- unique(df.mods$Ecov_sig)
obs.sigs  <- unique(df.mods$obs_sig)
l <- which(row.names(em_fits[[1]][[1]]$sdrep)=="Ecov_beta")
design <- expand.grid(obs.sigs,ecov_sigs,betas)

par(mfrow=c(3,4),mar=c(2,3,2,2),oma=c(2,2,2,2))
for(i in 1:nrow(design)){
  whic <- which(df.mods$Ecov_sig==design[i,2] & df.mods$obs_sig==design[i,1] & df.mods$beta==design[i,3])
  hist(unlist(lapply(whic, function(y) lapply(1:nsim, function(x) cor(exp(em_fits[[y]][[x]]$parList$log_NAA[,1]),
                                                                      exp(em_fits[[y]][[x]]$env$data$log_NAA[,1]))))),
       main=paste0("obs_sig=",obs.sigs[j],"  Ecov_sig=",ecov_sigs[i]), xlab='', xlim=c(0,1),breaks=15)
}
mtext(outer=TRUE,'Correlation Coefficient',side=1,line=0.5)


ecov_sigs <- unique(df.mods$Ecov_sig)
obs.sigs  <- unique(df.mods$obs_sig)
l <- which(row.names(em_fits[[1]][[1]]$sdrep)=="Ecov_beta")
design <- expand.grid(obs.sigs,ecov_sigs)

par(mfrow=c(2,2),mar=c(2,3,2,2),oma=c(2,2,2,2))
for(i in 1:nrow(design)){
  whic <- which(df.mods$Ecov_sig==design[i,2] & df.mods$obs_sig==design[i,1])
  hist(unlist(lapply(whic, function(y) lapply(1:nsim, function(x) cor(exp(em_fits[[y]][[x]]$parList$log_NAA[,1]),
                                                                      exp(em_fits[[y]][[x]]$env$data$log_NAA[,1]))))),
       main=paste0("obs_sig=",obs.sigs[j],"  Ecov_sig=",ecov_sigs[i]), xlab='', xlim=c(-1,1),breaks=20)
}
mtext(outer=TRUE,'Correlation Coefficient',side=1,line=0.5)




###############################################################################
## SSB ########################################################################
###############################################################################
betas    <- unique(df.mods$beta)
ecov_sigs <- unique(df.mods$Ecov_sig)
obs.sigs  <- unique(df.mods$obs_sig)
design <- expand.grid(obs.sigs,ecov_sigs,betas)

par(mfrow=c(3,4),mar=c(2,3,2,2),oma=c(2,2,2,2))
for(i in 1:nrow(design)){
  whic <- which(df.mods$Ecov_sig==design[i,2] & df.mods$obs_sig==design[i,1] & df.mods$beta==design[i,3])
  hist(unlist(lapply(whic, function(y) lapply(1:nsim, function(x){ 
    Waa <- em_fits[[y]][[x]]$env$data$waa[1,2:40,] * exp(em_fits[[y]][[x]]$parList$log_NAA) * em_fits[[y]][[x]]$env$data$mature[2:40,]
    return(cor(rowSums(Waa), sim_input[[y]][[x]]$data$SSB[-1]))
  } ))),
       xlab='', xlim=c(0.9,1),breaks=10, cex=0.1,main='')
  mtext(paste0("obs_sig=",design[i,1],"  Ecov_sig=",design[i,2],"  beta=",design[i,3]),cex=0.4)
}
mtext(outer=TRUE,'Correlation Coefficient',side=1,line=0.5)


ecov_sigs <- unique(df.mods$Ecov_sig)
obs.sigs  <- unique(df.mods$obs_sig)
design <- expand.grid(obs.sigs,ecov_sigs)

par(mfrow=c(2,2),mar=c(2,3,2,2),oma=c(2,2,2,2))
for(i in 1:nrow(design)){
  whic <- which(df.mods$Ecov_sig==design[i,2] & df.mods$obs_sig==design[i,1])
  hist(unlist(lapply(whic, function(y) lapply(1:nsim, function(x){ 
    Waa <- em_fits[[y]][[x]]$env$data$waa[1,2:40,] * exp(em_fits[[y]][[x]]$parList$log_NAA) * em_fits[[y]][[x]]$env$data$mature[2:40,]
    return(cor(rowSums(Waa), sim_input[[y]][[x]]$data$SSB[-1]))
  } ))),
  xlab='', xlim=c(0.9,1),breaks=10, cex=0.1,main='')
  mtext(paste0("obs_sig=",design[i,1],"  Ecov_sig=",design[i,2],"  beta=",design[i,3]),cex=0.4)
}
mtext(outer=TRUE,'Correlation Coefficient',side=1,line=0.5)







