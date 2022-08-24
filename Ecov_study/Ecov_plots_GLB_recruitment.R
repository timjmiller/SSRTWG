
#################################################
##--LOAD OMs, EMs, FITS--########################
#################################################
dir <- "~/dropbox/Working/state_space_assessments/SSRTWG/Ecov_study/results"
  
sim_input <- readRDS(file.path(dir, "om_sim_data_GLB_recruitment.RDS"))
em_input  <- readRDS(file.path(dir, "em_input_GLB_recruitment.RDS"))
em_fits   <- readRDS(file.path(dir, "em_fits_GLB_recruitment.RDS"))
df.mods   <- readRDS(file.path(dir, "om_sim_inputs_GLB_recruitment.RDS"))

yrs <- em_fits[[1]][[1]]$years


##--EXAMPLE PLOTS--############################
j <- 10
i <- 1

sim <- sim_input[[j]][[i]]$data$Ecov_re
est <- em_fits[[j]][[i]]$parList$Ecov_re

plot(-999,xlim=range(yrs), ylim=c(-1,1))
for(i in 1:2){
  sim <- sim_input[[j]][[i]]$data$Ecov_re
  est <- em_fits[[j]][[i]]$parList$Ecov_re
  
  lines(yrs,  sim, lty=i)
  lines(yrs,est,col='red', lty=i)
}

lines(yrs,  sim)
lines(yrs,est,col='red')



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



##--check which models did not converge--#########################
lapply(1:n.mods, function(y) unlist(lapply(1:nsim, function(x) check_convergence(em_fits[[y]][[x]],ret=TRUE)$convergence)))


##--PLOT ESTIMATED BETAS WITH SIMULATED VALUE--_####################
betas    <- unique(df.mods$beta)
obs.sigs <- unique(df.mods$obs_sig)
p <- which(row.names(em_fits[[1]][[1]]$sdrep)=="Ecov_beta")

pdf('beta_ecov_hists.pdf',height=4)
par(mfrow=c(3,2))
for(i in 1:length(betas)){
  for(j in 1:length(obs.sigs)){
    whic <- which(df.mods$beta == betas[i] & df.mods$obs_sig == obs.sigs[j])
    #hist(unlist(lapply(whic, function(y) lapply(1:nsim, function(x) em_fits[[y]][[x]]$parList$Ecov_beta[1,1,1,1]))),
    hist(unlist(lapply(whic, function(y) lapply(1:nsim, function(x) em_fits[[y]][[x]]$sdrep[p,1]))),
              main=paste0("obs_sig=",obs.sigs[j],"  beta=",beta[i]))
    abline(v=betas[i])
  }
}
dev.off()


betas    <- unique(df.mods$beta)
ecov_sigs <- unique(df.mods$Ecov_sig)
obs.sigs  <- unique(df.mods$obs_sig)
p <- which(row.names(em_fits[[1]][[1]]$sdrep)=="Ecov_beta")
par(mfrow=c(2,2))
for(i in 1:length(ecov_sigs)){
  for(j in 1:length(obs.sigs)){
    whic <- which(df.mods$Ecov_sig == ecov_sigs[i] & df.mods$obs_sig == obs.sigs[j])
    #hist(unlist(lapply(whic, function(y) lapply(1:nsim, function(x) em_fits[[y]][[x]]$parList$Ecov_beta[1,1,1,1]))),
    hist(unlist(lapply(whic, function(y) lapply(1:nsim, function(x) em_fits[[y]][[x]]$sdrep[p,1]))),
         main=paste0("obs_sig=",obs.sigs[j],"  Ecov_sig=",ecov_sigs[i]), xlab='')
    abline(v=betas[i],lwd=2)
  }
}


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








