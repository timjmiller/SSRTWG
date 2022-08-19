

yrs <- em_fits[[15]][[1]]$years

##_-PLOT SIMULATED AND ESTIMATION ECOV REs--##############################
pdf('ecov_re_diff.pdf',width=9)
par(mfrow=c(2,2))
for(i in 1:4){
  sim <- sim_input[[15]][[i]]$data$Ecov_re
  est <- em_fits[[15]][[i]]$parList$Ecov_re
  
  plot(yrs,  sim,xlim=c(1980,2021),ylim=range(c(sim,est)),ylab='')
  points(yrs,est,col='red')
}
dev.off()

##--PLOT ESTIMATED BETAS WITH SIMULATED VALUE--_####################
pdf('beta_ecov_hists.pdf',height=4)
ms <- c(1,15)
par(mfrow=c(1,2)) 
for(i in 1:length(ms)){
    hist(unlist(lapply(1:nsim, function(x){
    betas <- em_fits[[ms[i]]][[x]]$parList$Ecov_beta[1,1,1,1]
    return(betas)
})), main='',xlab='')
abline(v=df.mods[ms[i],]$beta,lwd=2)
}
dev.off()

