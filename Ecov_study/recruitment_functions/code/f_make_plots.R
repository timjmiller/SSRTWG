make_plots <- function(lls,om_inputs,ncols,nrows,ylims,xlims){
  
  pdf(file.path(here(),'Ecov_study','recruitment_functions','plots','S_R_om_sims.pdf'),height=8,width=8)
  par(mfrow=c(ncols,nrows),mar=rep(0,4),oma=c(5,5,2,2))
  #for(ll in 1:length(om_inputs)){
  for(ll in lls){
    om  <- fit_wham(om_inputs[[ll]], do.fit = FALSE, MakeADFun.silent = TRUE)
    lapply(1:(ncol*nrow),function(p) {
      sim <- om$simulate(complete=TRUE)
      ssb_sim <- sim$SSB
      r_sim   <- sim$NAA[,1]
      
      plot(ssb0,r0,type='l',
           xlim=xlims,
           ylim=ylims,
           #xlim=c(0,max(ssb0)),
           #ylim=c(0,max(c(r0,r_sim)))
           xaxt='n',yaxt='n')
      points(sim$SSB,sim$NAA[,1],pch=19,cex=0.8)
      if(p%in%seq(1,ncol*nrow,nrow)) axis(side=2,las=2)
      if(p%in%seq(nrow*ncol - ncol + 1, nrow*ncol)) axis(side=1)
    })
    nms <- lapply(2:10,function(j) paste0(colnames(df.oms[1,])[j],' = ',df.oms[ll,j]))
    mtext(paste(nms,collapse='  '),outer=TRUE,line=0,cex=0.6)
    mtext('SSB',outer=TRUE,side=1,line=3)
    mtext('Recruits',outer=TRUE,side=2,line=3.5)
  }
  dev.off()
  
  ##--PLOT ECOV EFFECTS--###################################################
  pdf(file.path(here(),'Ecov_study','recruitment_functions','plots','S_R_om_sims_residuals.pdf'),height=8,width=8)
  par(mfrow=c(ncols,nrows),mar=rep(0,4),oma=c(5,5,2,5))
  #for(ll in 1:length(om_inputs)){
  for(ll in lls){
    om  <- fit_wham(om_inputs[[ll]], do.fit = FALSE, MakeADFun.silent = TRUE)
    lapply(1:(ncol*nrow),function(p) {
      sim <- om$simulate(complete=TRUE)
      ssb_sim <- sim$SSB
      r_sim   <- sim$NAA[,1]
      ecov    <- sim$Ecov_obs
      r_res   <- r_sim/(a*ssb_sim)/(1+b*ssb_sim)
      
      plot(exp(-ecov),r_sim,xaxt='n',yaxt='n',pch=19,cex=0.8,log='y')
      if(p%in%seq(1,ncol*nrow,nrow)) axis(side=2,las=2,cex.axis=0.5)
      par(new=TRUE)
      plot(exp(-ecov),r_res,xaxt='n',yaxt='n',col='red',cex=0.8,log='y')
      if(p%in%seq(ncol,nrow*ncol,nrow)) axis(side=4,las=2,col='red',cex.axis=0.5)
      
      #cor(ecov,r_sim)
      #cor(ecov,r_res)
      
      
      if(p%in%seq(nrow*ncol - ncol + 1, nrow*ncol)) axis(side=1)
    })
    nms <- lapply(2:10,function(j) paste0(colnames(df.oms[1,])[j],' = ',df.oms[ll,j]))
    mtext(paste(nms,collapse='  '),outer=TRUE,line=0,cex=0.6)
    mtext(side=1,'Ecov',outer=TRUE,line=2.5)
    mtext(side=2,'Recruits',outer=TRUE,line=2.5)
    mtext(side=4,'Recruitmena Residual',outer=TRUE,line=2.5)
  }
  
  dev.off()
  
  ##--PLOT ECOV--#############################################
  pdf(file.path(here(),'Ecov_study','recruitment_functions','plots','S_R_om_sims_ecov.pdf'),height=8,width=8)
  par(mfrow=c(ncols,nrows),mar=rep(0,4),oma=c(5,5,3,3))
  #for(ll in 1:length(om_inputs)){
  for(ll in lls){
    om  <- fit_wham(om_inputs[[ll]], do.fit = FALSE, MakeADFun.silent = TRUE)
    lapply(1:(ncol*nrow),function(p) {
      sim  <- om$simulate(complete=TRUE)
      ecov <- sim$Ecov_obs
      year <- sim$Ecov_year
      
      plot(year,exp(-ecov),xaxt='n',yaxt='n',pch=19,cex=0.8,log='y')
      if(p%in%seq(1,ncol*nrow,nrow)) axis(side=2,las=2,cex.axis=0.5)
      #par(new=TRUE)
      #plot(ecov,r_res,xaxt='n',yaxt='n',col='red',cex=0.8,log='y')
      #if(p%in%seq(ncol,nrow*ncol,nrow)) axis(side=4,las=2,col='red',cex.axis=0.5)
      
      #cor(ecov,r_sim)
      #cor(ecov,r_res)
      
      
      if(p%in%seq(nrow*ncol - ncol + 1, nrow*ncol)) axis(side=1)
    })
    nms <- lapply(2:10,function(j) paste0(colnames(df.oms[1,])[j],' = ',df.oms[ll,j]))
    mtext(paste(nms,collapse='  '),outer=TRUE,line=0,cex=0.6)
    mtext(side=2,'Ecov',outer=TRUE,line=2.5)
  }
  dev.off()
}
