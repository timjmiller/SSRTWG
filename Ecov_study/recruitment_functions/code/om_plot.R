library(here)
library(wham) #make sure to use the right version of wham
source(file.path(here(), "Ecov_study", "mortality", "code", "sim_management.R"))
verify_version()

file_loc  <- 'Ecov_study/recruitment_functions/om_inputs_06_13_2023'
om_inputs <- readRDS(file.path(here(),file_loc,'om_inputs.rds'))
df.oms    <- readRDS(file.path(here(),file_loc,'df.oms.rds'))

#####################################################################
##--Recruitment parameter set in simulation--#########
#GLB: check this
a <- 5.955694e-01 
b <- 2.404283e-05

ncol <- 7
nrow <- 7

ylims <- c(0,1E5)
xlims <- c(0,5E5)

ssb0 <- seq(0,xlims[2],length.out=1000)
r0   <- (a*ssb0)/(1+b*ssb0)

#Subset stocks
lls <- which(df.oms$Ecov_how%in%c(1,2,4) & 
               df.oms$Fhist=='H-MSY' & 
               df.oms$Ecov_obs_sig==min(df.oms$Ecov_obs_sig) &
               df.oms$Ecov_effect==max(df.oms$Ecov_effect) & 
               df.oms$obs_error=='L' & 
               df.oms$R_sig==max(df.oms$R_sig) &
               df.oms$NAA_cor==min(df.oms$NAA_cor) & 
               df.oms$Ecov_re_cor==min(df.oms$Ecov_re_cor))

#PLOT
set.seed(1234)

pdf(file.path(here(),'Ecov_study','recruitment_functions','plots','S_R_om_sims.pdf'),height=10,width=10)
par(mfrow=c(ncol,nrow),mar=rep(0,4),oma=c(2,2,2,2))
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
}
dev.off()

##--PLOT ECOV EFFECTS--###################################################
set.seed(1234)


pdf(file.path(here(),'Ecov_study','recruitment_functions','plots','S_R_om_sims_residuals.pdf'),height=10,width=10)
par(mfrow=c(ncol,nrow),mar=rep(0,4),oma=c(rep(3,4)))
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
}
dev.off()



##--PLOT ECOV--#############################################
par(mfrow=c(ncol,nrow),mar=rep(0,4),oma=c(rep(3,4)))
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
}




