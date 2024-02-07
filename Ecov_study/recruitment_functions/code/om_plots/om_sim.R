library(here)
library(wham)

om_inputs <- readRDS(file.path(here(),"Ecov_study", "recruitment_functions", "inputs", "om_inputs.RDS"))
df.oms    <- readRDS(file.path(here(),"Ecov_study", "recruitment_functions", "inputs", "df.oms.RDS"))

om_sim <- function(om_input, nsim){
  om  <- fit_wham(om_input, do.fit = FALSE, MakeADFun.silent = TRUE)
  lapply(1:nsim, function(x) return(om$simulate(complete=TRUE)))
}

##--SR functions no ecov--###############################
is <- c(
  which(df.oms$Fhist=="H-MSY" & df.oms$R_sig==0.1 & df.oms$Ecov_how==0)[2],
  which(df.oms$Fhist=="MSY"   & df.oms$R_sig==0.1 & df.oms$Ecov_how==0)[2],
  which(df.oms$Fhist=="H-MSY" & df.oms$R_sig==1.0 & df.oms$Ecov_how==0)[2],
  which(df.oms$Fhist=="MSY"   & df.oms$R_sig==1.0 & df.oms$Ecov_how==0)[2])

labs <- c(expression("F = H-MSY, "*sigma['R']*' = 0.1'),
          expression("F = MSY, "*sigma['R']*' = 0.1'),
          expression("F = H-MSY, "*sigma['R']*' = 1.0'),
          expression("F = MSY, "*sigma['R']*' = 1.0'))

sims <- list()
for(i in 1:length(is)){
  sims[[i]] <-  om_sim(om_inputs[[is[i]]],nsim=n_sim)
}


pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','stock_recruit_no_ecov.pdf'),width=7,height=5)
par(mfrow=c(2,2),mar=c(2,2,2,2),oma=c(2,2,2,2),cex.axis=0.8)
for(i in 1:length(is)){
  R   <- sims[[i]][[1]]$NAA[,1]
  SSB <- sims[[i]][[1]]$SSB  
  if(i==1|i==2) plot(SSB,R,xlim=c(0,2.5E5),ylim=c(0,2.5E4),pch=20)
  if(i==3|i==4) plot(SSB,R,xlim=c(0,2.5E5),ylim=c(0,2.5E5),pch=20)
  mtext(labs[i])
}
mtext(outer=TRUE,side=1,'Spawning Stock Biomass',line=0.5)
mtext(outer=TRUE,side=2,'Recruitment',line=0.5)
dev.off()


##--SR functions with ecov--###############################
is <- c(
  which(df.oms$R_sig==0.1 & df.oms$Ecov_how==1 & df.oms$Ecov_re_cor==0.2 & df.oms$Ecov_effect==1.0)[1],
  which(df.oms$R_sig==0.1 & df.oms$Ecov_how==1 & df.oms$Ecov_re_cor==0.8 & df.oms$Ecov_effect==1.0)[1])
  
labs <- c(expression(sigma['R']*' = 0.1, '*rho["E"]*' = 0.2'),
          expression(sigma['R']*' = 0.1, '*rho["E"]*' = 0.8'))

sims <- list()
for(i in 1:length(is)){
  sims[[i]] <-  om_sim(om_inputs[[is[i]]],nsim=n_sim)
}


pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','recruit_ecov.pdf'),width=9,height=4)
par(mfrow=c(1,2),mar=c(2,2,2,0),oma=c(2,2,2,4),cex.axis=0.8)
for(i in 1:length(is)){
  R   <- sims[[i]][[1]]$NAA[,1]
  E   <- sims[[i]][[1]]$Ecov_x
  SSB <- sims[[i]][[1]]$SSB  
  years <- sims[[i]][[1]]$Ecov_year
  #plot(SSB,R,xlim=c(0,5E5),ylim=c(0,2E5),pch=20)
  plot(years,R,ylim=c(0,5E4),pch=20,xlim=c(1980,2022),yaxt='n')
  if(i==1) axis(side=2)
  par(new=TRUE)
  plot(years,E,yaxt='n',xaxt='n',type='l',col='red',bty='n',ylim=c(-0.5,0.5),xlim=c(1980,2022))
  mtext(labs[i])
  if(i==2) axis(side=4,col='red')
}
mtext(outer=TRUE,side=4,'Ecov',line=2.5,col='red')
mtext(outer=TRUE,side=2,'Recruitment',line=0.5)
dev.off()



