library(here)
library(tidyverse)
library(wham)
library(rpart)
library(rpart.plot)


#dir <- file.path(here::here(), "Ecov_study","recruitment_functions", "results_pile2")
res.path <- 'results'  # directory where simulation runs are (beta unstandardized)
res.dir <- 'results'   # 'results'     'results_beta_fix'   # results folder where AIC dataframes are
plot.dir <- 'plots'    # 'plots_lizruns'  'plots_beta_fix'   
plot.suffix <- '' 

## specify bad.grad.label and bad.se.value (these are the thresholds set in convergence_summaries.R to determine convergence) 
bad.grad.value <- 1E-6 #tim used 1E-6 (the abs of the exponent will be used for output suffix; ex: filename_grad_6.png)
bad.grad.label <- as.numeric(strsplit(as.character(bad.grad.value), split="-")[[1]][2])
bad.se.value <- 100 #tim used 100 (this value will be used for output suffix; ex: filename_se_100.png)

df.oms    <- readRDS(file.path(here::here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS"))
df.ems    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))
conv.runs <- readRDS(file.path(here::here(),'Ecov_study','recruitment_functions',res.dir, paste0("conv.runs_grad_", bad.grad.label, "_SE_", bad.se.value, plot.suffix, ".RDS")  ) )

conv.runs <- conv.runs %>%
  rename(Sim=sim) %>%
  mutate(ok.run =(bad.opt+ bad.conv+ bad.sdrep+ bad.grad+ bad.se.big) ) %>%
  replace_na(list(ok.run=1)) %>%
  relocate(OM, EM, Sim, ok.run, bad.opt, bad.conv, bad.sdrep, bad.grad, bad.se.big) %>%
  filter(ok.run==0)  # drop unconverged runs (there shouldn't be any at this point)

nyears <- 10  # number of projection years

df.recr=df.ssb=df.catch         <- matrix(NA, nrow=nrow(conv.runs)*nyears, ncol=12)
colnames(df.recr)=colnames(df.ssb)=colnames(df.catch) <- c('OM','EM','Sim','Year','re.cont.ecov','re.avg.ecov','re.use.ecov','rmse.cont.ecov','rmse.avg.ecov','rmse.use.ecov','ecov_slope','ssb_cv')

for(irun in 1:nrow(conv.runs)){
  print(irun)
  dat <- try(readRDS(file.path(here::here(),'Ecov_study','recruitment_functions',res.path, paste0("om", conv.runs$OM[irun], '/','sim',conv.runs$Sim[irun],'_','em', conv.runs$EM[irun],'.RDS') ) ) )
  if(class(dat)!='try-error'){

  ssb_cv     <- sd(dat$truth$SSB)/mean(dat$truth$SSB)
  ecov_slope <- summary(lm(dat$truth$Ecov_x ~ seq(1,40)))$coefficients[2,1]

  df.recr[((irun-1)*nyears+1):(irun*nyears),1] <- conv.runs$OM[irun]
  df.recr[((irun-1)*nyears+1):(irun*nyears),2] <- conv.runs$EM[irun]
  df.recr[((irun-1)*nyears+1):(irun*nyears),3] <- conv.runs$Sim[irun]
  df.recr[((irun-1)*nyears+1):(irun*nyears),4] <- 1:10
  df.recr[((irun-1)*nyears+1):(irun*nyears),5] <- (dat$proj$cont.ecov$rep$NAA[31:40,1] - dat$truth$NAA[31:40,1])/dat$truth$NAA[31:40,1]
  df.recr[((irun-1)*nyears+1):(irun*nyears),6] <- (dat$proj$avg.ecov$rep$NAA[31:40,1]  - dat$truth$NAA[31:40,1])/dat$truth$NAA[31:40,1]
  df.recr[((irun-1)*nyears+1):(irun*nyears),7] <- (dat$proj$use.ecov$rep$NAA[31:40,1]  - dat$truth$NAA[31:40,1])/dat$truth$NAA[31:40,1]  
  df.recr[((irun-1)*nyears+1):(irun*nyears),8] <- sqrt((dat$proj$cont.ecov$rep$NAA[31:40,1] - dat$truth$NAA[31:40,1])^2)
  df.recr[((irun-1)*nyears+1):(irun*nyears),9] <- sqrt((dat$proj$avg.ecov$rep$NAA[31:40,1]  - dat$truth$NAA[31:40,1])^2)
  df.recr[((irun-1)*nyears+1):(irun*nyears),10]<- sqrt((dat$proj$use.ecov$rep$NAA[31:40,1]  - dat$truth$NAA[31:40,1])^2)  
  df.recr[((irun-1)*nyears+1):(irun*nyears),11]<- ecov_slope
  df.recr[((irun-1)*nyears+1):(irun*nyears),12]<- ssb_cv  

  df.ssb[((irun-1)*nyears+1):(irun*nyears),1] <- conv.runs$OM[irun]
  df.ssb[((irun-1)*nyears+1):(irun*nyears),2] <- conv.runs$EM[irun]
  df.ssb[((irun-1)*nyears+1):(irun*nyears),3] <- conv.runs$Sim[irun]
  df.ssb[((irun-1)*nyears+1):(irun*nyears),4] <- 1:10
  df.ssb[((irun-1)*nyears+1):(irun*nyears),5] <- (dat$proj$cont.ecov$rep$SSB[31:40] - dat$truth$SSB[31:40])/dat$truth$SSB[31:40]
  df.ssb[((irun-1)*nyears+1):(irun*nyears),6] <- (dat$proj$avg.ecov$rep$SSB[31:40]  - dat$truth$SSB[31:40])/dat$truth$SSB[31:40]
  df.ssb[((irun-1)*nyears+1):(irun*nyears),7] <- (dat$proj$use.ecov$rep$SSB[31:40]  - dat$truth$SSB[31:40])/dat$truth$SSB[31:40]  
  df.ssb[((irun-1)*nyears+1):(irun*nyears),8] <- sqrt((dat$proj$cont.ecov$rep$SSB[31:40] - dat$truth$SSB[31:40])^2)
  df.ssb[((irun-1)*nyears+1):(irun*nyears),9] <- sqrt((dat$proj$avg.ecov$rep$SSB[31:40]  - dat$truth$SSB[31:40])^2)
  df.ssb[((irun-1)*nyears+1):(irun*nyears),10]<- sqrt((dat$proj$use.ecov$rep$SSB[31:40]  - dat$truth$SSB[31:40])^2)  
  df.ssb[((irun-1)*nyears+1):(irun*nyears),11]<- ecov_slope
  df.ssb[((irun-1)*nyears+1):(irun*nyears),12]<- ssb_cv  

  df.catch[((irun-1)*nyears+1):(irun*nyears),1] <- conv.runs$OM[irun]
  df.catch[((irun-1)*nyears+1):(irun*nyears),2] <- conv.runs$EM[irun]
  df.catch[((irun-1)*nyears+1):(irun*nyears),3] <- conv.runs$Sim[irun]
  df.catch[((irun-1)*nyears+1):(irun*nyears),4] <- 1:10
  df.catch[((irun-1)*nyears+1):(irun*nyears),5] <- (dat$proj$cont.ecov$rep$pred_catch[31:40] - dat$truth$pred_catch[31:40])/dat$truth$pred_catch[31:40]
  df.catch[((irun-1)*nyears+1):(irun*nyears),6] <- (dat$proj$avg.ecov$rep$pred_catch[31:40]  - dat$truth$pred_catch[31:40])/dat$truth$pred_catch[31:40]
  df.catch[((irun-1)*nyears+1):(irun*nyears),7] <- (dat$proj$use.ecov$rep$pred_catch[31:40]  - dat$truth$pred_catch[31:40])/dat$truth$pred_catch[31:40]  
  df.catch[((irun-1)*nyears+1):(irun*nyears),8] <- sqrt((dat$proj$cont.ecov$rep$pred_catch[31:40] - dat$truth$pred_catch[31:40])^2)
  df.catch[((irun-1)*nyears+1):(irun*nyears),9] <- sqrt((dat$proj$avg.ecov$rep$pred_catch[31:40]  - dat$truth$pred_catch[31:40])^2)
  df.catch[((irun-1)*nyears+1):(irun*nyears),10]<- sqrt((dat$proj$use.ecov$rep$pred_catch[31:40]  - dat$truth$pred_catch[31:40])^2)  
  df.catch[((irun-1)*nyears+1):(irun*nyears),11]<- ecov_slope
  df.catch[((irun-1)*nyears+1):(irun*nyears),12]<- ssb_cv  
  }
}

saveRDS(df.recr, file=file.path(here::here(),'Ecov_study','recruitment_functions',res.dir,'df.recr.proj.RDS'))
saveRDS(df.ssb,  file=file.path(here::here(),'Ecov_study','recruitment_functions',res.dir,'df.ssb.proj.RDS'))
saveRDS(df.catch,file=file.path(here::here(),'Ecov_study','recruitment_functions',res.dir,'df.catch.proj.RDS'))
