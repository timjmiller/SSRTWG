library(here)
library(tidyverse)
library(wham)
library(rpart)
library(rpart.plot)



#dir <- file.path(here::here(), "Ecov_study","recruitment_functions", "results_pile2")
res.path <- 'E:/results_beta_fix'  # directory where simulation runs are (beta unstandardized)
res.dir <- 'results_beta_fix'   # 'results'     'results_beta_fix'   # results folder where AIC dataframes are
plot.dir <- 'plots_beta_fix'    # 'plots_lizruns'  'plots_beta_fix'   
plot.suffix <- '_beta_fix' 

## specify bad.grad.label and bad.se.value (these are the thresholds set in convergence_summaries.R to determine convergence) 
bad.grad.value <- 1E-6 #tim used 1E-6 (the abs of the exponent will be used for output suffix; ex: filename_grad_6.png)
bad.grad.label <- as.numeric(strsplit(as.character(bad.grad.value), split="-")[[1]][2])
bad.se.value <- 100 #tim used 100 (this value will be used for output suffix; ex: filename_se_100.png)

df.oms          <- readRDS(file.path(here::here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS"))
df.ems    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))
conv.runs <- readRDS(file.path(here::here(),'Ecov_study','recruitment_functions',res.dir, paste0("conv.runs_grad_", bad.grad.label, "_SE_", bad.se.value, plot.suffix, ".RDS")  ) )

conv.runs <- conv.runs %>%
  rename(Sim=sim) %>%
  mutate(ok.run =(bad.opt+ bad.conv+ bad.sdrep+ bad.grad+ bad.se.big) ) %>%
  replace_na(list(ok.run=1)) %>%
  relocate(OM, EM, Sim, ok.run, bad.opt, bad.conv, bad.sdrep, bad.grad, bad.se.big) %>%
  filter(ok.run==0)  # drop unconverged runs (there shouldn't be any at this point)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# n_oms <- nrow(df.oms)
# n_ems <- nrow(df.ems)
# n_sims <- 100
# nyears <- 40
# nyrs.proj <- 10




# get bias (relative error) between peel 10 and truth ====
####  skipping this as it is redundant and time sucking #############################################
# t1 <- Sys.time()
# re.recr.proj <- lapply(1:n_oms,function(om){  #produces list of length n_oms with list of length n_ems, each with 3 projections of length 10 years X n_sims
#   print(paste0("om ",om))
#   lapply(1:n_ems, function(em){
#     sapply(1:n_sims, function(sim){
#       #print(paste0("om ",om," em ",em," sim ",sim))
#       re <- rep(NA, nyrs.proj)
#       if (file.exists( file.path(res.path, paste0("om", om, '/','sim',sim,'_','em',em,'.RDS') ) )  ) {
#       dat <- try(readRDS(file.path(res.path, paste0("om", om, '/','sim',sim,'_','em',em,'.RDS') )  ) )#,
#       # silent=TRUE)   )
#       if(class(dat)!='try-error'){
#         re.cont <- (dat$proj$cont.ecov$rep$NAA[31:40,1] - dat$truth$NAA[31:40,1])/ dat$truth$NAA[31:40,1]   
#         re.avg <- (dat$proj$avg.ecov$rep$NAA[31:40,1] - dat$truth$NAA[31:40,1])/ dat$truth$NAA[31:40,1]   
#         re.use <- (dat$proj$use.ecov$rep$NAA[31:40,1] - dat$truth$NAA[31:40,1])/ dat$truth$NAA[31:40,1]   
#         re <- list(re.cont=re.cont, re.avg=re.avg, re.use=re.use)
#       } #try-error
#       re
#     } #file-exists
#     } #sim function
#     ) #sapply
#   })
# })
# t2 <- Sys.time()
# t2-t1   # Time difference of 1.649097 hours
# saveRDS(re.recr.proj, file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0('re.recr.proj.', plot.suffix, '.RDS')) )
# 
# t3 <- Sys.time()
# re.ssb.proj <- lapply(1:n_oms,function(om){  #produces list of length n_oms with list of length n_ems, each with 3 projections of length 10 years X n_sims
#   print(paste0("om ",om))
#   lapply(1:n_ems, function(em){
#     sapply(1:n_sims, function(sim){
#       #print(paste0("om ",om," em ",em," sim ",sim))
#       re <- rep(NA, nyrs.proj)
#       if (file.exists( file.path(res.path, paste0("om", om, '/','sim',sim,'_','em',em,'.RDS') ) )  ) {
#         dat <- try(readRDS(file.path(res.path, paste0("om", om, '/','sim',sim,'_','em',em,'.RDS') )  ) )#,
#         # silent=TRUE)   )
#         if(class(dat)!='try-error'){
#           re.cont <- (dat$proj$cont.ecov$rep$SSB[31:40] - dat$truth$SSB[31:40])/ dat$truth$SSB[31:40]   
#           re.avg <- (dat$proj$avg.ecov$rep$SSB[31:40] - dat$truth$SSB[31:40])/ dat$truth$SSB[31:40]   
#           re.use <- (dat$proj$use.ecov$rep$SSB[31:40] - dat$truth$SSB[31:40])/ dat$truth$SSB[31:40]   
#           re <- list(re.cont=re.cont, re.avg=re.avg, re.use=re.use)
#         } #try-error
#         re
#       } #file-exists
#     } #sim function
#     ) #sapply
#   })
# })
# t4 <- Sys.time()
# t4-t3    # Time difference of 1.625668 hours
# saveRDS(re.ssb.proj, file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0('re.ssb.proj.', plot.suffix, '.RDS')) )
# 
# 
# 
# t5 <- Sys.time()
# re.catch.proj <- lapply(1:n_oms,function(om){  #produces list of length n_oms with list of length n_ems, each with 3 projections of length 10 years X n_sims
#   print(paste0("om ",om))
#   lapply(1:n_ems, function(em){
#     sapply(1:n_sims, function(sim){
#       #print(paste0("om ",om," em ",em," sim ",sim))
#       re <- rep(NA, nyrs.proj)
#       if (file.exists( file.path(res.path, paste0("om", om, '/','sim',sim,'_','em',em,'.RDS') ) )  ) {
#         dat <- try(readRDS(file.path(res.path, paste0("om", om, '/','sim',sim,'_','em',em,'.RDS') )  ) )#,
#         # silent=TRUE)   )
#         if(class(dat)!='try-error'){
#           re.cont <- (dat$proj$cont.ecov$rep$pred_catch[31:40] - dat$truth$pred_catch[31:40])/ dat$truth$pred_catch[31:40]   
#           re.avg <- (dat$proj$avg.ecov$rep$pred_catch[31:40] - dat$truth$pred_catch[31:40])/ dat$truth$pred_catch[31:40]   
#           re.use <- (dat$proj$use.ecov$rep$pred_catch[31:40] - dat$truth$pred_catch[31:40])/ dat$truth$pred_catch[31:40]   
#           re <- list(re.cont=re.cont, re.avg=re.avg, re.use=re.use)
#         } #try-error
#         re
#       } #file-exists
#     } #sim function
#     ) #sapply
#   })
# })
# t6 <- Sys.time()
# t6-t5   # Time difference of  1.622353 hours
# saveRDS(re.catch.proj, file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0('re.catch.proj.', plot.suffix, '.RDS')) )
# 


n_oms <- nrow(df.oms)
n_ems <- nrow(df.ems)
n_sims <- 100
n_proj <- 3   # number of different projection scenarios
nyears <- 10  # number of projection years




# create objects to store relative error (re) ====
re.recr.df <- matrix(NA, nrow=n_oms*n_ems*nyears*n_sims, ncol=(5+n_proj-1 )   )# OM, EM, Sim, Year, RE, df.om colnames
colnames(re.recr.df) <- c('OM', 'EM', 'Sim', 'Year', paste0(rep('RE',n_proj), seq(1,n_proj)) )# 'Cont', 'Avg', 'Obs'
re.ssb.df <- matrix(NA, nrow=n_oms*n_ems*nyears*n_sims, ncol=(5+n_proj-1 )   )# OM, EM, Sim, Year, RE, df.om colnames
colnames(re.ssb.df) <- c('OM', 'EM', 'Sim', 'Year', paste0(rep('RE',n_proj), seq(1,n_proj)) )
re.catch.df <- matrix(NA, nrow=n_oms*n_ems*nyears*n_sims, ncol=(5+n_proj-1 )   )# OM, EM, Sim, Year, RE, df.om colnames
colnames(re.catch.df) <- c('OM', 'EM', 'Sim', 'Year', paste0(rep('RE',n_proj), seq(1,n_proj)) )


re.recr.sims <- matrix(NA, nrow=n_sims*nyears, ncol=n_proj)
re.ssb.sims <- matrix(NA, nrow=n_sims*nyears, ncol=n_proj)
re.catch.sims  <- matrix(NA, nrow=n_sims*nyears, ncol=n_proj)


# create objects to store log_sd ====
sd.recr.df <- matrix(NA, nrow=n_oms*n_ems*nyears*n_sims, ncol=(5+n_proj-1 )   )# OM, EM, Sim, Year, RE, df.om colnames
colnames(sd.recr.df) <- c('OM', 'EM', 'Sim', 'Year', paste0(rep('SD',n_proj), seq(1,n_proj)) )# 'Cont', 'Avg', 'Obs'
sd.ssb.df <- matrix(NA, nrow=n_oms*n_ems*nyears*n_sims, ncol=(5+n_proj-1 )   )# OM, EM, Sim, Year, RE, df.om colnames
colnames(sd.ssb.df) <- c('OM', 'EM', 'Sim', 'Year', paste0(rep('SD',n_proj), seq(1,n_proj)) )
sd.catch.df <- matrix(NA, nrow=n_oms*n_ems*nyears*n_sims, ncol=(5+n_proj-1 )   )# OM, EM, Sim, Year, RE, df.om colnames
colnames(sd.catch.df) <- c('OM', 'EM', 'Sim', 'Year', paste0(rep('SD',n_proj), seq(1,n_proj)) )


sd.recr.sims <- matrix(NA, nrow=n_sims*nyears, ncol=n_proj)
sd.ssb.sims <- matrix(NA, nrow=n_sims*nyears, ncol=n_proj)
sd.catch.sims  <- matrix(NA, nrow=n_sims*nyears, ncol=n_proj)



# create objects to store AIC and convergence ====
# (not done)


t1 <- Sys.time()
k <- 1
kdf <- 1
for(iom in 1:n_oms) {
  print(paste0("om ", iom ))
  om.conv.runs <- conv.runs[conv.runs$OM==iom,]
  
  for (jem in 1:n_ems){
    kk<-1
    sim.set <- unique(om.conv.runs$Sim[om.conv.runs$EM==jem])
    if (length(sim.set)>0)  {
    for (ksim in first(sim.set):last(sim.set) ) {
      # for (ksim in 1:n_sims) {
      if(file.exists(file.path(res.path, paste0("om", iom, '/','sim', ksim,'_','em', jem,'.RDS') ) )) {
      dat <- try(readRDS(file.path(res.path, paste0("om", iom, '/','sim', ksim,'_','em', jem,'.RDS') ) ) )
      if(class(dat)!='try-error'){
        if(length(names(dat$fit))>0  )  {
          # don't think i need these 2 if statements
          # if(length(re.recr.proj[[iom]][[jem]])==n_sims*n_proj) {
          #   if(length( unlist(re.recr.proj[[iom]][[jem]][1,ksim]) )>0) {
          # re.recr.sims[((ksim-1)*nyears+1):(ksim*nyears),1:n_proj] <- sapply(re.recr.proj[[iom]][[jem]][,ksim], function(x) unlist(x))
          re.recr.sims[((ksim-1)*nyears+1):(ksim*nyears),1:n_proj] <- cbind(re.cont=(dat$proj$cont.ecov$rep$NAA[31:40,1]  - dat$truth$NAA[31:40,1] )/ dat$truth$NAA[31:40,1] ,   re.avg =(dat$proj$avg.ecov$rep$NAA[31:40,1]  - dat$truth$NAA[31:40,1] )/ dat$truth$NAA[31:40,1] ,    re.use =(dat$proj$use.ecov$rep$NAA[31:40,1]  - dat$truth$NAA[31:40,1] )/ dat$truth$pred_catch[31:40]    )
          
          # re.ssb.sims[((ksim-1)*nyears+1):(ksim*nyears),1:n_proj] <- sapply(re.ssb.proj[[iom]][[jem]][,ksim], function(x) unlist(x))
          re.ssb.sims[((ksim-1)*nyears+1):(ksim*nyears),1:n_proj] <- cbind(re.cont=(dat$proj$cont.ecov$rep$SSB[31:40] - dat$truth$SSB[31:40])/ dat$truth$SSB[31:40],   re.avg =(dat$proj$avg.ecov$rep$SSB[31:40] - dat$truth$SSB[31:40])/ dat$truth$SSB[31:40],    re.use =(dat$proj$use.ecov$rep$SSB[31:40] - dat$truth$SSB[31:40])/ dat$truth$SSB[31:40]    )
            
          
          # re.catch.sims[((ksim-1)*nyears+1):(ksim*nyears),1:n_proj] <- sapply(re.catch.proj[[iom]][[jem]][,ksim], function(x) unlist(x))
          re.catch.sims[((ksim-1)*nyears+1):(ksim*nyears),1:n_proj] <- cbind(re.cont=(dat$proj$cont.ecov$rep$pred_catch[31:40] - dat$truth$pred_catch[31:40])/ dat$truth$pred_catch[31:40],   re.avg =(dat$proj$avg.ecov$rep$pred_catch[31:40] - dat$truth$pred_catch[31:40])/ dat$truth$pred_catch[31:40],    re.use =(dat$proj$use.ecov$rep$pred_catch[31:40] - dat$truth$pred_catch[31:40])/ dat$truth$pred_catch[31:40]    )
          
          
          
          # grab log_sd stuff only has years 2:40 (first year is part of N1_pars) =====
          tmp.sdrep.rec <- summary(dat$proj$cont.ecov$sdrep$SE_par)
          tmp.ind.rec <- which(rownames(tmp.sdrep.rec) == "log_NAA")
          sd.recr.sims[((ksim-1)*nyears+1):(ksim*nyears),1:n_proj] <- cbind(sd.cont= tail(dat$proj$cont.ecov$sdrep$SE_par[[tmp.ind.rec]][,1], 10)  ,  sd.avg =  tail(dat$proj$avg.ecov$sdrep$SE_par[[tmp.ind.rec]][,1], 10) ,  sd.use =tail(dat$proj$use.ecov$sdrep$SE_par[[tmp.ind.rec]][,1], 10))
          
          
          tmp.sdrep.ssb <- summary(dat$proj$cont.ecov$sdrep$SE_rep)
          tmp.ind.ssb <- which(rownames(tmp.sdrep.ssb) == "log_SSB")
          sd.ssb.sims[((ksim-1)*nyears+1):(ksim*nyears),1:n_proj] <- cbind(sd.cont= tail(dat$proj$cont.ecov$sdrep$SE_rep[[tmp.ind.ssb]], 10)  ,  sd.avg =  tail(dat$proj$avg.ecov$sdrep$SE_rep[[tmp.ind.ssb]], 10) ,  sd.use =tail(dat$proj$use.ecov$sdrep$SE_rep[[tmp.ind.ssb]], 10) )
          
          
          tmp.sdrep.catch <- summary(dat$proj$cont.ecov$sdrep$SE_rep)
          tmp.ind.catch <- which(rownames(tmp.sdrep.catch) == "log_catch_proj")
          sd.catch.sims[((ksim-1)*nyears+1):(ksim*nyears),1:n_proj] <- cbind(sd.cont= tail(dat$proj$cont.ecov$sdrep$SE_rep[[tmp.ind.catch]], 10)  ,  sd.avg =  tail(dat$proj$avg.ecov$sdrep$SE_rep[[tmp.ind.catch]], 10) ,  sd.use =tail(dat$proj$use.ecov$sdrep$SE_rep[[tmp.ind.catch]], 10) )
          
          
          
          #   } # unlist(re.recr.proj)
          # }  # length(Re.recr.proj)
          # may not need next 5 lines =====
          # if(length(re.recr.proj[[iom]][[jem]])==n_sims) {
          #   if(length( unlist(re.recr.proj[[iom]][[jem]][ksim]) )>0) {
          #     re.recr.sims[((ksim-1)*nyears+1):(ksim*nyears),1:n_proj] <- sapply(re.recr.proj[[iom]][[jem]][[ksim]], function(x) unlist(x))
          #     re.ssb.sims[((ksim-1)*nyears+1):(ksim*nyears),1:n_proj] <- sapply(re.ssb.proj[[iom]][[jem]][[ksim]], function(x) unlist(x))
          #     re.catch.sims[((ksim-1)*nyears+1):(ksim*nyears),1:n_proj] <- sapply(re.catch.proj[[iom]][[jem]][[ksim]], function(x) unlist(x))
          #     
          #   }
          # }
          
            
        
        } # dat$fit>0
        
      }  # class dat
      kk <- kk+1
      kdf=kdf+1
      }  # file.exists
      
      # not necessary since i'm only reading in converged runs
      # if(!file.exists(file.path(res.path, paste0("om", iom, '/','sim', ksim,'_','em', jem,'.RDS') ) ) ) {
      #   re.recr.sims[((ksim-1)*nyears+1):(ksim*nyears),1:n_proj] <- matrix(NA, nrow=nyears, ncol=n_proj)  
      #   re.ssb.sims[((ksim-1)*nyears+1):(ksim*nyears),1:n_proj] <- matrix(NA, nrow=nyears, ncol=n_proj)  
      #   re.catch.sims[((ksim-1)*nyears+1):(ksim*nyears),1:n_proj] <- matrix(NA, nrow=nyears, ncol=n_proj)  
      # }
      
    } #ksim loop
    } #length(sim.set)
    
    # fill in relative error dataframes ====
    re.recr.df [((k-1)*nyears*n_sims+1):(n_sims*nyears*k),(5:(5+n_proj-1)) ] <- re.recr.sims
    re.recr.df[((k-1)*nyears*n_sims+1):(n_sims*nyears*k), 1:4] <- cbind(rep(iom, nyears*n_sims), rep(jem, nyears*n_sims), rep(seq(1,n_sims),each=nyears), rep(seq(1,nyears), n_sims))
    
    re.ssb.df [((k-1)*nyears*n_sims+1):(n_sims*nyears*k),(5:(5+n_proj-1))] <- re.ssb.sims
    re.ssb.df[((k-1)*nyears*n_sims+1):(n_sims*nyears*k), 1:4] <- cbind(rep(iom, nyears*n_sims), rep(jem, nyears*n_sims), rep(seq(1,n_sims),each=nyears), rep(seq(1,nyears), n_sims))

    re.catch.df [((k-1)*nyears*n_sims+1):(n_sims*nyears*k),(5:(5+n_proj-1)) ] <- re.catch.sims
    re.catch.df[((k-1)*nyears*n_sims+1):(n_sims*nyears*k), 1:4] <- cbind(rep(iom, nyears*n_sims), rep(jem, nyears*n_sims), rep(seq(1,n_sims),each=nyears), rep(seq(1,nyears), n_sims))
    
    
    # fill in log_sd dataframes ====
    sd.recr.df [((k-1)*nyears*n_sims+1):(n_sims*nyears*k),(5:(5+n_proj-1)) ] <- sd.recr.sims
    sd.recr.df[((k-1)*nyears*n_sims+1):(n_sims*nyears*k), 1:4] <- cbind(rep(iom, nyears*n_sims), rep(jem, nyears*n_sims), rep(seq(1,n_sims),each=nyears), rep(seq(1,nyears), n_sims))
    
    sd.ssb.df [((k-1)*nyears*n_sims+1):(n_sims*nyears*k),(5:(5+n_proj-1))] <- sd.ssb.sims
    sd.ssb.df[((k-1)*nyears*n_sims+1):(n_sims*nyears*k), 1:4] <- cbind(rep(iom, nyears*n_sims), rep(jem, nyears*n_sims), rep(seq(1,n_sims),each=nyears), rep(seq(1,nyears), n_sims))
    
    sd.catch.df [((k-1)*nyears*n_sims+1):(n_sims*nyears*k),(5:(5+n_proj-1)) ] <- sd.catch.sims
    sd.catch.df[((k-1)*nyears*n_sims+1):(n_sims*nyears*k), 1:4] <- cbind(rep(iom, nyears*n_sims), rep(jem, nyears*n_sims), rep(seq(1,n_sims),each=nyears), rep(seq(1,nyears), n_sims))
    
    
    k <- k+1
    
    # reset sims matrices to NA ====
    re.recr.sims <- matrix(NA, nrow=n_sims*nyears, ncol=n_proj)
    re.ssb.sims <- matrix(NA, nrow=n_sims*nyears, ncol=n_proj)
    re.catch.sims <- matrix(NA, nrow=n_sims*nyears, ncol=n_proj)
    sd.recr.sims <- matrix(NA, nrow=n_sims*nyears, ncol=n_proj)
    sd.ssb.sims <- matrix(NA, nrow=n_sims*nyears, ncol=n_proj)
    sd.catch.sims <- matrix(NA, nrow=n_sims*nyears, ncol=n_proj)
    
    
  } #jem loop
} #iom loop
t2 <- Sys.time()
t2-t1   # Time difference of 1.63261 hours

saveRDS(re.recr.df, file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "re.recr.proj.df", plot.suffix, ".RDS") ) )
saveRDS(re.ssb.df, file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0("re.ssb.proj.df", plot.suffix, ".RDS") ) )
saveRDS(re.catch.df, file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0("re.catch.proj.df", plot.suffix, ".RDS") ) )


saveRDS(sd.recr.df, file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "sd.recr.proj.df", plot.suffix, ".RDS") ) )
saveRDS(sd.ssb.df, file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0("sd.ssb.proj.df", plot.suffix, ".RDS") ) )
saveRDS(sd.catch.df, file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0("sd.catch.proj.df", plot.suffix, ".RDS") ) )

########################################################################
# read in RDS if already exist ==================
re.recr.df<- readRDS(file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "re.recr.proj.df", plot.suffix, ".RDS") ) )
re.ssb.df<- readRDS( file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0("re.ssb.proj.df", plot.suffix, ".RDS") ) )
re.catch.df<- readRDS( file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0("re.catch.proj.df", plot.suffix, ".RDS") ) )


sd.recr.df<- readRDS(file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "sd.recr.proj.df", plot.suffix, ".RDS") ) )
sd.ssb.df<- readRDS( file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0("sd.ssb.proj.df", plot.suffix, ".RDS") ) )
sd.catch.df<- readRDS(file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0("sd.catch.proj.df", plot.suffix, ".RDS") ) )
########################################################################


df.oms2 <- as_tibble(cbind(OM=seq(1,256), df.oms)) 
df.ems2 <- cbind(EM=seq(1,6), df.ems)
colnames(df.ems2) <- c("EM", "EM_ecov_how", "EM_r_mod") 
em_tib <- as_tibble(df.ems2) %>%
  mutate(SR=ifelse(EM_r_mod==2, 'Mean', 'BH')) %>%
  mutate(EM_mod = paste0(SR, "_", EM_ecov_how))


# recruitment ====
## relative error (re)
re.rec.tib <- as_tibble(re.recr.df)  %>%
  mutate_at(c('RE1', 'RE2', 'RE3'), as.numeric) %>%
  mutate_at(c('OM', 'EM', 'Sim', 'Year'), as.integer) %>%
  left_join(df.oms2) %>%
  left_join(em_tib) %>%
  mutate(mod.match=ifelse(EM_ecov_how==Ecov_how & EM_r_mod==recruit_mod, 1, 0)) %>%
  drop_na()
  # left_join(no.conv.runs) %>%                 #1=bad, 0=not bad
  # mutate(ok.run =(bad.opt+ bad.conv+ bad.sdrep+ bad.grad+ bad.se.big) ) %>%
  # replace_na(list(ok.run=1)) %>%
  # relocate(OM, EM, Sim, bad.opt, bad.conv, bad.sdrep, bad.grad, bad.se.big, ok.run) %>%
  # filter(ok.run==0)  # drop unconverged runs


re.rec.tib$R_sig <- factor(re.rec.tib$R_sig,labels = c("Rsig_0.1","Rsig_1.0"))
re.rec.tib$Ecov_effect <- factor(re.rec.tib$Ecov_effect,labels=c("Ecov_L","Ecov_H"))
re.rec.tib$Ecov_how    <- factor(re.rec.tib$Ecov_how,labels=c("Ecov_0", "Ecov_1","Ecov_2","Ecov_4"))
re.rec.tib$NAA_cor     <- factor(re.rec.tib$NAA_cor,labels=c("Rcor_L","Rcor_H"))
re.rec.tib$Fhist       <- factor(re.rec.tib$Fhist,labels=c("H-MSY","MSY") ) 
re.rec.tib$Ecov_re_cor <- factor(re.rec.tib$Ecov_re_cor,labels=c("EcovCor_L","EcovCor_H"))
re.rec.tib$obs_error   <- factor(re.rec.tib$obs_error,labels=c("ObsErr_H","ObsErr_L"))


## log_sd (sd)
sd.rec.tib <- as_tibble(sd.recr.df)  %>%
  mutate_at(c('SD1', 'SD2', 'SD3'), as.numeric) %>%
  mutate_at(c('OM', 'EM', 'Sim', 'Year'), as.integer) %>%
  left_join(df.oms2) %>%
  left_join(em_tib) %>%
  mutate(mod.match=ifelse(EM_ecov_how==Ecov_how & EM_r_mod==recruit_mod, 1, 0)) %>%
  drop_na()


sd.rec.tib$R_sig <- factor(sd.rec.tib$R_sig,labels = c("Rsig_0.1","Rsig_1.0"))
sd.rec.tib$Ecov_effect <- factor(sd.rec.tib$Ecov_effect,labels=c("Ecov_L","Ecov_H"))
sd.rec.tib$Ecov_how    <- factor(sd.rec.tib$Ecov_how,labels=c("Ecov_0", "Ecov_1","Ecov_2","Ecov_4"))
sd.rec.tib$NAA_cor     <- factor(sd.rec.tib$NAA_cor,labels=c("Rcor_L","Rcor_H"))
sd.rec.tib$Fhist       <- factor(sd.rec.tib$Fhist,labels=c("H-MSY","MSY") ) 
sd.rec.tib$Ecov_re_cor <- factor(sd.rec.tib$Ecov_re_cor,labels=c("EcovCor_L","EcovCor_H"))
sd.rec.tib$obs_error   <- factor(sd.rec.tib$obs_error,labels=c("ObsErr_H","ObsErr_L"))


# %%%%%%%%%%%%%%%%%%%%%%
# ssb ====
re.ssb.tib <- as_tibble(re.ssb.df)  %>%
  mutate_at(c('RE1', 'RE2', 'RE3'), as.numeric) %>%
  mutate_at(c('OM', 'EM', 'Sim', 'Year'), as.integer) %>%
  left_join(df.oms2) %>%
  left_join(em_tib) %>%
  mutate(mod.match=ifelse(EM_ecov_how==Ecov_how & EM_r_mod==recruit_mod, 1, 0)) %>%
  drop_na()


re.ssb.tib$R_sig <- factor(re.ssb.tib$R_sig,labels = c("Rsig_0.1","Rsig_1.0"))
re.ssb.tib$Ecov_effect <- factor(re.ssb.tib$Ecov_effect,labels=c("Ecov_L","Ecov_H"))
re.ssb.tib$Ecov_how    <- factor(re.ssb.tib$Ecov_how,labels=c("Ecov_0", "Ecov_1","Ecov_2","Ecov_4"))
re.ssb.tib$NAA_cor     <- factor(re.ssb.tib$NAA_cor,labels=c("Rcor_L","Rcor_H"))
re.ssb.tib$Fhist       <- factor(re.ssb.tib$Fhist,labels=c("H-MSY","MSY") ) 
re.ssb.tib$Ecov_re_cor <- factor(re.ssb.tib$Ecov_re_cor,labels=c("EcovCor_L","EcovCor_H"))
re.ssb.tib$obs_error   <- factor(re.ssb.tib$obs_error,labels=c("ObsErr_H","ObsErr_L"))



## log_sd (sd)
sd.ssb.tib <- as_tibble(sd.ssb.df)  %>%
  mutate_at(c('SD1', 'SD2', 'SD3'), as.numeric) %>%
  mutate_at(c('OM', 'EM', 'Sim', 'Year'), as.integer) %>%
  left_join(df.oms2) %>%
  left_join(em_tib) %>%
  mutate(mod.match=ifelse(EM_ecov_how==Ecov_how & EM_r_mod==recruit_mod, 1, 0)) %>%
  drop_na()


sd.ssb.tib$R_sig <- factor(sd.ssb.tib$R_sig,labels = c("Rsig_0.1","Rsig_1.0"))
sd.ssb.tib$Ecov_effect <- factor(sd.ssb.tib$Ecov_effect,labels=c("Ecov_L","Ecov_H"))
sd.ssb.tib$Ecov_how    <- factor(sd.ssb.tib$Ecov_how,labels=c("Ecov_0", "Ecov_1","Ecov_2","Ecov_4"))
sd.ssb.tib$NAA_cor     <- factor(sd.ssb.tib$NAA_cor,labels=c("Rcor_L","Rcor_H"))
sd.ssb.tib$Fhist       <- factor(sd.ssb.tib$Fhist,labels=c("H-MSY","MSY") ) 
sd.ssb.tib$Ecov_re_cor <- factor(sd.ssb.tib$Ecov_re_cor,labels=c("EcovCor_L","EcovCor_H"))
sd.ssb.tib$obs_error   <- factor(sd.ssb.tib$obs_error,labels=c("ObsErr_H","ObsErr_L"))



# catch ====
re.catch.tib <- as_tibble(re.catch.df)  %>%
  mutate_at(c('RE1', 'RE2', 'RE3'), as.numeric) %>%
  mutate_at(c('OM', 'EM', 'Sim', 'Year'), as.integer) %>%
  left_join(df.oms2) %>%
  left_join(em_tib) %>%
  mutate(mod.match=ifelse(EM_ecov_how==Ecov_how & EM_r_mod==recruit_mod, 1, 0)) %>%
  drop_na()


re.catch.tib$R_sig <- factor(re.catch.tib$R_sig,labels = c("Rsig_0.1","Rsig_1.0"))
re.catch.tib$Ecov_effect <- factor(re.catch.tib$Ecov_effect,labels=c("Ecov_L","Ecov_H"))
re.catch.tib$Ecov_how    <- factor(re.catch.tib$Ecov_how,labels=c("Ecov_0", "Ecov_1","Ecov_2","Ecov_4"))
re.catch.tib$NAA_cor     <- factor(re.catch.tib$NAA_cor,labels=c("Rcor_L","Rcor_H"))
re.catch.tib$Fhist       <- factor(re.catch.tib$Fhist,labels=c("H-MSY","MSY") ) 
re.catch.tib$Ecov_re_cor <- factor(re.catch.tib$Ecov_re_cor,labels=c("EcovCor_L","EcovCor_H"))
re.catch.tib$obs_error   <- factor(re.catch.tib$obs_error,labels=c("ObsErr_H","ObsErr_L"))

## log_sd (sd)
sd.catch.tib <- as_tibble(sd.catch.df)  %>%
  mutate_at(c('SD1', 'SD2', 'SD3'), as.numeric) %>%
  mutate_at(c('OM', 'EM', 'Sim', 'Year'), as.integer) %>%
  left_join(df.oms2) %>%
  left_join(em_tib) %>%
  mutate(mod.match=ifelse(EM_ecov_how==Ecov_how & EM_r_mod==recruit_mod, 1, 0)) %>%
  drop_na()


sd.catch.tib$R_sig <- factor(sd.catch.tib$R_sig,labels = c("Rsig_0.1","Rsig_1.0"))
sd.catch.tib$Ecov_effect <- factor(sd.catch.tib$Ecov_effect,labels=c("Ecov_L","Ecov_H"))
sd.catch.tib$Ecov_how    <- factor(sd.catch.tib$Ecov_how,labels=c("Ecov_0", "Ecov_1","Ecov_2","Ecov_4"))
sd.catch.tib$NAA_cor     <- factor(sd.catch.tib$NAA_cor,labels=c("Rcor_L","Rcor_H"))
sd.catch.tib$Fhist       <- factor(sd.catch.tib$Fhist,labels=c("H-MSY","MSY") ) 
sd.catch.tib$Ecov_re_cor <- factor(sd.catch.tib$Ecov_re_cor,labels=c("EcovCor_L","EcovCor_H"))
sd.catch.tib$obs_error   <- factor(sd.catch.tib$obs_error,labels=c("ObsErr_H","ObsErr_L"))


## regression trees by projection method (recr) ====
rf_RE1.rec   <- rpart(RE1   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=re.rec.tib, control=rpart.control(cp=0.01)) #  ====
imp.var <- rf_RE1.rec$frame[rf_RE1.rec$frame$var != '<leaf>',]
nodes_RE1.rec <- unique(imp.var[,1])
nodes_RE1.rec
# nothing    (restricting to EM 2,4,5,6  did not change this) data=re.rec.tib[re.rec.tib$EM %in% c(2,4,5,6),]; Year also not important; restricting to years 1-3 didn't matter   [re.rec.tib$Year<4,]



rf_RE2.rec   <- rpart(RE2   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=re.rec.tib, control=rpart.control(cp=0.01)) #  ====
imp.var <- rf_RE2.rec$frame[rf_RE2.rec$frame$var != '<leaf>',]
nodes_RE2.rec <- unique(imp.var[,1])
nodes_RE2.rec
# nothing    (restricting to EM 2,4,5,6  did not change this) data=re.rec.tib[re.rec.tib$EM %in% c(2,4,5,6),] Year also not important


rf_RE3.rec   <- rpart(RE3   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=re.rec.tib, control=rpart.control(cp=0.01)) #  ====
imp.var <- rf_RE3.rec$frame[rf_RE3.rec$frame$var != '<leaf>',]
nodes_RE3.rec <- unique(imp.var[,1])
nodes_RE3.rec
# nothing    (restricting to EM 2,4,5,6  did not change this) data=re.rec.tib[re.rec.tib$EM %in% c(2,4,5,6),] Year also not important


rf_SD1.rec   <- rpart(SD1   ~  R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=sd.rec.tib, control=rpart.control(cp=0.01)) #  ====
imp.var <- rf_SD1.rec$frame[rf_SD1.rec$frame$var != '<leaf>',]
nodes_SD1.rec <- unique(imp.var[,1])
nodes_SD1.rec
# nothing    (restricting to EM 2,4,5,6  did not change this) data=sd.rec.tib[sd.rec.tib$EM %in% c(2,4,5,6),]; Year also not important



rf_SD2.rec   <- rpart(SD2   ~  R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=sd.rec.tib, control=rpart.control(cp=0.01)) #  ====
imp.var <- rf_SD2.rec$frame[rf_SD2.rec$frame$var != '<leaf>',]
nodes_SD2.rec <- unique(imp.var[,1])
nodes_SD2.rec
# nothing    (restricting to EM 2,4,5,6  did not change this) data=sd.rec.tib[sd.rec.tib$EM %in% c(2,4,5,6),]; Year also not important


rf_SD3.rec   <- rpart(SD3   ~  R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=sd.rec.tib, control=rpart.control(cp=0.01)) #  ====
imp.var <- rf_SD3.rec$frame[rf_SD3.rec$frame$var != '<leaf>',]
nodes_SD3.rec <- unique(imp.var[,1])
nodes_SD3.rec
# nothing    (restricting to EM 2,4,5,6  did not change this) data=sd.rec.tib[sd.rec.tib$EM %in% c(2,4,5,6),]; Year also not important




## regression trees by projection method (ssb) ====
rf_RE1.ssb   <- rpart(RE1   ~  R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=re.ssb.tib, control=rpart.control(cp=0.01)) #  ====
imp.var <- rf_RE1.ssb$frame[rf_RE1.ssb$frame$var != '<leaf>',]
nodes_RE1.ssb <- unique(imp.var[,1])
nodes_RE1.ssb
# nothing    (restricting to EM 2,4,5,6  did not change this) data=re.ssb.tib[re.ssb.tib$EM %in% c(2,4,5,6),]; Year also not important



rf_RE2.ssb   <- rpart(RE2   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=re.ssb.tib, control=rpart.control(cp=0.01)) #  ====
imp.var <- rf_RE2.ssb$frame[rf_RE2.ssb$frame$var != '<leaf>',]
nodes_RE2.ssb <- unique(imp.var[,1])
nodes_RE2.ssb
# nothing    (restricting to EM 2,4,5,6  did not change this) data=re.ssb.tib[re.ssb.tib$EM %in% c(2,4,5,6),]


rf_RE3.ssb   <- rpart(RE3   ~  R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=re.ssb.tib, control=rpart.control(cp=0.01)) #  ====
imp.var <- rf_RE3.ssb$frame[rf_RE3.ssb$frame$var != '<leaf>',]
nodes_RE3.ssb <- unique(imp.var[,1])
nodes_RE3.ssb
# nothing    (restricting to EM 2,4,5,6  did not change this) data=re.ssb.tib[re.ssb.tib$EM %in% c(2,4,5,6),]



rf_SD1.ssb   <- rpart(SD1   ~  R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=sd.ssb.tib, control=rpart.control(cp=0.01)) #  ====
imp.var <- rf_SD1.ssb$frame[rf_SD1.ssb$frame$var != '<leaf>',]
nodes_SD1.ssb <- unique(imp.var[,1])
nodes_SD1.ssb
# nothing    (restricting to EM 2,4,5,6  did not change this) data=sd.ssb.tib[sd.ssb.tib$EM %in% c(2,4,5,6),]



rf_SD2.ssb   <- rpart(SD2   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=sd.ssb.tib, control=rpart.control(cp=0.01)) #  ====
imp.var <- rf_SD2.ssb$frame[rf_SD2.ssb$frame$var != '<leaf>',]
nodes_SD2.ssb <- unique(imp.var[,1])
nodes_SD2.ssb
# nothing    (restricting to EM 2,4,5,6  did not change this) data=sd.ssb.tib[sd.ssb.tib$EM %in% c(2,4,5,6),]


rf_SD3.ssb   <- rpart(SD3   ~  R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=sd.ssb.tib, control=rpart.control(cp=0.01)) #  ====
imp.var <- rf_SD3.ssb$frame[rf_SD3.ssb$frame$var != '<leaf>',]
nodes_SD3.ssb <- unique(imp.var[,1])
nodes_SD3.ssb
# nothing    (restricting to EM 2,4,5,6  did not change this) data=sd.ssb.tib[sd.ssb.tib$EM %in% c(2,4,5,6),]



## regression trees by projection method catch) ====
rf_RE1.catch   <- rpart(RE1   ~  R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=re.catch.tib, control=rpart.control(cp=0.01)) #  ====
imp.var <- rf_RE1.catch$frame[rf_RE1.catch$frame$var != '<leaf>',]
nodes_RE1.catch <- unique(imp.var[,1])
nodes_RE1.catch
# nothing    (restricting to EM 2,4,5,6  did not change this) data=re.catch.tib[re.catch.tib$EM %in% c(2,4,5,6),]



rf_RE2.catch   <- rpart(RE2   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=re.catch.tib, control=rpart.control(cp=0.01)) #  ====
imp.var <- rf_RE2.catch$frame[rf_RE2.catch$frame$var != '<leaf>',]
nodes_RE2.catch <- unique(imp.var[,1])
nodes_RE2.catch
# nothing    (restricting to EM 2,4,5,6  did not change this) data=re.catch.tib[re.catch.tib$EM %in% c(2,4,5,6),]


rf_RE3.catch   <- rpart(RE3   ~  R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=re.catch.tib, control=rpart.control(cp=0.01)) #  ====
imp.var <- rf_RE3.catch$frame[rf_RE3.catch$frame$var != '<leaf>',]
nodes_RE3.catch <- unique(imp.var[,1])
nodes_RE3.catch
# nothing    (restricting to EM 2,4,5,6  did not change this) data=re.catch.tib[re.catch.tib$EM %in% c(2,4,5,6),]



rf_SD1.catch   <- rpart(SD1   ~  R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=sd.catch.tib, control=rpart.control(cp=0.01)) #  ====
imp.var <- rf_SD1.catch$frame[rf_SD1.catch$frame$var != '<leaf>',]
nodes_SD1.catch <- unique(imp.var[,1])
nodes_SD1.catch
# nothing    (restricting to EM 2,4,5,6  did not change this) data=sd.catch.tib[sd.catch.tib$EM %in% c(2,4,5,6),]



rf_SD2.catch   <- rpart(SD2   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=sd.catch.tib, control=rpart.control(cp=0.01)) #  ====
imp.var <- rf_SD2.catch$frame[rf_SD2.catch$frame$var != '<leaf>',]
nodes_SD2.catch <- unique(imp.var[,1])
nodes_SD2.catch
# nothing    (restricting to EM 2,4,5,6  did not change this) data=sd.catch.tib[sd.catch.tib$EM %in% c(2,4,5,6),]


rf_SD3.catch   <- rpart(SD3   ~  R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=sd.catch.tib, control=rpart.control(cp=0.01)) #  ====
imp.var <- rf_SD3.catch$frame[rf_SD3.catch$frame$var != '<leaf>',]
nodes_SD3.catch <- unique(imp.var[,1])
nodes_SD3.catch
# nothing    (restricting to EM 2,4,5,6  did not change this) data=sd.catch.tib[sd.catch.tib$EM %in% c(2,4,5,6),]



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# summarize by important factors (none identified, just summarizing by R_sigma and Ecov_how) =====

## summarize  ALL Years (Relative Error) ====
### recruitment ====
re.rec.sum <- re.rec.tib %>%
  group_by( R_sig,  Ecov_how, EM, mod.match ) %>%
  summarise(RE1_mean=mean(RE1), RE1_median=median(RE1), RE1_var=var(RE1), across(RE1, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                  p97.5=~quantile(.,probs=0.975)) ),
            
            RE2_mean=mean(RE2), RE2_median=median(RE2), RE2_var=var(RE2), across(RE2, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                       p97.5=~quantile(.,probs=0.975)) ),
            
            RE3_mean=mean(RE3), RE3_median=median(RE3), RE3_var=var(RE3), across(RE3, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                       p97.5=~quantile(.,probs=0.975)) ),
  )%>%
  ungroup() %>%
  left_join(em_tib) %>%
  pivot_longer(cols=RE1_mean:RE3_p97.5, names_to=c("Projection", "Metric"), names_pattern="RE?(.*)_(.*)", values_to="Value") %>%
  relocate(Projection, Metric, Value) %>%
  pivot_wider(names_from=Metric,values_from=Value ) %>%
  relocate(Projection, median, var, p2.5, p97.5) %>%
  mutate(Proj.case=case_when(
    Projection=='1' ~ "Continue",
    Projection=='2' ~ "Avg5yr",
    Projection=='3' ~ "ObsEcov"
  ))

### ssb ====
re.ssb.sum <- re.ssb.tib %>%
  group_by( R_sig,  Ecov_how, EM, mod.match ) %>%
  summarise(RE1_mean=mean(RE1), RE1_median=median(RE1), RE1_var=var(RE1), across(RE1, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            RE2_mean=mean(RE2), RE2_median=median(RE2), RE2_var=var(RE2), across(RE2, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            RE3_mean=mean(RE3), RE3_median=median(RE3), RE3_var=var(RE3), across(RE3, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
  )%>%
  ungroup() %>%
  left_join(em_tib) %>%
  pivot_longer(cols=RE1_mean:RE3_p97.5, names_to=c("Projection", "Metric"), names_pattern="RE?(.*)_(.*)", values_to="Value") %>%
  relocate(Projection, Metric, Value) %>%
  pivot_wider(names_from=Metric,values_from=Value ) %>%
  relocate(Projection, median, var, p2.5, p97.5) %>%
  mutate(Proj.case=case_when(
    Projection=='1' ~ "Continue",
    Projection=='2' ~ "Avg5yr",
    Projection=='3' ~ "ObsEcov"
  ))


### catch ====
re.catch.sum <- re.catch.tib %>%
  group_by( R_sig,  Ecov_how, EM, mod.match ) %>%
  summarise(RE1_mean=mean(RE1), RE1_median=median(RE1), RE1_var=var(RE1), across(RE1, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            RE2_mean=mean(RE2), RE2_median=median(RE2), RE2_var=var(RE2), across(RE2, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            RE3_mean=mean(RE3), RE3_median=median(RE3), RE3_var=var(RE3), across(RE3, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
  )%>%
  ungroup() %>%
  left_join(em_tib) %>%
  pivot_longer(cols=RE1_mean:RE3_p97.5, names_to=c("Projection", "Metric"), names_pattern="RE?(.*)_(.*)", values_to="Value") %>%
  relocate(Projection, Metric, Value) %>%
  pivot_wider(names_from=Metric,values_from=Value ) %>%
  relocate(Projection, median, var, p2.5, p97.5) %>%
  mutate(Proj.case=case_when(
    Projection=='1' ~ "Continue",
    Projection=='2' ~ "Avg5yr",
    Projection=='3' ~ "ObsEcov"
  ))



# %%%%%%%%%%%%%%%%%%%
## summarize  ALL Years (log_sd) ====
### recruitment ====
sd.rec.sum <- sd.rec.tib %>%
  group_by( R_sig,  Ecov_how, EM, mod.match ) %>%
  summarise(SD1_mean=mean(SD1), SD1_median=median(SD1), SD1_var=var(SD1), across(SD1, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            SD2_mean=mean(SD2), SD2_median=median(SD2), SD2_var=var(SD2), across(SD2, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            SD3_mean=mean(SD3), SD3_median=median(SD3), SD3_var=var(SD3), across(SD3, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
  )%>%
  ungroup() %>%
  left_join(em_tib) %>%
  pivot_longer(cols=SD1_mean:SD3_p97.5, names_to=c("Projection", "Metric"), names_pattern="SD?(.*)_(.*)", values_to="Value") %>%
  relocate(Projection, Metric, Value) %>%
  pivot_wider(names_from=Metric,values_from=Value ) %>%
  relocate(Projection, median, var, p2.5, p97.5) %>%
  mutate(Proj.case=case_when(
    Projection=='1' ~ "Continue",
    Projection=='2' ~ "Avg5yr",
    Projection=='3' ~ "ObsEcov"
  ))

### ssb ====
sd.ssb.sum <- sd.ssb.tib %>%
  group_by( R_sig,  Ecov_how, EM, mod.match ) %>%
  summarise(SD1_mean=mean(SD1), SD1_median=median(SD1), SD1_var=var(SD1), across(SD1, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            SD2_mean=mean(SD2), SD2_median=median(SD2), SD2_var=var(SD2), across(SD2, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            SD3_mean=mean(SD3), SD3_median=median(SD3), SD3_var=var(SD3), across(SD3, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
  )%>%
  ungroup() %>%
  left_join(em_tib) %>%
  pivot_longer(cols=SD1_mean:SD3_p97.5, names_to=c("Projection", "Metric"), names_pattern="SD?(.*)_(.*)", values_to="Value") %>%
  relocate(Projection, Metric, Value) %>%
  pivot_wider(names_from=Metric,values_from=Value ) %>%
  relocate(Projection, median, var, p2.5, p97.5) %>%
  mutate(Proj.case=case_when(
    Projection=='1' ~ "Continue",
    Projection=='2' ~ "Avg5yr",
    Projection=='3' ~ "ObsEcov"
  ))


### catch ====
sd.catch.sum <- sd.catch.tib %>%
  group_by( R_sig,  Ecov_how, EM, mod.match ) %>%
  summarise(SD1_mean=mean(SD1), SD1_median=median(SD1), SD1_var=var(SD1), across(SD1, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            SD2_mean=mean(SD2), SD2_median=median(SD2), SD2_var=var(SD2), across(SD2, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            SD3_mean=mean(SD3), SD3_median=median(SD3), SD3_var=var(SD3), across(SD3, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
  )%>%
  ungroup() %>%
  left_join(em_tib) %>%
  pivot_longer(cols=SD1_mean:SD3_p97.5, names_to=c("Projection", "Metric"), names_pattern="SD?(.*)_(.*)", values_to="Value") %>%
  relocate(Projection, Metric, Value) %>%
  pivot_wider(names_from=Metric,values_from=Value ) %>%
  relocate(Projection, median, var, p2.5, p97.5) %>%
  mutate(Proj.case=case_when(
    Projection=='1' ~ "Continue",
    Projection=='2' ~ "Avg5yr",
    Projection=='3' ~ "ObsEcov"
  ))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PLOTS - ALL Years  ====

## RE recruitment ====
re.rec.proj.plot <- ggplot(re.rec.sum, aes(x=EM_mod, y=median, col=as.factor(mod.match) )) +
  facet_grid(R_sig + Ecov_how ~ Proj.case , scales="free_y"    ) +
  geom_point(size=2) +
  geom_hline(yintercept=0.0, col='#111111dd') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.4) )+
  ylab('RE(Recr) (2.5 - 50 - 97.5 percentiles)') +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  ggtitle('RE (across all 10 yrs of projection), 3 different Ecov assumptions in projections')
ggsave(re.rec.proj.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('re.rec.proj.plot', plot.suffix, '.png') ),  height=7, width=12)


## same plot, narrower y-lim
re.rec.proj.plot.ylim <- ggplot(re.rec.sum, aes(x=EM_mod, y=median, col=as.factor(mod.match) )) +
  facet_grid(R_sig + Ecov_how ~ Proj.case , scales="free_y"    ) +
  geom_point(size=2) +
  geom_hline(yintercept=0.0, col='#111111dd') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.4) )+
  ylab('RE(Recr) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-0.25,0.25)) +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  ggtitle('RE (across all 10 yrs of projection), 3 different Ecov assumptions in projections')
ggsave(re.rec.proj.plot.ylim, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('re.rec.proj.plot.ylim', plot.suffix, '.png') ),  height=7, width=12)


## SD recruitment 
## RE recruitment ====
sd.rec.proj.plot <- ggplot(sd.rec.sum, aes(x=EM_mod, y=median, col=as.factor(mod.match) )) +
  facet_grid(R_sig + Ecov_how ~ Proj.case , scales="free_y"    ) +
  geom_point(size=2) +
  geom_hline(yintercept=0.0, col='#111111dd') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.4) )+
  ylab('SD(Recr) (2.5 - 50 - 97.5 percentiles)') +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  ggtitle('CV (across all 10 yrs of projection), 3 different Ecov assumptions in projections')
ggsave(sd.rec.proj.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('sd.rec.proj.plot', plot.suffix, '.png') ),  height=7, width=12)





## RE ssb ====
re.ssb.proj.plot <- ggplot(re.ssb.sum, aes(x=EM_mod, y=median, col=as.factor(mod.match) )) +
  facet_grid(R_sig + Ecov_how ~ Proj.case , scales="free_y"    ) +
  geom_point(size=2) +
  geom_hline(yintercept=0.0, col='#111111dd') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.4) )+
  ylab('RE(SSB) (2.5 - 50 - 97.5 percentiles)') +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  ggtitle('RE (across all 10 yrs of projection), 3 different Ecov assumptions in projections')
ggsave(re.ssb.proj.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('re.ssb.proj.plot', plot.suffix, '.png') ),  height=7, width=12)


## same plot, narrower y-lim
re.ssb.proj.plot.ylim <- ggplot(re.ssb.sum, aes(x=EM_mod, y=median, col=as.factor(mod.match) )) +
  facet_grid(R_sig + Ecov_how ~ Proj.case , scales="free_y"    ) +
  geom_point(size=2) +
  geom_hline(yintercept=0.0, col='#111111dd') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.4) )+
  ylab('RE(SSB) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-0.25,0.25)) +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  ggtitle('RE (across all 10 yrs of projection), 3 different Ecov assumptions in projections')
ggsave(re.ssb.proj.plot.ylim, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('re.ssb.proj.plot.ylim', plot.suffix, '.png') ),  height=7, width=12)



## SD ssb ====
sd.ssb.proj.plot <- ggplot(sd.ssb.sum, aes(x=EM_mod, y=median, col=as.factor(mod.match) )) +
  facet_grid(R_sig + Ecov_how ~ Proj.case , scales="free_y"    ) +
  geom_point(size=2) +
  geom_hline(yintercept=0.0, col='#111111dd') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.4) )+
  ylab('CV(SSB) (2.5 - 50 - 97.5 percentiles)') +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  ggtitle('CV (across all 10 yrs of projection), 3 different Ecov assumptions in projections')
ggsave(sd.ssb.proj.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('sd.ssb.proj.plot', plot.suffix, '.png') ),  height=7, width=12)





## RE catch ====
re.catch.proj.plot <- ggplot(re.catch.sum, aes(x=EM_mod, y=median, col=as.factor(mod.match) )) +
  facet_grid(R_sig + Ecov_how ~ Proj.case , scales="free_y"    ) +
  geom_point(size=2) +
  geom_hline(yintercept=0.0, col='#111111dd') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.4) )+
  ylab('RE(Catch) (2.5 - 50 - 97.5 percentiles)') +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  ggtitle('RE (across all 10 yrs of projection), 3 different Ecov assumptions in projections')
ggsave(re.catch.proj.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('re.catch.proj.plot', plot.suffix, '.png') ),  height=7, width=12)


## same plot, narrower y-lim
re.catch.proj.plot.ylim <- ggplot(re.catch.sum, aes(x=EM_mod, y=median, col=as.factor(mod.match) )) +
  facet_grid(R_sig + Ecov_how ~ Proj.case , scales="free_y"    ) +
  geom_point(size=2) +
  geom_hline(yintercept=0.0, col='#111111dd') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.4) )+
  ylab('RE(Catch) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-0.25,0.25)) +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  ggtitle('RE (across all 10 yrs of projection), 3 different Ecov assumptions in projections')
ggsave(re.catch.proj.plot.ylim, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('re.catch.proj.plot.ylim', plot.suffix, '.png') ),  height=7, width=12)


## SD catch ====
sd.catch.proj.plot <- ggplot(sd.catch.sum, aes(x=EM_mod, y=median, col=as.factor(mod.match) )) +
  facet_grid(R_sig + Ecov_how ~ Proj.case , scales="free_y"    ) +
  geom_point(size=2) +
  geom_hline(yintercept=0.0, col='#111111dd') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.4) )+
  ylab('CV(Catch) (2.5 - 50 - 97.5 percentiles)') +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  ggtitle('CV (across all 10 yrs of projection), 3 different Ecov assumptions in projections')
ggsave(sd.catch.proj.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('sd.catch.proj.plot', plot.suffix, '.png') ),  height=7, width=12)




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# summarize by important factors (none identified, just summarizing by R_sigma and Ecov_how) =====

## summarize  FIRST 3 Years (Relative Error) ====
### recruitment ====
re.rec.sum3yr <- re.rec.tib %>%
  filter(Year<4) %>%
  group_by( R_sig,  Ecov_how, EM, mod.match ) %>%
  summarise(RE1_mean=mean(RE1), RE1_median=median(RE1), RE1_var=var(RE1), across(RE1, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            RE2_mean=mean(RE2), RE2_median=median(RE2), RE2_var=var(RE2), across(RE2, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            RE3_mean=mean(RE3), RE3_median=median(RE3), RE3_var=var(RE3), across(RE3, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
  )%>%
  ungroup() %>%
  left_join(em_tib) %>%
  pivot_longer(cols=RE1_mean:RE3_p97.5, names_to=c("Projection", "Metric"), names_pattern="RE?(.*)_(.*)", values_to="Value") %>%
  relocate(Projection, Metric, Value) %>%
  pivot_wider(names_from=Metric,values_from=Value ) %>%
  relocate(Projection, median, var, p2.5, p97.5) %>%
  mutate(Proj.case=case_when(
    Projection=='1' ~ "Continue",
    Projection=='2' ~ "Avg5yr",
    Projection=='3' ~ "ObsEcov"
  ))

### ssb ====
re.ssb.sum3yr <- re.ssb.tib %>%
  filter(Year<4) %>%
  group_by( R_sig,  Ecov_how, EM, mod.match ) %>%
  summarise(RE1_mean=mean(RE1), RE1_median=median(RE1), RE1_var=var(RE1), across(RE1, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            RE2_mean=mean(RE2), RE2_median=median(RE2), RE2_var=var(RE2), across(RE2, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            RE3_mean=mean(RE3), RE3_median=median(RE3), RE3_var=var(RE3), across(RE3, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
  )%>%
  ungroup() %>%
  left_join(em_tib) %>%
  pivot_longer(cols=RE1_mean:RE3_p97.5, names_to=c("Projection", "Metric"), names_pattern="RE?(.*)_(.*)", values_to="Value") %>%
  relocate(Projection, Metric, Value) %>%
  pivot_wider(names_from=Metric,values_from=Value ) %>%
  relocate(Projection, median, var, p2.5, p97.5) %>%
  mutate(Proj.case=case_when(
    Projection=='1' ~ "Continue",
    Projection=='2' ~ "Avg5yr",
    Projection=='3' ~ "ObsEcov"
  ))


### catch ====
re.catch.sum3yr <- re.catch.tib %>%
  filter(Year<4) %>%
  group_by( R_sig,  Ecov_how, EM, mod.match ) %>%
  summarise(RE1_mean=mean(RE1), RE1_median=median(RE1), RE1_var=var(RE1), across(RE1, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            RE2_mean=mean(RE2), RE2_median=median(RE2), RE2_var=var(RE2), across(RE2, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            RE3_mean=mean(RE3), RE3_median=median(RE3), RE3_var=var(RE3), across(RE3, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
  )%>%
  ungroup() %>%
  left_join(em_tib) %>%
  pivot_longer(cols=RE1_mean:RE3_p97.5, names_to=c("Projection", "Metric"), names_pattern="RE?(.*)_(.*)", values_to="Value") %>%
  relocate(Projection, Metric, Value) %>%
  pivot_wider(names_from=Metric,values_from=Value ) %>%
  relocate(Projection, median, var, p2.5, p97.5) %>%
  mutate(Proj.case=case_when(
    Projection=='1' ~ "Continue",
    Projection=='2' ~ "Avg5yr",
    Projection=='3' ~ "ObsEcov"
  ))



# %%%%%%%%%%%%%%%%%%%
## summarize  FIRST 3 Years (log_sd) ====
### recruitment ====
sd.rec.sum3yr <- sd.rec.tib %>%
  filter(Year<4) %>%
  group_by( R_sig,  Ecov_how, EM, mod.match ) %>%
  summarise(SD1_mean=mean(SD1), SD1_median=median(SD1), SD1_var=var(SD1), across(SD1, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            SD2_mean=mean(SD2), SD2_median=median(SD2), SD2_var=var(SD2), across(SD2, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            SD3_mean=mean(SD3), SD3_median=median(SD3), SD3_var=var(SD3), across(SD3, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
  )%>%
  ungroup() %>%
  left_join(em_tib) %>%
  pivot_longer(cols=SD1_mean:SD3_p97.5, names_to=c("Projection", "Metric"), names_pattern="SD?(.*)_(.*)", values_to="Value") %>%
  relocate(Projection, Metric, Value) %>%
  pivot_wider(names_from=Metric,values_from=Value ) %>%
  relocate(Projection, median, var, p2.5, p97.5) %>%
  mutate(Proj.case=case_when(
    Projection=='1' ~ "Continue",
    Projection=='2' ~ "Avg5yr",
    Projection=='3' ~ "ObsEcov"
  ))

### ssb ====
sd.ssb.sum3yr <- sd.ssb.tib %>%
  filter(Year<4) %>%
  group_by( R_sig,  Ecov_how, EM, mod.match ) %>%
  summarise(SD1_mean=mean(SD1), SD1_median=median(SD1), SD1_var=var(SD1), across(SD1, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            SD2_mean=mean(SD2), SD2_median=median(SD2), SD2_var=var(SD2), across(SD2, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            SD3_mean=mean(SD3), SD3_median=median(SD3), SD3_var=var(SD3), across(SD3, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
  )%>%
  ungroup() %>%
  left_join(em_tib) %>%
  pivot_longer(cols=SD1_mean:SD3_p97.5, names_to=c("Projection", "Metric"), names_pattern="SD?(.*)_(.*)", values_to="Value") %>%
  relocate(Projection, Metric, Value) %>%
  pivot_wider(names_from=Metric,values_from=Value ) %>%
  relocate(Projection, median, var, p2.5, p97.5) %>%
  mutate(Proj.case=case_when(
    Projection=='1' ~ "Continue",
    Projection=='2' ~ "Avg5yr",
    Projection=='3' ~ "ObsEcov"
  ))


### catch ====
sd.catch.sum3yr <- sd.catch.tib %>%
  filter(Year<4) %>%
  group_by( R_sig,  Ecov_how, EM, mod.match ) %>%
  summarise(SD1_mean=mean(SD1), SD1_median=median(SD1), SD1_var=var(SD1), across(SD1, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            SD2_mean=mean(SD2), SD2_median=median(SD2), SD2_var=var(SD2), across(SD2, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
            
            SD3_mean=mean(SD3), SD3_median=median(SD3), SD3_var=var(SD3), across(SD3, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                           p97.5=~quantile(.,probs=0.975)) ),
  )%>%
  ungroup() %>%
  left_join(em_tib) %>%
  pivot_longer(cols=SD1_mean:SD3_p97.5, names_to=c("Projection", "Metric"), names_pattern="SD?(.*)_(.*)", values_to="Value") %>%
  relocate(Projection, Metric, Value) %>%
  pivot_wider(names_from=Metric,values_from=Value ) %>%
  relocate(Projection, median, var, p2.5, p97.5) %>%
  mutate(Proj.case=case_when(
    Projection=='1' ~ "Continue",
    Projection=='2' ~ "Avg5yr",
    Projection=='3' ~ "ObsEcov"
  ))




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# PLOTS - Just first 3 Years  ====

## RE recruitment ====
re.rec.proj.plot3yr <- ggplot(re.rec.sum3yr, aes(x=EM_mod, y=median, col=as.factor(mod.match) )) +
  facet_grid(R_sig + Ecov_how ~ Proj.case , scales="free_y"    ) +
  geom_point(size=2) +
  geom_hline(yintercept=0.0, col='#111111dd') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.4) )+
  ylab('RE(Recr) (2.5 - 50 - 97.5 percentiles)') +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  ggtitle('RE (across first 3 yrs of projection), 3 different Ecov assumptions in projections')
ggsave(re.rec.proj.plot3yr, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('re.rec.proj.plot3yr', plot.suffix, '.png') ),  height=7, width=12)


## same plot, narrower y-lim
re.rec.proj.plot3yr.ylim <- ggplot(re.rec.sum3yr, aes(x=EM_mod, y=median, col=as.factor(mod.match) )) +
  facet_grid(R_sig + Ecov_how ~ Proj.case , scales="free_y"    ) +
  geom_point(size=2) +
  geom_hline(yintercept=0.0, col='#111111dd') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.4) )+
  ylab('RE(Recr) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-0.25,0.25)) +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  ggtitle('RE (across first 3 yrs of projection), 3 different Ecov assumptions in projections')
ggsave(re.rec.proj.plot3yr.ylim, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('re.rec.proj.plot3yr.ylim', plot.suffix, '.png') ),  height=7, width=12)


## SD recruitment 
## RE recruitment ====
sd.rec.proj.plot3yr <- ggplot(sd.rec.sum3yr, aes(x=EM_mod, y=median, col=as.factor(mod.match) )) +
  facet_grid(R_sig + Ecov_how ~ Proj.case , scales="free_y"    ) +
  geom_point(size=2) +
  geom_hline(yintercept=0.0, col='#111111dd') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.4) )+
  ylab('SD(Recr) (2.5 - 50 - 97.5 percentiles)') +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  ggtitle('CV (across first 3 yrs of projection), 3 different Ecov assumptions in projections')
ggsave(sd.rec.proj.plot3yr, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('sd.rec.proj.plot3yr', plot.suffix, '.png') ),  height=7, width=12)





## RE ssb ====
re.ssb.proj.plot3yr <- ggplot(re.ssb.sum3yr, aes(x=EM_mod, y=median, col=as.factor(mod.match) )) +
  facet_grid(R_sig + Ecov_how ~ Proj.case , scales="free_y"    ) +
  geom_point(size=2) +
  geom_hline(yintercept=0.0, col='#111111dd') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.4) )+
  ylab('RE(SSB) (2.5 - 50 - 97.5 percentiles)') +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  ggtitle('RE (across first 3 yrs of projection), 3 different Ecov assumptions in projections')
ggsave(re.ssb.proj.plot3yr, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('re.ssb.proj.plot3yr', plot.suffix, '.png') ),  height=7, width=12)


## same plot, narrower y-lim
re.ssb.proj.plot3yr.ylim <- ggplot(re.ssb.sum3yr, aes(x=EM_mod, y=median, col=as.factor(mod.match) )) +
  facet_grid(R_sig + Ecov_how ~ Proj.case , scales="free_y"    ) +
  geom_point(size=2) +
  geom_hline(yintercept=0.0, col='#111111dd') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.4) )+
  ylab('RE(SSB) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-0.25,0.25)) +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  ggtitle('RE (across first 3 yrs of projection), 3 different Ecov assumptions in projections')
ggsave(re.ssb.proj.plot3yr.ylim, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('re.ssb.proj.plot3yr.ylim', plot.suffix, '.png') ),  height=7, width=12)



## SD ssb ====
sd.ssb.proj.plot3yr <- ggplot(sd.ssb.sum3yr, aes(x=EM_mod, y=median, col=as.factor(mod.match) )) +
  facet_grid(R_sig + Ecov_how ~ Proj.case , scales="free_y"    ) +
  geom_point(size=2) +
  geom_hline(yintercept=0.0, col='#111111dd') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.4) )+
  ylab('CV(SSB) (2.5 - 50 - 97.5 percentiles)') +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  ggtitle('CV (across first 3 yrs of projection), 3 different Ecov assumptions in projections')
ggsave(sd.ssb.proj.plot3yr, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('sd.ssb.proj.plot3yr', plot.suffix, '.png') ),  height=7, width=12)





## RE catch ====
re.catch.proj.plot3yr <- ggplot(re.catch.sum3yr, aes(x=EM_mod, y=median, col=as.factor(mod.match) )) +
  facet_grid(R_sig + Ecov_how ~ Proj.case , scales="free_y"    ) +
  geom_point(size=2) +
  geom_hline(yintercept=0.0, col='#111111dd') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.4) )+
  ylab('RE(Catch) (2.5 - 50 - 97.5 percentiles)') +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  ggtitle('RE (acrossfirst 3 yrs of projection), 3 different Ecov assumptions in projections')
ggsave(re.catch.proj.plot3yr, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('re.catch.proj.plot3yr', plot.suffix, '.png') ),  height=7, width=12)


## same plot, narrower y-lim
re.catch.proj.plot3yr.ylim <- ggplot(re.catch.sum3yr, aes(x=EM_mod, y=median, col=as.factor(mod.match) )) +
  facet_grid(R_sig + Ecov_how ~ Proj.case , scales="free_y"    ) +
  geom_point(size=2) +
  geom_hline(yintercept=0.0, col='#111111dd') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.4) )+
  ylab('RE(Catch) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-0.25,0.25)) +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  ggtitle('RE (across first 3 yrs of projection), 3 different Ecov assumptions in projections')
ggsave(re.catch.proj.plot3yr.ylim, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('re.catch.proj.plot3yr.ylim', plot.suffix, '.png') ),  height=7, width=12)


## SD catch ====
sd.catch.proj.plot3yr <- ggplot(sd.catch.sum3yr, aes(x=EM_mod, y=median, col=as.factor(mod.match) )) +
  facet_grid(R_sig + Ecov_how ~ Proj.case , scales="free_y"    ) +
  geom_point(size=2) +
  geom_hline(yintercept=0.0, col='#111111dd') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.4) )+
  ylab('CV(Catch) (2.5 - 50 - 97.5 percentiles)') +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  ggtitle('CV (across first 3 yrs of projection), 3 different Ecov assumptions in projections')
ggsave(sd.catch.proj.plot3yr, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('sd.catch.proj.plot3yr', plot.suffix, '.png') ),  height=7, width=12)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# not updated below here

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#=================== Look at some runs ====
dir <- file.path(here::here(), "Ecov_study","recruitment_functions", "results_pile2")

sim=2
om.set <- c(65, 67, 105, 107, 114, 115)

proj.out <- c()

for (iom in 1:length(om.set))  {
  
  for (jem in 1:6) {
    test <- readRDS(file.path(dir, paste0("om", om.set[iom], '/','sim',sim,'_','em', jem,'.RDS') ) )
    
    rec.test <- cbind(OM=rep(om.set[iom], 40), EM=rep(jem, 40), Sim=rep(sim,40), Par=rep('Rec',40),
      Year=seq(1,40), Cont=test$proj$cont.ecov$rep$NAA[,1], Avg=test$proj$avg.ecov$rep$NAA[,1],
                      Use=test$proj$use.ecov$rep$NAA[,1], Truth=test$truth$NAA[,1])

    ssb.test <- cbind(OM=rep(om.set[iom], 40), EM=rep(jem, 40), Sim=rep(sim,40), Par=rep('SSB',40),
                      Year=seq(1,40), Cont=test$proj$cont.ecov$rep$SSB, Avg=test$proj$avg.ecov$rep$SSB,
                      Use=test$proj$use.ecov$rep$SSB, Truth=test$truth$SSB)
    
    catch.test <- cbind(OM=rep(om.set[iom], 40), EM=rep(jem, 40), Sim=rep(sim,40), Par=rep('Catch',40),
                      Year=seq(1,40), Cont=test$proj$cont.ecov$rep$pred_catch, Avg=test$proj$avg.ecov$rep$pred_catch,
                      Use=test$proj$use.ecov$rep$pred_catch, Truth=test$truth$pred_catch)
    
    proj.out <- rbind(proj.out, rec.test, ssb.test, catch.test)
    
  }
  
}

df.oms2 <- as_tibble(cbind(OM=seq(1,256), df.oms)) 
df.ems2 <- cbind(EM=seq(1,6), df.ems)
colnames(df.ems2) <- c("EM", "EM_ecov_how", "EM_r_mod")
em_tib <- as_tibble(df.ems2) %>%
  mutate(SR=ifelse(EM_r_mod==2, 'Mean', 'BH')) %>%
  mutate(EM_mod = paste0(SR, "_", EM_ecov_how))


proj.tib <- as_tibble(proj.out) %>%
  mutate_at(c('OM', 'EM', 'Sim', 'Year'), as.integer) %>%
  mutate_at(c('Cont', 'Avg', 'Use', 'Truth'), as.numeric) %>%
  pivot_longer(cols=c('Cont', 'Avg', 'Use'), names_to='Projection', values_to='Value') %>%
  left_join(df.oms2)  %>%
  left_join(em_tib)

proj.tib$R_sig <- factor(proj.tib$R_sig,labels = c("Rsig_0.1","Rsig_1.0"))
proj.tib$Ecov_effect <- factor(proj.tib$Ecov_effect,labels=c("Ecov_L","Ecov_H"))
proj.tib$Ecov_how    <- factor(proj.tib$Ecov_how,labels=c("0", "1","2","4"))
#proj.tib$NAA_cor     <- factor(proj.tib$NAA_cor,labels=c("L","H"))
proj.tib$Fhist       <- factor(proj.tib$Fhist,labels=c("H-MSY","MSY") ) 
proj.tib$Ecov_re_cor <- factor(proj.tib$Ecov_re_cor,labels=c("L","H"))
proj.tib$NAA_cor     <- factor(proj.tib$NAA_cor,labels=c("Rcor_L","Rcor_H"))
proj.tib$Fhist       <- factor(proj.tib$Fhist,labels=c("H-MSY") ) 
proj.tib$Ecov_re_cor <- factor(proj.tib$Ecov_re_cor,labels=c("Ecov_REcor_L","Ecov_REcor_H"))
proj.tib$obs_error <- factor(proj.tib$obs_error,labels=c("ObsErr_H","ObsErr_L"))



rec.proj_lowSigR <- ggplot(proj.tib[proj.tib$Par=='Rec' & proj.tib$R_sig=="Rsig_0.1",], aes(x=Year, y=Value, col=Projection, shape=as.factor(EM) )) +
  facet_grid(EM + NAA_cor ~ Ecov_effect + Ecov_re_cor    ) +
  geom_line() +
  geom_point() +
  geom_line(data=proj.tib[proj.tib$Par=='Rec' & proj.tib$R_sig=="Rsig_0.1",], aes(x=Year, y=Truth), col='black') +
  ggtitle('Rsigma=0.1') +
  ylab('Recr') +
  geom_vline(xintercept = 31, linetype='dashed', col='grey45', linewidth=1.5)  +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  # scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  guides(shape=guide_legend(title='EM')) 
#  ggtitle('RE all years')
ggsave(rec.proj_lowSigR, filename=file.path(here(),'Ecov_study','recruitment_functions','plots','rec.proj_lowSigR.plot.png'),  height=7, width=12)

rec.proj_hiSigR <- ggplot(proj.tib[proj.tib$Par=='Rec' & proj.tib$R_sig=="Rsig_1.0",], aes(x=Year, y=Value, col=Projection, shape=as.factor(EM) )) +
  facet_grid(EM + NAA_cor ~ Ecov_effect + Ecov_re_cor    ) +
  geom_line() +
  geom_point() +
  geom_line(data=proj.tib[proj.tib$Par=='Rec' & proj.tib$R_sig=="Rsig_1.0",], aes(x=Year, y=Truth), col='black') +
  ggtitle('Rsigma=1.0') +
  ylab('Recr') +
  geom_vline(xintercept = 31, linetype='dashed', col='grey45', linewidth=1.5)  +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  # scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  guides(shape=guide_legend(title='EM')) 
#  ggtitle('RE all years')
ggsave(rec.proj_hiSigR, filename=file.path(here(),'Ecov_study','recruitment_functions','plots','rec.proj_hiSigR.plot.png'),  height=7, width=12)




ssb.proj_lowSigR <- ggplot(proj.tib[proj.tib$Par=='SSB' & proj.tib$R_sig=="Rsig_0.1",], aes(x=Year, y=Value, col=Projection, shape=as.factor(EM) )) +
  facet_grid(EM + NAA_cor ~ Ecov_effect + Ecov_re_cor    ) +
  geom_line() +
  geom_point() +
  geom_line(data=proj.tib[proj.tib$Par=='SSB' & proj.tib$R_sig=="Rsig_0.1",], aes(x=Year, y=Truth), col='black') +
  ggtitle('Rsigma=0.1') +
  ylab('SSB') +
  geom_vline(xintercept = 31, linetype='dashed', col='grey45', linewidth=1.5)  +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  # scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  guides(shape=guide_legend(title='EM')) 
#  ggtitle('RE all years')
ggsave(ssb.proj_lowSigR, filename=file.path(here(),'Ecov_study','recruitment_functions','plots','ssb.proj_lowSigR.plot.png'),  height=7, width=12)

ssb.proj_hiSigR <- ggplot(proj.tib[proj.tib$Par=='SSB' & proj.tib$R_sig=="Rsig_1.0",], aes(x=Year, y=Value, col=Projection, shape=as.factor(EM) )) +
  facet_grid(EM + NAA_cor ~ Ecov_effect + Ecov_re_cor    ) +
  geom_line() +
  geom_point() +
  geom_line(data=proj.tib[proj.tib$Par=='SSB' & proj.tib$R_sig=="Rsig_1.0",], aes(x=Year, y=Truth), col='black') +
  ggtitle('Rsigma=1.0') +
  ylab('Recr') +
  geom_vline(xintercept = 31, linetype='dashed', col='grey45', linewidth=1.5)  +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  # scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  guides(shape=guide_legend(title='EM')) 
#  ggtitle('RE all years')
ggsave(ssb.proj_hiSigR, filename=file.path(here(),'Ecov_study','recruitment_functions','plots','ssb.proj_hiSigR.plot.png'),  height=7, width=12)



catch.proj_lowSigR <- ggplot(proj.tib[proj.tib$Par=='Catch' & proj.tib$R_sig=="Rsig_0.1",], aes(x=Year, y=Value, col=Projection, shape=as.factor(EM) )) +
  facet_grid(EM + NAA_cor ~ Ecov_effect + Ecov_re_cor    ) +
  geom_line() +
  geom_point() +
  geom_line(data=proj.tib[proj.tib$Par=='Catch' & proj.tib$R_sig=="Rsig_0.1",], 
            aes(x=Year, y=Truth), col='black') +
  ggtitle('Rsigma=0.1') +
  ylab('Catch') +
  geom_vline(xintercept = 31, linetype='dashed', col='grey45', linewidth=1.5)  +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  # scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  guides(shape=guide_legend(title='EM')) 
#  ggtitle('RE all years')
ggsave(catch.proj_lowSigR, filename=file.path(here(),'Ecov_study','recruitment_functions','plots','catch.proj_lowSigR.plot.png'),  height=7, width=12)

catch.proj_hiSigR <- ggplot(proj.tib[proj.tib$Par=='Catch' & proj.tib$R_sig=="Rsig_1.0",], aes(x=Year, y=Value, col=Projection, shape=as.factor(EM) )) +
  facet_grid(EM + NAA_cor ~ Ecov_effect + Ecov_re_cor    ) +
  geom_line() +
  geom_point() +
  geom_line(data=proj.tib[proj.tib$Par=='Catch' & proj.tib$R_sig=="Rsig_1.0",], aes(x=Year, y=Truth), col='black') +
  ggtitle('Rsigma=1.0') +
  ylab('Catch') +
  geom_vline(xintercept = 31, linetype='dashed', col='grey45', linewidth=1.5)  +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  # scale_color_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(col=guide_legend(title=NULL)) +
  guides(shape=guide_legend(title='EM')) 
#  ggtitle('RE all years')
ggsave(catch.proj_hiSigR, filename=file.path(here(),'Ecov_study','recruitment_functions','plots','catch.proj_hiSigR.plot.png'),  height=7, width=12)













