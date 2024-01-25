library(here)
library(tidyverse)
library(rpart)


# res.dir <- file.path(here::here(), "Ecov_study","recruitment_functions", "results")
res.path <- 'E:/results_beta_fix'  # directory where simulation runs are (beta unstandardized)
res.dir <- 'results_beta_fix'   # 'results'     'results_beta_fix'   # results folder where AIC dataframes are
plot.dir <- 'plots_beta_fix'    # 'plots_lizruns'  'plots_beta_fix'   
plot.suffix <- '_beta_fix'      # '_beta_fix'   '' 

## specify bad.grad.label and bad.se.value (these are the thresholds set in convergence_summaries.R to determine convergence) 
bad.grad.value <- 1E-6 #tim used 1E-6 (the abs of the exponent will be used for output suffix; ex: filename_grad_6.png)
bad.grad.label <- as.numeric(strsplit(as.character(bad.grad.value), split="-")[[1]][2])
bad.se.value <- 100 #tim used 100 (this value will be used for output suffix; ex: filename_se_100.png)


df.oms          <- readRDS(file.path(here::here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS"))
df.ems    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))
#no.conv.runs <- readRDS(file.path(here::here(),'Ecov_study','recruitment_functions','plots_lizruns', "no.conv.runs.RDS") )
conv.runs <- readRDS(file.path(here::here(),'Ecov_study','recruitment_functions',res.dir, paste0("conv.runs_grad_", bad.grad.label, "_SE_", bad.se.value, plot.suffix, ".RDS")  ) )

# conv.runs <- no.conv.runs %>%
#   rename(Sim=sim) %>%
#   mutate(ok.run =(bad.opt+ bad.conv+ bad.sdrep+ bad.grad+ bad.se.big) ) %>%
#   replace_na(list(ok.run=1)) %>%
#   relocate(OM, EM, Sim, ok.run, bad.opt, bad.conv, bad.sdrep, bad.grad, bad.se.big) %>%
#   filter(ok.run==0)  # drop unconverged runs

n_oms <- nrow(df.oms)
n_ems <- nrow(df.ems)
n_sims <- 100
nyears <- 40

#n.conv.runs <- nrow(conv.runs)


t1 <- Sys.time()
re.recr <- lapply(1:n_oms,function(om){  #produces list of length n_oms with list of length n_ems, each with matrix nyearsXn_sims
  print(paste0("om ",om ))
  lapply(1:n_ems, function(em){
    sapply(1:n_sims, function(sim){
      # print(paste0("om ",om," em ",em," sim ",sim))
      re <- rep(NA, nyears)
      dat <- try(readRDS(file.path(res.path, paste0("om", om, '/','sim',sim,'_','em',em,'.RDS') ) ) )#,
      # silent=TRUE)   )
      if(class(dat)!='try-error'){
        #dat$truth$NAA[,1] - dat$fit$rep$NAA[,1]  
        re <- (dat$fit$rep$NAA[,1] - dat$truth$NAA[,1])/ dat$truth$NAA[,1]   
      }
      re
    })
  })
})
t2 <- Sys.time()
t2-t1   #Time difference of 17.66439 mins (lizruns); Time difference of 2.353031 hours (gregruns)

saveRDS(re.recr,file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0('re.recr', plot.suffix, '.RDS') ) )
re.recr <- readRDS(file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0('re.recr', plot.suffix, '.RDS') ) )

t1 <- Sys.time()
re.ssb <- lapply(1:n_oms,function(om){  #produces list of length n_oms with list of length n_ems, each with matrix nyearsXn_sims ***except sometimes the matrix is a list instead of a matrix!
  print(paste0("om ",om ))
  lapply(1:n_ems, function(em){
    sapply(1:n_sims, function(sim){
      #print(paste0("om ",om," em ",em," sim ",sim))
      
      dat <- try(readRDS(file.path(res.path, paste0("om", om, '/','sim',sim,'_','em',em,'.RDS') ) ) )
      if(class(dat)!='try-error'){
        (dat$fit$rep$SSB - dat$truth$SSB)/ dat$truth$SSB   
      }
    })
  })
})
t2 <- Sys.time()
t2-t1  # Time difference of 1.764071 hours
saveRDS(re.ssb, file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0('re.ssb', plot.suffix, '.RDS') ) )
re.ssb <- readRDS( file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0('re.ssb', plot.suffix, '.RDS') ) )

t1 <- Sys.time()
re.fbar <- lapply(1:n_oms,function(om){  #produces list of length n_oms with list of length n_ems, each with matrix nyearsXn_sims
  print(paste0("om ",om ))
  lapply(1:n_ems, function(em){
    sapply(1:n_sims, function(sim){
      #print(paste0("om ",om," em ",em," sim ",sim))
      
      dat <- try(readRDS(file.path(res.path, paste0("om", om, '/','sim',sim,'_','em',em,'.RDS') ) ) )#,
      # silent=TRUE)   )
      if(class(dat)!='try-error'){
        (dat$fit$rep$Fbar - dat$truth$Fbar)/ dat$truth$Fbar   
      }
    })
  })
})
t2 <- Sys.time()
t2-t1 # Time difference of 1.781398 hours
saveRDS(re.fbar, file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0('re.fbar', plot.suffix, '.RDS') ))
re.fbar <- readRDS(file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0('re.fbar', plot.suffix, '.RDS') ) )


#ggplot(dplyr::bind_rows(ff[[1]], .id="EM"), aes(x=............) )

#summarize convergence by OM#, then drill into factors (did any combo consistently bomb?)
# 
# re.recr.all.yrs.med <- lapply(1:n_oms,function(om) {
#   lapply(1:n_ems, function(em)  {
#     sapply(1:n_sims, function(sim) {
#       med <- NA
#       dat <- try(readRDS(file.path(res.dir, paste0("om", om, '/','sim',sim,'_','em',em,'.RDS') ) ) )
#       if(class(dat)!='try-error'){
#       med <- quantile(re.recr[[om]][[em]][,sim], probs=0.5)
#       }
#       med
#     } )
#   } )
# }  )

nyears<-dim(re.recr[[1]][[1]])[1]

t1 <- Sys.time()
re.recr.all.yrs <- matrix(NA, nrow=(nrow(df.oms)*nrow(df.ems)), ncol=(7+2+ncol(df.oms)))
re.recr.df <- matrix(NA, nrow=n_oms*n_ems*nyears*n_sims, ncol=(5 )   )# OM, EM, Sim, Year, RE, df.om colnames
colnames(re.recr.df) <- c('OM', 'EM', 'Sim', 'Year', 'RE')
re.ssb.df <- matrix(NA, nrow=n_oms*n_ems*nyears*n_sims, ncol=(5 )   )# OM, EM, Sim, Year, RE, df.om colnames
colnames(re.ssb.df) <- c('OM', 'EM', 'Sim', 'Year', 'RE')
re.fbar.df <- matrix(NA, nrow=n_oms*n_ems*nyears*n_sims, ncol=(5 )   )# OM, EM, Sim, Year, RE, df.om colnames
colnames(re.fbar.df) <- c('OM', 'EM', 'Sim', 'Year', 'RE')


re.recr.sims <- matrix(NA, nrow=n_sims*nyears, ncol=1)
re.ssb.sims <- matrix(NA, nrow=n_sims*nyears, ncol=1)
re.fbar.sims  <- matrix(NA, nrow=n_sims*nyears, ncol=1)
colnames(re.recr.all.yrs) <- c('OM', 'EM', paste0(rep('P',7), c('025','05', '25','50','75','95', '975')), colnames(df.oms) )
k <- 1
kdf <- 1
  for(iom in 1:n_oms) {
    print(paste0("om ",iom ))
    for (jem in 1:n_ems){
        kk<-1
      for (ksim in 1:n_sims) {
        dat <- try(readRDS(file.path(res.path, paste0("om", iom, '/','sim',ksim,'_','em',jem,'.RDS') ) ) )
        if(class(dat)!='try-error'){
          if(length(re.recr[[iom]][[jem]])==n_sims*nyears)  re.recr.sims[((ksim-1)*nyears+1):(ksim*nyears),1] <- re.recr[[iom]][[jem]][,ksim]
          if(length(re.ssb[[iom]][[jem]])==n_sims*nyears) re.ssb.sims[((ksim-1)*nyears+1):(ksim*nyears),1] <- re.ssb[[iom]][[jem]][,ksim]
          if(length(re.fbar[[iom]][[jem]])==n_sims*nyears)  re.fbar.sims[((ksim-1)*nyears+1):(ksim*nyears),1] <- re.fbar[[iom]][[jem]][,ksim]
            
           
          
          if(length(re.recr[[iom]][[jem]])==n_sims) {
            if(length(re.recr[[iom]][[jem]][[ksim]])>0) {
              re.recr.sims[((ksim-1)*nyears+1):(ksim*nyears),1] <- re.recr[[iom]][[jem]][[ksim]]
            }
          }
          if(length(re.ssb[[iom]][[jem]])==n_sims) {
            if(length(re.ssb[[iom]][[jem]][[ksim]])>0) {
              re.ssb.sims[((ksim-1)*nyears+1):(ksim*nyears),1] <- re.ssb[[iom]][[jem]][[ksim]]
            }
          }
              if(length(re.fbar[[iom]][[jem]])==n_sims) {
                if(length(re.fbar[[iom]][[jem]][[ksim]])>0) {
              re.fbar.sims[((ksim-1)*nyears+1):(ksim*nyears),1] <- re.fbar[[iom]][[jem]][[ksim]]
                  
            }
          }
          
        }
          kk <- kk+1
          kdf=kdf+1
        
      } #ksim loop
        re.recr.all.yrs[k, 3:9] <- quantile(re.recr.sims, probs=c(0.025,0.05, 0.25,0.5,0.75,0.95, 0.975), na.rm=TRUE )
       re.recr.all.yrs[k, 1] <- iom
       re.recr.all.yrs[k, 2] <- jem 
       #re.recr.all.yrs[k, 10:21] <- df.oms[k,1:12]
       
        re.recr.df [((k-1)*nyears*n_sims+1):(n_sims*nyears*k),5 ] <- re.recr.sims
        re.recr.df[((k-1)*nyears*n_sims+1):(n_sims*nyears*k), 1:4] <- cbind(rep(iom, nyears*n_sims), rep(jem, nyears*n_sims), rep(seq(1,n_sims),each=nyears), rep(seq(1,nyears), n_sims))

        re.ssb.df [((k-1)*nyears*n_sims+1):(n_sims*nyears*k),5 ] <- re.ssb.sims
        re.ssb.df[((k-1)*nyears*n_sims+1):(n_sims*nyears*k), 1:4] <- cbind(rep(iom, nyears*n_sims), rep(jem, nyears*n_sims), rep(seq(1,n_sims),each=nyears), rep(seq(1,nyears), n_sims))

        re.fbar.df [((k-1)*nyears*n_sims+1):(n_sims*nyears*k),5 ] <- re.fbar.sims
        re.fbar.df[((k-1)*nyears*n_sims+1):(n_sims*nyears*k), 1:4] <- cbind(rep(iom, nyears*n_sims), rep(jem, nyears*n_sims), rep(seq(1,n_sims),each=nyears), rep(seq(1,nyears), n_sims))
        
        
         k <- k+1
        re.recr.sims <- matrix(NA, nrow=n_sims*nyears, ncol=1)
        re.ssb.sims <- matrix(NA, nrow=n_sims*nyears, ncol=1)
        re.fbar.sims <- matrix(NA, nrow=n_sims*nyears, ncol=1)
    } #jem loop
  } #iom loop
t2 <- Sys.time()
t2-t1  # Time difference of 1.623663 hours

saveRDS(re.recr.df, file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "re.recr", plot.suffix, ".df.RDS") ) )
saveRDS(re.ssb.df, file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "re.ssb", plot.suffix, ".df.RDS") ) )
saveRDS(re.fbar.df, file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "re.fbar", plot.suffix, ".df.RDS") ) )

# 

####################################################################################
## Read RDS dataframes (if already run above) ====
####################################################################################

re.recr.df <- readRDS( file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "re.recr", plot.suffix, ".df.RDS") ) )
re.ssb.df <- readRDS(  file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "re.ssb", plot.suffix, ".df.RDS") )  )
re.fbar.df <- readRDS( file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "re.fbar", plot.suffix, ".df.RDS") )  )



df.oms2 <- as_tibble(cbind(OM=seq(1,256), df.oms)) 
df.ems2 <- cbind(EM=seq(1,6), df.ems)
colnames(df.ems2) <- c("EM", "EM_ecov_how", "EM_r_mod") 
em_tib <- as_tibble(df.ems2) %>%
  mutate(SR=ifelse(EM_r_mod==2, 'Mean', 'BH')) %>%
  mutate(EM_mod = paste0(SR, "_", EM_ecov_how))


conv.runs.ok <- conv.runs %>%
  mutate(re.ok = paste0(OM, '.', sim, '.', EM)) %>%
  select(re.ok)

re.rec.tib <- as_tibble(re.recr.df) %>%
  left_join(df.oms2) %>%
  left_join(em_tib) %>%
  mutate(mod.match=ifelse(EM_ecov_how==Ecov_how & EM_r_mod==recruit_mod, 1, 0)) %>%
  mutate(re.ok= paste0(OM, '.', Sim, '.', EM) ) %>%
  relocate(OM, EM, Sim, Year, re.ok, RE) %>%
  filter(re.ok %in% conv.runs.ok$re.ok) 

re.ssb.tib <- as_tibble(re.ssb.df) %>%
  left_join(df.oms2) %>%
  left_join(em_tib) %>%
  mutate(mod.match=ifelse(EM_ecov_how==Ecov_how & EM_r_mod==recruit_mod, 1, 0)) %>%
  mutate(re.ok= paste0(OM, '.', Sim, '.', EM) ) %>%
  relocate(OM, EM, Sim, Year, re.ok, RE) %>%
  filter(re.ok %in% conv.runs.ok$re.ok) 



re.fbar.tib <- as_tibble(re.fbar.df) %>%
  left_join(df.oms2) %>%
  left_join(em_tib) %>%
  mutate(mod.match=ifelse(EM_ecov_how==Ecov_how & EM_r_mod==recruit_mod, 1, 0)) %>%
  mutate(re.ok= paste0(OM, '.', Sim, '.', EM) ) %>%
  relocate(OM, EM, Sim, Year, re.ok, RE) %>%
  filter(re.ok %in% conv.runs.ok$re.ok) 



re.rec.tib$R_sig <- factor(re.rec.tib$R_sig,labels = c("Rsig_0.1","Rsig_1.0"))
re.rec.tib$Ecov_effect <- factor(re.rec.tib$Ecov_effect,labels=c("Ecov_L","Ecov_H"))
re.rec.tib$Ecov_how    <- factor(re.rec.tib$Ecov_how,labels=c("Ecov_0", "Ecov_1","Ecov_2","Ecov_4"))
re.rec.tib$NAA_cor     <- factor(re.rec.tib$NAA_cor,labels=c("L","H"))
re.rec.tib$Fhist       <- factor(re.rec.tib$Fhist,labels=c("H-MSY","MSY") ) 
re.rec.tib$Ecov_re_cor <- factor(re.rec.tib$Ecov_re_cor,labels=c("L","H"))


re.ssb.tib$R_sig <- factor(re.ssb.tib$R_sig,labels = c("Rsig_0.1","Rsig_1.0"))
re.ssb.tib$Ecov_effect <- factor(re.ssb.tib$Ecov_effect,labels=c("Ecov_L","Ecov_H"))
re.ssb.tib$Ecov_how    <- factor(re.ssb.tib$Ecov_how,labels=c("Ecov_0", "Ecov_1","Ecov_2","Ecov_4"))
re.ssb.tib$NAA_cor     <- factor(re.ssb.tib$NAA_cor,labels=c("L","H"))
re.ssb.tib$Fhist       <- factor(re.ssb.tib$Fhist,labels=c("H-MSY","MSY") ) 
re.ssb.tib$Ecov_re_cor <- factor(re.ssb.tib$Ecov_re_cor,labels=c("L","H"))


re.fbar.tib$R_sig <- factor(re.fbar.tib$R_sig,labels = c("Rsig_0.1","Rsig_1.0"))
re.fbar.tib$Ecov_effect <- factor(re.fbar.tib$Ecov_effect,labels=c("Ecov_L","Ecov_H"))
re.fbar.tib$Ecov_how    <- factor(re.fbar.tib$Ecov_how,labels=c("Ecov_0", "Ecov_1","Ecov_2","Ecov_4"))
re.fbar.tib$NAA_cor     <- factor(re.fbar.tib$NAA_cor,labels=c("L","H"))
re.fbar.tib$Fhist       <- factor(re.fbar.tib$Fhist,labels=c("H-MSY","MSY") ) 
re.fbar.tib$Ecov_re_cor <- factor(re.fbar.tib$Ecov_re_cor,labels=c("L","H"))


####################################################################################
#  regression trees to find nodes for bias  ====
####################################################################################

# recr ====
  
rf_recr_all.yrs <- rpart(RE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how , data=re.rec.tib, control=rpart.control(cp=0.01))
imp.var <- rf_recr_all.yrs$frame[rf_recr_all.yrs$frame$var != '<leaf>',]
nodes_recr_all.yrs <- unique(imp.var[,1])
# nothing


rf_recr_last10.yrs <- rpart(RE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how , data=re.rec.tib[re.rec.tib$Year>30,], control=rpart.control(cp=0.01))
imp.var <- rf_recr_last10.yrs$frame[rf_recr_last10.yrs$frame$var != '<leaf>',]
nodes_recr_last10.yrs <- unique(imp.var[,1])
# "R_sig"

rf_recr_final.yrs <- rpart(RE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how , data=re.rec.tib[re.rec.tib$Year>39,], control=rpart.control(cp=0.01))
imp.var <- rf_recr_final.yrs$frame[rf_recr_final.yrs$frame$var != '<leaf>',]
nodes_recr_final.yrs <- unique(imp.var[,1])
# nothing


# ssb ====

rf_ssb_all.yrs <- rpart(RE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how , data=re.rec.tib, control=rpart.control(cp=0.01))
imp.var <- rf_ssb_all.yrs$frame[rf_ssb_all.yrs$frame$var != '<leaf>',]
nodes_ssb_all.yrs <- unique(imp.var[,1])
# nothing


rf_ssb_last10.yrs <- rpart(RE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how , data=re.rec.tib[re.rec.tib$Year>30,], control=rpart.control(cp=0.01))
imp.var <- rf_ssb_last10.yrs$frame[rf_ssb_last10.yrs$frame$var != '<leaf>',]
nodes_ssb_last10.yrs <- unique(imp.var[,1])
# "R_sig"

rf_ssb_final.yrs <- rpart(RE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how , data=re.rec.tib[re.rec.tib$Year>39,], control=rpart.control(cp=0.01))
imp.var <- rf_ssb_final.yrs$frame[rf_ssb_final.yrs$frame$var != '<leaf>',]
nodes_ssb_final.yrs <- unique(imp.var[,1])
# nothing


# fbar ====

rf_fbar_all.yrs <- rpart(RE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how , data=re.rec.tib, control=rpart.control(cp=0.01))
imp.var <- rf_fbar_all.yrs$frame[rf_fbar_all.yrs$frame$var != '<leaf>',]
nodes_fbar_all.yrs <- unique(imp.var[,1])
# nothing


rf_fbar_last10.yrs <- rpart(RE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how , data=re.rec.tib[re.rec.tib$Year>30,], control=rpart.control(cp=0.01))
imp.var <- rf_fbar_last10.yrs$frame[rf_fbar_last10.yrs$frame$var != '<leaf>',]
nodes_fbar_last10.yrs <- unique(imp.var[,1])
# "R_sig"

rf_fbar_final.yrs <- rpart(RE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how , data=re.rec.tib[re.rec.tib$Year>39,], control=rpart.control(cp=0.01))
imp.var <- rf_fbar_final.yrs$frame[rf_fbar_final.yrs$frame$var != '<leaf>',]
nodes_fbar_final.yrs <- unique(imp.var[,1])
# "R_sig"


#################################################
# summarize  all years  ====

re.rec.sum <- re.rec.tib %>%
  group_by( R_sig, Ecov_how, EM, mod.match ) %>%
  # group_by( Fhist, R_sig, Ecov_effect, Ecov_how, EM, mod.match ) %>%
  summarise(median.re=median(RE), var.re=var(RE), across(RE, list(p2.5=~quantile(.,probs=0.0275 ),
                                                             p97.5=~quantile(.,probs=0.975)) )
            )%>%
  ungroup() %>%
  left_join(em_tib)

re.ssb.sum <- re.ssb.tib %>%
  group_by(  R_sig, Ecov_how, EM, mod.match ) %>%
  # group_by( Fhist, R_sig, Ecov_effect, Ecov_how, EM, mod.match ) %>%
  summarise(median.re=median(RE), var.re=var(RE), across(RE, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                  p97.5=~quantile(.,probs=0.975)) )
  )%>%
  ungroup() %>%
  left_join(em_tib)

re.fbar.sum <- re.fbar.tib %>%
  group_by( R_sig,  Ecov_how, EM, mod.match ) %>%
  # group_by( Fhist, R_sig, Ecov_effect, Ecov_how, EM, mod.match ) %>%
  summarise(median.re=median(RE), var.re=var(RE), across(RE, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                  p97.5=~quantile(.,probs=0.975)) )
  )%>%
  ungroup() %>%
  left_join(em_tib)

#################################################
# summarize  last 10 years  ====


re.rec.sum.last.10 <- re.rec.tib %>%
  filter(Year>30) %>%
  group_by(  R_sig,  Ecov_how, EM, mod.match ) %>%
  # group_by( Fhist, R_sig, Ecov_effect, Ecov_how, EM, mod.match ) %>%
  summarise(median.re=median(RE), var.re=var(RE), across(RE, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                  p97.5=~quantile(.,probs=0.975)) )
  )%>%
  ungroup() %>%
  left_join(em_tib)

re.ssb.sum.last.10 <- re.ssb.tib %>%
  filter(Year>30) %>%
  group_by(  R_sig,  Ecov_how, EM, mod.match ) %>%
  # group_by( Fhist, R_sig, Ecov_effect, Ecov_how, EM, mod.match ) %>%
  summarise(median.re=median(RE), var.re=var(RE), across(RE, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                  p97.5=~quantile(.,probs=0.975)) )
  )%>%
  ungroup() %>%
  left_join(em_tib)

re.fbar.sum.last.10 <- re.fbar.tib %>%
  filter(Year>30) %>%
  group_by(  R_sig,  Ecov_how, EM, mod.match ) %>%
  # group_by( Fhist, R_sig, Ecov_effect, Ecov_how, EM, mod.match ) %>%
  summarise(median.re=median(RE), var.re=var(RE), across(RE, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                  p97.5=~quantile(.,probs=0.975)) )
  )%>%
  ungroup() %>%
  left_join(em_tib)


#################################################
# summarize final year  ====


re.rec.sum.last.yr <- re.rec.tib %>%
  filter(Year==40) %>%
  group_by( R_sig, Ecov_how, EM, mod.match ) %>%
  # group_by( Fhist, R_sig, Ecov_effect, Ecov_how, EM, mod.match ) %>%
  summarise(median.re=median(RE), var.re=var(RE), across(RE, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                  p97.5=~quantile(.,probs=0.975)) )
  )%>%
  ungroup() %>%
  left_join(em_tib)


re.ssb.sum.last.yr <- re.ssb.tib %>%
  filter(Year==40) %>%
  group_by(  R_sig,  Ecov_how, EM, mod.match ) %>%
  # group_by( Fhist, R_sig, Ecov_effect, Ecov_how, EM, mod.match ) %>%
  summarise(median.re=median(RE), var.re=var(RE), across(RE, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                  p97.5=~quantile(.,probs=0.975)) )
  )%>%
  ungroup() %>%
  left_join(em_tib)

re.fbar.sum.last.yr <- re.fbar.tib %>%
  filter(Year==40) %>%
  group_by(  R_sig,  Ecov_how, EM, mod.match ) %>%
  # group_by( Fhist, R_sig, Ecov_effect, Ecov_how, EM, mod.match ) %>%
  summarise(median.re=median(RE), var.re=var(RE), across(RE, list(p2.5=~quantile(.,probs=0.0275 ),
                                                                  p97.5=~quantile(.,probs=0.975)) )
  )%>%
  ungroup() %>%
  left_join(em_tib)




#################################################
# RE Plots  ====

re.rec.all.yrs.plot <- ggplot(re.rec.sum, aes(x=EM_mod, y=median.re, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~    R_sig  ) +
  # facet_grid(R_sig+Ecov_how ~ Fhist+Ecov_effect     ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE_p2.5), ymax=(RE_p97.5), width=.5) )+
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
  ggtitle('RE all years')
ggsave(re.rec.all.yrs.plot, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('re.rec.all.yrs.plot', plot.suffix, '.png') ),  height=7, width=12)


re.ssb.all.yrs.plot <- ggplot(re.ssb.sum, aes(x=EM_mod, y=median.re, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~    R_sig  ) +
  # facet_grid(R_sig+Ecov_how ~ Fhist+Ecov_effect     ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE_p2.5), ymax=(RE_p97.5), width=.5) )+
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
  ggtitle('RE all years')
ggsave(re.ssb.all.yrs.plot, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('re.ssb.all.yrs.plot', plot.suffix, '.png') ),  height=7, width=12)


re.fbar.all.yrs.plot <- ggplot(re.fbar.sum, aes(x=EM_mod, y=median.re, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~    R_sig  ) +
  # facet_grid(R_sig+Ecov_how ~ Fhist+Ecov_effect     ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE_p2.5), ymax=(RE_p97.5), width=.5) )+
  ylab('RE(Fbar) (2.5 - 50 - 97.5 percentiles)') +
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
  ggtitle('RE all years')
ggsave(re.fbar.all.yrs.plot, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('re.fbar.all.yrs.plot', plot.suffix,'.png') ),  height=7, width=12)

# last 10 years ====
re.rec.last.10.plot <- ggplot(re.rec.sum.last.10, aes(x=EM_mod, y=median.re, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~    R_sig  ) +
  # facet_grid(R_sig+Ecov_how ~ Fhist+Ecov_effect     ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE_p2.5), ymax=(RE_p97.5), width=.5) )+
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
  ggtitle('RE last 10 years')
ggsave(re.rec.last.10.plot, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('re.rec.last.10.plot', plot.suffix, '.png') ),  height=7, width=12)



re.ssb.last.10.plot <- ggplot(re.ssb.sum.last.10, aes(x=EM_mod, y=median.re, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~    R_sig  ) +
  # facet_grid(R_sig+Ecov_how ~ Fhist+Ecov_effect     ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE_p2.5), ymax=(RE_p97.5), width=.5) )+
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
  ggtitle('RE last 10 years')
ggsave(re.ssb.last.10.plot, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('re.ssb.last.10.plot', plot.suffix, '.png')),  height=7, width=12)



re.fbar.last.10.plot <- ggplot(re.fbar.sum.last.10, aes(x=EM_mod, y=median.re, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~    R_sig  ) +
  # facet_grid(R_sig+Ecov_how ~ Fhist+Ecov_effect     ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE_p2.5), ymax=(RE_p97.5), width=.5) )+
  ylab('RE(Fbar) (2.5 - 50 - 97.5 percentiles)') +
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
  ggtitle('RE last 10 years')
ggsave(re.fbar.last.10.plot, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir, paste0('re.fbar.last.10.plot', plot.suffix, '.png') ),  height=7, width=12)



# last year ====
re.rec.last.yr.plot <- ggplot(re.rec.sum.last.yr, aes(x=EM_mod, y=median.re, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~    R_sig  ) +
  # facet_grid(R_sig+Ecov_how ~ Fhist+Ecov_effect     ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE_p2.5), ymax=(RE_p97.5), width=.5) )+
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
  ggtitle('RE final year')
ggsave(re.rec.last.yr.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('re.rec.last.yr.plot', plot.suffix, '.png') ),  height=7, width=12)


re.ssb.last.yr.plot <- ggplot(re.ssb.sum.last.yr, aes(x=EM_mod, y=median.re, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~    R_sig  ) +
  # facet_grid(R_sig+Ecov_how ~ Fhist+Ecov_effect     ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE_p2.5), ymax=(RE_p97.5), width=.5) )+
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
  ggtitle('RE final year')
ggsave(re.ssb.last.yr.plot, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('re.ssb.last.yr.plot', plot.suffix, '.png') ),  height=7, width=12)


re.fbar.last.yr.plot <- ggplot(re.fbar.sum.last.yr, aes(x=EM_mod, y=median.re, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~    R_sig  ) +
  # facet_grid(R_sig+Ecov_how ~ Fhist+Ecov_effect     ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE_p2.5), ymax=(RE_p97.5), width=.5) )+
  ylab('RE(Fbar) (2.5 - 50 - 97.5 percentiles)') +
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
  ggtitle('RE final year')
ggsave(re.fbar.last.yr.plot, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir, paste0('re.fbar.last.yr.plot', plot.suffix, '.png') ),  height=7, width=12)






############################################################################################################

# rho transformation
# all rho pars use this -1 + 2/(1 + exp(-k*x)) where k = 1 for Ecov and 2 for everything else

#OM w/ SRR = ; OM w/o SRR = 
#facet is EM?
# compare EM=OM vs EM=min(AIC) -- bias?


# estimability of sigmaR and rho_y?

# ref pts? (F40 and/or Fmsy; SSB40 and/or SSBmsy)



