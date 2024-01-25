library(here)
library(tidyverse)
library(wham)
library(rpart)
library(rpart.plot)


source(file.path(here::here(), "Ecov_study","recruitment_functions", "code", "analysis", "functions_for_analysis.R" ) )


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

n_oms <- nrow(df.oms)
n_ems <- nrow(df.ems)
n_sims <- 100
nyears <- 40

n.conv.runs <- nrow(conv.runs)
ten.pct <- round(n.conv.runs/10, 0)


# folders <- list.files(dir, pattern='om')
# folders <- folders[-(folders %in% 'om1')]
# om.set <- sapply(folders, function(x) substr(x, 3, (nchar(x))) )
# om.set <- sort(as.numeric(om.set))
# n_om.set <- length(om.set)
# n_om.set*n_ems*n_sims
# nyrs.proj <- 10
# 
# keep.good <-as_tibble(expand.grid(OM = om.set,
#                          EM=seq(1,n_ems), 
#                          Sim=seq(1,n_sims)  ,
#                                     stringsAsFactors = FALSE) ) %>%
#   left_join(conv.runs)  %>%
#   drop_na(bad.sum)   #this drops the rows that are not in conv.runs
# 
# bad.iters <-as_tibble(expand.grid(OM = om.set,
#                                   EM=seq(1,n_ems), 
#                                   Sim=seq(1,n_sims)  ,
#                                   stringsAsFactors = FALSE) ) %>%
#   left_join(conv.runs)  %>%
#   filter(is.na(bad.sum)  )  #this keeps only bad runs



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# get mohn's rho for SSB, Fbar, Recr + NAA ====
t1 <- Sys.time()
rho.peels <- lapply(1:n_oms,function(om){ 
  print(paste0("om ",om))
  lapply(1:n_ems, function(em){
    sapply(1:n_sims, function(sim){
      #print(paste0("om ",om," em ",em," sim ",sim))
      mrho <- rep(NA, 12)  # rho for all ages (10) + SSB and F
      if (file.exists(file.path(res.path, paste0("om", om, '/','sim',sim,'_','em',em,'.RDS') ))) {
        dat <- try(readRDS(file.path(res.path, paste0("om", om, '/','sim',sim,'_','em',em,'.RDS') ) ) )#,
        # silent=TRUE)   )
        if(class(dat)!='try-error'){
          mrho <- c(mohns_rho_set_peel(model=dat, npeels=7, ny=40, na=10), OM=om, EM=em, SIM=sim) 
          mrho
         
        } #try-error
      } #file-exists
      
    } #sim function
 
    ) #sapply across sims
      
  }) # lapply across em
          
  
}) # lapply across om
t2 <- Sys.time()
t2-t1
# Time difference of 1.6924 hours (gregs runs)
saveRDS(rho.peels, file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0('rho.peels', plot.suffix, '.RDS') ) )
# rho.peels <- readRDS(file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0('rho.peels', plot.suffix, '.RDS') ) )

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# get rho for random effects on recruitment ====

t1 <- Sys.time()
randeff.peels <- lapply(1:n_oms,function(om){ 
  print(paste0("om ",om))
  lapply(1:n_ems, function(em){
    sapply(1:n_sims, function(sim){
      mrho <- NA  # rho for recruitment
      if (file.exists(file.path(res.path, paste0("om", om, '/','sim',sim,'_','em',em,'.RDS') ))) {
        dat <- try(readRDS(file.path(res.path, paste0("om", om, '/','sim',sim,'_','em',em,'.RDS') ) ) )#,
        # silent=TRUE)   )
        if(class(dat)!='try-error'){
          mrho <- c(mohns_rho_randeff_peel(model=dat, npeels=7, ny=40), OM=om, EM=em, SIM=sim) 
          mrho
          
        } #try-error
      } #file-exists
      
    } #sim function
    
    ) #sapply across sims
    
  }) # lapply across em
  
  
}) # lapply across om
t2 <- Sys.time()
t2-t1
# Time difference of 1.668574 hours (gregs runs)
saveRDS(randeff.peels, file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0('randeff.peels', plot.suffix, '.RDS') ) )



# df.names <- row.names(rho.peels[[1]][[1]])
# rho.df <-(sapply(rho.peels, function(x)  unlist(x)) )
# rownames(rho.df) <- rep(df.names, nrow(rho.df)/length(df.names) )
# rownames(rho.df) <- c(paste0(df.names,rep())

# rho.tib <- as_tibble(t(rho.df))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# reformat the mohn's rho values for SSB, Fbar and Rec ====
rho.df2 <- c()
t1 <- Sys.time()
k=1
for (iom in 1:n_oms) {
   tmp <- rho.peels[[k]]
   
  for (jem in 1:n_ems) {
       tmp2 <- (tmp[[jem]])
       dummy <- matrix(NA, nrow=15, ncol=n_sims)

      if(mode(tmp2)=="list" & length(tmp2)==10 )  {
         
         for (ksim in 1:n_sims) {
           if(!is.null(tmp2[[ksim]]) )  dummy[ ,ksim] <- tmp2[[ksim]]
         }
         
      }
       if(mode(tmp2)!="list") {
       if( mode(tmp2)=="numeric" & dim(tmp2)[1]==15 )  {
         
             dummy <- tmp2
       }
       }
         
       
       rho.df2 <- rbind(rho.df2, t(dummy))
       
    
  } # jem loop
   k=k+1
}

t2 <- Sys.time()
t2-t1    # Time difference of 14.54442 secs



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# reformat the mohn's rho values for random effects on recruitment ====
rho.randeff.df2 <- c()
tmp<-c()
tmp2<-c()
dummy<-c()
t1 <- Sys.time()
k=1
for (iom in 1:n_oms) {
  tmp <- randeff.peels[[k]]
  
  for (jem in 1:n_ems) {
    tmp2 <- (tmp[[jem]])
    dummy <- matrix(NA, nrow=4, ncol=n_sims)
    
    if(mode(tmp2)=="list" & length(tmp2)==n_sims )  {
      
      for (ksim in 1:n_sims) {
        if(!is.null(tmp2[[ksim]]) )  dummy[ ,ksim] <- tmp2[[ksim]]
      }
      
    }
    if(mode(tmp2)!="list") {
      if( mode(tmp2)=="numeric" & dim(tmp2)[1]==4 )  {
        
        dummy <- tmp2
      }
    }
    
    
    rho.randeff.df2 <- rbind(rho.randeff.df2, t(dummy))
    
    
  } # jem loop
  k=k+1
}

t2 <- Sys.time()
t2-t1    # Time difference of 10.38174 secs


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  analyze mohn's rho for SSB, Fbar, Recruitment ====

df.oms2 <- as_tibble(cbind(OM=seq(1,256), df.oms)) 
df.ems2 <- cbind(EM=seq(1,6), df.ems)
colnames(df.ems2) <- c("EM", "EM_ecov_how", "EM_r_mod")
em_tib <- as_tibble(df.ems2) %>%
  mutate(SR=ifelse(EM_r_mod==2, 'Mean', 'BH')) %>%
  mutate(EM_mod = paste0(SR, "_", EM_ecov_how))

rho.tib <- as_tibble(rho.df2) %>%
  pivot_longer(cols=SSB:N10, names_to='Quantity', values_to="Rho") %>%
  left_join(df.oms2) %>%
  left_join(em_tib) %>%
  drop_na()   %>%
  mutate(mod.match=ifelse(EM_ecov_how==Ecov_how & EM_r_mod==recruit_mod, 1, 0)) %>%
  filter(Quantity %in% c('SSB', 'Fbar', 'R')) %>%
  filter(Rho<10)

saveRDS(rho.tib, file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0('rho.SSB.Fbar.Recr', plot.suffix, '.RDS') ) )

# look at regression trees for explanatory factors ====
library(rpart.plot)

rho.tib$R_sig <- factor(rho.tib$R_sig,labels = c("Rsig_0.1","Rsig_1.0"))
rho.tib$Ecov_effect <- factor(rho.tib$Ecov_effect,labels=c("Ecov_L","Ecov_H"))
rho.tib$Ecov_how    <- factor(rho.tib$Ecov_how,labels=c("Ecov_0", "Ecov_1","Ecov_2","Ecov_4"))
rho.tib$NAA_cor     <- factor(rho.tib$NAA_cor,labels=c("Rcor_L","Rcor_H"))
rho.tib$Fhist       <- factor(rho.tib$Fhist,labels=c("H-MSY","MSY") ) 
rho.tib$Ecov_re_cor <- factor(rho.tib$Ecov_re_cor,labels=c("EcovCor_L","EcovCor_H"))
rho.tib$obs_error   <- factor(rho.tib$obs_error,labels=c("ObsErr_L","ObsErr_H"))
rho.tib$mod.match <- factor(rho.tib$mod.match, labels=c('1', '0'))

# regression tree for mohn's rho ======================================
rf_SSB   <- rpart(Rho   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=rho.tib[rho.tib$Quantity=='SSB',], control=rpart.control(cp=0.01))
imp.var <- rf_SSB$frame[rf_SSB$frame$var != '<leaf>',]
nodes_rf_SSB <- unique(imp.var[,1])
nodes_rf_SSB
# nothing  (including mod.match made no diff)


rf_Fbar   <- rpart(Rho   ~   R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=rho.tib[rho.tib$Quantity=='Fbar',], control=rpart.control(cp=0.01))
imp.var <- rf_Fbar$frame[rf_Fbar$frame$var != '<leaf>',]
nodes_rf_Fbar <- unique(imp.var[,1])
nodes_rf_Fbar
#  "obs_error"   (including mod.match made no diff)


rf_R <- rpart(Rho ~  R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=rho.tib[rho.tib$Quantity=='R',], control=rpart.control(cp=0.01))
imp.var <- rf_R$frame[rf_R$frame$var != '<leaf>',]
nodes_rf_R <- unique(imp.var[,1])
nodes_rf_R
#   "R_sig"     "obs_error"   (including mod.match made no diff)



pdf(file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('reg_tree_Rho', plot.suffix, '.pdf') ),
    height=4,width=7)
par(mfrow=c(1,3), oma=c(5,0,0,0))
prp(rf_SSB,yesno=FALSE,type=4,clip.right.labs=TRUE)
mtext('a) SSB',adj=0,line=2.5)
prp(rf_Fbar,yesno=FALSE,type=4,clip.right.labs=TRUE)
mtext('b) Fbar',adj=0,line=2.5)
title(sub= "Mohn's rho",
      adj=0.5, outer=TRUE)
prp(rf_R,yesno=FALSE,type=4,clip.right.labs=TRUE)
mtext('c) Recr',adj=0,line=2.5)


dev.off()



# now summarize by important factors ====
rho.tib.sum <- rho.tib %>%
  filter(Quantity %in% c('SSB', 'Fbar', 'R' ) ) %>%
  group_by( EM, Ecov_how, mod.match, Quantity, R_sig, Fhist, obs_error) %>%
  summarise(mean.rho=mean(Rho, na.rm=T), var.rho=var(Rho, na.rm=T), med.rho=median(Rho, na.rm=T),
            p2.5=quantile(Rho, prob=0.025), p97.5=quantile(Rho, prob=0.975)) %>%
  ungroup()  %>%
  left_join(em_tib)
  #drop_na() %>%


# now plot distribution summary of rho values  ====

ssb.rho.plot <- ggplot(rho.tib.sum[rho.tib.sum$Quantity=="SSB",], aes(x=EM_mod, y=med.rho, col=as.factor(mod.match)  )) +
  facet_grid(R_sig+Ecov_how ~ Fhist+obs_error     ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0, col='grey35') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.5) )+
  ylab('Mohn`s Rho SSB (2.5 - 50 - 97.5 percentiles)') +
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
  guides(col=guide_legend(title=NULL))
ggsave(ssb.rho.plot, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('ssb.rho', plot.suffix, '.png') ),  height=7, width=14)

recr.rho.plot <- ggplot(rho.tib.sum[rho.tib.sum$Quantity=="R",], aes(x=EM_mod, y=med.rho, col=as.factor(mod.match)  )) +
  facet_grid(R_sig+Ecov_how ~ Fhist+obs_error     ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0, col='grey35') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.5) )+
  ylab('Mohn`s Rho Recr (2.5 - 50 - 97.5 percentiles)') +
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
  guides(col=guide_legend(title=NULL))
ggsave(recr.rho.plot, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('recr.rho', plot.suffix, '.png') ),  height=7, width=14)


fbar.rho.plot <- ggplot(rho.tib.sum[rho.tib.sum$Quantity=="Fbar",], aes(x=EM_mod, y=med.rho, col=as.factor(mod.match)  )) +
  facet_grid(R_sig+Ecov_how ~ Fhist+obs_error     ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0, col='grey35') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.5) )+
  ylab('Mohn`s Rho Fbar (2.5 - 50 - 97.5 percentiles)') +
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
  guides(col=guide_legend(title=NULL))
ggsave(fbar.rho.plot, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('fbar.rho', plot.suffix, '.png') ),  height=7, width=14)






# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#  analyze mohn's rho Recruitment random effects  ====

df.oms2 <- as_tibble(cbind(OM=seq(1,256), df.oms)) 
df.ems2 <- cbind(EM=seq(1,6), df.ems)
colnames(df.ems2) <- c("EM", "EM_ecov_how", "EM_r_mod")
em_tib <- as_tibble(df.ems2) %>%
  mutate(SR=ifelse(EM_r_mod==2, 'Mean', 'BH')) %>%
  mutate(EM_mod = paste0(SR, "_", EM_ecov_how))

rho.randeff.tib <- as_tibble(rho.randeff.df2) %>%
  left_join(df.oms2) %>%
  left_join(em_tib) %>%
  drop_na()   %>%
  mutate(mod.match=ifelse(EM_ecov_how==Ecov_how & EM_r_mod==recruit_mod, 1, 0)) %>%
  filter(abs(R_dev)<10)

saveRDS(rho.randeff.tib, file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0('rho.Recr_devs', plot.suffix, '.RDS') ) )

# look at regression trees for explanatory factors ====
library(rpart.plot)

rho.randeff.tib$R_sig <- factor(rho.randeff.tib$R_sig,labels = c("Rsig_0.1","Rsig_1.0"))
rho.randeff.tib$Ecov_effect <- factor(rho.randeff.tib$Ecov_effect,labels=c("Ecov_L","Ecov_H"))
rho.randeff.tib$Ecov_how    <- factor(rho.randeff.tib$Ecov_how,labels=c("Ecov_0", "Ecov_1","Ecov_2","Ecov_4"))
rho.randeff.tib$NAA_cor     <- factor(rho.randeff.tib$NAA_cor,labels=c("Rcor_L","Rcor_H"))
rho.randeff.tib$Fhist       <- factor(rho.randeff.tib$Fhist,labels=c("H-MSY","MSY") ) 
rho.randeff.tib$Ecov_re_cor <- factor(rho.randeff.tib$Ecov_re_cor,labels=c("EcovCor_L","EcovCor_H"))
rho.randeff.tib$obs_error   <- factor(rho.randeff.tib$obs_error,labels=c("ObsErr_L","ObsErr_H"))
rho.randeff.tib$mod.match <- factor(rho.randeff.tib$mod.match, labels=c('1', '0'))

# regression tree  ======================================

rf_Rdev <- rpart(R_dev ~  R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=rho.randeff.tib, control=rpart.control(cp=0.01))
imp.var <- rf_Rdev$frame[rf_Rdev$frame$var != '<leaf>',]
nodes_rf_Rdev <- unique(imp.var[,1])
nodes_rf_Rdev
#   "R_sig"     (including mod.match made no diff)



pdf(file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('reg_tree_Rdev_retro', plot.suffix, '.pdf') ),
    height=4,width=7)
par(mfrow=c(1,3), oma=c(5,0,0,0))
prp(rf_Rdev,yesno=FALSE,type=4,clip.right.labs=TRUE)
mtext('R_dev',adj=0,line=2.5)
title(sub= "Mohn's rho on recruitment random effects",
      adj=0.5, outer=TRUE)

dev.off()



# now summarize by important factors  ====
rho.randeff.tib.sum <- rho.randeff.tib %>%
  group_by( EM, Ecov_how, mod.match, R_sig, Fhist) %>%
  summarise(mean.rho=mean(R_dev, na.rm=T), var(R_dev, na.rm=T), med.rho=median(R_dev, na.rm=T),
            p2.5=quantile(R_dev, prob=0.025), p97.5=quantile(R_dev, prob=0.975)) %>%
  ungroup()  %>%
  left_join(em_tib)
#drop_na() %>%


# now plot distribution summary of rho values  ====

recr.dev.rho.plot <- ggplot(rho.randeff.tib.sum, aes(x=EM_mod, y=med.rho, col=as.factor(mod.match)  )) +
  facet_grid(Ecov_how ~  R_sig  + Fhist     ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0, col='grey35') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.5) )+
  ylab('Mohn`s Rho for Recruitment random effects (2.5 - 50 - 97.5 percentiles)') +
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
  guides(col=guide_legend(title=NULL))
ggsave(recr.dev.rho.plot, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('recr.dev.rho', plot.suffix, '.png') ),  height=7, width=14)



# same plot but narrower ylim
recr.dev.rho.plot.ylim <- ggplot(rho.randeff.tib.sum, aes(x=EM_mod, y=med.rho, col=as.factor(mod.match)  )) +
  facet_grid(Ecov_how ~  R_sig  + Fhist     ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0, col='grey35') +
  geom_errorbar(aes(ymin=(p2.5), ymax=(p97.5), width=.5) )+
  ylab('Mohn`s Rho for Recruitment random effects (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-1,1)) +
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
  guides(col=guide_legend(title=NULL))
ggsave(recr.dev.rho.plot.ylim, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('recr.dev.rho.ylim', plot.suffix, '.png') ),  height=7, width=14)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











