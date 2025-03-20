# summarize convergence from AIC_weight.RDS, which is created from 
## dataframe.R with function call to get_aic_convergence_info

library('here')
library('tidyverse')
library('xtable')
library(rpart)

res.dir     <- 'results'  # 'results'     'results_beta_fix'
plot.dir    <- 'plots'    # 'plots_lizruns'  'plots_beta_fix'   
plot.suffix <- ''         # '_beta_fix'   '' 

n.sims <- 100

AIC_all     <- readRDS(file.path(here(),'Ecov_study','recruitment_functions',res.dir,paste0('AIC', plot.suffix, '.rds') ) )
AIC_weight  <- readRDS(file.path(here(),'Ecov_study','recruitment_functions',res.dir,paste0('AIC_weight', plot.suffix, '.rds') ) )
df.oms      <- readRDS(file.path(here::here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS"))
df.ems      <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))
om.fails.df <- readRDS(file.path(here::here("Ecov_study", 'recruitment_functions', res.dir,'om.fails.df.RDS')))


bad.grad.value <- 1E-6 #tim used 1E-6 (the abs of the exponent will be used for output suffix; ex: filename_grad_6.png)
bad.grad.label <- as.numeric(strsplit(as.character(bad.grad.value), split="-")[[1]][2])
bad.se.value   <- 100 #tim used 100 (this value will be used for output suffix; ex: filename_se_100.png)

# analyze AIC_all (info just saved for best model for given OM)  ====
## specify cut-off for non-converged runs for AIC_all ====
# convergence IF opt == 1; conv ==0; sdrep==0; max_grad<max.grad.threshold; SE_par_max< not too big;
bad.opt    <- which(as.numeric(AIC_all$opt) != 1)
bad.conv   <- which(as.numeric(AIC_all$conv) !=0)
bad.sdrep  <- which(as.numeric(AIC_all$sdrep)>0)
bad.grad   <- which(as.numeric(AIC_all$max_grad)>bad.grad.value) 
bad.se.big <- which(as.numeric(AIC_all$SE_par_max)>bad.se.value ) #tim used 100, but threshold depends on scale of parameter (could refine this category to be parameter specific or try to specify threshold on relative scale? "CV-like"?)
bad.se.na  <- which( is.na(AIC_all$SE_par_max))

bad.runs        <- c(bad.opt, bad.conv, bad.sdrep, bad.grad, bad.se.big, bad.se.na)
bad.runs.unique <- unique(c(bad.opt, bad.conv, bad.sdrep, bad.grad, bad.se.big, bad.se.na))
#length(bad.runs)
#length(unique(bad.runs))
#length(bad.opt)
#length(bad.conv)
#length(bad.sdrep)
#length(bad.grad)
#length(bad.se.big)


#characterize the bad.runs
AIC.bad <- as_tibble(AIC_all[bad.runs.unique,]) %>%
  mutate(EM_SR_name = ifelse(EM_SR==3, "BH", "Mean")) %>%
  mutate(EM_mod = paste0(EM_SR_name, "_", EM_ecov)) %>%
  relocate(EM_mod, Model) %>%
  mutate(R_sig = paste0('R_sig_', R_sig), Ecov_effect=paste0('Ecov_B_', Ecov_effect))

AIC_best_conv.runs <- AIC_all[-bad.runs.unique,]
saveRDS(AIC_best_conv.runs, file =  file.path(here(),'Ecov_study','recruitment_functions', res.dir, paste0("AIC_best_conv.runs_grad_", bad.grad.label, "_SE_", bad.se.value, plot.suffix, ".RDS")  ) )

pct.converge      <- nrow(AIC_best_conv.runs)/nrow(AIC_all)
pct.fail.converge <- 1-pct.converge
n.bad.runs        <- length(unique(bad.runs) )


bad.mods.plot <- ggplot(AIC.bad, aes(x=EM_mod)) +
  facet_grid(Ecov_effect  ~ Fhist +R_sig ) +
  geom_bar(fill='#2099aa99') +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=12)) +
  theme(axis.text.x = element_text(size = 8,angle=90))   + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  ylab('Number failed convg. checks (EM in OM-Sim with lowest AIC)') +
  labs(subtitle=paste0(100*round(pct.fail.converge,3), '% of lowest AIC model failed 1 or more convergence checks; max(abs(gradient)) > ', bad.grad.label, ' and/or par_SE > ', bad.se.value ))
ggsave(bad.mods.plot, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir, paste0("bad.mods_lowestAIC.plot_grad_",bad.grad.label, "_SE_", bad.se.value, plot.suffix,".png") ),  height=7, width=12)

# analyze AIC_weight (info for all OM - EM - Sim)  ====
## specify cut-off for non-converged runs for AIC_weight (same as for AIC_all) ====
# convergence IF opt == 1; conv ==0; sdrep==0; max_grad<max.grad.threshold; SE_par_max< not too big;
## first, identify bad runs
bad.opt    <- which(as.numeric(AIC_weight$opt) != 1)  # 0
bad.conv   <- which(as.numeric(AIC_weight$conv) !=0)  #2178
bad.sdrep  <- which(as.numeric(AIC_weight$sdrep)>0)   #5703
bad.grad   <- which(as.numeric(AIC_weight$max_grad)>bad.grad.value) # 11,154             #tim used 1E-6
bad.se.big <- which(as.numeric(AIC_weight$SE_par_max)>bad.se.value ) # 3==>33,325       100==>19,613        #tim used 100
bad.se.na  <- which( is.na(AIC_weight$SE_par_max))  #1936

bad.runs   <- c(bad.opt, bad.conv, bad.sdrep, bad.grad, bad.se.big, bad.se.na)
length(bad.runs)
length(unique(bad.runs))
bad.runs.unique <- unique(c(bad.opt, bad.conv, bad.sdrep, bad.grad, bad.se.big, bad.se.na))


AIC_weight_conv.runs <- AIC_weight[-bad.runs.unique,]
pct.converge         <-nrow(AIC_weight_conv.runs)/nrow(AIC_weight)
pct.fail.convg       <- 1-pct.converge
n.bad.runs           <- length(unique(bad.runs) )
saveRDS(AIC_weight_conv.runs, file =  file.path(here(),'Ecov_study','recruitment_functions', res.dir, paste0("AIC_weight_conv.runs_grad_", bad.grad.label, "_SE_", bad.se.value, plot.suffix,".RDS")  ) )


conv.summary <- as_tibble(AIC_weight) %>%
  mutate(bad.opt=ifelse(as.numeric(AIC_weight$opt) != 1, 1,0),     #1=bad, 0=not bad
         bad.conv = ifelse(as.numeric(AIC_weight$conv) !=0, 1,0),       #1=bad, 0=not bad
         bad.sdrep = ifelse(as.numeric(AIC_weight$sdrep)>0, 1,0),     #1=bad, 0=not bad
         bad.grad = ifelse(as.numeric(AIC_weight$max_grad)>bad.grad.value, 1,0),   #1=bad, 0=not bad
         bad.se.big = ifelse(as.numeric(AIC_weight$SE_par_max)>bad.se.value , 1,0),   #1=bad, 0=not bad
         bad.se.na = ifelse(is.na(AIC_weight$SE_par_max), 1, 0) 
  ) %>%
  replace_na(list(bad.opt=99, bad.conv=99, bad.sdrep=99, bad.grad=99, bad.se.big=99, bad.se.na=99)) %>%
  select(OM, EM, sim, bad.opt, bad.conv, bad.sdrep, bad.grad, bad.se.big, bad.se.na)  %>%
  mutate(bad.grad.value=bad.grad.value, bad.se.value=bad.se.value) %>%
  arrange(OM, sim, EM) %>%
  relocate(OM, sim, EM, bad.grad.value, bad.se.value)

saveRDS(conv.summary, file.path(here(),'Ecov_study','recruitment_functions', res.dir, paste0("conv.summary_grad_", bad.grad.label, "_SE_", bad.se.value, plot.suffix, ".RDS")  ) )


# isolate OM - EM - SIM that did not converge
nonconv.runs <- conv.summary %>%
  mutate(bad.sum=rowSums(across(c(bad.opt, bad.conv, bad.sdrep, bad.grad, bad.se.big, bad.se.na)))) %>%
  filter(bad.sum!=0) %>%  # bad.sum>0 had 1s or NAs (replaced by 99) for convergence criteria
  relocate(OM, EM, sim, bad.sum)
saveRDS(nonconv.runs, file.path(here(),'Ecov_study','recruitment_functions', res.dir, paste0("nonconv.runs_grad_", bad.grad.label, "_SE_", bad.se.value, plot.suffix, ".RDS") ) )


# identify the runs that did converge
conv.runs <- conv.summary %>%
  mutate(bad.sum=rowSums(across(c(bad.opt, bad.conv, bad.sdrep, bad.grad, bad.se.big, bad.se.na)))) %>%
  filter(bad.sum==0) %>%  # if bad.sum==0, all conv checks passed
  relocate(OM, EM, sim, bad.sum)
saveRDS(conv.runs, file.path(here(),'Ecov_study','recruitment_functions', res.dir, paste0("conv.runs_grad_", bad.grad.label, "_SE_", bad.se.value, plot.suffix, ".RDS") ) )



df.oms2 <- as_tibble(cbind(OM=seq(1,nrow(df.oms)), df.oms)) 
df.ems2 <- cbind(EM=seq(1,nrow(df.ems)), df.ems)
colnames(df.ems2) <- c("EM", "EM_ecov_how", "EM_r_mod")
em_tib <- as_tibble(df.ems2) %>%
  mutate(SR=ifelse(EM_r_mod==2, 'Mean', 'BH')) %>%
  mutate(EM_mod = paste0(SR, "_", EM_ecov_how))

non.conv.run.info <- nonconv.runs %>%
  left_join(as_tibble(df.oms2)) %>%
  left_join(em_tib) %>%
  mutate(R_sig = paste0('R_sig_', R_sig), Ecov_effect=paste0('Ecov_B_', Ecov_effect), 
         Ecov_re_sig=paste0('Ecov_re_sig_', Ecov_re_sig))


pct.converge<-nrow(conv.runs)/nrow(AIC_weight)
pct.fail.convg <- 1-pct.converge
n.bad.runs <- length((bad.runs.unique) )


bad.mods_all.plot <- ggplot(non.conv.run.info, aes(x=EM_mod)) +
  facet_grid(Ecov_effect ~ Fhist +R_sig   ) +
  geom_bar(fill='#2099aa99') +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=12)) +
  theme(axis.text.x = element_text(size = 8, angle=90))   + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  ylab('Number failed convg. checks') +
  labs(subtitle=paste0(100*round(pct.fail.convg,3), '% of all OM-EM-Sims failed 1 or more convergence checks; max(abs(gradient)) > ', bad.grad.value, ' and/or par_SE > ', bad.se.value  ))
ggsave(bad.mods_all.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0("bad.mods_all_grad_", bad.grad.label, "_SE_", bad.se.value, plot.suffix, ".plot.png" ) ),  height=7, width=12)


###################################################################
# characterize crashes ====  


om_fail_tib <- as_tibble(om.fails.df) %>%
  rename(EM=em.fails, Iter=iter.fails) %>% 
  left_join(df.oms2) %>%
  left_join(em_tib) %>%
  mutate(mod.match=ifelse(EM_ecov_how==Ecov_how & EM_r_mod==recruit_mod, 1, 0)) %>%
  mutate(fail.iter.em = paste0(N.fail, '_', EM)) %>%
  relocate(OM, N.fail, Iter, EM, EM_mod)

sim_tib <- as_tibble(cbind(OM=rep(seq(1, nrow(df.oms)), each=n.sims*nrow(df.ems)) , Iter=rep(seq(1,n.sims),nrow(df.oms)*nrow(df.ems)) , EM=rep(seq(1,6), nrow(df.oms)*n.sims) ) ) %>%
  left_join(df.oms2)
  

om_all_tib <- sim_tib %>%
  full_join(om_fail_tib) %>%
  mutate(Crash=ifelse(is.na(N.fail), 0, 1)) %>%
  relocate(OM, Iter, EM, Crash)

om_fail_tib$R_sig <- factor(om_fail_tib$R_sig,labels = c("Rsig_0.1","Rsig_1.0"))
om_fail_tib$Ecov_how    <- factor(om_fail_tib$Ecov_how,labels=c("Ecov_0", "Ecov_1","Ecov_2","Ecov_4"))
om_fail_tib$NAA_cor     <- factor(om_fail_tib$NAA_cor,labels=c("Rcor_L","Rcor_H"))
om_fail_tib$Fhist       <- factor(om_fail_tib$Fhist,labels=c("H-MSY","MSY") )
om_fail_tib$Ecov_effect <- factor(om_fail_tib$Ecov_effect,labels=c("Ecov_L", "Ecov_H"))
om_fail_tib$Ecov_re_cor <- factor(om_fail_tib$Ecov_re_cor,labels=c("EcovCor_L","EcovCor_H"))

om_all_tib$R_sig <- factor(om_all_tib$R_sig,labels = c("Rsig_0.1","Rsig_1.0"))
om_all_tib$Ecov_how    <- factor(om_all_tib$Ecov_how,labels=c("Ecov_0", "Ecov_1","Ecov_2","Ecov_4"))
om_all_tib$NAA_cor     <- factor(om_all_tib$NAA_cor,labels=c("Rcor_L","Rcor_H"))
om_all_tib$Fhist       <- factor(om_all_tib$Fhist,labels=c("H-MSY","MSY") )
om_all_tib$Ecov_effect <- factor(om_all_tib$Ecov_effect,labels=c("Ecov_L", "Ecov_H"))
om_all_tib$Ecov_re_cor <- factor(om_all_tib$Ecov_re_cor,labels=c("EcovCor_L","EcovCor_H"))

# regression tree for Crashes ======================================
rf_crash   <- rpart(Crash   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + EM, data=om_all_tib, control=rpart.control(cp=0.01))
imp.var <- rf_crash$frame[rf_crash$frame$var != '<leaf>',]
nodes_crash <- unique(imp.var[,1])
# "EM"  

crashed.models.plot <- ggplot(om_fail_tib, aes(x=EM_mod, fill=as.factor(mod.match) )) +
  facet_grid(R_sig  + Ecov_how ~ Fhist  + Ecov_effect  ) +
  geom_bar() +
  ylab('N crashes') +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=11)) +
  theme(axis.text.x = element_text(size = 9))   + 
  theme(axis.text.y = element_text(size = 11))   +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  theme(legend.position="bottom") +
  theme(legend.text = element_text(colour="black", size = 14, face = "bold")) +
  scale_fill_discrete(labels = c("OM != EM", "OM = EM")) +
  guides(fill=guide_legend(title=NULL))
ggsave(crashed.models.plot, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('crashed.models.plot', plot.suffix, '.png') ),  height=7, width=12)



# Characterize parameters that were problemmatic across all OM-EM-Sims ====
# unique(AIC_bad_SE_par$SE_par_max_name)
# unique(AIC_bad_SE_par$max_grad_name)

bad_grad_par <- as_tibble(AIC_weight[bad.grad,]) %>%
  replace_na(list(max_grad_name = '0_Model.crash')) %>%
  group_by( max_grad_name) %>%
    summarise(ntimes.bad.grad=n()) %>%
  rename(Par=max_grad_name)

bad_SE_par <- as_tibble(AIC_weight[bad.se.big,]) %>%
  group_by( SE_par_max_name) %>%
  summarise(ntimes.bad.se=n()) %>% 
  rename(Par=SE_par_max_name)




bad_par_table <- bad_grad_par %>%
  full_join(bad_SE_par) %>%
  arrange(Par) %>%
  rename(N_Bad.Grad=ntimes.bad.grad, N_Big.SE=ntimes.bad.se) 
write.csv(bad_par_table, file=file.path(here(),'Ecov_study','recruitment_functions','tables', paste0("Bad.param.summary.table_grad_", bad.grad.label, "_SE_", bad.se.value, plot.suffix, ".csv")  ), row.names = FALSE )


bad_par_long <- bad_par_table %>%
  mutate(tot.bad.grad=length(bad.grad), tot.big.se=length(bad.se.big)) %>%
  mutate(Pct_Bad.Grad = N_Bad.Grad/tot.bad.grad, Pct_Big.SE=N_Big.SE/tot.big.se ) %>%
  select(-c(tot.bad.grad, tot.big.se) ) %>%
  pivot_longer(cols=-Par,  names_to=c(".value", 'Condition'), names_sep="_" )

bad_par_N_barplot <- ggplot(bad_par_long, aes(x=N, y=Par)) +
  facet_grid(~Condition, labeller = as_labeller(c(Bad.Grad=paste0('max(abs(gradient))> ', bad.grad.value),
                                                  Big.SE = paste0('max(SE_par)> ', bad.se.value)
                                                  ) 
                                                )
             ) +
  geom_col(fill="#2099aa99") +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=12)) +
  theme(axis.text.x = element_text(size = 10))   + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  ylab('Parameter') +
  xlab('Count') +
  labs(subtitle=paste0('Number of times each parameter had the largest gradient or the biggest SE' ))
ggsave(bad_par_N_barplot , filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir, paste0("bad_par_N_barplot_grad_", bad.grad.label, "_SE_", bad.se.value, plot.suffix, ".png")  ),  height=7, width=12)


bad_par_pct_barplot <- ggplot(bad_par_long, aes(x=Pct, y=Par)) +
  facet_grid(~Condition, labeller = as_labeller(c(Bad.Grad=paste0('max(abs(gradient))> ', bad.grad.value),
                                                  Big.SE = paste0('max(SE_par)> ', bad.se.value)
  ) 
  )
  ) +
  geom_col(fill="#2099aa99") +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=12)) +
  theme(axis.text.x = element_text(size = 10))   + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  ylab('Parameter') +
  xlab('Proportion') +
labs(subtitle=paste0('Proportion of times each parameter had the largest gradient or the biggest SE' ))
ggsave(bad_par_pct_barplot, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir, paste0("bad_par_pct_barplot_grad_", bad.grad.label, "_SE_", bad.se.value, plot.suffix, ".png" ) ),  height=7, width=12)



##################################################################################################
###  make plot for 2nd largest SE parameter (since first largest was dominated by mean_rec_pars)   ====  



bad_SE2_par <- as_tibble(AIC_weight[bad.se.big,]) %>%
  group_by( SE_par_max2_name) %>%
  summarise(ntimes.bad.se2=n()) %>% 
  rename(Par=SE_par_max2_name)




bad_par2_table <- bad_grad_par %>%
  full_join(bad_SE2_par) %>%
  arrange(Par) %>%
  rename(N_Bad.Grad=ntimes.bad.grad, N_Big.SE2=ntimes.bad.se2) 
write.csv(bad_par2_tib, file=file.path(here(),'Ecov_study','recruitment_functions','tables', paste0("Bad.param2.summary.table_grad_", bad.grad.label, "_SE_", bad.se.value, plot.suffix, ".csv")  ), row.names = FALSE )


bad_par2_long <- bad_par2_table %>%
  mutate(tot.bad.grad=length(bad.grad), tot.big.se2=length(bad.se.big)) %>%
  mutate(Pct_Bad.Grad = N_Bad.Grad/tot.bad.grad, Pct_Big.SE2=N_Big.SE2/tot.big.se2 ) %>%
  select(-c(tot.bad.grad, tot.big.se2) ) %>%
  pivot_longer(cols=-Par,  names_to=c(".value", 'Condition'), names_sep="_" )

bad_par2_N_barplot <- ggplot(bad_par2_long, aes(x=N, y=Par)) +
  facet_grid(~Condition, labeller = as_labeller(c(Bad.Grad=paste0('max(abs(gradient))> ', bad.grad.value),
                                                  Big.SE2 = paste0('2nd_max(SE_par)> ', bad.se.value)
  ) 
  )
  ) +
  geom_col(fill="#2099aa99") +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=12)) +
  theme(axis.text.x = element_text(size = 10))   + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  ylab('Parameter') +
  xlab('Count') +
  labs(subtitle=paste0('Number of times each parameter had the largest gradient or the biggest SE' ))
ggsave(bad_par2_N_barplot , filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir, paste0("bad_par_2nd_N_barplot_grad_", bad.grad.label, "_SE_", bad.se.value, plot.suffix, ".png")  ),  height=7, width=12)




