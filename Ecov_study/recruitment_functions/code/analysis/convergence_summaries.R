# summarize convergence from AIC_weight.RDS, which is created from 
## dataframe.R with function call to get_aic_convergence_info

library("here")
library('tidyverse')


AIC_all <- readRDS(file.path(here(),'Ecov_study','recruitment_functions','results','AIC.rds'))
AIC_weight <- readRDS(file.path(here(),'Ecov_study','recruitment_functions','results','AIC_weight.rds'))
df.oms          <- readRDS(file.path(here::here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS"))
df.ems    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))



# analyze AIC_all (info just saved for best model for given OM)  ====
## specify cut-off for non-converged runs for AIC_all ====
# convergence IF opt == 1; conv ==0; sdrep==0; max_grad<max.grad.threshold; SE_par_max< not too big;
bad.opt <- which(AIC_all$opt != 1)
bad.conv <- which(AIC_all$conv !=0)
bad.sdrep <- which(AIC_all$sdrep>0)
bad.grad <- which(AIC_all$max_grad>1E-6) #tim used 1E-6
bad.se.big <- which(AIC_all$SE_par_max>3 ) #tim used 100, but threshold depends on scale of parameter (could refine this category to be parameter specific or try to specify threshold on relative scale? "CV-like?)
bad.se.na <- which( is.na(AIC_all$SE_par_max))

bad.runs <- c(bad.opt, bad.conv, bad.sdrep, bad.grad, bad.se.big, bad.se.na)
length(bad.runs)
length(unique(bad.runs))
bad.runs.unique <- unique(c(bad.opt, bad.conv, bad.sdrep, bad.grad, bad.se.big, bad.se.na))

sim.conv <- rep(0, nrow(AIC_all))
sim.conv[bad.runs.unique] <- 1

#characterize the bad.runs
AIC.bad <- as_tibble(AIC_all[bad.runs,]) %>%
  mutate(EM_SR_name = ifelse(EM_SR==3, "BH", "Mean")) %>%
  mutate(EM_mod = paste0(EM_SR_name, "_", EM_ecov)) %>%
  relocate(EM_mod, Model) %>%
  mutate(R_sig = paste0('R_sig_', R_sig), Ecov_effect=paste0('Ecov_B_', Ecov_effect))

AIC <- AIC_all[-bad.runs,]
pct.converge<-nrow(AIC)/nrow(AIC_all)
pct.fail.convg <- 1-pct.converge
n.bad.runs <- length(unique(bad.runs) )


bad.mods.plot <- ggplot(AIC.bad, aes(x=EM_mod)) +
  facet_grid(Ecov_effect  ~ Fhist +R_sig ) +
  geom_bar(fill='#2099aa99') +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=12)) +
  theme(axis.text.x = element_text(size = 10))   + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  ylab('Number failed convg. checks (OM-EM-Sim with lowest AIC)') +
  labs(subtitle=paste0(100*round(pct.fail.convg,3), '% of lowest AIC model failed 1 or more convergence checks' ))
ggsave(bad.mods.plot, filename=file.path(here(),'Ecov_study','recruitment_functions','plots', "bad.mods_lowestAIC.plot.png"),  height=7, width=12)





# analyze AIC_weight (info for all OM - EM - Sim)  ====
## specify cut-off for non-converged runs for AIC_weight (same as for AIC_all) ====
# convergence IF opt == 1; conv ==0; sdrep==0; max_grad<max.grad.threshold; SE_par_max< not too big;
## first, identify bad runs
bad.opt <- which(as.numeric(AIC_weight$opt) != 1)  # 0
bad.conv <- which(as.numeric(AIC_weight$conv) !=0)  #2178
bad.sdrep <- which(as.numeric(AIC_weight$sdrep)>0)   #5703
bad.grad <- which(as.numeric(AIC_weight$max_grad)>1E-6) # 11,154             #tim used 1E-6
bad.se.big <- which(as.numeric(AIC_weight$SE_par_max)>100 ) # 3==>33,325       100==>19,613        #tim used 100
bad.se.na <- which( is.na(AIC_weight$SE_par_max))  #1936

bad.runs <- c(bad.opt, bad.conv, bad.sdrep, bad.grad, bad.se.big, bad.se.na)
length(bad.runs)
length(unique(bad.runs))
bad.runs.unique <- unique(c(bad.opt, bad.conv, bad.sdrep, bad.grad, bad.se.big, bad.se.na))


conv.summary <- as_tibble(AIC_weight) %>%
  mutate(bad.opt=ifelse(as.numeric(AIC_weight$opt) != 1, 1,0),     #1=bad, 0=not bad
         bad.conv = ifelse(as.numeric(AIC_weight$conv) !=0, 1,0),       #1=bad, 0=not bad
         bad.sdrep = ifelse(as.numeric(AIC_weight$sdrep)>0, 1,0),     #1=bad, 0=not bad
         bad.grad = ifelse(as.numeric(AIC_weight$max_grad)>1E-6, 1,0),   #1=bad, 0=not bad
         bad.se.big = ifelse(as.numeric(AIC_weight$SE_par_max)>100 , 1,0),   #1=bad, 0=not bad
         bad.se.na = ifelse(is.na(AIC_weight$SE_par_max), 1, 0) 
  ) %>%
  replace_na(list(bad.opt=99, bad.conv=99, bad.sdrep=99, bad.grad=99, bad.se.big=99, bad.se.na=99)) %>%
  select(OM, EM, sim, bad.opt, bad.conv, bad.sdrep, bad.grad, bad.se.big, bad.se.na)  %>%
  arrange(OM, sim, EM)

saveRDS(conv.summary, file.path(here(),'Ecov_study','recruitment_functions','results', "conv.summary.RDS") )


# isolate OM - EM - SIM that did not converge
nonconv.runs <- conv.summary %>%
  mutate(bad.sum=rowSums(across(c(bad.opt, bad.conv, bad.sdrep, bad.grad, bad.se.big, bad.se.na)))) %>%
  filter(bad.sum!=0) %>%  # bad.sum>0 had 1s or NAs (replaced by 99) for convergence criteria
  relocate(OM, EM, sim, bad.sum)
saveRDS(nonconv.runs, file.path(here(),'Ecov_study','recruitment_functions','results', "nonconv.runs.RDS") )


# identify the runs that did converge
conv.runs <- conv.summary %>%
  mutate(bad.sum=rowSums(across(c(bad.opt, bad.conv, bad.sdrep, bad.grad, bad.se.big, bad.se.na)))) %>%
  filter(bad.sum==0) %>%  # if bad.sum==0, all conv checks passed
  relocate(OM, EM, sim, bad.sum)
saveRDS(conv.runs, file.path(here(),'Ecov_study','recruitment_functions','results', "conv.runs.RDS") )



df.oms2 <- as_tibble(cbind(OM=seq(1,256), df.oms)) 
df.ems2 <- cbind(EM=seq(1,6), df.ems)
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
n.bad.runs <- length(unique(bad.runs.unique) )


bad.mods_all.plot <- ggplot(non.conv.run.info, aes(x=EM_mod)) +
  facet_grid(Ecov_effect + Ecov_re_sig ~ Fhist +R_sig   ) +
  geom_bar(fill='#2099aa99') +
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=12)) +
  theme(axis.text.x = element_text(size = 10))   + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 13))   + 
  theme(axis.title.y = element_text(size = 13))   +
  ylab('Number failed convg. checks (OM-EM-Sim with lowest AIC)') +
  labs(subtitle=paste0(100*round(pct.fail.convg,3), '% of all OM-EM-Sims failed 1 or more convergence checks' ))
ggsave(bad.mods_all.plot, filename=file.path(here(),'Ecov_study','recruitment_functions','plots', "bad.mods_all.plot.png"),  height=7, width=12)



# Characterize parameters that were problemmatic across all OM-EM-Sims ====
unique(AIC_bad_SE_par$SE_par_max_name)
unique(AIC_bad_SE_par$max_grad_name)

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
write.csv(bad_par_tib, file=file.path(here(),'Ecov_study','recruitment_functions','tables', "Bad.param.summary.table.csv"), row.names = FALSE )


bad_par_long <- bad_par_table %>%
  mutate(tot.bad.grad=length(bad.grad), tot.big.se=length(bad.se.big)) %>%
  mutate(Pct_Bad.Grad = N_Bad.Grad/tot.bad.grad, Pct_Big.SE=N_Big.SE/tot.big.se ) %>%
  select(-c(tot.bad.grad, tot.big.se) ) %>%
  pivot_longer(cols=-Par,  names_to=c(".value", 'Condition'), names_sep="_" )

bad_par_N_barplot <- ggplot(bad_par_long, aes(x=N, y=Par)) +
  facet_grid(~Condition) +
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
ggsave(bad_par_N_barplot , filename=file.path(here(),'Ecov_study','recruitment_functions','plots', "bad_par_N_barplot .png"),  height=7, width=12)


bad_par_pct_barplot <- ggplot(bad_par_long, aes(x=Pct, y=Par)) +
  facet_grid(~Condition) +
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
ggsave(bad_par_pct_barplot, filename=file.path(here(),'Ecov_study','recruitment_functions','plots', "bad_par_pct_barplot .png"),  height=7, width=12)
  
