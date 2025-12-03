# summarize convergence from AIC_weight.RDS, which is created from 
## dataframe.R with function call to get_aic_convergence_info
# ** modifying bar plot from frequencies to proportions
# ** liz brooks dec 2025

library('here')
library('tidyverse')
# library('xtable')


rds.dir     <- 'rds_output'
# res.dir     <- 'results'  #  
plot.dir    <- 'plots'    #  
plot.suffix <- ''         # 


AIC_all1     <- as_tibble(readRDS(file.path(here(),'Ecov_study','recruitment_functions',rds.dir,paste0('AIC', plot.suffix, '.rds') ) ))
AIC_weight1  <- as_tibble(readRDS(file.path(here(),'Ecov_study','recruitment_functions',rds.dir,paste0('AIC_weight', plot.suffix, '.rds') ) ) )
df.oms      <- as_tibble(readRDS(file.path(here::here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS")) )
df.ems      <- as_tibble(readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS")) )
om.fails.df <- as_tibble(readRDS(file.path(here::here("Ecov_study", 'recruitment_functions', rds.dir,'om.fails.df.RDS'))) )


bad.grad.value <- 1E-6 #tim used 1E-6 (the abs of the exponent will be used for output suffix; ex: filename_grad_6.png)
bad.grad.label <- as.numeric(strsplit(as.character(bad.grad.value), split="-")[[1]][2])
bad.se.value   <- 100 #tim used 100 (this value will be used for output suffix; ex: filename_se_100.png)

# data-manipulation  ====
## specify cut-off for non-converged runs (see case_when statement) ====
# convergence IF opt == 1 & conv ==0 & sdrep==0 & max_grad<max.grad.threshold & SE_par_max< not too big;

#data frame of AIC for lowest model across all OM sims / factor
AIC_all <- AIC_all1 %>%
  mutate(EM_SR_name = ifelse(EM_SR==3, "BH", "Mean")) %>%
  mutate(EM_mod = paste0(EM_SR_name, "_", EM_ecov)) %>%
  relocate(EM_mod, Model) %>%
  mutate(R_sig = paste0('R_sig_', R_sig), Ecov_effect=paste0('Ecov_B_', Ecov_effect)) %>%
  mutate(bad.run=case_when(
    as.numeric(opt) != 1 ~ 1,
    as.numeric(conv) !=0 ~ 1,
    as.numeric(sdrep)>0  ~ 1,
    as.numeric(max_grad)>bad.grad.value ~ 1,
    as.numeric(SE_par_max)>bad.se.value ~ 1,
    is.na(SE_par_max)  ~ 1,
    TRUE               ~ 0   #this applies to all remaining cases (converged runs)
  )
  )  %>%
  mutate(bad.run = as.logical(bad.run)) %>%
  mutate(conv.run = as.integer(!bad.run) ) %>%
  mutate(N=1)

#data frame of AIC and AIC_weight for all models across all OM sims / factor
AIC_weight <- AIC_weight1 %>%
  mutate(EM_SR_name = ifelse(EM_r_mod==3, "BH", "Mean")) %>%
  mutate(EM_mod = paste0(EM_SR_name, "_", EM_ecov_how)) %>%
  relocate(EM_mod, Model) %>%
  mutate(R_sig = paste0('R_sig_', R_sig), Ecov_effect=paste0('Ecov_B_', Ecov_effect)) %>%
  mutate(bad.run=case_when(
    as.numeric(opt) != 1 ~ 1,
    as.numeric(conv) !=0 ~ 1,
    as.numeric(sdrep)>0  ~ 1,
    as.numeric(max_grad)>bad.grad.value ~ 1,
    as.numeric(SE_par_max)>bad.se.value ~ 1,
    is.na(SE_par_max)  ~ 1,
    TRUE               ~ 0   #this applies to all remaining cases (converged runs)
  )
  )  %>%
  mutate(bad.run = as.logical(bad.run)) %>%
  mutate(conv.run = !bad.run) %>%
  mutate(N=1)


# calculate % that don't converge  ====

pct.fail.converge.aic.best <- sum(AIC_all$bad.run) / nrow(AIC_all)              # 0.08378472
pct.fail.converge.aic.weight <- sum(AIC_weight$bad.run) / nrow(AIC_weight)              # 0.2650035



# plots for the model with lowest AIC ====
# proportion plot (misleading because different N for each EM)
bad.mods_lowestAIC.prop.plot <-   ggplot(AIC_all, aes(x=EM_mod, y=N, fill=bad.run)) +
  facet_grid(Ecov_effect  ~ Fhist +R_sig ) +
  geom_col(position="fill" ) +
    scale_fill_manual(values=c('TRUE'='#C94C16', 'FALSE'=NA))+
  theme_light()  +
  theme(strip.background =element_rect(fill="white", color="grey65"))+
  theme(strip.text = element_text(colour = 'black', size=10)) +
  theme(axis.text.x = element_text(size = 10, angle=90))   + 
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.title.x = element_text(size = 12))   + 
  theme(axis.title.y = element_text(size = 11))   +
    theme(legend.position = "none") +
    ylab('Proportion failed convergence checks (when EM had lowest AIC)') +
  xlab('EM')  +
  labs(subtitle=paste0(100*round(pct.fail.converge.aic.best,3), '% of lowest AIC model failed 1 or more convergence checks, max(abs(gradient)) > ', 
                       bad.grad.label, ' and/or SE(parameter) > ', bad.se.value ))

  ggsave(bad.mods_lowestAIC.prop.plot, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir, paste0("bad.mods_lowestAIC.prop.plot_grad_",bad.grad.label, "_SE_", bad.se.value, plot.suffix,".pdf") ),  
         units="in", height=7, width=12, device = "pdf", dpi=300)


# number in each EM category differs, so plot as frequency of converged vs non-converged
## note labeling switch for legend (legend title is converged and categories for Yes or No) 
  freq.good.bad.mods_lowestAIC.plot <-   ggplot(AIC_all, aes(x=EM_mod, y=N, fill=bad.run)) +
    facet_grid(Ecov_effect  ~ Fhist +R_sig ) +
    geom_bar(position = "stack", stat = "identity" ) +
    scale_fill_manual(values=c('TRUE'='#C94C16' , 'FALSE'='#2798F5'), labels=c('Yes', 'No')  )+
    theme_light()  +
    theme(strip.background =element_rect(fill="white", color="grey65"))+
    theme(strip.text = element_text(colour = 'black', size=10)) +
    theme(axis.text.x = element_text(size = 10, angle=90))   + 
    theme(axis.text.y = element_text(size = 10)) +
    theme(axis.title.x = element_text(size = 12))   + 
    theme(axis.title.y = element_text(size = 11))   +
    labs(fill = "Converged") +
    theme(legend.position = "bottom") +
    ylab('N of Sims where EM had lowest AIC') +
    xlab('EM') +
    labs(subtitle=paste0(100*round(pct.fail.converge.aic.best,3), '% of lowest AIC model failed 1 or more convergence checks, max(abs(gradient)) > ', 
                         bad.grad.label, ' and/or SE(parameter) > ', bad.se.value ))
  
  ggsave(freq.good.bad.mods_lowestAIC.plot, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir, paste0("bad.mods_lowestAIC.freq.plot_grad_",bad.grad.label, "_SE_", bad.se.value, plot.suffix,".pdf") ),  
         units="in", height=7, width=12, device = "pdf", dpi=300)
  
  
# plots for the all EM models (not just lowest) ====
# proportion plot ( N is 2400 for each EM)  
  bad.mods_allAIC.prop.plot <-   ggplot(AIC_weight, aes(x=EM_mod, y=N, fill=bad.run)) +
    facet_grid(Ecov_effect  ~ Fhist +R_sig ) +
    geom_col(position="fill" ) +
    scale_fill_manual(values=c('TRUE'='#C94C16', 'FALSE'=NA))+
    theme_light()  +
    theme(strip.background =element_rect(fill="white", color="grey65"))+
    theme(strip.text = element_text(colour = 'black', size=10)) +
    theme(axis.text.x = element_text(size = 10, angle=90))   + 
    theme(axis.text.y = element_text(size = 10)) +
    theme(axis.title.x = element_text(size = 12))   + 
    theme(axis.title.y = element_text(size = 11))   +
    theme(legend.position = "none") +
    ylab('Proportion failed convergence checks') +
    xlab('EM') +
    labs(subtitle=paste0(100*round(pct.fail.converge.aic.weight,3), '% of all OM-EM-Simulations failed 1 or more convergence checks, max(abs(gradient)) > ', bad.grad.value, ' and/or SE(parameter) > ', bad.se.value  ))
  
  ggsave(bad.mods_allAIC.prop.plot, filename=file.path(here(),'Ecov_study','recruitment_functions',plot.dir, paste0("bad.mods_allAIC.prop.plot_grad_",bad.grad.label, "_SE_", bad.se.value, plot.suffix,".pdf") ),  
         units="in", height=7, width=12, device = "pdf", dpi=300)
  
