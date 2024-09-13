# analysis for factors influencing RE ====
par_mat_tib$R_sig <- factor(par_mat_tib$R_sig,labels = c("Rsig_0.1","Rsig_0.5","Rsig_1.0"))
par_mat_tib$Ecov_effect <- factor(par_mat_tib$Ecov_effect,labels=c("Ecov_L","Ecov_H"))
par_mat_tib$Ecov_how    <- factor(par_mat_tib$Ecov_how,labels=c("Ecov_0", "Ecov_1","Ecov_2","Ecov_4"))
par_mat_tib$NAA_cor     <- factor(par_mat_tib$NAA_cor,labels=c("Rcor_L","Rcor_H"))
par_mat_tib$Fhist       <- factor(par_mat_tib$Fhist,labels=c("H-MSY","MSY") ) 
par_mat_tib$Ecov_re_cor <- factor(par_mat_tib$Ecov_re_cor,labels=c("EcovCor_L","EcovCor_H"))
par_mat_tib$obs_error   <- factor(par_mat_tib$obs_error,labels=c("ObsErr_H","ObsErr_L"))

# regression tree for correct form (or SR or Ecov) ======================================
rf_RE.SRa   <- rpart(RE.mean_rec1   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=par_mat_tib[par_mat_tib$EM>2,], control=rpart.control(cp=0.01)) # *** only EM>2 ====
imp.var <- rf_RE.SRa$frame[rf_RE.SRa$frame$var != '<leaf>',]
nodes_RE.SRa <- unique(imp.var[,1])
nodes_RE.SRa
# nothing   # still nothing (liz 384)

rf_RE.SRb   <- rpart(RE.mean_rec2   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=par_mat_tib[par_mat_tib$EM>2,], control=rpart.control(cp=0.01)) # *** only EM>2 ====
imp.var <- rf_RE.SRb$frame[rf_RE.SRb$frame$var != '<leaf>',]
nodes_RE.SRb <- unique(imp.var[,1])
nodes_RE.SRb
# nothing   # still nothing (liz 384)

rf_RE.NAA_sigma   <- rpart(RE.NAA_sigma   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=par_mat_tib, control=rpart.control(cp=0.01)) 
imp.var <- rf_RE.NAA_sigma$frame[rf_RE.NAA_sigma$frame$var != '<leaf>',]
nodes_RE.NAA_sigma <- unique(imp.var[,1])
nodes_RE.NAA_sigma
# "NAA_cor"     # NAA_cor (liz 384)


rf_RE.NAA_rho   <- rpart(RE.NAA_rho   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=par_mat_tib, control=rpart.control(cp=0.01))
imp.var <- rf_RE.NAA_rho$frame[rf_RE.NAA_rho$frame$var != '<leaf>',]
nodes_RE.NAA_rho <- unique(imp.var[,1])
nodes_RE.NAA_rho
# "Fhist"   "R_sig"   "NAA_cor"   #  nothing (liz 384)

rf_RE.Ecov_process_pars1   <- rpart(RE.Ecov_process_pars1   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=par_mat_tib, control=rpart.control(cp=0.01))
imp.var <- rf_RE.Ecov_process_pars1$frame[rf_RE.Ecov_process_pars1$frame$var != '<leaf>',]
nodes_RE.Ecov_process_pars1 <- unique(imp.var[,1])
nodes_RE.Ecov_process_pars1
# nothing  # still nothing (liz 384)

rf_RE.Ecov_process_pars2   <- rpart(RE.Ecov_process_pars2   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=par_mat_tib, control=rpart.control(cp=0.01))
imp.var <- rf_RE.Ecov_process_pars2$frame[rf_RE.Ecov_process_pars2$frame$var != '<leaf>',]
nodes_RE.Ecov_process_pars2 <- unique(imp.var[,1])
nodes_RE.Ecov_process_pars2
# "Ecov_re_cor"  # "Ecov_re_cor" (liz 384)

rf_RE.Ecov_process_pars3   <- rpart(RE.Ecov_process_pars3   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + obs_error, data=par_mat_tib, control=rpart.control(cp=0.01))
imp.var <- rf_RE.Ecov_process_pars3$frame[rf_RE.Ecov_process_pars3$frame$var != '<leaf>',]
nodes_RE.Ecov_process_pars3 <- unique(imp.var[,1])
nodes_RE.Ecov_process_pars3
# nothing  # "Ecov_re_cor" (liz 384)


# plot the regression trees that showed some factor sensitivity ====
pdf(file.path(here(),'Ecov_study','recruitment_functions',plot.dir,paste0('reg_tree_NAA_RE', plot.suffix, '.pdf') ),
    height=4,width=7)
par(mfrow=c(1,2), oma=c(5,0,0,0))
prp(rf_RE.NAA_sigma,yesno=FALSE,type=4,clip.right.labs=TRUE)
mtext('a) RE(R_sigma)',adj=0,line=2.5, cex=0.9)
prp(rf_RE.NAA_rho,yesno=FALSE,type=4,clip.right.labs=TRUE)
mtext('b) RE(R_rho)',adj=0,line=2.5, cex=0.9)
# title(sub=paste0(100*round(aic.best.pct.conv,2), "% of runs converged (", aic.best.n.bad, " out of ", aic.best.n, " failed)" ),   adj=0.5, outer=TRUE)

dev.off()

# summarize the RE that we want to look at ======

par_mean_rec1_sum <- par_mat_tib %>%
  left_join(em_tib) %>%
  select(RE.mean_rec1, NAA_cor, obs_error, R_sig, Ecov_re_cor, Fhist, EM, EM_ecov_how, EM_r_mod, Ecov_how, recruit_mod )  %>%
  filter(EM>2) %>%
  mutate(mod.match=ifelse(EM_ecov_how==as.numeric(substr(Ecov_how, 6,6)) & EM_r_mod==recruit_mod, 1, 0)) %>%
  group_by( R_sig,  Fhist,  Ecov_how, EM, mod.match ) %>%
  summarise(mean.RE.mean_rec1=mean(RE.mean_rec1), var.RE.mean_rec1=var(RE.mean_rec1), 
            median.RE.mean_rec1=median(RE.mean_rec1), across(RE.mean_rec1,                                                                                                       list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                                                                                                                      p97.5=~quantile(.,probs=0.975)))
  )%>%
  ungroup() %>%
  left_join(em_tib)

#  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par_mean_rec2_sum <- par_mat_tib %>%
  left_join(em_tib) %>%
  select(RE.mean_rec2, NAA_cor, obs_error, R_sig, Ecov_re_cor, Fhist, EM, EM_ecov_how, EM_r_mod, Ecov_how, recruit_mod )  %>%
  filter(EM>2) %>%
  mutate(mod.match=ifelse(EM_ecov_how==as.numeric(substr(Ecov_how, 6,6)) & EM_r_mod==recruit_mod, 1, 0)) %>%
  group_by( R_sig,  Fhist,  Ecov_how, EM, mod.match ) %>%
  summarise(mean.RE.mean_rec2=mean(RE.mean_rec2, na.rm=T), var.daic=var(RE.mean_rec2, na.rm=T), ## calculating for all EM even though only has meaning for EM>2
            median.RE.mean_rec2=median(RE.mean_rec2, na.rm=T), across(RE.mean_rec2, 
                                                                      list(p2.5=~quantile(.,probs=0.0275 , na.rm=T),
                                                                           p97.5=~quantile(.,probs=0.975, na.rm=T))) 
  )%>%
  ungroup() %>%
  left_join(em_tib)

#  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par_R_sigma_sum <- par_mat_tib %>%
  left_join(em_tib) %>%
  select(RE.NAA_sigma, NAA_cor, obs_error, R_sig, Ecov_re_cor, Fhist, EM, EM_ecov_how, EM_r_mod, Ecov_how, recruit_mod )  %>%
  mutate(mod.match=ifelse(EM_ecov_how==as.numeric(substr(Ecov_how, 6,6)) & EM_r_mod==recruit_mod, 1, 0)) %>%
  group_by( R_sig,  Fhist, NAA_cor, Ecov_how, EM, mod.match ) %>%
  summarise(mean.RE.NAA_sigma=mean(RE.NAA_sigma), var.RE.NAA_sigma=var(RE.NAA_sigma), 
            median.RE.NAA_sigma=median(RE.NAA_sigma), across(RE.NAA_sigma,                                                                                                       list(p2.5=~quantile(.,probs=0.0275 ),                                                                                p97.5=~quantile(.,probs=0.975))) 
  )%>%
  ungroup() %>%
  left_join(em_tib)

#  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par_R_cor_sum <- par_mat_tib %>%
  left_join(em_tib) %>%
  select(RE.NAA_rho, NAA_cor, obs_error, R_sig, Ecov_re_cor, Fhist, EM, EM_ecov_how, EM_r_mod, Ecov_how, recruit_mod )  %>%
  mutate(mod.match=ifelse(EM_ecov_how==as.numeric(substr(Ecov_how, 6,6)) & EM_r_mod==recruit_mod, 1, 0)) %>%
  group_by( R_sig,  Fhist, NAA_cor, Ecov_how, EM, mod.match ) %>%
  summarise(mean.RE.NAA_rho=mean(RE.NAA_rho), var.RE.NAA_rho=var(RE.NAA_rho), 
            median.RE.NAA_rho=median(RE.NAA_rho), across(RE.NAA_rho,                                                                                                       list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                                                                                                                p97.5=~quantile(.,probs=0.975))) 
  )%>%
  ungroup() %>%
  left_join(em_tib)


#  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par_Ecov1_sum <- par_mat_tib %>%
  left_join(em_tib) %>%
  select(RE.Ecov_process_pars1, Ecov_re_cor, obs_error, R_sig, Ecov_re_cor, Fhist, EM, EM_ecov_how, EM_r_mod, Ecov_how, recruit_mod )  %>%
  mutate(mod.match=ifelse(EM_ecov_how==as.numeric(substr(Ecov_how, 6,6)) & EM_r_mod==recruit_mod, 1, 0)) %>%
  group_by(   Fhist, Ecov_re_cor, Ecov_how, EM, mod.match ) %>%
  summarise( mean.RE.Ecov_process_pars1=mean(RE.Ecov_process_pars1), var.RE.Ecov_process_pars1=var(RE.Ecov_process_pars1), 
             median.RE.Ecov_process_pars1=median(RE.Ecov_process_pars1), across(RE.Ecov_process_pars1,                                                                                                       list(p2.5=~quantile(.,probs=0.0275 ) ,
                                                                                                                                                                                                                  p97.5=~quantile(.,probs=0.975)))
  )%>%
  ungroup() %>%
  left_join(em_tib)


#  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par_Ecov2_sum <- par_mat_tib %>%
  left_join(em_tib) %>%
  select(RE.Ecov_process_pars2, Ecov_re_cor, obs_error, R_sig, Ecov_re_cor, Fhist, EM, EM_ecov_how, EM_r_mod, Ecov_how, recruit_mod )  %>%
  mutate(mod.match=ifelse(EM_ecov_how==as.numeric(substr(Ecov_how, 6,6)) & EM_r_mod==recruit_mod, 1, 0)) %>%
  group_by(   Fhist, Ecov_re_cor, Ecov_how, EM, mod.match ) %>%
  summarise( mean.RE.Ecov_process_pars2=mean(RE.Ecov_process_pars2), var.RE.Ecov_process_pars=var(RE.Ecov_process_pars2), 
             median.RE.Ecov_process_pars2=median(RE.Ecov_process_pars2), across(RE.Ecov_process_pars2,                                                                                                       list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                                                                                                                                                  p97.5=~quantile(.,probs=0.975)))
  )%>%
  ungroup() %>%
  left_join(em_tib)


#  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par_Ecov3_sum <- par_mat_tib %>%
  left_join(em_tib) %>%
  select(RE.Ecov_process_pars3, Ecov_re_cor, obs_error, R_sig, Ecov_re_cor, Fhist, EM, EM_ecov_how, EM_r_mod, Ecov_how, recruit_mod )  %>%
  mutate(mod.match=ifelse(EM_ecov_how==as.numeric(substr(Ecov_how, 6,6)) & EM_r_mod==recruit_mod, 1, 0)) %>%
  group_by(  Fhist, Ecov_re_cor, Ecov_how, EM, mod.match ) %>%
  summarise(mean.RE.Ecov_process_pars3=mean(RE.Ecov_process_pars3), var.RE.Ecov_process_pars3=var(RE.Ecov_process_pars3), 
            median.RE.Ecov_process_pars3=median(RE.Ecov_process_pars3), across(RE.Ecov_process_pars3,                                                                                                       list(p2.5=~quantile(.,probs=0.0275 ),
                                                                                                                                                                                                                 p97.5=~quantile(.,probs=0.975)))
  )%>%
  ungroup() %>%
  left_join(em_tib)
#################################################
# RE Plots  ====

RE.mean_rec1.plot <- ggplot(par_mean_rec1_sum, aes(x=EM_mod, y=median.RE.mean_rec1, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist + R_sig   ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.mean_rec1_p2.5), ymax=(RE.mean_rec1_p97.5), width=.5) )+
  ylab('RE(Recr_par_a) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-1,5)) +
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
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in Beverton-Holt parameter a')
ggsave(RE.mean_rec1.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.mean_rec1.plot', plot.suffix, '.png') ),  height=7, width=12)

# same plot, with narrower ylim
RE.mean_rec1.plot.ylim <- ggplot(par_mean_rec1_sum, aes(x=EM_mod, y=median.RE.mean_rec1, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist + R_sig   ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.mean_rec1_p2.5), ymax=(RE.mean_rec1_p97.5), width=.5) )+
  ylab('RE(Recr_par_a) (2.5 - 50 - 97.5 percentiles)') +
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
  guides(col=guide_legend(title=NULL)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in Beverton-Holt parameter a')
ggsave(RE.mean_rec1.plot.ylim, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.mean_rec1.plot.ylim', plot.suffix, '.png') ),  height=7, width=12)


RE.mean_rec2.plot <- ggplot(par_mean_rec2_sum, aes(x=EM_mod, y=median.RE.mean_rec2, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist + R_sig   ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.mean_rec2_p2.5), ymax=(RE.mean_rec2_p97.5), width=.5) )+
  ylab('RE(Recr_par_b) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-1,5)) +
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
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in Beverton-Holt parameter b')
ggsave(RE.mean_rec2.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.mean_rec2.plot', plot.suffix, '.png') ),  height=7, width=12)

# same plot, with narrower ylim
RE.mean_rec2.plot.ylim <- ggplot(par_mean_rec2_sum, aes(x=EM_mod, y=median.RE.mean_rec2, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist + R_sig   ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.mean_rec2_p2.5), ymax=(RE.mean_rec2_p97.5), width=.5) )+
  ylab('RE(Recr_par_b) (2.5 - 50 - 97.5 percentiles)') +
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
  guides(col=guide_legend(title=NULL)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in Beverton-Holt parameter b')
ggsave(RE.mean_rec2.plot.ylim, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.mean_rec2.plot.ylim', plot.suffix, '.png') ),  height=7, width=12)


RE.R_sigma.plot <- ggplot(par_R_sigma_sum, aes(x=EM_mod, y=median.RE.NAA_sigma, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist + R_sig + NAA_cor  ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.NAA_sigma_p2.5), ymax=(RE.NAA_sigma_p97.5), width=.5) )+
  ylab('RE(R_sigma) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-1,3)) +
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
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in R_sigma')
ggsave(RE.R_sigma.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.R_sigma.plot', plot.suffix, '.png') ),  height=8, width=14)

# same plot, with narrower ylim
RE.R_sigma.plot.ylim <- ggplot(par_R_sigma_sum, aes(x=EM_mod, y=median.RE.NAA_sigma, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist + R_sig + NAA_cor  ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.NAA_sigma_p2.5), ymax=(RE.NAA_sigma_p97.5), width=.5) )+
  ylab('RE(R_sigma) (2.5 - 50 - 97.5 percentiles)') +
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
  guides(col=guide_legend(title=NULL)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in R_sigma')
ggsave(RE.R_sigma.plot.ylim, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.R_sigma.plot.ylim', plot.suffix, '.png') ),  height=8, width=14)



RE.R_cor.plot <- ggplot(par_R_cor_sum, aes(x=EM_mod, y=median.RE.NAA_rho, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist + R_sig + NAA_cor  ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.NAA_rho_p2.5), ymax=(RE.NAA_rho_p97.5), width=.5) )+
  ylab('RE(R_sigma) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-1,3)) +
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
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in R_rho')
ggsave(RE.R_cor.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.R_cor.plot', plot.suffix, '.png') ),  height=8, width=14)



# mean ecov is very precisely estimated; probably because  (?):
# unique(df.oms$Ecov_obs_sig)  is a fixed input (not estimated)
# [1] 0.1

RE.Ecov.mean.plot <- ggplot(par_Ecov1_sum, aes(x=EM_mod, y=median.RE.Ecov_process_pars1, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist +  Ecov_re_cor  ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.Ecov_process_pars1_p2.5), ymax=(RE.Ecov_process_pars1_p97.5), width=.5) )+
  ylab('RE(Ecov mean) (2.5 - 50 - 97.5 percentiles)') +
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
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in Ecov mean')
ggsave(RE.Ecov.mean.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.Ecov.mean.plot', plot.suffix, '.png') ),  height=8, width=14)




RE.Ecov.sigma.plot <- ggplot(par_Ecov2_sum, aes(x=EM_mod, y=median.RE.Ecov_process_pars2, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist +  Ecov_re_cor  ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.Ecov_process_pars2_p2.5), ymax=(RE.Ecov_process_pars2_p97.5), width=.5) )+
  ylab('RE(Ecov_sigma) (2.5 - 50 - 97.5 percentiles)') +
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
  guides(col=guide_legend(title=NULL)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in Ecov process sigma')
ggsave(RE.Ecov.sigma.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.Ecov.sigma.plot', plot.suffix, '.png') ),  height=8, width=14)



RE.Ecov.rho.plot <- ggplot(par_Ecov3_sum, aes(x=EM_mod, y=median.RE.Ecov_process_pars3, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist +  Ecov_re_cor  ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.Ecov_process_pars3_p2.5), ymax=(RE.Ecov_process_pars3_p97.5), width=.5) )+
  ylab('RE(Ecov_cor) (2.5 - 50 - 97.5 percentiles)') +
  coord_cartesian(ylim=c(-5,5)) +
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
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in Ecov process autocorrelation')
ggsave(RE.Ecov.rho.plot, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.Ecov.rho.plot', plot.suffix, '.png') ),  height=8, width=14)


# same plot, with narrower ylim
RE.Ecov.rho.plot.ylim <- ggplot(par_Ecov3_sum, aes(x=EM_mod, y=median.RE.Ecov_process_pars3, col=as.factor(mod.match) )) +
  facet_grid(Ecov_how ~ Fhist +  Ecov_re_cor  ) +
  geom_point(size=2.5) +
  geom_hline(yintercept=0.0, col='#111111dd', linewidth=0.5) +
  geom_errorbar(aes(ymin=(RE.Ecov_process_pars3_p2.5), ymax=(RE.Ecov_process_pars3_p97.5), width=.5) )+
  ylab('RE(Ecov_cor) (2.5 - 50 - 97.5 percentiles)') +
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
  guides(col=guide_legend(title=NULL)) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('RE in Ecov process autocorrelation')
ggsave(RE.Ecov.rho.plot.ylim, filename=file.path(here(),'Ecov_study','recruitment_functions', plot.dir, paste0('RE.Ecov.rho.plot.ylim', plot.suffix, '.png') ),  height=8, width=14)

# i think i need a plot for beta?
