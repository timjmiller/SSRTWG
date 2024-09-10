recr.df <- readRDS( file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "error.recr", plot.suffix, ".df.RDS") ) )
ssb.df  <- readRDS( file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "error.ssb",  plot.suffix, ".df.RDS") )  )
fbar.df <- readRDS( file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "error.fbar", plot.suffix, ".df.RDS") )  )



df.oms2 <- as_tibble(cbind(OM=seq(1,nrow(df.oms)), df.oms)) 
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


re.rec.tib$R_sig <- factor(re.rec.tib$R_sig,labels = c("Rsig_0.1","Rsig_0.5","Rsig_1.0"))
re.rec.tib$Ecov_effect <- factor(re.rec.tib$Ecov_effect,labels=c("Ecov_L","Ecov_H"))
re.rec.tib$Ecov_how    <- factor(re.rec.tib$Ecov_how,labels=c("Ecov_0", "Ecov_1","Ecov_2","Ecov_4"))
re.rec.tib$NAA_cor     <- factor(re.rec.tib$NAA_cor,labels=c("Rcor_L","Rcor_H"))
re.rec.tib$Fhist       <- factor(re.rec.tib$Fhist,labels=c("H-MSY","MSY") ) 
re.rec.tib$Ecov_re_cor <- factor(re.rec.tib$Ecov_re_cor,labels=c("Ecor_L","Ecor_H"))
re.rec.tib$obs_error <- factor(re.rec.tib$obs_error,labels=c("Obs_H","Obs_L"))


re.ssb.tib$R_sig <- factor(re.ssb.tib$R_sig,labels = c("Rsig_0.1","Rsig_0.5","Rsig_1.0"))
re.ssb.tib$Ecov_effect <- factor(re.ssb.tib$Ecov_effect,labels=c("Ecov_L","Ecov_H"))
re.ssb.tib$Ecov_how    <- factor(re.ssb.tib$Ecov_how,labels=c("Ecov_0", "Ecov_1","Ecov_2","Ecov_4"))
re.ssb.tib$NAA_cor     <- factor(re.ssb.tib$NAA_cor,labels=c("Rcor_L","Rcor_H"))
re.ssb.tib$Fhist       <- factor(re.ssb.tib$Fhist,labels=c("H-MSY","MSY") ) 
re.ssb.tib$Ecov_re_cor <- factor(re.ssb.tib$Ecov_re_cor,labels=c("Ecor_L","Ecor_H"))
re.ssb.tib$obs_error <- factor(re.ssb.tib$obs_error,labels=c("Obs_H","Obs_L"))


re.fbar.tib$R_sig <- factor(re.fbar.tib$R_sig,labels = c("Rsig_0.1","Rsig_0.5","Rsig_1.0"))
re.fbar.tib$Ecov_effect <- factor(re.fbar.tib$Ecov_effect,labels=c("Ecov_L","Ecov_H"))
re.fbar.tib$Ecov_how    <- factor(re.fbar.tib$Ecov_how,labels=c("Ecov_0", "Ecov_1","Ecov_2","Ecov_4"))
re.fbar.tib$NAA_cor     <- factor(re.fbar.tib$NAA_cor,labels=c("Rcor_L","Rcor_H"))
re.fbar.tib$Fhist       <- factor(re.fbar.tib$Fhist,labels=c("H-MSY","MSY") ) 
re.fbar.tib$Ecov_re_cor <- factor(re.fbar.tib$Ecov_re_cor,labels=c("Ecor_L","Ecor_H"))
re.fbar.tib$obs_error <- factor(re.fbar.tib$obs_error,labels=c("Obs_H","Obs_L"))


####################################################################################
#  regression trees to find nodes for bias  ====
####################################################################################

# recr ====

rf_recr_all.yrs <- rpart(RE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how +obs_error, data=re.rec.tib, control=rpart.control(cp=0.01))
imp.var <- rf_recr_all.yrs$frame[rf_recr_all.yrs$frame$var != '<leaf>',]
nodes_recr_all.yrs <- unique(imp.var[,1])
nodes_recr_all.yrs
# nothing  # "obs_error" (liz 384)


rf_recr_last10.yrs <- rpart(RE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how +obs_error, data=re.rec.tib[re.rec.tib$Year>30,], control=rpart.control(cp=0.01))
imp.var <- rf_recr_last10.yrs$frame[rf_recr_last10.yrs$frame$var != '<leaf>',]
nodes_recr_last10.yrs <- unique(imp.var[,1])
nodes_recr_last10.yrs
# "R_sig"  # "obs_error" "R_sig" (liz 384)

rf_recr_final.yrs <- rpart(RE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how +obs_error, data=re.rec.tib[re.rec.tib$Year>39,], control=rpart.control(cp=0.01))
imp.var <- rf_recr_final.yrs$frame[rf_recr_final.yrs$frame$var != '<leaf>',]
nodes_recr_final.yrs <- unique(imp.var[,1])
nodes_recr_final.yrs
# nothing   # "obs_error" "R_sig" (liz 384)


# ssb ====

rf_ssb_all.yrs <- rpart(RE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how +obs_error, data=re.ssb.tib, control=rpart.control(cp=0.01))
imp.var <- rf_ssb_all.yrs$frame[rf_ssb_all.yrs$frame$var != '<leaf>',]
nodes_ssb_all.yrs <- unique(imp.var[,1])
nodes_ssb_all.yrs
# nothing  # (liz 384)


rf_ssb_last10.yrs <- rpart(RE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how +obs_error, data=re.ssb.tib[re.rec.tib$Year>30,], control=rpart.control(cp=0.01))
imp.var <- rf_ssb_last10.yrs$frame[rf_ssb_last10.yrs$frame$var != '<leaf>',]
nodes_ssb_last10.yrs <- unique(imp.var[,1])
nodes_ssb_last10.yrs
# nothing   # nothing (liz 384)

rf_ssb_final.yrs <- rpart(RE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how +obs_error, data=re.ssb.tib[re.rec.tib$Year>39,], control=rpart.control(cp=0.01))
imp.var <- rf_ssb_final.yrs$frame[rf_ssb_final.yrs$frame$var != '<leaf>',]
nodes_ssb_final.yrs <- unique(imp.var[,1])
nodes_ssb_final.yrs
# nothing  #  nothing  (liz 384)


# fbar ====

rf_fbar_all.yrs <- rpart(RE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how +obs_error, data=re.fbar.tib, control=rpart.control(cp=0.01))
imp.var <- rf_fbar_all.yrs$frame[rf_fbar_all.yrs$frame$var != '<leaf>',]
nodes_fbar_all.yrs <- unique(imp.var[,1])
nodes_fbar_all.yrs
# nothing   #  nothing  (liz 384)


rf_fbar_last10.yrs <- rpart(RE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how +obs_error, data=re.fbar.tib[re.rec.tib$Year>30,], control=rpart.control(cp=0.01))
imp.var <- rf_fbar_last10.yrs$frame[rf_fbar_last10.yrs$frame$var != '<leaf>',]
nodes_fbar_last10.yrs <- unique(imp.var[,1])
nodes_fbar_last10.yrs
# nothing   #  nothing  (liz 384)

rf_fbar_final.yrs <- rpart(RE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how +obs_error, data=re.fbar.tib[re.rec.tib$Year>39,], control=rpart.control(cp=0.01))
imp.var <- rf_fbar_final.yrs$frame[rf_fbar_final.yrs$frame$var != '<leaf>',]
nodes_fbar_final.yrs <- unique(imp.var[,1])
nodes_fbar_final.yrs
# nothing   #  nothing  (liz 384)


#################################################
# summarize  all years  ====

re.rec.sum <- re.rec.tib %>%
  group_by( obs_error, R_sig, Ecov_how, EM, mod.match ) %>%
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
  group_by( obs_error, R_sig,  Ecov_how, EM, mod.match ) %>%
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
  group_by(obs_error, R_sig, Ecov_how, EM, mod.match ) %>%
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
  facet_grid(Ecov_how ~  obs_error +  R_sig  ) +
  # facet_grid(Ecov_how ~    R_sig  ) +
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
  facet_grid(Ecov_how ~ obs_error +   R_sig  ) +
  # facet_grid(Ecov_how ~    R_sig  ) +
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
  facet_grid(Ecov_how ~ obs_error +  R_sig  ) +
  # facet_grid(Ecov_how ~    R_sig  ) +
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



