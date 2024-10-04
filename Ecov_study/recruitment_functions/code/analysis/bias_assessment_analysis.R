library(here)
library(tidyverse)
library(rpart)
library(rpart.plot)

res.path    <- file.path(here::here(),"Ecov_study", "recruitment_functions", "results")  # directory where simulation 
res.dir     <- 'results'   # 'results'     'results_beta_fix'   # results folder where AIC dataframes are
plot.dir    <- 'plots'    # 'plots_lizruns'  'plots_beta_fix'  
table.dir   <- 'tables'   # 'results'     'results_beta_fix'
plot.suffix <- ''      # '_beta_fix'   '' 

df.oms    <- readRDS(file.path(here::here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS"))
df.ems    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))

bad.grad.value <- 1E-6 #tim used 1E-6 (the abs of the exponent will be used for output suffix; ex: filename_grad_6.png)
bad.grad.label <- as.numeric(strsplit(as.character(bad.grad.value), split="-")[[1]][2])
bad.se.value   <- 100 #tim used 100 (this value will be used for output suffix; ex: filename_se_100.png)


recr.df <- as.data.frame(readRDS( file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "error.recr", plot.suffix, ".df.RDS") ) )) %>%
  mutate(ecov_how=as.factor(ecov_how)) 

ssb.df  <- as.data.frame(readRDS( file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "error.ssb",  plot.suffix, ".df.RDS") )  )) %>%
  mutate(ecov_how=as.factor(ecov_how)) 

fbar.df <- as.data.frame(readRDS( file.path( here::here(),'Ecov_study','recruitment_functions',res.dir , paste0( "error.fbar", plot.suffix, ".df.RDS") )  )) %>%
  mutate(ecov_how=as.factor(ecov_how)) 

####################################################################################
#  regression trees to find nodes for bias  ====
####################################################################################
# recr ====
rf_recr_re_all   <- rpart(RE   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + ecov_how, 
                        data=recr.df, control=rpart.control(cp=0.001))
rf_recr_rmse_all <- rpart(RMSE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + ecov_how, 
                        data=recr.df, control=rpart.control(cp=0.001))


rf_recr_re_ten   <- rpart(RE   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + ecov_how, 
                          data=recr.df, control=rpart.control(cp=0.01))
rf_recr_rmse_ten <- rpart(RMSE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + ecov_how, 
                          data=error.recr.df, control=rpart.control(cp=0.01))
rf_recr_re_all   <- rpart(RE   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + ecov_how, 
                          data=error.recr.df, control=rpart.control(cp=0.01))
rf_recr_rmse_all <- rpart(RMSE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + ecov_how, 
                          data=error.recr.df, control=rpart.control(cp=0.01))

# ssb ====
rf_ssb_re_all   <- rpart(RE   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + ecov_how, 
                        data=ssb.df, control=rpart.control(cp=0.01))
rf_ssb_rmse_all <- rpart(RMSE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + ecov_how, 
                        data=recr.df, control=rpart.control(cp=0.01))
rf_recr_re_ten   <- rpart(RE   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + ecov_how, 
                          data=recr.df, control=rpart.control(cp=0.01))
rf_recr_rmse_ten <- rpart(RMSE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how, 
                          data=error.recr.df, control=rpart.control(cp=0.01))
rf_recr_re_all   <- rpart(RE   ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how, 
                          data=error.recr.df, control=rpart.control(cp=0.01))
rf_recr_rmse_all <- rpart(RMSE ~ R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how, 
                          data=error.recr.df, control=rpart.control(cp=0.01))


lm_recr_re_all <- lm(RE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + 
                      EM_ecov_how + ecov_how + r_mod,
                      data=rrecr.df)

lm_recr_rmse_all <- lm(RMSE ~ ssb_cv + ecov_slope + obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + 
                      EM_ecov_how + ecov_how + r_mod,
                      data=rrecr.df)


nn <- nrow(recr.df)
ii <- sample(1:nn,size=10000)

rrecr.df <- recr.df[ii,]





glm_SR   <- glm(correct_SR   ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                data=AIC_best, family='binomial')
glm_form <- glm(correct_form ~ obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                data=AIC_best, family='binomial')

glmm_ecov <- glmer(correct_ecov ~  (1|sim) + obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                   data=AIC_best, family='binomial')
glmm_SR   <- glmer(correct_SR   ~  (1|sim) + obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                   data=AIC_best, family='binomial')
glmm_form <- glmer(correct_form ~  (1|sim) + obs_error + R_sig + Fhist + NAA_cor + Ecov_re_cor + Ecov_effect + Ecov_how + ecov_slope + ssb_cv, 
                   data=AIC_best, family='binomial')

FITS_AIC <- list(glm_ecov,glm_SR,glm_form,
                 glmm_ecov,glmm_SR,glmm_form)
saveRDS(FITS_AIC,file=file.path(here::here(),'Ecov_study','recruitment_functions','results','FITS_AIC.rds'))




pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','trees_error.pdf'),height=5,width=8)
par(mfrow=c(3,2), oma=c(5,0,0,0))

prp(rf_recr_re_all,yesno=FALSE,type=4,clip.right.labs=TRUE)
  mtext('a) recr_re_all',adj=0, line=-1, cex=0.9)
prp(rf_recr_rmse_all,yesno=FALSE,type=4,clip.right.labs=TRUE)
  mtext(expression('b) recr_rmse_all'),adj=0,line=-1, cex=0.9)
  #title(sub=paste0(100*round(aic.best.pct.conv,2), "% of runs converged (", aic.best.n.bad, " out of ", aic.best.n, " failed)" ),   adj=0.5, outer=TRUE)
#prp(tree_form,yesno=FALSE,type=4,clip.right.labs=TRUE)
#  mtext(expression('c) SRR and'~'E'['cov']~'Y/N'),adj=0,line=-1, cex=0.9)
dev.off()


##--MAKE PLOTS--###########################
labels <- c(#expression(italic('intercept')),
            expression(sigma['obs']~'= H'),
            expression(sigma['r']~'= 0.3'),
            expression(sigma['r']~'= 0.5'),
            expression(sigma['r']~'= 1.0'),
            expression(italic('F')~'= L-H'),
            expression(italic('F')~'= H-MSY'),
            expression(rho['NAA']~' = H'), 
            expression(rho['E'['cov']]~' = H'), 
            expression(beta['E'['cov']]~' = H'), 
            expression(italic('f(E'['cov']*')')~'= 1'), 
            expression(italic('f(E'['cov']*')')~'= 2'), 
            expression(Delta*'E'['cov']*' = M'), 
            expression(Delta*'E'['cov']*' = H'), 
            expression(italic('CV'['SSB']~'= M')),
            expression(italic('CV'['SSB']~'= H')))

pdf(file.path(here::here(),'Ecov_study','recruitment_functions','plots','effects.pdf'),height=4.5,width=4.5)
par(mfrow=c(3,1),mar=c(1,2,1,2),oma=c(6,2,2,2),cex.axis=0.9)
ylims <- c(-5,5)
plotlm(glm_SR,add=FALSE,ylim=ylims,labels=rep(NA,length(labels)))
plotlm(glmm_SR,add=TRUE,ylim=ylims,labels=rep(NA,length(labels)))
  mtext('a) SR Y/N',adj=0.0)
plotlm(glm_ecov,add=FALSE,ylim=ylims,labels=rep(NA,length(labels)))
plotlm(glmm_ecov,add=TRUE,ylim=ylims,labels=rep(NA,length(labels)))
  mtext(expression('b) E'['cov']~'Y/N'),adj=0.0)
plotlm(glm_form,add=FALSE,ylim=ylims,labels=rep(NA,length(labels)))
plotlm(glmm_form,add=TRUE,ylim=ylims,labels=labels)
  mtext(expression('c) SR & E'['cov']~'Y/N'),adj=0.0)
mtext(outer=TRUE,expression(Delta*'log['~italic('p/(1-p)')*']'),side=2)
dev.off()



























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



