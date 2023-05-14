## Cole made this custom file to quickly look at results for
## R&D. Note the script to run things has modified code w/ what
## to return

library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
theme_set(theme_bw())
get_maxgrad <- function(fit){
  if(!fit$model$optimized){
    return(NULL)
  }
  maxgrad <- max(abs(fit$fit$final_gradient))
  return(maxgrad)
}
get_ts <- function(fits){
  ff <- function(fit){
    if(!fit$model$optimized){
     ##  print(fit$model[1,1:3])
      return(NULL)
    }
    ssb <- data.frame(par='SSB',
                      year=1:40,
                      est=fit$fit$rep$SSB,
                      truth=fit$truth$SSB)
   ## if(fit$model$sdreport) ssb$se <- fit$sdrep$SE_rep$log_SSB
    recruits <- data.frame(par='recruits',
                      year=1:40,
                      est=fit$fit$rep$NAA[,1],
                      truth=fit$truth$NAA[,1])
   ## if(fit$model$sdreport) ssb$se <- fit$sdrep$SE_rep$log_SSB
    ## this fails w/ cbind for some reason??
    f <- data.frame(par='F',
                    year=1:40,
                    est=fit$fit$rep$F,
                    truth=fit$truth$F)
    ## if(fit$model$sdreport) f$se <- fit$sdrep$SE_rep$F
    ## this fails w/ cbind for some reason??
    ecov <- data.frame(par='Ecov_out',
                       year=1:40,
                       ## take first one -- all the same
                       est=fit$fit$rep$Ecov_out[,1,1],
                       truth=fit$truth$Ecov_out[,1,1])
    ts <- bind_rows(ssb,f, recruits, ecov) %>% bind_cols(fit$model) %>%
      mutate(rel_error=(est-truth)/truth, abs_error=est-truth,
      sim=as.factor(im), maxgrad=get_maxgrad(fit))
    return(ts)
  }
  lapply(fits, function(i) ff(i)) %>% bind_rows() %>% add_labels
}
get_selex <- function(fits){
  ff <- function(fit){
    if(!fit$model$optimized){
     ##  print(fit$model)
      return(NULL)
    }
    s1 <- data.frame(par='selex_fishery',
                     age=1:10,
                     est=fit$fit$rep$selAA[[1]][1,],
                     truth=fit$truth$selAA[[1]][1,])
    s2 <- data.frame(par='selex_survey1',
                     age=1:10,
                     est=fit$fit$rep$selAA[[2]][1,],
                     truth=fit$truth$selAA[[2]][1,])
    s3 <- data.frame(par='selex_survey2',
                     age=1:50, # not really age but whatever
                     est=fit$fit$rep$selAL[[3]][1,],
                     truth=fit$truth$selAL[[3]][1,])
    selex <- bind_rows(s1,s2,s3) %>% bind_cols(fit$model) %>%
      mutate(rel_error=(est-truth)/truth, abs_error=est-truth,
             sim=as.factor(im), maxgrad=get_maxgrad(fit))
    return(selex)
  }
  lapply(fits, function(i) ff(i)) %>% bind_rows()%>% add_labels
}
## get_waa <- function(fits){
##   ff <- function(fit){
##     if(!fit$model$optimized){
##       ## print(fit$model)
##       return(NULL)
##     }
##     waa <- data.frame(par='waa',
##                       age=1:10,
##                       est=fit$fit$rep$pred_waa[5,40,],
##                       truth=fit$truth$pred_waa[5,40,]) %>% bind_cols(fit$model) %>%
##       mutate(rel_error=(est-truth)/truth, abs_error=est-truth,
##              sim=as.factor(im),  maxgrad=get_maxgrad(fit))
##     return(waa)
##   }
##   lapply(fits, function(i) ff(i)) %>% bind_rows() %>% add_labels
## }
get_laa <- function(fits){
  ff <- function(fit){
    if(!fit$model$optimized){
      ## print(fit$model)
      return(NULL)
    }
    laa <- list()
    for(year in seq_len(nrow(fit$fit$rep$LAA))){
      laa[[year]] <- data.frame(par='laa', age=1:10, year=year,
                                est=fit$fit$rep$LAA[year,],
                                truth=fit$truth$LAA[year,]) %>%
        bind_cols(fit$model)
    }
    laa <- laa %>% bind_rows() %>%
        mutate(rel_error=(est-truth)/truth, abs_error=est-truth,
               sim=as.factor(im),  maxgrad=get_maxgrad(fit))
    return(laa)
  }
  lapply(fits, function(i) ff(i)) %>% bind_rows() %>% add_labels
}
get_waa <- function(fits){
  ff <- function(fit){
    if(!fit$model$optimized){
      ## print(fit$model)
      return(NULL)
    }
    waa <- list()
    for(year in seq_len(nrow(fit$fit$rep$pred_waa[5,,]))){
      waa[[year]] <- data.frame(par='waa', age=1:10, year=year,
                                est=fit$fit$rep$pred_waa[5,year,],
                                truth=fit$truth$pred_waa[5,year,]) %>%
        bind_cols(fit$model)
    }
    waa <- waa %>% bind_rows() %>%
        mutate(rel_error=(est-truth)/truth, abs_error=est-truth,
               sim=as.factor(im),  maxgrad=get_maxgrad(fit))
    return(waa)
  }
  lapply(fits, function(i) ff(i)) %>% bind_rows() %>% add_labels
}
get_pars <- function(fits){
  ff <- function(fit){
    if(!fit$model$optimized){
      ## print(fit$model)
      return(NULL)
    }
    pars <- merge(fit$ompars, fit$empars, by='par2') %>%
      filter(grepl(x=par.y, "Ecov|growth_a|SD_par|mean_rec_pars|logit_q|log_F1|log_N1_pars|logit_selpars|catch_paa_pars|index_paa_pars"))
    pars <- pars %>% select(par=par.x, par2, truth=value.x, est=value.y) %>%
      bind_cols(fit$model) %>%
      mutate(rel_error=(est-truth)/truth, abs_error=est-truth,
             sim=as.factor(im),  maxgrad=get_maxgrad(fit))
    return(pars)
  }
  lapply(fits, function(i) ff(i)) %>% bind_rows() %>% add_labels
}
get_growth <- function(fits){
  ff <- function(fit){
    if(!fit$model$optimized){
     ##  print(fit$model)
      return(NULL)
    }
    ## !!will need to update this if the growth estimation changes!!
    p1 <- data.frame(par=c('k', 'L1', 'Linf'),
                     est=exp(fit$fit$rep$growth_a[1:3,1]),
                     truth=exp(fit$truth$growth_a[1:3,1]))
    ## p2 <- data.frame(par=c('SDold'),
    ##                  est=fit$empars$value[fit$empars$par=='SDgrowth_par'],
    ##                  truth=fit$ompars$value[fit$ompars$par=='SDgrowth_par'][2])
    ## p2 <- data.frame(par=c('SD1', 'SD2'),
    ##                  est=fit$fit$rep$expSD,
    ##                  truth=fit$truth$expSD)
    p2 <- NULL
    growth <- bind_rows(p1,p2) %>% bind_cols(fit$model) %>%
      mutate(rel_error=(est-truth)/truth, abs_error=est-truth,
      sim=as.factor(im),  maxgrad=get_maxgrad(fit))
    return(growth)
  }
  warning("need to fix growth pars")
  lapply(fits, function(i) ff(i)) %>% bind_rows() %>% add_labels
}
add_labels <- function(df){
  emlabs <- c('OM: Ecov_beta=0.0', 'OM: Ecov_beta=0.5')
  omlabs <- c('EM: Ecov on L1', 'EM: Constant',
              'EM: AR(1) on L1', 'EM: LAA constant',
              'EM: LAA 2d-AR(1)')
  omf <- function(om) factor(om, 1:2, emlabs)
  emf <- function(em) factor(em, 1:5, omlabs)
  mutate(df, omf=omf(om), emf=emf(em))
}

fits <- list.files('../results', pattern='RDS', recursive=1,
                      full.names=TRUE) %>% lapply(readRDS)
saveRDS(fits, '../results/fits.RDS')

fits <- readRDS("../results/fits.RDS")

## check convergence stats
models <- lapply(fits, function(x) x$model) %>% bind_rows %>% add_labels
group_by(models, omf, emf) %>%
  summarize(pct.converged=mean(optimized), n.converged=sum(optimized), n.run=length(optimized))

## Quick exploration via plots

ts <- get_ts(fits)
ggplot(filter(ts,par=='SSB'), aes(year, truth, group=sim, color=abs(maxgrad)>1)) +
  geom_line() + facet_grid(omf~emf, scales='free') + labs(y='SSB truth')
ggplot(filter(ts,par=='Ecov_out'), aes(year, truth, group=sim)) +
  geom_line() + facet_grid(omf~emf, scales='free') + labs(y='Ecov_out truth')
ggplot(filter(ts,par=='Ecov_out'), aes(year, est, group=sim)) +
  geom_line() + facet_grid(omf~emf, scales='free') + labs(y='Ecov_out est')
g <- ggplot(filter(ts, par!='Ecov_out'), aes(year, rel_error, group=sim, color=abs(maxgrad)>1)) +
  geom_hline(yintercept=0, col=2, lwd=1) +
  geom_line(alpha=.5) + coord_cartesian(ylim=c(-1,1))+
  facet_grid(par~omf+emf, scales='free')
ggsave("../plots/relerror_ts_by_year.png", g, width=15 , height=5)
g <- ggplot(filter(ts, par!='Ecov_out' & abs(maxgrad)<1), aes(emf, rel_error)) +
  geom_hline(yintercept=0, col=2, lwd=1) +
  geom_violin() +# coord_cartesian(ylim=c(-1,1))+
  facet_grid(par~omf, scales='free') +
  theme(axis.text.x = element_text(angle = 90)) + labs(x=NULL)
ggsave("../plots/relerror_ts.png", g, width=7 , height=5)

# pars <- get_pars(fits) %>% filter(abs(maxgrad)<1)
g <- ggplot(pars, aes(par2, rel_error)) +
  geom_hline(yintercept=0, col=2, lwd=1) + #ylim(-3,3)+
  geom_violin() +
    geom_jitter(width=.3, height=0, alpha=.3) +
  #geom_jitter(width=.1, height=0, alpha=.5, aes(color=abs(maxgrad)>1)) +
  coord_cartesian(ylim=c(-1,1))+
  facet_grid(omf+emf~par, scales='free') +
  theme(axis.text.x = element_text(angle = 90)) + labs(x=NULL) +
  theme(panel.spacing = unit(0, "cm"))
g
ggsave("../plots/relerror_pars.png", g, width=10 , height=10)

ecov <- filter(pars, grepl(x=par, pattern='Ecov'))
g <- ggplot(ecov, aes(emf, rel_error)) +
  geom_hline(yintercept=0, col=2, lwd=1) + #ylim(-3,3)+
  geom_violin() +
  geom_jitter(width=.3, height=0, alpha=.3) +
  coord_cartesian(ylim=c(-1,1))+
  facet_grid(omf~par2, scales='free')+
  theme(axis.text.x = element_text(angle = 90)) + labs(x=NULL)
## g <- ggplot(ecov, aes(par2, rel_error)) +
##   geom_hline(yintercept=0, col=2, lwd=1) + #ylim(-3,3)+
##   geom_violin() +
##   coord_cartesian(ylim=c(-1,1))+
##   facet_grid(omf+emf~par, scales='free')
ggsave("../plots/relerror_pars_ecov.png", g, width=6 , height=4)

growthpars <- filter(pars, grepl(x=par, pattern='growth_a'))
g <- ggplot(growthpars, aes(emf, rel_error)) +
  geom_hline(yintercept=0, col=2, lwd=1) + #ylim(-3,3)+
  geom_violin() +
  geom_jitter(width=.3, height=0, alpha=.3) +
  coord_cartesian(ylim=c(-1,1))+
  facet_grid(omf~par2, scales='free')+
  theme(axis.text.x = element_text(angle = 90)) + labs(x=NULL)
ggsave("../plots/relerror_pars_growth.png", g, width=5 , height=4)
## test <- filter(ts, par=='SSB' & sim==1 & em==2)
## plot(test$year, test$truth)
## lines(test$year, test$est)


selex <- get_selex(fits) %>% filter(abs(maxgrad)<1)
ggplot(selex, aes(age, est, group=sim)) +
  geom_line() + facet_grid(par~emf+om)
g <- ggplot(selex, aes(age, abs_error, group=sim)) + geom_line(alpha=.5) +
  geom_hline(yintercept=0, col=2, lwd=1) +
  facet_grid(par~omf+emf)
ggsave("../plots/abserror_selex.png", g, width=12 , height=5)

growth <- get_growth(fits) %>% filter(abs(maxgrad)<1)
## ggplot(growth, aes(par, abs_error)) + geom_violin() +
##   facet_grid(.~em_growth_est) + geom_hline(yintercept=0, col=2)
g <- ggplot(growth, aes(emf, rel_error)) +
  geom_hline(yintercept=0, col=2) +
  geom_violin() +
  geom_jitter(width=.3, height=0, alpha=.3) +
  facet_grid(par~omf) +
  coord_cartesian(ylim=c(-1,1))+
  theme(axis.text.x = element_text(angle = 90)) + labs(x=NULL)
ggsave("../plots/abserror_growth.png", g, width=7 , height=5)


## ## need to get this by year too
## waa <- get_waa(fits) %>% filter(abs(maxgrad)<1)
## g <- ggplot(waa, aes(age, rel_error, group=sim)) + geom_line() +
##   facet_grid(omf~emf)+
##   geom_hline(yintercept=0, col=2) + labs(y='rel_error WAA')
## g <- ggplot(waa, aes(age, abs_error, group=sim)) + geom_line() +
##   facet_grid(omf~emf)+
##   geom_hline(yintercept=0, col=2) + labs(y='abs_error WAA')
## ggsave("../plots/abserror_waa.png", g, width=7 , height=5)

waa <- get_waa(fits) %>% filter(abs(maxgrad)<1)
tmp <- filter(waa, year==40) # terminal year
g <- ggplot(tmp, aes(age, rel_error, group=sim)) + geom_line() +
  facet_grid(omf~emf)+
  geom_hline(yintercept=0, col=2) + labs(y='rel_error WAA')
g <- ggplot(tmp, aes(age, abs_error, group=sim)) + geom_line() +
  facet_grid(omf~emf)+
  geom_hline(yintercept=0, col=2) + labs(y='abs_error WAA')
ggsave("../plots/abserror_waa_terminal_year.png", g, width=7 , height=5)
## by year and age, mre= median relative error
tmp <- waa %>% group_by(omf, emf, age, year) %>% summarize(n=n(), mre=median(rel_error))
g <- ggplot(tmp, aes(year, age, size=abs(mre), color=mre>0)) +
  geom_point(alpha=.5) +
  facet_grid(omf~emf)
ggsave("../plots/waa_year_by_age.png", g, width=10 , height=5)

laa <- get_laa(fits) %>% filter(abs(maxgrad)<1)
tmp <- filter(laa, year==40) # terminal year
g <- ggplot(tmp, aes(age, rel_error, group=sim)) + geom_line() +
  facet_grid(omf~emf)+
  geom_hline(yintercept=0, col=2) + labs(y='rel_error LAA')
g <- ggplot(tmp, aes(age, abs_error, group=sim)) + geom_line() +
  facet_grid(omf~emf)+
  geom_hline(yintercept=0, col=2) + labs(y='abs_error LAA')
ggsave("../plots/abserror_laa_terminal_year.png", g, width=7 , height=5)
## by year and age, mre= median relative error
tmp <- laa %>% group_by(omf, emf, age, year) %>% summarize(n=n(), mre=median(rel_error))
g <- ggplot(tmp, aes(year, age, size=abs(mre), color=mre>0)) +
  geom_point(alpha=.5) +
  facet_grid(omf~emf)
ggsave("../plots/laa_year_by_age.png", g, width=10 , height=5)
