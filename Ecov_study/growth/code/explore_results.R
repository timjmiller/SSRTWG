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
      print(fit$model[1,1:3])
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
  lapply(fits, function(i) ff(i)) %>% bind_rows()
}
get_selex <- function(fits){
  ff <- function(fit){
    if(!fit$model$optimized){
      print(fit$model)
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
  lapply(fits, function(i) ff(i)) %>% bind_rows()
}
get_waa <- function(fits){
  ff <- function(fit){
    if(!fit$model$optimized){
      print(fit$model)
      return(NULL)
    }
    waa <- data.frame(par='waa',
                      age=1:10,
                      est=fit$fit$rep$pred_waa[5,1,],
                      truth=fit$truth$pred_waa[5,1,]) %>% bind_cols(fit$model) %>%
      mutate(rel_error=(est-truth)/truth, abs_error=est-truth,
             sim=as.factor(im),  maxgrad=get_maxgrad(fit))
    return(waa)
  }
  lapply(fits, function(i) ff(i)) %>% bind_rows()
}
get_pars <- function(fits){
  ff <- function(fit){
    if(!fit$model$optimized){
      print(fit$model)
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
  lapply(fits, function(i) ff(i)) %>% bind_rows()
}
get_growth <- function(fits){
  ff <- function(fit){
    if(!fit$model$optimized){
      print(fit$model)
      return(NULL)
    }
    ## !!will need to update this if the growth estimation changes!!
    p1 <- data.frame(par=c('logk', 'logLinf'),
                     est=fit$fit$rep$growth_a[1:2,1],
                     truth=fit$truth$growth_a[1:2,1])
    p2 <- data.frame(par=c('SDold'),
                     est=fit$empars$value[fit$empars$par=='SDgrowth_par'],
                     truth=fit$ompars$value[fit$ompars$par=='SDgrowth_par'][2])
    ## p2 <- data.frame(par=c('SD1', 'SD2'),
    ##                  est=fit$fit$rep$expSD,
    ##                  truth=fit$truth$expSD)
    growth <- bind_rows(p1,p2) %>% bind_cols(fit$model) %>%
      mutate(rel_error=(est-truth)/truth, abs_error=est-truth,
      sim=as.factor(im),  maxgrad=get_maxgrad(fit))
    return(growth)
  }
  lapply(fits, function(i) ff(i)) %>% bind_rows()
}


fits <- list.files('../results', pattern='RDS', recursive=1,
                      full.names=TRUE) %>%
  lapply(function(i) readRDS(i))
## check convergence stats
models <- lapply(fits, function(x) x$model) %>% bind_rows
group_by(models, em, om) %>%
  summarize(pct.converged=mean(optimized), n.converged=sum(optimized))
## table(models$optimized,models$em, models$om)
## table(models$sdreport,models$em)

## fit <- fits[[2]]

ts <- get_ts(fits) %>% filter(maxgrad<.1)
ggplot(filter(ts,par=='SSB'), aes(year, truth, group=sim, color=factor(om_Ecov_effect))) +
  geom_line() + facet_grid(om~em, scales='free') + labs(y='SSB truth')
ggplot(filter(ts,par=='Ecov_out'), aes(year, truth, group=sim, color=factor(om_Ecov_effect))) +
  geom_line() + facet_grid(om~em, scales='free') + labs(y='Ecov_out truth')
ggplot(filter(ts,par=='Ecov_out'), aes(year, est, group=sim, color=factor(om_Ecov_effect))) +
  geom_line() + facet_grid(om~em, scales='free') + labs(y='Ecov_out est')
ggplot(filter(ts, par!='Ecov_out'), aes(year, rel_error, group=sim)) + geom_line(alpha=.5) +# ylim(-3,3)+
  geom_hline(yintercept=0, col=2, lwd=1) +
  facet_grid(par~om_Ecov_effect+em_Ecov_est, scales='free')

pars <- get_pars(fits)
ggplot(pars, aes(par2, rel_error)) + geom_violin() +
  geom_hline(yintercept=0, col=2, lwd=1) + #ylim(-3,3)+
  facet_grid(om_Ecov_effect+em_Ecov_est~par, scales='free')
ecov <- filter(pars, grepl(x=par, pattern='Ecov'))
ggplot(ecov, aes(par2, abs_error)) + geom_violin() +
  geom_hline(yintercept=0, col=2, lwd=1) + #ylim(-3,3)+
  facet_grid(om_Ecov_effect+em_Ecov_est~par, scales='free')

## test <- filter(ts, par=='SSB' & sim==1 & em==2)
## plot(test$year, test$truth)
## lines(test$year, test$est)


selex <- get_selex(fits)
ggplot(selex, aes(age, est, group=sim)) +
  geom_line() + facet_grid(par~em_Ecov_est+om_Ecov_effect)
ggplot(selex, aes(age, abs_error, group=sim)) + geom_line(alpha=.5) +
  geom_hline(yintercept=0, col=2, lwd=1) + facet_grid(par~em_Ecov_est+om_Ecov_effect)

growth <- get_growth(fits)
ggplot(growth, aes(par, abs_error)) + geom_violin() +
  facet_grid(.~em_growth_est) + geom_hline(yintercept=0, col=2)
ggplot(growth, aes(par, rel_error)) + geom_violin() +
  facet_grid(.~em_growth_est) + geom_hline(yintercept=0, col=2)


## need to get this by year too
waa <- get_waa(fits)
ggplot(waa, aes(age, rel_error, group=sim)) + geom_line() +
  facet_grid(par~em_Ecov_est+om_Ecov_effect)+
  geom_hline(yintercept=0, col=2) + labs(y='rel_error WAA')
ggplot(waa, aes(age, abs_error, group=sim)) + geom_line() +
  facet_grid(par~em_Ecov_est+om_Ecov_effect)+
  geom_hline(yintercept=0, col=2) + labs(y='abs_error WAA')

