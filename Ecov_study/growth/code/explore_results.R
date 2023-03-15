library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
theme_set(theme_bw())
get_ts <- function(fits){
  ff <- function(fit){
    if(!fit$model$optimized){
      print(fit$model)
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
    ts <- bind_rows(ssb,f, recruits) %>% bind_cols(fit$model) %>%
      mutate(rel_error=(est-truth)/truth, abs_error=est-truth, sim=as.factor(im))
    return(ts)
  }
  lapply(fits, function(i) ff(i)) %>% bind_rows()
}
get_selex <- function(fit){
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
      mutate(rel_error=(est-truth)/truth, abs_error=est-truth, sim=as.factor(im))
    return(selex)
  }
  lapply(fits, function(i) ff(i)) %>% bind_rows()
}
get_waa <- function(fit){
  ff <- function(fit){
    if(!fit$model$optimized){
      print(fit$model)
      return(NULL)
    }
    waa <- data.frame(par='waa',
                      age=1:10,
                      est=fit$fit$rep$pred_waa[5,1,],
                      truth=fit$truth$pred_waa[5,1,]) %>% bind_cols(fit$model) %>%
      mutate(rel_error=(est-truth)/truth, abs_error=est-truth, sim=as.factor(im))
    return(waa)
  }
  lapply(fits, function(i) ff(i)) %>% bind_rows()
}
get_pars <- function(fit){
  ff <- function(fit){
    if(!fit$model$optimized){
      print(fit$model)
      return(NULL)
    }
    pars <- merge(fit$ompars, fit$empars, by='par2') %>%
      filter(grepl(x=par.y, "growth_a|SD_par|mean_rec_pars|logit_q|log_F1|log_N1_pars|logit_selpars|catch_paa_pars|index_paa_pars"))
    pars <- pars %>% select(par=par.x, par2, truth=value.x, est=value.y) %>%
      bind_cols(fit$model) %>%
      mutate(rel_error=(est-truth)/truth, abs_error=est-truth, sim=as.factor(im))
    return(pars)
  }
  lapply(fits, function(i) ff(i)) %>% bind_rows()
}
get_growth <- function(fit){
  print(fit$model)
  ff <- function(fit){
    if(!fit$model$optimized){
      print(fit$model)
      return(NULL)
    }
    p1 <- data.frame(par=c('logk', 'logLinf'),
                     est=fit$fit$rep$growth_a[1:2,1],
                     truth=fit$truth$growth_a[1:2,1])
    p2 <- data.frame(par=c('SD1', 'SD2'),
                     est=fit$fit$rep$expSD,
                     truth=fit$truth$expSD)
    growth <- bind_rows(p1,p2) %>% bind_cols(fit$model) %>%
      mutate(rel_error=(est-truth)/truth, abs_error=est-truth, sim=as.factor(im))
    return(growth)
  }
  lapply(fits, function(i) ff(i)) %>% bind_rows()
}


fits <- list.files('../results', pattern='RDS', recursive=1,
                      full.names=TRUE) %>%
  lapply(function(i) readRDS(i))

## check convergence stats
models <- lapply(fits, function(x) x$model) %>% bind_rows
table(models$optimized,models$em)
table(models$sdreport,models$em)

## fit <- fits[[2]]

ts <- get_ts(fits)
ggplot(ts, aes(year, truth, group=sim, color=em_growth_est)) +
  geom_line() + facet_wrap('par', scales='free')
ggplot(ts, aes(year, rel_error, group=sim)) + geom_line(alpha=.5) +
  geom_hline(yintercept=0, col=2, lwd=1) + facet_grid(par~em_growth_est)

pars <- get_pars(fits)
ggplot(pars, aes(par2, rel_error)) + geom_violin() +
  geom_hline(yintercept=0, col=2, lwd=1) +
  facet_grid(em_growth_est~par, scales='free')

## test <- filter(ts, par=='SSB' & sim==1 & em==2)
## plot(test$year, test$truth)
## lines(test$year, test$est)


selex <- get_selex(fits)
ggplot(selex, aes(age, est, group=sim, color=em_growth_est)) +
  geom_line() + facet_grid(par~em_growth_est)
ggplot(selex, aes(age, abs_error, group=sim)) + geom_line(alpha=.5) +
  geom_hline(yintercept=0, col=2, lwd=1) + facet_grid(par~em_growth_est)

growth <- get_growth(fits)
ggplot(growth, aes(par, abs_error)) + geom_violin() +
  facet_grid(.~em_growth_est) + geom_hline(yintercept=0, col=2)
ggplot(growth, aes(par, rel_error)) + geom_violin() +
  facet_grid(.~em_growth_est) + geom_hline(yintercept=0, col=2)

waa <- get_waa(fits)
ggplot(waa, aes(age, rel_error, group=sim)) + geom_line() +
  facet_grid(.~em_growth_est) + geom_hline(yintercept=0, col=2)
ggplot(waa, aes(age, abs_error, group=sim)) + geom_line() +
  facet_grid(.~em_growth_est) + geom_hline(yintercept=0, col=2)

