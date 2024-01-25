# test allowing estimation of FAA freely, 'WHAM-as-SAM'
# add option selectivity$model = 'none'

# selectivity$model
#   'none', WHAM-as-SAM estimation of FAA/QAA matrix freely, intended for use with random effects (warning if re = none or ar1_y). Fixes full F/q at upper bound, allows model to estimate max selected age (can differ by year)
#   'age-specific', estimate sel-at-age and F as fixed effects, should fix one sel-at-age at 1 (warning if no ages are fixed). Can add random effects.
#   'logistic', 'double-logistic', 'decreasing-logistic' as before. Can add random effects.

# selectivity$re
#   "none", warn if model=none selectivity for all ages and years will = 1 (probably don't want this unless all ages fixed, eg YOY index)
#   "iid",
#   "ar1", warn if model='age-specific' and no ages are fixed
#   "ar1_y", warn if model=none selectivity for all ages will be the same (probably don't want this), but could work well with model = 'age-specific'
#   "2dar1"

#   model = none, re = iid or 2dar1: fix F1 and Fdevs at upper bound, estimate FAA by age and year directly
#   model = none, re = ar1: estimate F1, Fdevs, selAA re by age (fix mean logit_selpars at 0.5)
#   ? model = age-specific, re = ar1_y: fix F1 and Fdevs, estimate selAA re by year, estimate selAA as fixed effects by age
#   model = age-specific, re = ar1: estimate F1, Fdevs, sel mean, sig, cor_age (should fix one age at 1)

# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="devel")
# devtools::install_github("brianstock/wham", dependencies=TRUE, ref="devel")
library(wham)
devtools::load_all("C://Users//a37789//Documents//wham", compile=FALSE, recompile=FALSE)
library(here)
library(tidyverse)
source(file.path("C://Users//a37789//Documents//noaa//SSRTWG","common_code", "set_selectivity.R")) # overwrite wham's set_selectivity() with ours

# test using wham ex4
wham.dir <- find.package("wham")
asap3 <- read_asap3_dat(file.path(wham.dir,"inst","extdata","ex1_SNEMAYT.dat"))
n.b = asap3$dat$n_indices + asap3$dat$n_fleets

# model = age-specific, re = ar1 (by age), no ages fixed 
# should trigger warning to fix one age at 1 or use model = none
input1 <- prepare_wham_input(asap3, recruit_model=2,
                            selectivity=list(model=rep("age-specific",3),
                                             initial_pars=list(c(0.1,0.5,0.5,1,1,1),
                                                               c(0.5,0.5,0.5,1,0.5,0.5),
                                                               c(0.5,0.5,1,1,1,1)), 
                                             fix_pars=list(4:6,4,3:6)),
                            NAA_re = list(sigma='rec+1',cor='iid'),
                            age_comp = "logistic-normal-miss0")
input = set_selectivity(input1, selectivity=list(model=rep("age-specific",3),
                                    initial_pars=list(c(0.1,0.5,0.5,1,1,1),
                                                      c(0.5,0.5,0.5,0.5,0.5,0.5),
                                                      c(0.5,0.5,1,1,1,1)), 
                                    fix_pars=list(4:6,NULL,3:6), 
                                    re=c('none','ar1','none')))
input$par$logit_selpars
# [1,] -2.197225    0    0  Inf  Inf  Inf   NA   NA   NA    NA    NA    NA
# [2,]  0.000000    0    0    0    0    0   NA   NA   NA    NA    NA    NA
# [3,]  0.000000    0  Inf  Inf  Inf  Inf   NA   NA   NA    NA    NA    NA
matrix(as.numeric(input$map$logit_selpars), nrow=n.b) # estimate mean of AR1
# [1,]    1    4    6   NA   NA   NA   NA   NA   NA    NA    NA    NA
# [2,]    2    2    2    2    2    2   NA   NA   NA    NA    NA    NA
# [3,]    3    5   NA   NA   NA   NA   NA   NA   NA    NA    NA    NA
input$par$sel_repars
# [1,] -2.302585    0    0
# [2,] -2.302585    0    0
# [3,] -2.302585    0    0
matrix(as.numeric(input$map$sel_repars), nrow=n.b) # estimate 2 AR1 pars, sigma and correlation
# [1,]   NA   NA   NA
# [2,]    1    2   NA
# [3,]   NA   NA   NA
matrix(as.numeric(input$map$selpars_re), nrow=input$data$n_years_model) # one re for each age, constant in time
# [1,]    1    2    3    4    5    6
# [2,]    1    2    3    4    5    6
# [3,]    1    2    3    4    5    6
# ...
# [42,]    1    2    3    4    5    6
# [43,]    1    2    3    4    5    6
# [44,]    1    2    3    4    5    6
unique(input$map$selpars_re) # one re for each age
# [1] 1 2 3 4 5 6
# Levels: 1 2 3 4 5 6
mod2 <- fit_wham(input, do.sdrep=T, do.retro=F, do.osa=F, do.proj=F)  

round(t(sapply(mod2$rep$selAA, function(x) x[1,]/max(x[1,]))),4)
# [1,] 0.0139 0.4310 0.8860    1 1.0000 1.0000
# [2,] 0.0162 0.2778 0.6751    1 0.9899 0.8403
# [3,] 0.2404 0.8167 1.0000    1 1.0000 1.0000
round(t(sapply(mod2$rep$selAA, function(x) x[1,])),4)
# [1,] 0.0139 0.4310 0.8860 1.000 1.0000 1.0000
# [2,] 0.0153 0.2628 0.6386 0.946 0.9364 0.7949
# [3,] 0.2404 0.8167 1.0000 1.000 1.0000 1.0000
round(mod2$parList$logit_selpars,5)
# [1,] -4.26523 -0.27784  2.05023      Inf      Inf      Inf   NA   NA   NA    NA    NA    NA
# [2,] -0.02421 -0.02421 -0.02421 -0.02421 -0.02421 -0.02421   NA   NA   NA    NA    NA    NA
# [3,] -1.15031  1.49384      Inf      Inf      Inf      Inf   NA   NA   NA    NA    NA    NA
mod2$parList$logit_q
 # -7.641498 -8.633531

# ------------------------------------------------------------------------------
# now find and fix age with highest sel at 1
# model = age-specific, re = ar1 (by age), one age fixed --> should NOT trigger warning
which(mod2$rep$selAA[[2]][1,] == max(mod2$rep$selAA[[2]][1,]))
# 4
input = set_selectivity(input1, selectivity=list(model=rep("age-specific",3),
                                                 initial_pars=list(c(0.1,0.5,0.5,1,1,1),
                                                                   c(0.5,0.5,0.5,1,0.5,0.5),
                                                                   c(0.5,0.5,1,1,1,1)), 
                                                 fix_pars=list(4:6,4,3:6),
                                                 re=c('none','ar1','none')))
input$par$logit_selpars
matrix(as.numeric(input$map$logit_selpars), nrow=n.b)
unique(input$map$selpars_re)
# [1] 1 2 3 4 5
# Levels: 1 2 3 4 5
mod3 <- fit_wham(input, do.sdrep=T, do.retro=F, do.osa=F, do.proj=F)  
round(t(sapply(mod3$rep$selAA, function(x) x[1,])),4)
# [1,] 0.0140 0.4341 0.8921    1 1.0000 1.0000
# [2,] 0.0154 0.2639 0.6356    1 0.9132 0.7934
# [3,] 0.2420 0.8208 1.0000    1 1.0000 1.0000
round(t(sapply(mod2$rep$selAA, function(x) x[1,]/max(x[1,]))),4)
# [1,] 0.0139 0.4310 0.8860    1 1.0000 1.0000
# [2,] 0.0162 0.2778 0.6751    1 0.9899 0.8403
# [3,] 0.2404 0.8167 1.0000    1 1.0000 1.0000
round(mod3$parList$logit_selpars,5)
# [1,] -4.25412 -0.26515  2.11220  Inf      Inf      Inf   NA   NA   NA    NA    NA    NA
# [2,] -0.46399 -0.46399 -0.46399  Inf -0.46399 -0.46399   NA   NA   NA    NA    NA    NA
# [3,] -1.14151  1.52177      Inf  Inf      Inf      Inf   NA   NA   NA    NA    NA    NA
mod3$parList$logit_q
# -7.647629 -8.635871
mod2$parList$logit_q
# -7.641498 -8.633531

# ------------------------------------------------------------------------------
# fleet model = none, re = ar1 (by age)
input = set_selectivity(input1, selectivity=list(model=c("none","age-specific","age-specific"),
                                                 initial_pars=list(rep(.01,6),
                                                                   c(0.5,0.5,0.5,1,0.5,0.5),
                                                                   c(0.5,0.5,1,1,1,1)), 
                                                 fix_pars=list(NULL,4,3:6),
                                                 re=c('2dar1','none','none')))
input$par$logit_selpars
matrix(as.numeric(input$map$logit_selpars), nrow=n.b)
input$par$sel_repars
matrix(as.numeric(input$map$sel_repars), nrow=n.b)
matrix(as.numeric(input$map$selpars_re), nrow=input$data$n_years_model)

mod <- fit_wham(input, do.sdrep=T, do.retro=F, do.osa=F, do.proj=F)  
round(t(sapply(mod3$rep$selAA, function(x) x[1,])),4)


# index model = none, re = ar1 (by age)
input = set_selectivity(input1, selectivity=list(model=c("age-specific","none","age-specific"),
                                                 initial_pars=list(c(0.1,0.5,0.5,1,1,1),
                                                                   c(0.5,0.5,0.5,0.5,0.5,0.5),
                                                                   c(0.5,0.5,1,1,1,1)), 
                                                 fix_pars=list(4:6,NULL,3:6),
                                                 re=c('none','ar1','none')))
input$par$logit_selpars
matrix(as.numeric(input$map$logit_selpars), nrow=n.b)
unique(input$map$selpars_re)
# [1] 1 2 3 4 5
# Levels: 1 2 3 4 5
mod3 <- fit_wham(input, do.sdrep=T, do.retro=F, do.osa=F, do.proj=F)  
round(t(sapply(mod3$rep$selAA, function(x) x[1,])),4)
# [1,] 0.0140 0.4341 0.8921    1 1.0000 1.0000
# [2,] 0.0154 0.2639 0.6356    1 0.9132 0.7934
# [3,] 0.2420 0.8208 1.0000    1 1.0000 1.0000
round(t(sapply(mod2$rep$selAA, function(x) x[1,]/max(x[1,]))),4)
# [1,] 0.0139 0.4310 0.8860    1 1.0000 1.0000
# [2,] 0.0162 0.2778 0.6751    1 0.9899 0.8403
# [3,] 0.2404 0.8167 1.0000    1 1.0000 1.0000
round(mod3$parList$logit_selpars,5)
# [1,] -4.25412 -0.26515  2.11220  Inf      Inf      Inf   NA   NA   NA    NA    NA    NA
# [2,] -0.46399 -0.46399 -0.46399  Inf -0.46399 -0.46399   NA   NA   NA    NA    NA    NA
# [3,] -1.14151  1.52177      Inf  Inf      Inf      Inf   NA   NA   NA    NA    NA    NA
mod3$parList$logit_q
# -7.647629 -8.635871
mod2$parList$logit_q
# -7.641498 -8.633531

