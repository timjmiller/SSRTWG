# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="devel")
library(wham)
devtools::load_all("C://Users//a37789//Documents//wham", compile=FALSE, recompile=FALSE)
library(here)
library(tidyverse)
source(file.path("C://Users//a37789//Documents//noaa//SSRTWG","common_code", "set_selectivity.R")) # overwrite wham's set_selectivity() with ours
source(file.path("C://Users//a37789//Documents//noaa//SSRTWG","common_code", "make_basic_info.R"))

groundfish_info <- make_basic_info() # 10 ages
n.b <- groundfish_info$n_fleets + groundfish_info$n_indices
M = list(initial_means = rep(0.2, length(groundfish_info$ages)))
NAA_re = list(N1_pars = exp(10)*exp(-(0:(length(groundfish_info$ages)-1))*M$initial_means[1]))

# logistic, all estimated
selectivity = list(model = rep("logistic", n.b))
input = prepare_wham_input(basic_info = groundfish_info, selectivity = selectivity, M = M, NAA_re = NAA_re)

# selpars matrix, rows are blocks
#   columns     parameters
#    1:n.ages    age-specific
#    next 2      logistic
#    last 4      double-logistic
input$par$logit_selpars
matrix(as.numeric(input$map$logit_selpars), nrow=n.b)

# age-specific, all estimated
selectivity = list(model = rep("age-specific", n.b))
input = prepare_wham_input(basic_info = groundfish_info, selectivity = selectivity, M = M, NAA_re = NAA_re)
input$par$logit_selpars
matrix(as.numeric(input$map$logit_selpars), nrow=n.b)

# age-specific, fix flat-topped in fleet (block 1) ages 6-10
#   initialize all estimated ages at 0.5 and 6:10 at 1
selectivity = list(model = rep("age-specific", n.b),
                    fix_pars = list(6:10,NULL,NULL),
                    initial_pars = list(c(rep(0.5,5),rep(1,5)),
                                        rep(0.5,10),
                                        rep(0.5,10)))
input = prepare_wham_input(basic_info = groundfish_info, selectivity = selectivity, M = M, NAA_re = NAA_re)
input = set_selectivity(input, selectivity) # overwrite selectivity with SSRTWG set_selectivity function
input$par$logit_selpars
matrix(as.numeric(input$map$logit_selpars), nrow=n.b)

# age-specific, fix flat-topped in fleet (block 1) ages 6-10 as above,
#   but also have ages 3-5 share same (estimated) selectivity
selectivity = list(model = rep("age-specific", n.b),
                    initial_pars = list(c(rep(0.5,5),rep(1,5)),
                                        rep(0.5,10),
                                        rep(0.5,10)),
                    map_pars = list(c(1,2,3,3,3,rep(NA,5)),
                                    4:13,
                                    14:23))
input = prepare_wham_input(basic_info = groundfish_info, selectivity = selectivity, M = M, NAA_re = NAA_re)
input = set_selectivity(input, selectivity) # overwrite selectivity with SSRTWG set_selectivity function
input$par$logit_selpars
matrix(as.numeric(input$map$logit_selpars), nrow=n.b)

# logistic, have the 2 indices share same estimated selectivity
selectivity = list(model = rep("logistic", n.b),
                    map_pars = list(1:2, 3:4, 3:4))
input = prepare_wham_input(basic_info = groundfish_info, selectivity = selectivity, M = M, NAA_re = NAA_re)
input = set_selectivity(input, selectivity) # overwrite selectivity with SSRTWG set_selectivity function
input$par$logit_selpars
matrix(as.numeric(input$map$logit_selpars), nrow=n.b)

# AR1_age random effects, fleet and indices
selectivity = list(model = rep("age-specific", n.b),
                   re = rep("ar1", n.b))
input = prepare_wham_input(basic_info = groundfish_info, selectivity = selectivity, M = M, NAA_re = NAA_re)
input = set_selectivity(input, selectivity) # overwrite selectivity with SSRTWG set_selectivity function
input$par$logit_selpars
matrix(as.numeric(input$map$logit_selpars), nrow=n.b)

# AR1_age random effects, fleet only
selectivity = list(model = rep("age-specific", n.b),
                   re = c("ar1","none","none"))
input = prepare_wham_input(basic_info = groundfish_info, selectivity = selectivity, M = M, NAA_re = NAA_re)
input = set_selectivity(input, selectivity) # overwrite selectivity with SSRTWG set_selectivity function
input$par$logit_selpars
matrix(as.numeric(input$map$logit_selpars), nrow=n.b)

# AR1_age random effects, index 1 only
selectivity = list(model = rep("age-specific", n.b),
                   re = c("none","ar1","none"))
input = prepare_wham_input(basic_info = groundfish_info, selectivity = selectivity, M = M, NAA_re = NAA_re)
input = set_selectivity(input, selectivity) # overwrite selectivity with SSRTWG set_selectivity function
input$par$logit_selpars
matrix(as.numeric(input$map$logit_selpars), nrow=n.b)
mod <- fit_wham(input, do.sdrep=T, do.retro=F, do.osa=F, do.proj=F)  

# ---------------------------------------------------------------------------
# test using map_pars instead of fix_pars with wham ex4
wham.dir <- find.package("wham")
asap3 <- read_asap3_dat(file.path(wham.dir,"inst","extdata","ex1_SNEMAYT.dat"))
input <- prepare_wham_input(asap3, recruit_model=2,
                    selectivity=list(model=rep("age-specific",3),
                                    initial_pars=list(c(0.1,0.5,0.5,1,1,1),
                                                       c(0.5,0.5,0.5,1,0.5,0.5),
                                                       c(0.5,0.5,1,1,1,1)), 
                                    fix_pars=list(4:6,4,3:6)),
                    NAA_re = list(sigma='rec+1',cor='iid'),
                    age_comp = "logistic-normal-miss0")
input$par$logit_selpars
# [1,] -2.197225    0    0  Inf  Inf  Inf   NA   NA   NA    NA    NA    NA
# [2,]  0.000000    0    0  Inf    0    0   NA   NA   NA    NA    NA    NA
# [3,]  0.000000    0  Inf  Inf  Inf  Inf   NA   NA   NA    NA    NA    NA
matrix(as.numeric(input$map$logit_selpars), nrow=n.b)
# [1,]    1    4    7   NA   NA   NA   NA   NA   NA    NA    NA    NA
# [2,]    2    5    8   NA    9   10   NA   NA   NA    NA    NA    NA
# [3,]    3    6   NA   NA   NA   NA   NA   NA   NA    NA    NA    NA
mod <- fit_wham(input, do.sdrep=T, do.retro=F, do.osa=F, do.proj=F)  
round(t(sapply(mod$rep$selAA, function(x) x[1,])),4)
# [1,] 0.0136 0.4290 0.8968    1 1.0000 1.0000
# [2,] 0.0147 0.2537 0.6089    1 0.9799 0.7888
# [3,] 0.2333 0.7965 1.0000    1 1.0000 1.0000
mod$parList$logit_q
 # -7.617794 -8.624950

# now use map_pars instead of fix_pars
input <- prepare_wham_input(asap3, recruit_model=2,
                            selectivity=list(model=rep("age-specific",3),
                                             initial_pars=list(c(0.1,0.5,0.5,1,1,1),
                                                               c(0.5,0.5,0.5,1,0.5,0.5),
                                                               c(0.5,0.5,1,1,1,1)), 
                                             map_pars=list(c(1:3,NA,NA,NA),
                                                           c(1:3,NA,5:6),
                                                           c(1,2,NA,NA,NA,NA))),
                            NAA_re = list(sigma='rec+1',cor='iid'),
                            age_comp = "logistic-normal-miss0")
input$par$logit_selpars
# [1,] -2.197225    0    0   10   10   10   NA   NA   NA    NA    NA    NA
# [2,]  0.000000    0    0   10    0    0   NA   NA   NA    NA    NA    NA
# [3,]  0.000000    0   10   10   10   10   NA   NA   NA    NA    NA    NA
matrix(as.numeric(input$map$logit_selpars), nrow=n.b)
# [1,]    1    4    7   10   13   16   19   22   NA    NA    NA    NA
# [2,]    2    5    8   11   14   17   20   23   NA    NA    NA    NA
# [3,]    3    6    9   12   15   18   21   24   NA    NA    NA    NA
input = set_selectivity(input, selectivity=list(model=rep("age-specific",3),
                                               initial_pars=list(c(0.1,0.5,0.5,1,1,1),
                                                                 c(0.5,0.5,0.5,1,0.5,0.5),
                                                                 c(0.5,0.5,1,1,1,1)), 
                                               map_pars=list(c(1:3,NA,NA,NA),
                                                             c(4:6,NA,7:8),
                                                             c(9:10,NA,NA,NA,NA)))) # overwrite selectivity with SSRTWG set_selectivity function
input$par$logit_selpars
# [1,] -2.197225    0    0  Inf  Inf  Inf   NA   NA   NA    NA    NA    NA
# [2,]  0.000000    0    0  Inf    0    0   NA   NA   NA    NA    NA    NA
# [3,]  0.000000    0  Inf  Inf  Inf  Inf   NA   NA   NA    NA    NA    NA
matrix(as.numeric(input$map$logit_selpars), nrow=n.b)
# [1,]    1    2    3   NA   NA   NA   NA   NA   NA    NA    NA    NA
# [2,]    4    5    6   NA    7    8   NA   NA   NA    NA    NA    NA
# [3,]    9   10   NA   NA   NA   NA   NA   NA   NA    NA    NA    NA
mod2 <- fit_wham(input, do.sdrep=T, do.retro=F, do.osa=F, do.proj=F)  
round(t(sapply(mod2$rep$selAA, function(x) x[1,])),4)
# [1,] 0.0136 0.4290 0.8968    1 1.0000 1.0000
# [2,] 0.0147 0.2537 0.6089    1 0.9799 0.7888
# [3,] 0.2333 0.7965 1.0000    1 1.0000 1.0000
all.equal(mod2$opt$objective,mod$opt$objective)
# TRUE

# ------------------------------------------------------------
# test AR1 fix
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
matrix(as.numeric(input$map$logit_selpars), nrow=n.b)
unique(input$map$selpars_re)
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

# find and fix age with highest sel at 1
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
round(mod3$parList$logit_selpars,5)
# [1,] -4.25412 -0.26515  2.11220  Inf      Inf      Inf   NA   NA   NA    NA    NA    NA
# [2,] -0.46399 -0.46399 -0.46399  Inf -0.46399 -0.46399   NA   NA   NA    NA    NA    NA
# [3,] -1.14151  1.52177      Inf  Inf      Inf      Inf   NA   NA   NA    NA    NA    NA
mod3$parList$logit_q
# -7.647629 -8.635871


