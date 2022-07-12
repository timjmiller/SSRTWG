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

# test using wham ex4
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
mod <- fit_wham(input, do.sdrep=T, do.retro=F, do.osa=F, do.proj=F)  
round(t(sapply(mod$rep$selAA, function(x) x[1,])),4)
# [1,] 0.0136 0.4290 0.8968    1 1.0000 1.0000
# [2,] 0.0147 0.2537 0.6089    1 0.9799 0.7888
# [3,] 0.2333 0.7965 1.0000    1 1.0000 1.0000
mod$parList$logit_q
 # -7.617794 -8.624950

input = set_selectivity(input, selectivity=list(model=rep("age-specific",3),
                                    initial_pars=list(c(0.1,0.5,0.5,1,1,1),
                                                      c(0.5,0.5,0.5,0.5,0.5,0.5),
                                                      c(0.5,0.5,1,1,1,1)), 
                                    fix_pars=list(4:6,NULL,3:6),
                                    re=c('none','ar1','none')))
input$par$logit_selpars
matrix(as.numeric(input$map$logit_selpars), nrow=n.b)
mod2 <- fit_wham(input, do.sdrep=T, do.retro=F, do.osa=F, do.proj=F)  

# catchability is now estimated much higher (4 on logit) and all selAA are really low (-9 on logit)
round(t(sapply(mod2$rep$selAA, function(x) x[1,]/max(x[1,]))),4)
# [1,] 0.0141 0.4320 0.8960    1 1.0000 1.0000
# [2,] 0.0156 0.2564 0.6118    1 0.9796 0.7912
# [3,] 0.2430 0.8089 1.0000    1 1.0000 1.0000
t(sapply(mod$rep$selAA, function(x) x[1,]))
#              [,1]         [,2]         [,3]         [,4]         [,5]
# [1,] 1.414048e-02 0.4320158652 0.8960452940 1.0000000000 1.0000000000
# [2,] 7.781320e-06 0.0001276102 0.0003044413 0.0004976551 0.0004874906
# [3,] 2.429613e-01 0.8089107753 1.0000000000 1.0000000000 1.0000000000
#              [,6]
# [1,] 1.0000000000
# [2,] 0.0003937389
# [3,] 1.0000000000
mod2$parList$logit_selpars
mod2$parList$logit_q
 # 3.995047 -8.630140

