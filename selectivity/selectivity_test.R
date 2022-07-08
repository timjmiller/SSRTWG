# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="devel")
library(wham)
devtools::load_all("C://Users//a37789//Documents//wham", compile=FALSE, recompile=FALSE)
library(here)
library(tidyverse)
source(file.path("C://Users//a37789//Documents//noaa//SSRTWG","common_code", "set_selectivity.R")) # overwrite wham's set_selectivity() with ours
source(file.path("C://Users//a37789//Documents//noaa//SSRTWG","common_code", "make_basic_info.R"))

groundfish_info <- make_basic_info() # 10 ages
M = list(initial_means = rep(0.2, length(groundfish_info$ages)))
NAA_re = list(N1_pars = exp(10)*exp(-(0:(length(groundfish_info$ages)-1))*M$initial_means[1]))

# logistic, all estimated
selectivity = list(model = c(rep("logistic", groundfish_info$n_fleets),rep("logistic", groundfish_info$n_indices)))
input = prepare_wham_input(basic_info = groundfish_info, selectivity = selectivity, M = M, NAA_re = NAA_re)

# selpars matrix, rows are blocks
#   columns     parameters
#    1:n.ages    age-specific
#    next 2      logistic
#    last 4      double-logistic
input$par$logit_selpars
matrix(as.numeric(input$map$logit_selpars), nrow=(groundfish_info$n_fleets + groundfish_info$n_indices))

# age-specific, all estimated
selectivity = list(model = c(rep("age-specific", groundfish_info$n_fleets),rep("age-specific", groundfish_info$n_indices)))
input = prepare_wham_input(basic_info = groundfish_info, selectivity = selectivity, M = M, NAA_re = NAA_re)
input$par$logit_selpars
matrix(as.numeric(input$map$logit_selpars), nrow=(groundfish_info$n_fleets + groundfish_info$n_indices))

# age-specific, fix flat-topped in fleet (block 1) ages 6-10
#   initialize all estimated ages at 0.5 and 6:10 at 1
selectivity = list(model = c(rep("age-specific", groundfish_info$n_fleets),rep("age-specific", groundfish_info$n_indices)),
                    fix_pars = list(6:10,NULL,NULL),
                    initial_pars = list(c(rep(0.5,5),rep(1,5)),
                                        rep(0.5,10),
                                        rep(0.5,10)))
input = prepare_wham_input(basic_info = groundfish_info, selectivity = selectivity, M = M, NAA_re = NAA_re)
input = set_selectivity(input, selectivity) # overwrite selectivity with SSRTWG set_selectivity function
input$par$logit_selpars
matrix(as.numeric(input$map$logit_selpars), nrow=(groundfish_info$n_fleets + groundfish_info$n_indices))

# age-specific, fix flat-topped in fleet (block 1) ages 6-10 as above,
#   but also have ages 3-5 share same (estimated) selectivity
selectivity = list(model = c(rep("age-specific", groundfish_info$n_fleets),rep("age-specific", groundfish_info$n_indices)),
                    initial_pars = list(c(rep(0.5,5),rep(1,5)),
                                        rep(0.5,10),
                                        rep(0.5,10)),
                    map_pars = list(c(1,2,3,3,3,rep(NA,5)),
                                    4:13,
                                    14:23))
input = prepare_wham_input(basic_info = groundfish_info, selectivity = selectivity, M = M, NAA_re = NAA_re)
input = set_selectivity(input, selectivity) # overwrite selectivity with SSRTWG set_selectivity function
input$par$logit_selpars
matrix(as.numeric(input$map$logit_selpars), nrow=(groundfish_info$n_fleets + groundfish_info$n_indices))

# logistic, have the 2 indices share same estimated selectivity
selectivity = list(model = c(rep("logistic", groundfish_info$n_fleets),rep("logistic", groundfish_info$n_indices)),
                    map_pars = list(1:2, 3:4, 3:4))
input = prepare_wham_input(basic_info = groundfish_info, selectivity = selectivity, M = M, NAA_re = NAA_re)
input = set_selectivity(input, selectivity) # overwrite selectivity with SSRTWG set_selectivity function
input$par$logit_selpars
matrix(as.numeric(input$map$logit_selpars), nrow=(groundfish_info$n_fleets + groundfish_info$n_indices))

