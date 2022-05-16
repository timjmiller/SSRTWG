set_N1 = function()

source(file.path(here(), "common_code", "make_basic_info.R"))
source(file.path(here(), "common_code", "get_FMSY.R"))
groundfish_info <- make_basic_info()

gf_selectivity = list(
  model = c(rep("logistic", groundfish_info$n_fleets),rep("logistic", groundfish_info$n_indices)),
  initial_pars = rep(list(c(5,1)), groundfish_info$n_fleets + groundfish_info$n_indices)) #fleet, index

gf_M = list(initial_means = rep(0.2, length(groundfish_info$ages)))

gf_NAA_re = list(
  N1_pars = exp(10)*exp(-(0:(length(groundfish_info$ages)-1))*gf_M$initial_means[1]),
  sigma = "rec", #random about mean
  cor="iid", #random effects are independent
  use_steepness = 1,
  #recruit_model = 2, #random effects with a constant mean
  recruit_model = 3, #B-H
  recruit_pars = c(0.75,exp(10))
)


gf_input <- prepare_wham_input(basic_info = groundfish_info, selectivity = gf_selectivity, NAA_re = gf_NAA_re, M= gf_M)
gf_input$data$bias_correct_oe[] = 0
gf_input$data$bias_correct_pe[] = 0
gf_input$data$FMSY_init[] = F40
gf_input$data$FMSY_init[] = 0.6
temp <- fit_wham(gf_input, do.fit = FALSE, MakeADFun.silent = FALSE)
F40 = exp(temp$rep$log_FXSPR[1])


get_FMSY(steepness = 0.70, NAA_re = gf_NAA_re, basic_info = groundfish_info, selectivity = gf_selectivity, 
         M = gf_M, F40)

#will get steepness such that Fmsy = F40
steepness = get_steepness(steepness_ini =0.75)
gf_NAA_re$recruit_pars[1] = steepness
groundfish_info$F[] = F40

gf_input <- prepare_wham_input(basic_info = groundfish_info, selectivity = gf_selectivity, NAA_re = gf_NAA_re, M= gf_M)
gf_input$data$FMSY_init[] = F40
temp <- fit_wham(gf_input, do.fit = FALSE, MakeADFun.silent = FALSE)
Jan1NAA_per_recruit = wham:::get_SPR(F = F40,M = temp$rep$MAA[1,],temp$rep$selAA[[1]][1,],
  mat= rep(1,gf_input$data$n_ages), waa = rep(1,gf_input$data$n_ages), fracyrssb = 0,at.age = TRUE)
gf_NAA_re$N1_pars = exp(temp$rep$log_R_MSY[1]) * Jan1NAA_per_recruit

gf_input <- prepare_wham_input(basic_info = groundfish_info, selectivity = gf_selectivity, NAA_re = gf_NAA_re, M= gf_M)
gf_input$par$log_NAA[,1] = temp$rep$log_R_MSY[1]
gf_input$par$log_NAA_sigma[] = -100
gf_input$data$FMSY_init[] = F40
temp <- fit_wham(gf_input, do.fit = FALSE, MakeADFun.silent = FALSE)

sim = temp$simulate(complete=TRUE)

#input already has basic_info, M, NAA_re, selectivity specified

  
  input = 
  
}

gf_NAA_re$recruit_pars[1] = 0.6928989
gf_input <- prepare_wham_input(basic_info = groundfish_info, selectivity = gf_selectivity, NAA_re = gf_NAA_re, M= gf_M)
