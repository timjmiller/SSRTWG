make_om <- function(Fhist = "Fmsy", N1_state = "Fmsy",
                    selectivity = NULL, M = NULL, NAA_re = NULL,
                    catchability = NULL, growth, LW,
                    ecov = NULL, age_comp = "logistic-normal-miss0",
                    len_comp='multinomial',
                    brp_year = 1,
                    eq_F_init = 0.3, om_input = TRUE,
                    max_mult_Fmsy = 2.5, min_mult_Fmsy = 1,
                    F_change_time = 0.5, df.oms = NULL) {

  basic_info <- make_basic_info()
  basic_info$fracyr_indices[,1] = 0.25
  basic_info$fracyr_indices[,2] = 0.75

  ny <- length(basic_info$years)
  na <- 10
  basic_info$lengths <- seq(1,100, by=2)
  nlbins <- length(basic_info$lengths)
  basic_info$n_lengths <- nlbins
  ## Turn on marginal lengths comps for survey two and turn off the
  ## age comps
  basic_info$use_index_paa <- cbind(rep(1,ny), rep(0,ny))
  basic_info$index_pal <- array(data=1/nlbins, dim=c(2,40, nlbins))
  basic_info$use_index_pal <- cbind(rep(0,ny), rep(1,ny))
  ## for now assuming multinomial N=100
  basic_info$index_NeffL <- cbind(rep(1,ny), rep(100,ny))

  #overfishing_mult = 2.5 #multiplier for Fmsy for overfishing
  input <-
    wham::prepare_wham_input(basic_info = basic_info, growth=growth,
                       LW=LW, len_comp=len_comp,
                       selectivity = selectivity, NAA_re = NAA_re, M= M, ecov = ecov,
                       age_comp = age_comp, catchability = catchability)
  input$data$FXSPR_init[] = eq_F_init
  input$data$FMSY_init[] = eq_F_init
  # You will need to add a conditional for the next two lines if you want to play with AR1 and RW, the indices will change:
  input$par$Ecov_process_pars[2,] = df.oms$Ecov_re_sig # This is cond sd for the AR1 Ecov process. double check this
  input$par$Ecov_process_pars[3,] = df.oms$Ecov_re_cor # This is phi for the AR1 Ecov process. double check this
  input$par$Ecov_beta[5,1,1,] = df.oms$Ecov_effect # Effect. the index will vary if the number of fleets change
  #if you want to change the %Spawning Potential
  #input$data$percentSPR = 35 #for example
  temp <- fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE)
  #print(temp$rep$MAA)
  #print(temp$rep$selAA[[1]])
  print(exp(temp$rep$log_FXSPR))
  print(exp(temp$rep$log_FMSY))
#  stop()

  F40 = exp(temp$rep$log_FXSPR[brp_year])
  if(is.na(F40) || F40 < 1e-5 || F40 > 10) stop("search for F40 failed, try changing eq_F_init. \n")
  cat(paste0("F40 = ", round(F40,3), ".\n"))
  Fmsy = exp(temp$rep$log_FMSY[brp_year])
  cat(paste0("Fmsy = ", round(Fmsy,3), ".\n"))


  #will get steepness such that Fmsy = F40
  #print(input$data$FMSY_init)
  if(temp$input$data$use_steepness==1){
    cat("Finding steepness such that Fmsy = F40.\n")
    steepness = get_steepness(steepness_ini = 0.75, input = input, NAA_re = NAA_re, F40 = F40,
      brp_year = brp_year)
    NAA_re$recruit_pars[1] = steepness$par
    input = set_NAA(input, NAA_re)
    #input$data$FMSY_init[] = F40
    temp <- fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE)
    Fmsy = exp(temp$rep$log_FMSY[brp_year])
    if(is.na(Fmsy) || Fmsy < 1e-5 || Fmsy > 10) stop("search for Fmsy failed, try changing eq_F_init. \n")
  }
  if(round(Fmsy,3) != round(F40,3)){
    cat(paste0("NOTE: Fmsy (", round(Fmsy,3), ") is not equal to F40 (", round(F40,3), ")\n"))
  }
  #set up initial numbers at age according to equilibrium conditions
  #N1_state == "Fmsy"
  eq_F_N1 = exp(temp$rep$log_FMSY[brp_year]) # = F40
  if(N1_state == "overfished"){
    eq_F_N1 = eq_F_N1 * max_mult_Fmsy
  }
  if(N1_state == "unfished"){
    eq_F_N1 = 0
  }
  Jan1_NAA_per_recruit = wham:::get_SPR(F = eq_F_N1, M = temp$rep$MAA[brp_year,],
    sel = temp$rep$selAA[[1]][brp_year,], mat= rep(1,input$data$n_ages),
    waa = rep(1,input$data$n_ages), fracyrssb = 0,at.age = TRUE)
  Jan1_SSB_per_recruit = wham:::get_SPR(F = eq_F_N1, M = temp$rep$MAA[brp_year,],
    sel = temp$rep$selAA[[1]][brp_year,], mat= input$data$mature[brp_year,],
    waa = input$data$waa[input$data$waa_pointer_ssb,brp_year,], fracyrssb = input$data$fracyr_SSB[brp_year])
  SR_a = exp(temp$rep$log_SR_a[brp_year])
  SR_b = exp(temp$rep$log_SR_b[brp_year])
  #beverton-holt!!!!
  if(input$data$recruit_model < 3) stop("NAA_re$recruit_model must be >= 3 for B-H or Ricker recruitment")
  if(input$data$recruit_model == 3) {
    eq_R_N1 = (SR_a - 1/Jan1_SSB_per_recruit) / SR_b;
  } else {
    cat(paste0("FYI: you are using a Ricker model for this operating model.\n"))
    eq_R_N1 = (log(SR_a) + log(Jan1_SSB_per_recruit))/(SR_b * Jan1_SSB_per_recruit)
  }
  NAA_re$N1_pars = eq_R_N1 * Jan1_NAA_per_recruit
  input = set_NAA(input, NAA_re)
  ## input = set_ecov(input, ecov)
  ## input$par$Ecov_pro -- set this manually
  input = set_q(input, catchability)
  ## what is this doing and why do I need it?
  ##  input = set_selectivity(input, selectivity)
  input = set_M(input, M)
  #set F relative to Fmsy. This function is in get_FMSY.R
  input = set_F_scenario(input, Fhist, Fmsy = Fmsy, max_mult = max_mult_Fmsy, min_mult= min_mult_Fmsy,
    change_time = F_change_time)
  if(om_input) return(input)
  else return(fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE))
}
