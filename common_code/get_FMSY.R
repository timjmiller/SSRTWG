get_FMSY = function(steepness = 0.75, input, NAA_re, brp_year){
  NAA_re$recruit_pars[1] = steepness
  #print(input$data$FMSY_init)
  input = set_NAA(input, NAA_re)
  #print(paste0("get_FMSY: Fmsy_init = ", input$data$FMSY_init[brp_year]))
  #stop()
  #input <- prepare_wham_input(basic_info = basic_info, selectivity = selectivity, NAA_re = NAA_re, M= M)
  temp <- fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE)
  return(exp(temp$rep$log_FMSY[brp_year]))
}

get_steepness = function(steepness_ini, input, NAA_re, F40, brp_year){
  #find steepness that gives Fmsy = F40
  obj = function(steepness, input, NAA_re, F40, brp_year){ 
  #print(input$data$FMSY_init)
    FMSY = get_FMSY(steepness = steepness,input = input, NAA_re = NAA_re, brp_year = brp_year)
    print(paste0("F40: ", round(F40,3), ", current steepness: ", round(steepness,3), 
      ", current FMSY: ", round(FMSY,3)))
    return((FMSY - F40)^2) #minimize squared difference
  }
  out = nlminb(steepness_ini, obj, input = input, NAA_re = NAA_re, F40 = F40, brp_year= brp_year)
  return(out)
}

set_F_scenario = function(input, Fhist, Fmsy, max_mult = 1, min_mult = 1, change_time = 0.5){
  nby <- input$data$n_years_model
  year_change <- floor(nby * change_time)
  if(!Fhist %in% c("Fmsy","H-L","H")) {
    stop("Fhist must be 'Fmsy'','H-L'', or 'H'. Edit set_F_scenario to allow other options.")
  }

  if(Fhist == "Fmsy"){
    cat("OM will have F=Fmsy for all years.\n")
    input$par$log_F1[] = log(Fmsy)
    input$par$F_devs[] = 0
  }
  if(Fhist == "H"){
    cat("OM will have F= max_mult * Fmsy for all years.\n")
    input$par$log_F1[] = log(max_mult * Fmsy)
    input$par$F_devs[] = 0
  }
  if(Fhist == "H-L"){
    cat("OM will have F decrease abruptly from max_mult x Fmsy to min_mult * Fmsy at ", change_time, 
      " x n_model_years = ", year_change, ".\n")
    input$par$log_F1[] = log(max_mult * Fmsy)
    input$par$F_devs[] = 0
    input$par$F_devs[year_change-1,] = log(min_mult) - log(max_mult)
  }
  if(Fhist == "L-H"){
    cat("OM will have F increase abruptly from min_mult x Fmsy to max_mult * Fmsy at ", change_time, 
      " x n_model_years = ", year_change, ".\n")
    input$par$log_F1[] = log(min_mult * Fmsy)
    input$par$F_devs[] = 0
    input$par$F_devs[year_change-1,] = log(max_mult) - log(min_mult)
  }
  return(input)
}
