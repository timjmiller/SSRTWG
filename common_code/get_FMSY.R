get_FMSY = function(steepness = 0.75, NAA_re, basic_info, selectivity, M, F40 = exp(temp$rep$log_FXSPR[1])){
  NAA_re$recruit_pars[1] = steepness
  input <- prepare_wham_input(basic_info = basic_info, selectivity = selectivity, NAA_re = NAA_re, M= M)
  input$data$FMSY_init[] = F40
  temp <- fit_wham(input, do.fit = FALSE, MakeADFun.silent = FALSE)
  return(exp(temp$rep$log_FMSY[1]))
}

get_steepness = function(steepness_ini){
  #find steepness that gives Fmsy = F40
  obj = function(steepness){ 
    FMSY = get_FMSY(steepness,NAA_re = gf_NAA_re, basic_info = groundfish_info, selectivity = gf_selectivity, 
         M = gf_M, F40=F40)
    print(c(FMSY,F40))
    return((FMSY - F40)^2)
  }
  out = nlminb(steepness_ini, obj)
  return(out)
}
