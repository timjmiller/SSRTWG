# get the base input values for IBMWG analyses
# these will be modified by other functions,
# but this gets the defaults

# 2020/06/25 - recruitment parameters need fixing
# 2020/07/01 - TJM modded so BH can take steepness as input, agreed on OM vals now complete

get_base_input <- function(n_selblocks = 1, Fhist = 1, Fmsy_scale = 2.5, scaa = FALSE, nprojyrs = 40) {
  input <- list()
  input$nf <- nf <- 1 #number of fleets (we only need 1)
  input$ni <- ni <- 2 #number of indices
  input$na <- na <- 10 #number of ages
  input$modyears <- 1970:2019 #50 years
  #input$maturity <- 1/(1 + exp(-1*(1:input$na - 5))) #maturity at age
  input$maturity <- c(0.04, 0.25, 0.60, 0.77, 0.85, 0.92, 1, 1, 1, 1)
  input$fracyr_spawn <- 0 #when spawning occurs within annual time step
  #L = LVB_pars[1]*(1-exp(-LVB_pars[2]*(1:na - 0))) #make up LVB growth
  #W = LW_pars[1]*L^LW_pars[2] #make up L-W function to give WAA
  W <- c(0.15, 0.5, 0.9, 1.4, 2.0, 2.6, 3.2, 4.1, 5.9, 9.0)
  input$waa_catch <- t(matrix(W, na, nf)) #WAA for each fleet
  input$waa_indices <- t(matrix(W, na, ni)) #WAA for each index
  input$waa_totcatch = input$waa_ssb = input$waa_jan1 = rbind(W) #WAA for total catch, SSB, Jan 1 pop
  input$catch_cv <- rep(0.1, nf) #CVs for aggregate catch for each fleet
  input$catch_Neff <- rep(200, nf) #Effective sample size for age comp for each fleet
  input$index_cv <- c(0.3,0.4) #rep(0.3, ni) #CVs for aggregate indices
  input$index_Neff <- rep(100, ni) #Effectin sample size for age compe for each index
  #input$fracyr_indices = (1:ni)/(ni+1) #when index observations occur within annual time step
  input$fracyr_indices <- c(0.25,0.75)
  input$sel_model_fleets = rep(2, nf) #logistic selectivity for each fleet
  #input$sel_model_indices = 2
  input$sel_model_indices = rep(2,ni) #logistic selectivity for each index

  sel.list=list(model=rep("logistic",ni+nf), re=rep("none",ni+nf), 
    initial_pars= list(
      c(3.57,1), #fishery (see factorial pruning)
      c(1.8, 1/6), #survey 1 (see factorial pruning)
      c(1.2, 1/5.5))) #survey 2 (see factorial pruning)
  
  equilib.sel = 1/(1+exp(-sel.list$initial_pars[[1]][2]*(1:na - sel.list$initial_pars[[1]][1])))
  equilib.sel = equilib.sel/max(equilib.sel)
  input$sel.list = sel.list
  input$NAA.list = list(sigma='rec',cor='ar1_y')
  input$proj.list = list(n.yrs=nprojyrs, use.last.F=TRUE, use.avg.F=FALSE, use.FXSPR=FALSE,
                                                proj.F=NULL, proj.catch=NULL, avg.yrs=NULL,
                                                cont.ecov=TRUE, use.last.ecov=FALSE, avg.ecov.yrs=NULL, proj.ecov=NULL, cont.Mre=NULL)


  #input$q = 0.3
  #input$q = (1:ni)/(ni+1) #catchability for each index
  input$q <- c(0.0002, 0.0001)
  
  #input$F = matrix(rep(0.2/nf,length(input$modyears)), length(input$modyears), nf) #Annual Full F for each fleet
  input$F <- matrix(rep(c(seq(0.1,0.7,length.out = 1995-1971+1),
                      seq(0.7,0.25,length.out = 2010-1995+1)[-1],
                      rep(0.25,2020-2010)),nf),
                    ncol = nf)
  input$M = rep(0.2, na) #Nat. Mort at age
  #input$N1 = exp(10)*exp(-(0:(na-1))*input$M[1]) #Initial numbers at age
  #recruit_model = 2 #random devations around mean. 3 = BH (mean_rec_pars = a,b R = aS/(1 + bS)), 4 = Ricker (mean_rec_pars = a,b)
  input$recruit_model = 3 #Beverton-Holt
  #input$use_steepness = 0 #don't use steepness
  input$use_steepness = 1 #mean_rec_pars = h,R0
  #input$mean_rec_pars = numeric(c(0,1,2,2)[recruit_model])
  h <- 0.75 #a = 4*0.7/(0.3*25.8) #h=0.7, phi0=25.8
  R0 <- 10000 #b = (a - 1/25.8)/exp(10) #R0 = exp(10)
  #if(recruit_model == 2) input$mean_rec_pars[] = exp(10)
  #if(recruit_model == 3) 
  spr0 = wham:::get_SPR(0, M=input$M, sel = rep(1,na), mat=input$maturity, waassb=input$waa_catch, fracyrssb=input$fracyr_spawn, at.age = FALSE)
  a = 4*h/((1-h)*spr0)
  b = (a - 1/spr0)/R0
  eq.yield = function(log_F){
    spr = wham:::get_SPR(exp(log_F), M=input$M, sel = equilib.sel, mat=input$maturity, waassb=input$waa_catch, fracyrssb=input$fracyr_spawn, at.age = FALSE)
    ypr = wham:::get_YPR(exp(log_F), M=input$M, sel = equilib.sel, waacatch=input$waa_catch, at.age = FALSE)
    R_F = (a - 1/sum(spr))/b
    Y_F = R_F * ypr
    #print(Y_F)
    return(-Y_F)
  }
  Fmsy = exp(nlminb(log(0.2), eq.yield)$par)
  #print(paste0("Fmsy = ", Fmsy))
  
  #Fhist == 1 means overfishing for first half of base period, then Fmsy for second half
  prop.F = rep(Fmsy_scale, length(input$modyears)/2)
  prop.F[length(prop.F) + 1:(length(input$modyears)-length(prop.F))] = 1
  #Fhist == 2 means overfishing for whole base period
  if(Fhist == 2) prop.F = rep(Fmsy_scale,length(input$modyears))
  input$F = cbind(Fmsy * prop.F)
  #print(input$F)
  F1 = input$F[1]
  sprF1 = wham:::get_SPR(F1, M=input$M, equilib.sel, mat=input$maturity, waassb=input$waa_catch, fracyrssb=input$fracyr_spawn, at.age = FALSE)
  nprF1 = wham:::get_SPR(F1, M=input$M, equilib.sel, mat=rep(1,na), waassb=rep(1,na), fracyrssb=input$fracyr_spawn, at.age = TRUE)
  R_F1 = (a - 1/sum(sprF1))/b
  input$N1 <- R_F1*nprF1 #Initial numbers at age
  
  input$mean_rec_pars = c(h, R0)
  #if(recruit_model == 4) input$mean_rec_pars[2] = exp(-10)
  
  input$NAA_rho <- c(0, 0.4) #AR1 rho for age and year (Here: just AR1 with year for recruitment deviations)
  input$NAA_sigma <- 0.5*prod(sqrt(1-input$NAA_rho^2)) #recruitment sd, Marginal SD = 0.5. Need more values if full state-space abundance at age

  if(scaa){
    input$use_steepness = 0
    input$mean_rec_pars = input$mean_rec_pars[2]
    input$recruit_model = 2
    input$NAA.list = NULL
    #input$proj.list = NULL
  }

  om = prepare_wham_om_input(input, recruit_model = input$recruit_model, selectivity=input$sel.list, NAA_re = input$NAA.list, proj.opts = input$proj.list)

  if(n_selblocks == 2){ #two fishery selectivity blocks, second starting second half of base period.
    selpars2 = c(5,1)
    om$data$n_selblocks = 4
    om$data$selblock_models = rep(2,4)
    om$data$selblock_models_re = rep(1,4)
    om$data$selblock_pointer_fleets[26:50,1] = 4
    om$data$selblock_years[26:50,1] = 0
    om$data$selblock_years = cbind(om$data$selblock_years, rep(0:1, each = 25))
    om$data$n_years_selblocks[1] = 25
    om$data$n_years_selblocks = c(om$data$n_years_selblocks, 25)
    om$data$selpars_est = rbind(om$data$selpars_est, rep(c(0,1,0), c(10,2,4)))
    om$data$n_selpars = rep(2,4)
    om$data$n_selpars_est = rep(2,4)
    om$data$selpars_lower = om$data$selpars_lower[c(1:3,3),]
    om$data$selpars_upper = om$data$selpars_upper[c(1:3,3),]
    om$par$logit_selpars = om$par$logit_selpars[c(1:3,1),]
    om$par$logit_selpars[4,11] = log(selpars2[1]-0) - log(10-selpars2[1]) #= 0 logit(a50 = 5), slope = 1 like first selblock
    om$par$sel_repars = om$par$sel_repars[c(1:3,3),]
    om$map$logit_selpars = as.factor(rbind(matrix(as.integer(om$map$logit_selpars),3,16),rep(c(NA,7:8,NA),c(10,1,1,4))))
    om$map$sel_repars = as.factor(rep(NA, length(om$par$sel_repars)))

    #Fmsy, base period full F and initial N1 change based on recent selectivity
    #change equilibrium selectivity
    equilib.sel = 1/(1+exp(-selpars2[2]*(1:na - selpars2[1])))
    equilib.sel = equilib.sel/max(equilib.sel)
    Fmsy = exp(nlminb(log(0.2), eq.yield)$par)
    #print(paste0("Fmsy = ", Fmsy))
    prop.F = rep(Fmsy_scale, length(input$modyears)/2)
    prop.F[length(prop.F) + 1:(length(input$modyears)-length(prop.F))] = 1
    if(Fhist == 2) prop.F = rep(Fmsy_scale,length(input$modyears))
    F = Fmsy * prop.F
    #print(F)
    om$par$log_F1 = log(F[1])
    om$par$F_devs = cbind(diff(log(F)))
    sprF1 = wham:::get_SPR(F[1], M=input$M, equilib.sel, mat=input$maturity, waassb=input$waa_catch, fracyrssb=input$fracyr_spawn, at.age = FALSE)
    nprF1 = wham:::get_SPR(F[1], M=input$M, equilib.sel, mat=rep(1,na), waassb=rep(1,na), fracyrssb=input$fracyr_spawn, at.age = TRUE)
    R_F1 = (a - 1/sum(sprF1))/b
    om$par$log_N1_pars <- log(R_F1*nprF1) #Initial numbers at age    
  }
  
  om$data$Fbar_ages = 10 #not sure if this is needed anywhere, but I did need it for examining retro values.
  return(om)

  #Below intial NAA will be changed before simulations conducted
  #input$N1 <- R0*exp(-(0:(na-1))*input$M[1]) #Initial numbers at age
} #end get_base_input function
