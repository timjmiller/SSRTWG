library(wham) #make sure to use the right version of wham
library(tidyr)
library(dplyr)
library(here)
source(file.path(here(), "common_code", "make_basic_info.R"))
source(file.path(here(), "common_code", "set_NAA.R"))
source(file.path(here(), "common_code", "set_M.R"))
source(file.path(here(), "common_code", "set_q.R"))
source(file.path(here(), "common_code", "set_ecov.R"))
source(file.path(here(), "common_code", "set_selectivity.R"))
source(file.path(here(), "common_code", "set_simulation_options.R"))
source(file.path(here(), "common_code", "get_FMSY.R"))
source(file.path(here(), "Ecov_study", "mortality", "code", "make_om.R"))
source(file.path(here(), "Ecov_study", "mortality", "code", "sim_management.R"))
verify_version()


set_F_scenario = function(input, Fhist, Fmsy, max_mult = 1, min_mult = 1, change_time = 0.5){
  nby <- input$data$n_years_model
  year_change <- floor(nby * change_time)
  if(!Fhist %in% c("Fmsy","H-L","H","L-H")) {
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




naa_om_inputs <- readRDS(file.path(here::here(),"Project_0","inputs", "NAA_om_inputs.RDS"))

#SR parameters are the same for all naa_om models 
temp <- fit_wham(naa_om_inputs[[1]], do.fit = FALSE, MakeADFun.silent = TRUE)
SRab <- exp(c(temp$rep$log_SR_a[1], temp$rep$log_SR_b[1]))

#Operating model factors
#Constant factors
NAA_re       <- c("rec")  #recruitment random effects
Ecov_obs_sig <- c(0.1)    #ecov observation error
Ecov_re_sig  <- c(0.1)
obs_error    <- c("L","H")

#variable factors
R_sig       <- c(0.1,0.3,0.5,1.0)        #recruitment random effects sigma
Fhist       <- c("H-MSY","MSY","L-H")  #fishing history
NAA_cor     <- c(0.2,0.8)        #correlation of recruitment random effects
Ecov_re_cor <- c(0.2,0.8)        #correlation of ecov random effects
Ecov_effect <- c(0.1, 1.0)       #beta coefficients; need to modify according to functional form
Ecov_how    <- c(0,1,2)          #ecov-recruiment functional form
recruit_mod <- c(3)

df.oms <- expand.grid(NAA_re = NAA_re,
                      Ecov_obs_sig=Ecov_obs_sig, 
                      Ecov_re_sig=Ecov_re_sig, 
                      obs_error = obs_error,
                      R_sig = R_sig,
                      Fhist = Fhist,
                      NAA_cor = NAA_cor,
                      Ecov_re_cor=Ecov_re_cor, 
                      Ecov_effect = Ecov_effect, 
                      Ecov_how=Ecov_how,
                      stringsAsFactors = FALSE)

L_N_sigma   = c(L = 0.3, H = 1.5)   #logistic-normal age comp SDs for L/H observation error (both indices AND CATCH!!!????)
index_sigma = c(L = 0.1, H = 0.4) #(log) index SDs for L/H observation error 

n.mods       <- dim(df.oms)[1] 
df.oms$Model <- paste0("om_",1:n.mods)
df.oms       <- df.oms %>% select(Model, everything()) # moves Model to first col

##--BETA STANDARDIZATION--##
#source(file.path(here(),"Ecov_study","recruitment_functions","code","beta_standardization.r"))

saveRDS(df.oms, file.path(here(), "Ecov_study", "recruitment_functions", "inputs", "df.oms.RDS"))

gf_info = make_basic_info()

#selectivity is not changing
gf_selectivity <- list(model       =c(rep("logistic", gf_info$n_fleets),rep("logistic", gf_info$n_indices)),
                       initial_pars=rep(list(c(5,1)), gf_info$n_fleets + gf_info$n_indices)) #fleet, index

#M set is not changing
gf_M <- list(initial_means=rep(0.2, length(gf_info$ages)))

#NAA_re set up that can be changed for each OM scenario
gf_NAA_re = list(
  N1_pars = exp(10)*exp(-(0:(length(gf_info$ages)-1))*gf_M$initial_means[1]),
  sigma = "rec", 
  cor   ="ar1_y", 
  use_steepness = 0,
  recruit_model = 3, #B-H
  recruit_pars  = SRab) #defined above from naa_om_inputs

gf_ecov <- list(
  label = "Ecov",
  process_model = "ar1",
  lag = 0,
  mean = cbind(rep(0, length(gf_info$years))),
  year = gf_info$years,
  use_obs = cbind(rep(1, length(gf_info$years))),
  process_mean_vals = 0
)

beta_vals <- list(rep(list(matrix(0,1,length(gf_info$ages))), 4))

#make inputs for operating model (smaller objects to save, can recreate simulated data sets)
om_inputs = list()
for(i in 1:NROW(df.oms)){
  print(paste0("row ", i))
  NAA_re = gf_NAA_re
  NAA_re$sigma = "rec"
  NAA_re$sigma_vals = df.oms$R_sig[i]
  NAA_re$cor_vals = df.oms$NAA_cor[i]
  
  ecov_i = gf_ecov
  ecov_i$logsigma = cbind(rep(log(df.oms$Ecov_obs_sig[i]), length(ecov_i$year)))
  ecov_i$process_sig_vals = df.oms$Ecov_re_sig[i]
  ecov_i$process_cor_vals = df.oms$Ecov_re_cor[i]
  ecov_i$how = df.oms$Ecov_how[i]
  ecov_i$where = "recruit"
  beta_vals_i = beta_vals
  beta_vals_i[[1]][[1]][] <- df.oms$Ecov_effect[i]
  ecov_i$beta_vals = beta_vals_i
  
  Fhist. = "Fmsy"
  if(Fhist. == "Fmsy"){ 
	max_mult = 1
  	min_mult = 1 # fishing at Fmsy
	}

  if(df.oms$Fhist[i] == "H-MSY"){ 
	Fhist. = "H-L"
  	max_mult = 2.5 # fishing at 2.5 x Fmsy
	}

  if(df.oms$Fhist[i] == "L-H"){
	Fhist. = "L-H"
	min_mult = 0.25
	max_mult = 1
	}
 
  
  om_inputs[[i]] = make_om(Fhist       = Fhist., 
                           N1_state    = "overfished", 
                           selectivity = gf_selectivity, 
                           M           = gf_M, 
                           NAA_re      = NAA_re, 
                           ecov        = ecov_i, 
                           age_comp    = "logistic-normal-miss0", 
                           brp_year    = 1, 
                           eq_F_init   = 0.3, 
                           om_input    = TRUE, 
                           max_mult_Fmsy = max_mult, 
                           min_mult_Fmsy = min_mult)
  #turn off bias correction
  om_inputs[[i]] = set_simulation_options(om_inputs[[i]], 
                                          simulate_data = TRUE, 
                                          simulate_process = TRUE, 
                                          simulate_projection = TRUE, 
                                          bias_correct_pe = FALSE, 
                                          bias_correct_oe = FALSE)
  #set L-N SD parameters for catch and index age comp
  om_inputs[[i]]$par$catch_paa_pars[,1] = log(L_N_sigma[df.oms$obs_error[i]])
  om_inputs[[i]]$par$index_paa_pars[,1] = log(L_N_sigma[df.oms$obs_error[i]])
  om_inputs[[i]]$data$agg_catch_sigma[] = 0.1
  om_inputs[[i]]$data$agg_index_sigma[] = index_sigma[df.oms$obs_error[i]]
  #change Neff so scalar doesn't affect L-N SD
  om_inputs[[i]]$data$catch_Neff[] = 1
  om_inputs[[i]]$data$index_Neff[] = 1
}

