# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="devel")
if(file.exists("c:/Users/timothy.j.miller")) {
  library(wham, lib.loc = "c:/work/wham/old_packages/77bbd94")
} else library(wham) #make sure to use the right version of wham
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

naa_om_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "NAA_om_inputs.RDS"))
# temp = sapply(naa_om_inputs, function(x) {
#   temp = fit_wham(x, do.fit = FALSE, MakeADFun.silent = TRUE)
#   return(temp$rep$log_SR_a)
# })
#SR parameters are the same for all naa_om models 
temp = fit_wham(naa_om_inputs[[1]], do.fit = FALSE, MakeADFun.silent = TRUE)
SRab = exp(c(temp$rep$log_SR_a[1], temp$rep$log_SR_b[1]))

###############Estimating model inputs
#Estimating model factors
SR_model = c(2)
M_est = c(TRUE, FALSE)
re_config = c("rec","rec+1", "rec+M")
Ecov_est = c(TRUE,FALSE)

#create data.frame defining estimation models data.fram
df.ems <- expand.grid(M_est = M_est, re_config = re_config, Ecov_est = Ecov_est, stringsAsFactors = FALSE)
saveRDS(df.ems, file.path(here(),"Ecov_study", "mortality", "inputs", "df.ems.RDS"))

#same as naa_om_setup.R
gf_info = make_basic_info()

#same as naa_om_setup.R
#selectivity is not changing
gf_selectivity = list(
  model = c(rep("logistic", gf_info$n_fleets),rep("logistic", gf_info$n_indices)),
  initial_pars = rep(list(c(5,1)), gf_info$n_fleets + gf_info$n_indices)) #fleet, index

#different from naa_om_setup.R
#M set up that can be changed for each EM
#gf_M = list(initial_means = rep(0.2, length(gf_info$ages)))
gf_M = list(initial_means = rep(0.2, length(gf_info$ages)), model = "age-specific")

#same as naa_om_setup.R
#NAA_re set up that can be changed for each EM
gf_NAA_re = list(
  N1_pars = exp(10)*exp(-(0:(length(gf_info$ages)-1))*gf_M$initial_means[1]),
  sigma = "rec", #random about mean
  cor="iid", #random effects are independent
  use_steepness = 0,
  recruit_model = 2 #random effects with a constant mean
  #recruit_model = 3#, #B-H
  #recruit_pars = SRab
)

gf_ecov <- list(
  label = "Ecov",
  process_model = "ar1",
  logsigma = cbind(rep(log(0.1), length(gf_info$year))),
  lag = 0,
  mean = cbind(rep(0, length(gf_info$years))),
  year = gf_info$years,
  use_obs = cbind(rep(1, length(gf_info$years))),
  where = "none",
  how = 0
)

#make inputs for estimating model (smaller objects to save, can overwrinte data elements with simulated data)
em_inputs = list()
for(i in 1:NROW(df.ems)){
  print(paste0("row ", i))
  NAA_re_i = gf_NAA_re
  M_i = gf_M
  ecov_i = gf_ecov
  selectivity = gf_selectivity

  #NAA_re$sigma_vals = df.oms$R_sig[i] #don't start at true values?
  if(df.ems$re_config[i] == "rec+1") { #full random effects for NAA
    NAA_re_i$sigma = "rec+1"
    #NAA_re$sigma_vals[2] = df.oms$NAA_sig[i] #don't start at true values?
  }
  if(df.ems$re_config[i] == "rec+M") { #full random effects for NAA
    M_i$re = "ar1_y"
    #NAA_re$sigma_vals[2] = df.oms$NAA_sig[i] #don't start at true values?
  }
  if(df.ems$M_est[i]) { #estimate mean M
    M_i$est_ages = 1:length(gf_info$ages)
    #NAA_re$sigma_vals[2] = df.oms$NAA_sig[i] #don't start at true values?
  }
  if(df.ems$Ecov_est[i]){
    ecov_i$how = 1
    ecov_i$where = "M"
  }
  basic_info <- make_basic_info()
  basic_info$fracyr_indices[,1] = 0.25
  basic_info$fracyr_indices[,2] = 0.75

  em_inputs[[i]] <- prepare_wham_input(basic_info = basic_info, selectivity = selectivity, NAA_re = NAA_re_i, M= M_i, ecov = ecov_i,
    age_comp = "logistic-normal-miss0")  
  #have to do this after the input is made
  #source("c:/work/wham_master/wham/R/set_ecov.R")
  #em_input[[i]] <- set_ecov(em_inputs[[i]], ecov = ecov_i)
  em_inputs[[i]] = set_M(em_inputs[[i]], M = M_i) #this set_M will change parameter values. needed for iid M_re
  #let's estimate AR1 for all these models
  # if(df.ems$re_config[i] == "rec+M") { 
  #   em_inputs[[i]]$map$M_repars = factor(c(1,NA,NA)) #still need to fix rho which is set to 0 above
  # }
  if(df.ems$M_est[i]) { #estimate mean M constant across ages
    em_inputs[[i]]$map$M_a = factor(rep(1,length(em_inputs[[i]]$par$M_a)))
  }

  #turn off bias correction
  em_inputs[[i]] = set_simulation_options(em_inputs[[i]], simulate_data = TRUE, simulate_process = TRUE, simulate_projection = TRUE, 
    bias_correct_pe = FALSE, bias_correct_oe = FALSE)
  #change Neff so scalar doesn't affect L-N SD
  em_inputs[[i]]$data$catch_Neff[] = 1
  em_inputs[[i]]$data$index_Neff[] = 1
  em_inputs[[i]]$data$FXSPR_init[] = 0.3
  em_inputs[[i]]$data$FMSY_init[] = 0.3
}

saveRDS(em_inputs, file.path(here(),"Ecov_study", "mortality", "inputs", "em_inputs.RDS"))
