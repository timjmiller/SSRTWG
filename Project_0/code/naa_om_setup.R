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
source(file.path(here(), "common_code", "set_selectivity.R"))
source(file.path(here(), "common_code", "set_q.R"))
source(file.path(here(), "common_code", "set_simulation_options.R"))
source(file.path(here(), "common_code", "get_FMSY.R"))
source(file.path(here(), "Project_0", "code", "make_om.R"))
source(file.path(here(), "Project_0", "code", "sim_management.R"))
verify_version()

write.dir <- file.path(here(),"Project_0", "inputs")

#if(!exists("write.dir")) write.dir = getwd()
if(!dir.exists(write.dir)) dir.create(write.dir, recursive = T)
setwd(write.dir)

#Operating model factors
#NAA sigmas for each scenario
R_sig <- c(0.5,1.5)
NAA_sig <- c(NA,0.25,0.5)
#F time series
Fhist = c("H-MSY","MSY")
#how much observation error
obs_error = c("L", "H")
df.oms <- expand.grid(R_sig = R_sig, NAA_sig = NAA_sig, Fhist = Fhist, 
  obs_error = obs_error, stringsAsFactors = FALSE)

#logistic-normal age comp SDs for L/H observation error (both indices AND CATCH!!!????)
L_N_sigma = c(L = 0.3, H = 1.5)
#(log) index SDs for L/H observation error 
index_sigma = c(L = 0.1, H = 0.4)

n.mods <- dim(df.oms)[1] #24 operating model scenarios
df.oms$Model <- paste0("om_",1:n.mods)
df.oms <- df.oms %>% select(Model, everything()) # moves Model to first col
# look at model table
df.oms
saveRDS(df.oms, file.path(here(),"Project_0", "inputs", "df.oms.RDS"))

# naa.sim.jobs = expand.grid(om = df.oms$Model, sim = 1:500)
# naa.sim.jobs$wham_commit = NA
# naa.sim.jobs$member = NA
# naa.sim.jobs$member[naa.sim.jobs$sim %in% 1:5] = "TJM"
# saveRDS(naa.sim.jobs, file.path(here(),"Project_0", "inputs", "naa.sim.jobs.RDS"))



gf_info = make_basic_info()

#selectivity is not changing
gf_selectivity = list(
  model = c(rep("logistic", gf_info$n_fleets),rep("logistic", gf_info$n_indices)),
  initial_pars = rep(list(c(5,1)), gf_info$n_fleets + gf_info$n_indices)) #fleet, index

#M set is not changing
gf_M = list(initial_means = rep(0.2, length(gf_info$ages)))

#NAA_re set up that can be changed for each OM scenario
gf_NAA_re = list(
  N1_pars = exp(10)*exp(-(0:(length(gf_info$ages)-1))*gf_M$initial_means[1]),
  sigma = "rec", #random about mean
  cor="iid", #random effects are independent
  use_steepness = 1,
  #recruit_model = 2, #random effects with a constant mean
  recruit_model = 3, #B-H
  recruit_pars = c(0.75,exp(10))
)

#make inputs for operating model (smaller objects to save, can recreate simulated data sets)
om_inputs = list()
for(i in 1:NROW(df.oms)){
  print(paste0("row ", i))
  NAA_re = gf_NAA_re
  NAA_re$sigma_vals = df.oms$R_sig[i]
  if(!is.na(df.oms$NAA_sig[i])) { #full random effects for NAA
    NAA_re$sigma = "rec+1"
    NAA_re$sigma_vals[2] = df.oms$NAA_sig[i]
  }
  Fhist. = "Fmsy"
  if(df.oms$Fhist[i] == "H-MSY") Fhist. = "H-L"
  max_mult = 2.5 # fishing at 2.5 x Fmsy
  if(Fhist. == "Fmsy") max_mult = 1
  min_mult = 1 # fishing at Fmsy
  
  om_inputs[[i]] = make_om(Fhist = Fhist., N1_state = "overfished", selectivity = gf_selectivity, catchability = NULL,
    M = gf_M, NAA_re = NAA_re, age_comp = "logistic-normal-miss0", brp_year = 1, eq_F_init = 0.3, 
    om_input = TRUE, max_mult_Fmsy = max_mult, min_mult_Fmsy = min_mult)
  #turn off bias correction
  om_inputs[[i]] = set_simulation_options(om_inputs[[i]], simulate_data = TRUE, simulate_process = TRUE, simulate_projection = TRUE, 
    bias_correct_pe = FALSE, bias_correct_oe = FALSE)
  #set L-N SD parameters for catch and index age comp
  om_inputs[[i]]$par$catch_paa_pars[,1] = log(L_N_sigma[df.oms$obs_error[i]])
  om_inputs[[i]]$par$index_paa_pars[,1] = log(L_N_sigma[df.oms$obs_error[i]])
  om_inputs[[i]]$data$agg_catch_sigma[] = 0.1
  om_inputs[[i]]$data$agg_index_sigma[] = index_sigma[df.oms$obs_error[i]]
  #change Neff so scalar doesn't affect L-N SD
  om_inputs[[i]]$data$catch_Neff[] = 1
  om_inputs[[i]]$data$index_Neff[] = 1
}

saveRDS(om_inputs, file.path(here(),"Project_0","inputs", "NAA_om_inputs.RDS"))


#check equilibrium
om_msy = make_om(Fhist = "Fmsy", N1_state = "Fmsy", selectivity = gf_selectivity, 
  M = gf_M, NAA_re = gf_NAA_re, brp_year = 1, eq_F_init = 0.3)

input = om_msy#$input
input$par$log_NAA_sigma[] = -100 #no process error
temp <- fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE)
sim = temp$simulate(complete=TRUE)
sim$NAA #NAA should stay the same throughout time series

#start out overfished and overfishing the whole time series
om_overfished = make_om(Fhist = "H", N1_state = "overfished", selectivity = gf_selectivity, 
  M = gf_M, NAA_re = gf_NAA_re, brp_year = 1, eq_F_init = 0.3)
input = om_overfished#$input
input$par$log_NAA_sigma[] = -100 #no process error
temp <- fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE)
sim = temp$simulate(complete=TRUE)
sim$NAA #NAA should stay the same throughout time series
