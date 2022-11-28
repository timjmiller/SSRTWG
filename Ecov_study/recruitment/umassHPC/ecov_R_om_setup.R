# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="devel")
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
source(file.path(here(), "Project_0", "code", "make_om.R"))
source(file.path(here(), "Project_0", "code", "sim_management.R"))
verify_version()

write.dir <- file.path(here(),"Ecov_study", "recruitment", "inputs")

naa_om_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "NAA_om_inputs.RDS"))
#SR parameters are the same for all naa_om models 
temp = fit_wham(naa_om_inputs[[1]], do.fit = FALSE, MakeADFun.silent = TRUE)
SRab = exp(c(temp$rep$log_SR_a[1], temp$rep$log_SR_b[1]))

#if(!exists("write.dir")) write.dir = getwd()
if(!dir.exists(write.dir)) dir.create(write.dir, recursive = T)
setwd(write.dir)

#number of simulations for each scenario
#nsim = 1000
#nsim = 100
nsim=10

Ecov_where <- c("recruit")
Ecov_mean  <- 0
Ecov_sig   <- c("L" = 0.1, "H" = 0.5) #units?
ar1_y      <- c("L" = 0.5, "H" = 0.95)
beta       <- c("L" = 0.1 ,"H" = 1.0) #units?
obs_sig    <- c("L" = 1e-5,"H" = 0.25) #units?

df.e.r.oms       <- expand.grid(Ecov_sig=Ecov_sig, Ecov_phi=ar1_y, Ecov_mean = Ecov_mean, beta = beta, Ecov_where = Ecov_where, obs_sig = obs_sig)
n.mods           <- dim(df.e.r.oms)[1] #108 scenarios
df.e.r.oms$Model <- paste0("om_",1:n.mods)
df.e.r.oms       <- df.e.r.oms %>% select(Model, everything()) # moves Model to first col
df.e.r.oms       <- cbind(Recruitment = 2, NAA_re = "rec+1", df.e.r.oms) #recruit model not yet discussed by WG.
df.e.r.oms$nsim  <- rep(nsim,nrow(df.e.r.oms))

saveRDS(df.mods,file.path(write.dir, "om_sim_inputs_GLB_recruitment.RDS"))

gf_info = make_basic_info()

#selectivity is not changing
gf_selectivity = list(
  model = c(rep("logistic", gf_info$n_fleets),rep("logistic", gf_info$n_indices)),
  initial_pars = rep(list(c(5,1)), gf_info$n_fleets + gf_info$n_indices)) #fleet, index

#M set is not changing
#gf_M = list(model = "age-specific", 
#  initial_means = rep(0.2, length(gf_info$ages)),
#  re = "ar1_y"#, # This is needed to set up operating models with a single annual re, iid or ar1_y
#  #sigma_vals = 0.1,
#  #cor_vals = 0
#  )
gf_M <- list(initial_means=rep(0.2, length(gf_info$ages)))

#NAA_re set up that can be changed for each OM scenario
gf_NAA_re = list(
  N1_pars = exp(10)*exp(-(0:(length(gf_info$ages)-1))*gf_M$initial_means[1]),
  sigma = "rec", #random about mean
  cor="iid", #random effects are independent
  use_steepness = 0,
  #recruit_model = 2, #random effects with a constant mean
  recruit_model = 3, #B-H
  recruit_pars = SRab) #defined above from naa_om_inputs

#
gf_catchability = list(prior_sd = c(NA, 0.3))

#make inputs for operating model (smaller objects to save, can recreate simulated data sets)
om_inputs = list()
for(i in 1:NROW(df.M.oms)){
  print(paste0("row ", i))
  NAA_re = gf_NAA_re
  #NAA_re$sigma_vals = df.e.r.oms$obs_sig[i]
  
  om_inputs[[i]] = make_om(selectivity = gf_selectivity, catchability = gf_catchability,   #Fhist, N1_state, as default
    M = gf_M, NAA_re = NAA_re, age_comp = "logistic-normal-miss0", brp_year = 1, eq_F_init = 0.3, 
    om_input = TRUE)
  #turn off bias correction
  om_inputs[[i]] = set_simulation_options(om_inputs[[i]], simulate_data = TRUE, simulate_process = TRUE, simulate_projection = TRUE, 
    bias_correct_pe = FALSE, bias_correct_oe = FALSE)
  #set L-N SD parameters for catch and index age comp
  om_inputs[[i]]$par$catch_paa_pars[,1] = log(L_N_sigma[df.M.oms$obs_error[i]])
  om_inputs[[i]]$par$index_paa_pars[,1] = log(L_N_sigma[df.M.oms$obs_error[i]])
  om_inputs[[i]]$data$agg_catch_sigma[] = 0.1
  om_inputs[[i]]$data$agg_index_sigma[] = index_sigma[df.M.oms$obs_error[i]]
  #change Neff so scalar doesn't affect L-N SD
  om_inputs[[i]]$data$catch_Neff[] = 1
  om_inputs[[i]]$data$index_Neff[] = 1
}

saveRDS(om_inputs, file.path(here(),"Project_0","inputs", "M_om_inputs.RDS"))
