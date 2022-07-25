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
source(file.path(here(), "common_code", "set_simulation_options.R"))
source(file.path(here(), "common_code", "get_FMSY.R"))
source(file.path(here(), "Project_0", "code", "make_om.R"))
source(file.path(here(), "Project_0", "code", "sim_management.R"))
verify_version()

write.dir <- file.path(here(),"Project_0", "inputs")

naa_om_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "NAA_om_inputs.RDS"))
#SR parameters are the same for all operating models 
temp = fit_wham(naa_om_inputs[[1]], do.fit = FALSE, MakeADFun.silent = TRUE)
SRab = exp(c(temp$rep$log_SR_a[1], temp$rep$log_SR_b[1]))

#if(!exists("write.dir")) write.dir = getwd()
if(!dir.exists(write.dir)) dir.create(write.dir, recursive = T)
setwd(write.dir)

#Operating model factors
#NAA sigmas for each scenario
R_sig <- c(0.5)
#F time series

Sel_sig <- c(0.1,0.5)
Sel_cor <- c(0, 0.9)

Fhist = c("H-MSY","MSY")
#how much observation error
obs_error = c("L", "H")
df.Sel.oms <- expand.grid(R_sig = R_sig, Sel_sig = Sel_sig, Sel_cor = Sel_cor, Fhist = Fhist, 
  obs_error = obs_error, stringsAsFactors = FALSE)

#logistic-normal age comp SDs for L/H observation error (both indices AND CATCH!!!????)
L_N_sigma = c(L = 0.3, H = 1.5)
#(log) index SDs for L/H observation error 
index_sigma = c(L = 0.1, H = 0.4)

n.mods <- dim(df.Sel.oms)[1] #24 operating model scenarios
df.Sel.oms$Model <- paste0("om_",1:n.mods)
df.Sel.oms <- df.Sel.oms %>% select(Model, everything()) # moves Model to first col
# look at model table
df.Sel.oms
saveRDS(df.Sel.oms, file.path(here(),"Project_0", "inputs", "df.Sel.oms.RDS"))


gf_info = make_basic_info()

#selectivity is not changing
gf_selectivity = list(
  model = c(rep("logistic", gf_info$n_fleets),rep("logistic", gf_info$n_indices)),
  initial_pars = rep(list(c(5,1)), gf_info$n_fleets + gf_info$n_indices),
  re = c("2dar1", "none", "none")) #fleet, index

#M set is not changing
gf_M = list(model = "age-specific", 
  initial_means = rep(0.2, length(gf_info$ages))#,
  #re = "ar1_y"#, # This is needed to set up operating models with a single annual re, iid or ar1_y
  #sigma_vals = 0.1,
  #cor_vals = 0
  )

#NAA_re set up that can be changed for each OM scenario
gf_NAA_re = list(
  N1_pars = exp(10)*exp(-(0:(length(gf_info$ages)-1))*gf_M$initial_means[1]),
  sigma = "rec", #random about mean
  cor="iid", #random effects are independent
  use_steepness = 0,
  #recruit_model = 2, #random effects with a constant mean
  recruit_model = 3, #B-H
  recruit_pars = SRab) #defined above from naa_om_inputs

#make inputs for operating model (smaller objects to save, can recreate simulated data sets)
om_inputs = list()
for(i in 1:NROW(df.Sel.oms)){
  print(paste0("row ", i))
  NAA_re = gf_NAA_re
  NAA_re$sigma_vals = df.Sel.oms$R_sig[i]
  
  Fhist. = "Fmsy"
  max_mult = 2.5 # fishing at 2.5 x Fmsy
  min_mult = 1 # fishing at Fmsy
  
  #M
  Sel_i = gf_selectivity
  cond.sig = df.Sel.oms$Sel_sig[i] * sqrt(1-df.Sel.oms$Sel_cor[i]^2) #defining marginal variance, but wham estimates conditional var.
  Sel_i$initial_sigma = rep(cond.sig, gf_info$n_fleets + gf_info$n_indices)
  Sel_i$initial_cor = list(c(0,df.Sel.oms$Sel_cor[i]), 0, 0) #independent slope and a50, but perhaps correlation over time.
  #Sel_i$map_sigma = c(1,NA,NA) #doesn't matter for operating models
  #Sel_i$map_cor = list(1, NA, NA)

  if(df.Sel.oms$Fhist[i] == "H-MSY") Fhist. = "H-L"
  om_inputs[[i]] = make_om(Fhist = Fhist., N1_state = "overfished", selectivity = Sel_i, 
    M = gf_M, NAA_re = NAA_re, age_comp = "logistic-normal-miss0", brp_year = 1, eq_F_init = 0.3, 
    om_input = TRUE, max_mult_Fmsy = max_mult, min_mult_Fmsy = min_mult)
  #temp = fit_wham(om_inputs[[1]], do.fit = F)
  #temp$fn()
  #temp$report()$selpars[[1]]
  om_inputs[[i]]$map$selpars_re
  #matrix(temp$simulate()$selpars_re, ncol = 2)
  #turn off bias correction
  om_inputs[[i]] = set_simulation_options(om_inputs[[i]], simulate_data = TRUE, simulate_process = TRUE, simulate_projection = TRUE, 
    bias_correct_pe = FALSE, bias_correct_oe = FALSE)
  #set L-N SD parameters for catch and index age comp
  om_inputs[[i]]$par$catch_paa_pars[,1] = log(L_N_sigma[df.Sel.oms$obs_error[i]])
  om_inputs[[i]]$par$index_paa_pars[,1] = log(L_N_sigma[df.Sel.oms$obs_error[i]])
  om_inputs[[i]]$data$agg_catch_sigma[] = 0.1
  om_inputs[[i]]$data$agg_index_sigma[] = index_sigma[df.Sel.oms$obs_error[i]]
  #change Neff so scalar doesn't affect L-N SD
  om_inputs[[i]]$data$catch_Neff[] = 1
  om_inputs[[i]]$data$index_Neff[] = 1
}

saveRDS(om_inputs, file.path(here(),"Project_0","inputs", "Sel_om_inputs.RDS"))
