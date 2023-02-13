# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="devel")
library(wham)
library(tidyr)
library(dplyr)
library(here)
#library(doParallel)
#x <- detectCores()      
#registerDoParallel(x-1) #leave one core for other tasks
#writeLines(paste(x), "cores_detected.txt") #print how many cores were used   

source(file.path(here(), "common_code", "make_basic_info.R"))
source(file.path(here(), "common_code", "set_NAA.R"))
source(file.path(here(), "common_code", "set_M.R"))
source(file.path(here(), "common_code", "set_q.R"))
source(file.path(here(), "common_code", "set_ecov.R"))
source(file.path(here(), "common_code", "set_selectivity.R"))
source(file.path(here(), "common_code", "set_simulation_options.R"))
source(file.path(here(), "common_code", "get_FMSY.R"))
source(file.path(here(), "Ecov_study", "recruitment", "make_om.R"))
#source(file.path(here(), "Project_0", "code", "make_om.R"))

write.dir <- file.path(here(),"Ecov_study", "recruitment", "results") # create directory for analysis

if(!exists("write.dir")) write.dir = getwd()  #if we don't specify above, set as current wd
if(!dir.exists(write.dir)) dir.create(write.dir, recursive = T)  #if the write.dir directory doesn't exist, create it
setwd(write.dir)

naa_om_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "NAA_om_inputs.RDS"))
#SR parameters are the same for all naa_om models 
#temp = fit_wham(naa_om_inputs[[1]], do.fit = FALSE, MakeADFun.silent = TRUE)
#SRab = exp(c(temp$rep$log_SR_a[1], temp$rep$log_SR_b[1]))
SRab <- c(5.955694e-01, 2.404283e-05)

nsim = 100 #number of simulations for each scenario

#Operating model factors
#NAA sigmas for each scenario
R_sig <- c(0.5)
NAA_sig <- c(0.3)
#F time series

Ecov_obs_sig <- c(0.1, 0.5)
Ecov_re_sig  <- c(0.1,0.5)
Ecov_re_cor  <- c(0, 0.5) #GLB: higher level?
Ecov_effect  <- c(0, 0.25, 0.5)
Fhist        <- c("H-MSY","MSY")
#how much observation error
obs_error = c("L", "H")
NAA_R_re <- c("rec","rec+1")
df.oms <- expand.grid(NAA_R_re = NAA_R_re, 
                      Ecov_obs_sig=Ecov_obs_sig, Ecov_re_sig=Ecov_re_sig, Ecov_re_cor=Ecov_re_cor, Ecov_effect = Ecov_effect, 
                      Fhist = Fhist, obs_error = obs_error, stringsAsFactors = FALSE)

#logistic-normal age comp SDs for L/H observation error (both indices AND CATCH!!!????)
L_N_sigma = c(L = 0.3, H = 1.5)
#(log) index SDs for L/H observation error 
index_sigma = c(L = 0.1, H = 0.4)


n.mods <- dim(df.oms)[1] #288 operating model scenarios
df.oms$Model <- paste0("om_",1:n.mods)
df.oms <- df.oms %>% select(Model, everything()) # moves Model to first col
# look at model table
df.oms


gf_info = make_basic_info()

#selectivity is not changing
gf_selectivity = list(
  model = c(rep("logistic", gf_info$n_fleets),rep("logistic", gf_info$n_indices)),
  initial_pars = rep(list(c(5,1)), gf_info$n_fleets + gf_info$n_indices)) #fleet, index

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
  #sigma = "rec", #random about mean
  cor="iid", #random effects are independent
  use_steepness = 0,
  #recruit_model = 2, #random effects with a constant mean
  recruit_model = 3, #B-H
  recruit_pars = SRab) #defined above from naa_om_inputs

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


##--MAKE OMS--###############################
om_inputs = em_inputs <- list()

for(i in 1:NROW(df.oms)){
  print(paste0("row ", i))
  NAA_re = gf_NAA_re
  #M
  M_i = gf_M
  if(df.oms$NAA_R_re[i] == "rec"){
    NAA_re$sigma = "rec"
    NAA_re$sigma_vals = R_sig
  } 
  if(df.oms$NAA_R_re[i] == "rec+1"){
    NAA_re$sigma = "rec+1"
    NAA_re$sigma_vals = c(R_sig, NAA_sig)
  }
  ecov_i = gf_ecov
  ecov_i$logsigma = cbind(rep(log(df.oms$Ecov_obs_sig[i]), length(ecov_i$year)))
  ecov_i$process_sig_vals = df.oms$Ecov_re_sig[i]
  ecov_i$process_cor_vals = df.oms$Ecov_re_cor[i]
  if(df.oms$Ecov_effect[i] < 1e-7){
    ecov_i$how = 0
    ecov_i$where = "none"
  } else {
    ecov_i$how = 1
    ecov_i$where = "recruit"
    beta_vals_i = beta_vals
    beta_vals_i[[1]][[2]][] <- df.oms$Ecov_effect[i]
    ecov_i$beta_vals = beta_vals_i
  }
  
  Fhist. = "Fmsy"
  if(df.oms$Fhist[i] == "H-MSY") Fhist. = "H-L"
  max_mult = 2.5 # fishing at 2.5 x Fmsy
  if(Fhist. == "Fmsy") max_mult = 1
  min_mult = 1 # fishing at Fmsy
  
  om_inputs[[i]] = make_om(Fhist = Fhist., N1_state = "overfished", selectivity = gf_selectivity, 
                           M = M_i, NAA_re = NAA_re, ecov = ecov_i, age_comp = "logistic-normal-miss0", brp_year = 1, eq_F_init = 0.3, 
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

##################################################################
## MAKE SIMS AND EMS #############################################
##################################################################
sim_input=em_input <- list()

##--EM=OM--############
for(i in 1:n.mods){
print(paste0("model ", i))
  input <- om_inputs[[i]]
  om    <- fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE)
  
  sim_input[[i]] = lapply(1:nsim, function(x) {      #all RE and data are simulated
    sim = om$simulate(complete=TRUE)
    input_i <- input
    input_i$data <- sim
    return(input_i)
  })
  em_input[[i]] = lapply(1:nsim, function(x) {     #put in data simulated from operating model
    input_i = input
    input_i$data = sim_input[[i]][[x]]$data #put in simulated operating model data
    return(input_i)
  })
}

##--EM!=OM--#################
for(i in (n.mods+1):(2*n.mods)){
  sim_input[[i]] <- sim_input[[i-n.mods]]
  
  input <- om_inputs[[i-n.mods]]
  input[[1]]$Ecov_how=0
  em_input[[i]] <- lapply(1:nsim, function(x){
    input_i <- input
    input_i$data <- sim_input[[i]][[x]]$data
    return(input_i)
  })
}
  
saveRDS(em_input,  file.path(write.dir, "em_input.RDS"))


