# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="devel")
library(wham)
library(tidyr)
library(dplyr)
library(here)
files_to_source = list.files(file.path(here(), "common_code"), full.names=TRUE, pattern = "*.R")
sapply(files_to_source, source)
source(file.path(here(), "Project_0", "code", "make_om.R"))

###############Estimating model inputs
#Estimating model factors
SR_model = c(2,3)
M_est = c(FALSE,TRUE)
re_config = c("rec","rec+1", "M_re", "sel_re", "q_re")

#create data.frame defining estimation models data.fram
df.ems <- expand.grid(SR_model = SR_model, M_est = M_est, re_config = re_config, stringsAsFactors = FALSE)
saveRDS(df.ems, file.path(here(),"Project_0", "inputs", "df.ems.RDS"))

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
gf_M = list(initial_means = 0.2, model = "constant")

#same as naa_om_setup.R
#NAA_re set up that can be changed for each EM
gf_NAA_re = list(
  N1_pars = exp(10)*exp(-(0:(length(gf_info$ages)-1))*gf_M$initial_means[1]),
  sigma = "rec", #random about mean
  cor="iid", #random effects are independent
  use_steepness = 1,
  #recruit_model = 2, #random effects with a constant mean
  recruit_model = 3, #B-H
  recruit_pars = c(0.75,exp(10))
)

#make inputs for estimating model (smaller objects to save, can overwrinte data elements with simulated data)
em_inputs = list()
for(i in 1:NROW(df.ems)){
  print(paste0("row ", i))
  NAA_re = gf_NAA_re
  M = gf_M
  NAA_re$recruit_model = df.ems$SR_model[i] #set to estimate S-R or not
  if(df.ems$SR_model[i] == 2) NAA_re$recruit_pars = exp(10)
  #NAA_re$sigma_vals = df.oms$R_sig[i] #don't start at true values?
  if(df.ems$re_config[i] == "rec+1") { #full random effects for NAA
    NAA_re$sigma = "rec+1"
    #NAA_re$sigma_vals[2] = df.oms$NAA_sig[i] #don't start at true values?
  }
  if(df.ems$M_est[i]) { #full random effects for NAA
    M$est_ages = 1
    #NAA_re$sigma_vals[2] = df.oms$NAA_sig[i] #don't start at true values?
  }

  if(df.ems$re_config[i] == "M_re") { #full random effects for NAA
    M$re = "ar1_y"
    M$cor_vals = 0 #iid year effects
    #NAA_re$sigma_vals[2] = df.oms$NAA_sig[i] #don't start at true values?
  }
  catchability = NULL
  if(df.ems$re_config[i] == "q_re") { #full random effects for NAA
    catchability = list(re = c("iid","none"))
  }
  selectivity = gf_selectivity
  if(df.ems$re_config[i] == "sel_re") { #full random effects for NAA
    selectivity$re = c("iid", "none", "none")
    #NAA_re$sigma_vals[2] = df.oms$NAA_sig[i] #don't start at true values?
  }
  em_inputs[[i]] <- prepare_wham_input(basic_info = make_basic_info(), selectivity = selectivity, NAA_re = NAA_re, M= M,
    age_comp = "logistic-normal-miss0")  
  em_inputs[[i]] = set_M(em_inputs[[i]], M = M) #this set_M will change parameter values
  if(df.ems$re_config[i] == "M_re") {
    em_inputs[[i]]$map$M_repars = factor(c(1,NA,NA)) #still need to fix rho=0
  }
  input$data$FXSPR_init[] = 0.3
  input$data$FMSY_init[] = 0.3

  #turn off bias correction
  em_inputs[[i]] = set_simulation_options(em_inputs[[i]], simulate_data = TRUE, simulate_process = TRUE, simulate_projection = TRUE, 
    bias_correct_pe = FALSE, bias_correct_oe = FALSE)
  #change Neff so scalar doesn't affect L-N SD
  em_inputs[[i]]$data$catch_Neff[] = 1
  em_inputs[[i]]$data$index_Neff[] = 1
}

saveRDS(em_inputs, file.path(here(),"Project_0", "inputs", "em_inputs.RDS"))
