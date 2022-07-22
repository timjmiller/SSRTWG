# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="devel")
if(file.exists("c:/Users/timothy.j.miller")) {
  library(wham, lib.loc = "c:/work/wham/old_packages/77bbd94")
} else library(wham) #make sure to use the right version of wham
library(tidyr)
library(dplyr)
library(here)
files_to_source = list.files(file.path(here(), "common_code"), full.names=TRUE, pattern = "*.R")
sapply(files_to_source, source)
source(file.path(here(), "Project_0", "code", "make_om.R"))

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
gf_M = list(initial_means = rep(0.2, length(gf_info$ages)), model = "age-specific")

#same as naa_om_setup.R
#NAA_re set up that can be changed for each EM
gf_NAA_re = list(
  N1_pars = exp(10)*exp(-(0:(length(gf_info$ages)-1))*gf_M$initial_means[1]),
  sigma = "rec", #random about mean
  cor="iid", #random effects are independent
  use_steepness = 0
  #recruit_model = 2, #random effects with a constant mean
  #recruit_model = 3#, #B-H
  #recruit_pars = SRab
)

#make inputs for estimating model (smaller objects to save, can overwrinte data elements with simulated data)
em_inputs = list()
for(i in 1:NROW(df.ems)){
  print(paste0("row ", i))
  NAA_re = gf_NAA_re
  M = gf_M
  NAA_re$recruit_model = df.ems$SR_model[i] #set to estimate S-R or not
  if(df.ems$SR_model[i] == 3) NAA_re$recruit_pars = SRab
  if(df.ems$SR_model[i] == 2) NAA_re$recruit_pars = exp(10)
  #NAA_re$sigma_vals = df.oms$R_sig[i] #don't start at true values?
  if(df.ems$re_config[i] == "rec+1") { #full random effects for NAA
    NAA_re$sigma = "rec+1"
    #NAA_re$sigma_vals[2] = df.oms$NAA_sig[i] #don't start at true values?
  }
  if(df.ems$M_est[i]) { #estimate mean M
    M$est_ages = 1:length(gf_info$ages)
    #NAA_re$sigma_vals[2] = df.oms$NAA_sig[i] #don't start at true values?
  }

  if(df.ems$re_config[i] == "M_re") { #full random effects for NAA
    M$re = "ar1_y" #this will make it so that there is just one random effect, need to set and fix cor values for iid assumption
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
    catchability = catchability, age_comp = "logistic-normal-miss0")  
  #have to do this after the input is made
  em_inputs[[i]] = set_M(em_inputs[[i]], M = M) #this set_M will change parameter values. needed for iid M_re
  if(df.ems$re_config[i] == "M_re") { 
    em_inputs[[i]]$map$M_repars = factor(c(1,NA,NA)) #still need to fix rho which is set to 0 above
  }
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

#saveRDS(em_inputs, file.path(here(),"Project_0", "inputs", "em_inputs.RDS"))

#add a few more that are relevant to different operating models
#em_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "em_inputs.RDS"))
#df.ems = readRDS(file.path(here::here(),"Project_0","inputs", "df.ems.RDS"))


#df.ems = rbind(df.ems, cbind(SR_model = c(2:3,2:3), M_est = c(FALSE,FALSE,TRUE,TRUE), re_config = "rec+1"))
df.ems = rbind(df.ems, cbind.data.frame(SR_model = c(2:3,2:3), M_est = c(FALSE,FALSE,TRUE,TRUE), re_config = "M_re"))
df.ems = rbind(df.ems, cbind.data.frame(SR_model = c(2:3,2:3), M_est = c(FALSE,FALSE,TRUE,TRUE), re_config = "sel_re"))
df.ems = rbind(df.ems, cbind.data.frame(SR_model = c(2:3,2:3), M_est = c(FALSE,FALSE,TRUE,TRUE), re_config = "q_re"))

df.ems$M_re_cor = NA
df.ems$M_re_cor[9:12] = "iid"
df.ems$M_re_cor[21:24] = "ar1_y"
df.ems$sel_re_cor = NA
df.ems$sel_re_cor[13:16] = "iid"
df.ems$sel_re_cor[25:28] = "ar1_y"
df.ems$q_re_cor = NA
df.ems$q_re_cor[17:20] = "iid"
df.ems$q_re_cor[29:32] = "ar1"
df.ems
saveRDS(df.ems, file.path(here(),"Project_0", "inputs", "df.ems.RDS"))


for(i in 21:NROW(df.ems)){
  print(paste0("row ", i))
  NAA_re = gf_NAA_re
  M = gf_M
  NAA_re$recruit_model = df.ems$SR_model[i] #set to estimate S-R or not
  if(df.ems$SR_model[i] == 3) NAA_re$recruit_pars = SRab
  if(df.ems$SR_model[i] == 2) NAA_re$recruit_pars = exp(10)
  #if(df.ems$SR_model[i] == 2) NAA_re$recruit_pars = exp(10)
  #NAA_re$sigma_vals = df.oms$R_sig[i] #don't start at true values?
  if(df.ems$M_est[i]) { #estimate mean M
    M$est_ages = 1:length(gf_info$ages)
    #NAA_re$sigma_vals[2] = df.oms$NAA_sig[i] #don't start at true values?
  }

  if(df.ems$re_config[i] == "M_re") { #full random effects for NAA
    M$re = "ar1_y" #this is always used but cor may be 0 for iid option in rows 9-12 above
  }
  catchability = NULL
  if(df.ems$re_config[i] == "q_re") { #full random effects for NAA
    catchability = list(re = c(df.ems$q_re_cor[i],"none")) #iid or ar1
  }
  selectivity = gf_selectivity
  if(df.ems$re_config[i] == "sel_re") { #full random effects for NAA
    selectivity$re = c(df.ems$sel_re_cor[i], "none", "none") #iid or ar1_y
    #NAA_re$sigma_vals[2] = df.oms$NAA_sig[i] #don't start at true values?
  }
  em_inputs[[i]] <- prepare_wham_input(basic_info = make_basic_info(), selectivity = selectivity, NAA_re = NAA_re, M= M,
    catchability = catchability, age_comp = "logistic-normal-miss0")
  if(df.ems$M_est[i]) { #estimate mean M
    em_inputs[[i]]$map$M_a = factor(rep(1,length(em_inputs[[i]]$par$M_a)))
  }
  #below isn't necessary for rows 21-32.
  #em_inputs[[i]] = set_M(em_inputs[[i]], M = M) #this set_M will change parameter values
  #if(df.ems$re_config[i] == "M_re") {
  #}

  #turn off bias correction
  em_inputs[[i]] = set_simulation_options(em_inputs[[i]], simulate_data = TRUE, simulate_process = TRUE, simulate_projection = TRUE, 
    bias_correct_pe = FALSE, bias_correct_oe = FALSE)
  #change Neff so scalar doesn't affect L-N SD
  em_inputs[[i]]$data$catch_Neff[] = 1
  em_inputs[[i]]$data$index_Neff[] = 1
  em_inputs[[i]]$data$FXSPR_init[] = 0.3
  em_inputs[[i]]$data$FMSY_init[] = 0.3
}

saveRDS(em_inputs, file.path(here(),"Project_0", "inputs", "em_inputs.RDS"))
