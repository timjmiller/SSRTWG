library(wham) #make sure to use the right version of wham
library(tidyr)
library(dplyr)
library(here)
source(file.path(here(), "common_code", "make_basic_info.R"))
source(file.path(here(), "Ecov_study", "mortality", "code", "sim_management.R"))
verify_version()

naa_om_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "NAA_om_inputs.RDS"))

#SR parameters are the same for all naa_om models 
temp = fit_wham(naa_om_inputs[[1]], do.fit = FALSE, MakeADFun.silent = TRUE)
SRab = exp(c(temp$rep$log_SR_a[1], temp$rep$log_SR_b[1]))

##Estimating model inputs
#Estimating model factors
ecov_how <- c(0,1,2,4)
r_mod <- c("BH")
#create data.frame defining estimation models data.fram
df.ems <- expand.grid(ecov_how=ecov_how, r_mod=r_mod, stringsAsFactors = FALSE)
saveRDS(df.ems, file.path(here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))

#same as naa_om_setup.R
gf_info = make_basic_info()

#same as naa_om_setup.R
#selectivity is not changing
gf_selectivity <- list(model       =c(rep("logistic", gf_info$n_fleets),rep("logistic", gf_info$n_indices)),
                       initial_pars=rep(list(c(5,1)), gf_info$n_fleets + gf_info$n_indices)) #fleet, index

#M set is not changing
gf_M <- list(initial_means=rep(0.2, length(gf_info$ages)))

#same as naa_om_setup.R
#NAA_re set up that can be changed for each EM
gf_NAA_re <- list(
  N1_pars = exp(10)*exp(-(0:(length(gf_info$ages)-1))*gf_M$initial_means[1]),
  sigma   = "rec", #random about mean
  cor     = "ar1_y", #random effects are independent
  use_steepness = 0,
  recruit_model = 3, #BH
  recruit_pars  = SRab
)

gf_ecov <- list(
  label         = "Ecov",
  process_model = "ar1",
  logsigma      = cbind(rep(log(0.1), length(gf_info$year))),
  lag           = 0,
  mean          = cbind(rep(0, length(gf_info$years))),
  year          = gf_info$years,
  use_obs       = cbind(rep(1, length(gf_info$years))),
  how           = 0,
  where         = 'none'
)

em_inputs = list()
for(i in 1:NROW(df.ems)){
  print(paste0("row ", i))
  NAA_re_i    <- gf_NAA_re
  M_i         <- gf_M
  ecov_i      <- gf_ecov
  selectivity <- gf_selectivity
  
  if(df.ems$ecov_how[i]!=0){
    ecov_i$how   <- df.ems$ecov_how[i]
    ecov_i$where <- "recruit"
  }
  
  basic_info <- make_basic_info()
  basic_info$fracyr_indices[,1] = 0.25
  basic_info$fracyr_indices[,2] = 0.75
  
  em_inputs[[i]] <- prepare_wham_input(basic_info = basic_info,
                                       NAA_re = NAA_re_i,
                                       M= M_i,
                                       ecov = ecov_i, 
                                       selectivity = selectivity, 
                                       age_comp = "logistic-normal-miss0")  

  #turn off bias correction
  em_inputs[[i]] = set_simulation_options(em_inputs[[i]], 
                                          simulate_data=TRUE, 
                                          simulate_process=TRUE, 
                                          simulate_projection=TRUE, 
                                          bias_correct_pe=FALSE, 
                                          bias_correct_oe = FALSE)
  
  #change Neff so scalar doesn't affect L-N SD
  em_inputs[[i]]$data$catch_Neff[] = 1
  em_inputs[[i]]$data$index_Neff[] = 1
  em_inputs[[i]]$data$FXSPR_init[] = 0.3
  em_inputs[[i]]$data$FMSY_init[] = 0.3
}

saveRDS(em_inputs, file.path(here(),"Ecov_study", "recruitment_functions", "inputs", "em_inputs.RDS"))
