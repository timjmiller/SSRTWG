# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="devel")
# if(file.exists("c:/Users/timothy.j.miller")) {
#   library(wham, lib.loc = "c:/work/wham/old_packages/77bbd94")
# } else library(wham) #make sure to use the right version of wham
library(wham, lib.loc= "C:/Users/liz.brooks/Documents/R/win-library/4.1")
# mod$wham_version  #shows git commit # for a fitted mod
# mod <- readRDS("C:/Users/liz.brooks/R/ssrtwg_recruitment/M.ecov_m1.rds")
# mod$wham_version
# > wham_commit
# [1] "fa8e16a"    should be 77bbd94
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
source(file.path(here(), "Ecov_study", "recruitment", "code", "M_make_om.R"))  #liz
source(file.path(here(), "Ecov_study", "recruitment", "code", "M_sim_management.R"))  #liz

verify_version()  #  ==== !!!!
# Error in verify_version() : 
#   your wham commit:1.0.6.9000, commit: fa8e16ais not the required commit:77bbd94.
# Install the right commit using 
# library(devtools)
#  devtools::install_github('timjmiller/wham', dependencies=TRUE, 
#                           ref="77bbd94")

write.dir <- file.path(here(),"Ecov_study", "recruitment", "inputs")
do.write <- TRUE #liz added

naa_om_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "NAA_om_inputs.RDS"))
#SR parameters are the same for all naa_om models 
temp = fit_wham(naa_om_inputs[[1]], do.fit = FALSE, MakeADFun.silent = TRUE)
SRab = exp(c(temp$rep$log_SR_a[1], temp$rep$log_SR_b[1]))

#if(!exists("write.dir")) write.dir = getwd()
if(!dir.exists(write.dir)) dir.create(write.dir, recursive = T)
#setwd(write.dir)


#number of simulations for each scenario
#nsim = 1000
#nsim = 100

#Operating model factors
#NAA sigmas for each scenario
R_sig <- c( 0.5)  #c(0.02, 0.5)  #c(0.5, 1.0)
NAA_re_cor <- c(0.2)  #c(0, 0.5) 
NAA_re_sig <- c(0.2)

#F time series

M_sig <- 0  #0.3
M_cor <- 0

Ecov_obs_sig <- c(0.02, 0.5)  #c(0.1, 0.5)
Ecov_re_sig <- c(0.3)  #c(0.1,0.5)
Ecov_re_cor <- c( 0.5 ,0)
Ecov_effect <- c(0, 0.5, 3)  #c(0, 1.0, 1.5)  #c(0, 0.25, 0.5)
Ecov_how <- c(0,1,2,4)
Ecov_where <- "recruit"    # "recruit"   "M"  "q"
Fhist = c("H-MSY","MSY")
#how much observation error
obs_error = c("L", "H")
NAA_M_re <- c("rec","rec+1")     #c("rec","rec+1", "rec+M")

# df.oms <- expand.grid(NAA_M_re = NAA_M_re, NAA_cor=NAA_cor,  R_sig=R_sig, 
#                       Ecov_obs_sig=Ecov_obs_sig, Ecov_re_sig=Ecov_re_sig, Ecov_re_cor=Ecov_re_cor, Ecov_how=Ecov_how, Ecov_effect = Ecov_effect, 
#                      Fhist = Fhist, obs_error = obs_error, stringsAsFactors = FALSE)
df.oms <- expand.grid(Ecov_how=Ecov_how, Ecov_effect = Ecov_effect, R_sig=R_sig, NAA_M_re = NAA_M_re,  
 NAA_re_cor=NAA_re_cor,  Ecov_obs_sig=Ecov_obs_sig, Ecov_re_sig=Ecov_re_sig, Ecov_re_cor=Ecov_re_cor, 
  Fhist = Fhist, obs_error = obs_error, stringsAsFactors = FALSE)

df.oms$NAA_cor = ifelse(df.oms$NAA_re_cor==0, "iid", "ar1_y")
df.oms$NAA_re_sig = ifelse(df.oms$NAA_M_re=="rec", 0, NAA_re_sig)

#logistic-normal age comp SDs for L/H observation error (both indices AND CATCH!!!????)
L_N_sigma = c(L = 0.3, H = 1.5)
#(log) index SDs for L/H observation error 
index_sigma = c(L = 0.1, H = 0.4)

n.mods <- dim(df.oms)[1] #288 operating model scenarios
n.mods
df.oms$Model <- paste0("om_",1:n.mods)
df.oms <- df.oms %>% select(Model, everything()) # moves Model to first col
# look at model table
df.oms
#saveRDS(df.oms, file.path(here(),"Project_0", "inputs", "df.oms.RDS"))


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
# base_om = make_om(Fhist = "Fmsy", N1_state = "overfished", selectivity = gf_selectivity, 
#     M = gf_M, NAA_re = gf_NAA_re, age_comp = "logistic-normal-miss0", brp_year = 1, eq_F_init = 0.3, 
#     om_input = TRUE, max_mult_Fmsy = 1, min_mult_Fmsy = 1)

#make inputs for operating model (smaller objects to save, can recreate simulated data sets)
om_inputs = list()
#for(i in 1:NROW(df.oms)){
for(i in 1:24){
#for(i in 197:204){
  print(paste0("row ", i))
  #NAA_re = gf_NAA_re
  #M
  M_i = gf_M
  NAA_re_i = gf_NAA_re
  if(df.oms$NAA_M_re[i] == "rec"){
    NAA_re_i$sigma = "rec"
    NAA_re_i$sigma_vals = df.oms$R_sig[i] * sqrt(1-df.oms$NAA_re_cor[i]^2) 
    NAA_re_i$cor = df.oms$NAA_cor[i]   # testing iid and ar1_y
    NAA_re_i$cor_vals = c(df.oms$NAA_re_cor[i])
  } 
  if(df.oms$NAA_M_re[i] == "rec+1"){
    NAA_re_i$sigma = "rec+1"
    NAA_re_i$sigma_vals[1] = df.oms$R_sig[i] * sqrt(1-df.oms$N1_cor[i]^2)
    NAA_re_i$sigma_vals[2] = df.oms$NAA_sig[i] * sqrt(1- (df.oms$N1_cor[i]*df.oms$N2plus_cor[i])^2)
    NAA_re_i$cor = df.oms$NAA_cor[i]
    NAA_re_i$cor_vals = df.oms$NAA_re_cor[i]
    
  }
  if(df.oms$NAA_M_re[i] == "rec+M"){  #NOT planning on testing M re with NAA[1] re, but this should do it if desired
    NAA_re_i$sigma = "rec"
    NAA_re_i$sigma_vals = df.oms$R_sig[i]
    NAA_re_i$cor = df.oms$NAA_cor[i]   # testing iid and ar1_y
    NAA_re_i$cor_vals = c(df.oms$NAA_re_cor[i])
    M_i$re = "ar1_y"
    M_i$sigma_vals = M_sig * sqrt(1-M_cor^2) #defining marginal variance, but wham estimates conditional var.
    M_i$cor_vals = M_cor
  }

  ecov_i = gf_ecov
  ecov_i$logsigma = cbind(rep(log(df.oms$Ecov_obs_sig[i]), length(ecov_i$year)))
  ecov_i$process_sig_vals = df.oms$Ecov_re_sig[i]
  ecov_i$process_cor_vals = df.oms$Ecov_re_cor[i]
  if(Ecov_where=="recruit") beta.row=1
  if(Ecov_where=="Mortality") beta.row=2  #not sure what to do if where= catchability
  if(df.oms$Ecov_effect[i] < 1e-7){
    ecov_i$how =     0
    ecov_i$where = "none"
  } else {
    ecov_i$how = df.oms$Ecov_how[i]   #1
    ecov_i$where = Ecov_where  # "recruit"  #"M"
    beta_vals_i = beta_vals
    beta_vals_i[[1]][[beta.row]][] <- df.oms$Ecov_effect[i]
    ecov_i$beta_vals = beta_vals_i
  }


  Fhist. = "Fmsy"
  if(df.oms$Fhist[i] == "H-MSY") Fhist. = "H-L"
  max_mult = 2.5 # fishing at 2.5 x Fmsy
  if(Fhist. == "Fmsy") max_mult = 1
  min_mult = 1 # fishing at Fmsy
  
  om_inputs[[i]] = make_om(Fhist = Fhist., N1_state = "overfished", selectivity = gf_selectivity, 
    M = M_i, NAA_re = NAA_re_i, ecov = ecov_i, age_comp = "logistic-normal-miss0", brp_year = 1, eq_F_init = 0.3, 
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

#####  save om_inputs ====
if (do.write==TRUE)  saveRDS(om_inputs, file.path(here(), "Ecov_study", "recruitment", "inputs", "om_inputs.RDS"))
if (do.write==TRUE)  saveRDS(df.oms, file.path(here(), "Ecov_study", "recruitment", "inputs", "df.oms.RDS"))
if (do.write==TRUE)  write.csv(df.oms, file.path(here(), "Ecov_study", "recruitment", "inputs", "df.oms.csv"), row.names = FALSE)


#start out at MSY and continue
om_msy = make_om(Fhist = "Fmsy", N1_state = "overfished", selectivity = gf_selectivity, 
  M = gf_M, NAA_re = NAA_re, brp_year = 1, eq_F_init = 0.3, max_mult_Fmsy = 1)
input = om_msy#$input
input$par$log_NAA_sigma[] = -100 #no process error
temp1 <- fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE)
sim1 = temp$simulate(complete=TRUE)
sim1$NAA #NAA should stay the same throughout time series

#start out overfished and continue for 20 years
om_msy = make_om(Fhist = "H-L", N1_state = "overfished", selectivity = gf_selectivity, 
  M = gf_M, NAA_re = NAA_re, brp_year = 1, eq_F_init = 0.3, max_mult_Fmsy = 2.5)
input = om_msy#$input
input$par$log_NAA_sigma[] = -100 #no process error
temp2 <- fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE)
sim2 = temp$simulate(complete=TRUE)
sim2$NAA #NAA should stay the same throughout time series


# just using the seeds from mortality (liz)
# #I don't think we want to use the same (e.g. 1000) seeds for everything.
# set.seed(8675309)
# seeds = sample(x = (-1e9):(1e9), size = NROW(df.oms)*1000, replace = FALSE)
# seeds <- lapply(1:NROW(df.oms), function(x) seeds[(1:1000) + 1000*(x-1)])
# saveRDS(seeds, file.path(here(), "Ecov_study", "recruitment", "inputs","seeds.RDS"))
# seeds = readRDS(file.path(here::here(), "Ecov_study", "recruitment", "inputs","seeds.RDS"))



