# OM configuration 

library(tidyr)
library(dplyr)
library(here)
library(ggplot2)
library(wham)
source(file.path(here(), "common_code", "make_basic_info.R"))
source(file.path(here(), "common_code", "set_NAA.R"))
source(file.path(here(), "common_code", "set_M.R"))
source(file.path(here(), "common_code", "set_q.R"))
source(file.path(here(), "common_code", "set_ecov.R"))
source(file.path(here(), "common_code", "set_selectivity.R"))
source(file.path(here(), "common_code", "set_simulation_options.R"))
source(file.path(here(), "common_code", "get_FMSY.R"))
source(file.path(here(), "Ecov_study", "growth", "code", "make_om.R"))
source(file.path(here(), "Ecov_study", "growth", "code", "sim_management.R"))
#source("make_om.R")

write.dir <- file.path(here(),"Ecov_study", "growth", "inputs")

#if(!exists("write.dir")) write.dir = getwd()
if(!dir.exists(write.dir)) dir.create(write.dir, recursive = T)
#setwd(write.dir)


#number of simulations for each scenario
#nsim = 1000
#nsim = 100

#Operating model factors
#NAA sigmas for each scenario
NAA_sig <- c(0.3)
Ecov_obs_sig <- c(0.1)
Ecov_re_sig <- c(0.1) # This parameter is not being specified correctly in the OM. Make sure you use the right scale
# For example. 0 == (sigma = 1 internally)
Ecov_re_cor <- c(0) # This parameter is not being specified correctly in the OM. Make sure you use the right scale 
Ecov_effect <- c(0.5) # This parameter is not being specified correctly in the OM
Fhist = c("MSY")
#how much observation error
obs_error = c("L")
NAA_M_re <- c("rec")
growth_par = 1:3 # K, Linf, or L1
df.oms <- expand.grid(NAA_M_re = NAA_M_re,
  Ecov_obs_sig=Ecov_obs_sig, Ecov_re_sig=Ecov_re_sig, Ecov_re_cor=Ecov_re_cor, Ecov_effect = Ecov_effect,
  Fhist = Fhist, obs_error = obs_error, growth_par = growth_par,  stringsAsFactors = FALSE)

n.mods <- dim(df.oms)[1] 
df.oms$Model <- paste0("om_",1:n.mods)
df.oms <- df.oms %>% select(Model, everything()) # moves Model to first col
saveRDS(df.oms, file.path(here(),"Ecov_study", "growth", "inputs", "df.oms.RDS"))


gf_info = make_basic_info()
#selectivity is not changing
gf_selectivity = list(
  model = c("double-normal", "logistic"),
  initial_pars = list(c(4,-1,1.5,1.4,0,0.4), c(1,0.5)))
gf_M = list(model = "constant",
            initial_means = 0.2)
#NAA_re set up that can be changed for each OM scenario
gf_NAA_re = list(
  N1_pars = c(1e+05, 0),
  sigma = "rec", #random about mean
  cor = "iid", #random effects are independent
  recruit_model = 2,
  N1_model = 1) #defined above from naa_om_inputs
gf_ecov <- list(
  label = "Ecov",
  process_model = "ar1",
  lag = 0,
  mean = cbind(rep(0, length(gf_info$years))),
  year = gf_info$years,
  ages = list(1:10),
  use_obs = cbind(rep(1, length(gf_info$years))),
  where = list('growth'),
  where_subindex = 1 # Could be on K, Linf, or L1
)

### ------------------------------------------------------------
##  Growth configuration:
Linf <- 85
k <- 0.3
t0 <- 0
a_LW <- exp(-12.1)
b_LW <- 3.2
L_a <- Linf*(1-exp(-k*(1:10 - t0)))
W_a <- a_LW*L_a^b_LW
CV_a <- .1
gf_growth <- list(model='vB_classic', init_vals=c(k, Linf, L_a[1]),
                  SD_vals=c(CV_a*L_a[1], CV_a*L_a[10]))
gf_LW <- list(init_vals=c(a_LW, b_LW))

### --------------------------------------------------

beta_vals <- list(rep(list(matrix(0,1,length(gf_info$ages))), 7)) # R, M, Q, + 4



# Create OMs --------------------------------------------------------------
om_inputs = list()
for(i in 1:NROW(df.oms)){
  print(paste0("row ", i))
  
  # Recruitment information:
  NAA_re = gf_NAA_re
  # Natural mortality:
  M_i = gf_M

  ecov_i = gf_ecov
  ecov_i$logsigma = cbind(rep(log(df.oms$Ecov_obs_sig[i]), length(ecov_i$year)))
  if(df.oms$Ecov_effect[i] < 1e-7){
    ecov_i$how = 0
    ecov_i$where = "none"
  } else {
    ecov_i$how = 1
    ecov_i$where = "growth"
    ecov_i$where_subindex = df.oms$growth_par[i]
  }
  Fhist = "Fmsy"
  max_mult = 1.1 
  min_mult = 1 # fishing at Fmsy
  om_inputs[[i]] <-
    make_om(Fhist = Fhist, N1_state = "overfished", selectivity = gf_selectivity,
            M = M_i, NAA_re = NAA_re, ecov = ecov_i,
            growth=gf_growth, LW=gf_LW,
            age_comp = "logistic-normal-miss0", brp_year = 1, eq_F_init = 0.3,
            om_input = TRUE, max_mult_Fmsy = max_mult, min_mult_Fmsy = min_mult,
            df.oms = df.oms[i,]) 
  om_inputs[[i]] = set_simulation_options(om_inputs[[i]], simulate_data = TRUE, simulate_process = TRUE, simulate_projection = FALSE,
    bias_correct_pe = FALSE, bias_correct_oe = FALSE)
  #set L-N SD parameters for catch and index age comp
  om_inputs[[i]]$par$catch_paa_pars[,1] = log(0.3)
  om_inputs[[i]]$par$index_paa_pars[,1] = log(0.3)
  om_inputs[[i]]$data$agg_catch_sigma[] = 0.05
  om_inputs[[i]]$data$agg_index_sigma[] = 0.2
  om_inputs[[i]]$data$catch_Neff[] = 1
  om_inputs[[i]]$data$index_Neff[] = 1
}

saveRDS(om_inputs, file.path(here(), "Ecov_study", "growth", "inputs", "om_inputs.RDS"))

#I don't think we want to use the same (e.g. 1000) seeds for everything.
set.seed(8675309)
seeds = sample(x = (-1e9):(1e9), size = NROW(df.oms)*1000, replace = FALSE)
seeds <- lapply(1:NROW(df.oms), function(x) seeds[(1:1000) + 1000*(x-1)])
saveRDS(seeds, file.path(here(), "Ecov_study", "growth", "inputs","seeds.RDS"))
seeds = readRDS(file.path(here::here(), "Ecov_study", "growth", "inputs","seeds.RDS"))
