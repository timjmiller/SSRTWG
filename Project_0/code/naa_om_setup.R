# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="devel")
library(wham)
library(tidyr)
library(dplyr)
library(here)
source(file.path(here(), "common_code", "make_basic_info.R"))
source(file.path(here(), "common_code", "set_NAA.R"))
source(file.path(here(), "common_code", "get_FMSY.R"))
source(file.path(here(), "Project_0", "code", "make_om.R"))

write.dir <- file.path(here(),"Project_0", "NAA_om_scenarios")

if(!exists("write.dir")) write.dir = getwd()
if(!dir.exists(write.dir)) dir.create(write.dir, recursive = T)
setwd(write.dir)

#number of simulations for each scenario
#nsim = 1000
nsim = 100

#NAA sigmas for each scenario
R_sig <- c(0.5,1.5)
NAA_sig <- c(NA,0.25,0.5)
Fhist = c("H-MSY","MSY")
obs_error = c("L", "H")
df.oms <- expand.grid(R_sig = R_sig, NAA_sig = NAA_sig, Fhist = Fhist, 
  obs_error = obs_error, stringsAsFactors = FALSE)

SR_model = c(2,3)
M_est = c(FALSE,TRUE)
re_config = c("rec","rec+1", "M_re", "F_re", "q_re")

df.ems <- expand.grid(SR_model = SR_model, M_est = M_est, re_config = re_config, stringsAsFactors = FALSE)


# ar1_y <- c(0)
# ar1_a <- c(0)
L_N_sigma = c(L = 0.3, H = 1.5)
index_sigma = c(L = 0.1, H = 0.4)

n.mods <- dim(df.oms)[1] #48 scenarios
df.oms$Model <- paste0("om_",1:n.mods)
df.oms <- df.oms %>% select(Model, everything()) # moves Model to first col
# look at model table
df.oms

gf_selectivity = list(
  model = c(rep("logistic", groundfish_info$n_fleets),rep("logistic", groundfish_info$n_indices)),
  initial_pars = rep(list(c(5,1)), groundfish_info$n_fleets + groundfish_info$n_indices)) #fleet, index
gf_M = list(initial_means = rep(0.2, length(groundfish_info$ages)))

gf_NAA_re = list(
  N1_pars = exp(10)*exp(-(0:(length(groundfish_info$ages)-1))*gf_M$initial_means[1]),
  sigma = "rec", #random about mean
  cor="iid", #random effects are independent
  use_steepness = 1,
  #recruit_model = 2, #random effects with a constant mean
  recruit_model = 3, #B-H
  recruit_pars = c(0.75,exp(10))
)

#make inputs for operating model (smaller objects to save)
om_inputs = list()
for(i in 1:NROW(df.oms)){
  print(paste0("row ", i))
  NAA_re = gf_NAA_re
  NAA_re$sigma_vals = df.oms$R_sig[i]
  if(!is.na(df.oms$NAA_sig[i])) {
    NAA_re$sigma = "rec+1"
    NAA_re$sigma_vals[2] = df.oms$NAA_sig[i]
  }
  # obs_sigmas = c(agg_index_sigma= agg_index_sigma[df.oms$obs_error[i]], 
  #   L_N_sigma = L_N_sigma[df.oms$obs_error[i]])
  Fhist. = "Fmsy"
  max_mult = 2.5
  min_mult = 1
  if(df.oms$Fhist == "H-Fmsy") Fhist. = "H-L"
  om_inputs[[i]] = make_om(Fhist = Fhist., N1_state = "overfished", selectivity = gf_selectivity, 
    M = gf_M, NAA_re = NAA_re, age_comp = "logistic-normal-miss0", brp_year = 1, eq_F_init = 0.3, 
    om_input = TRUE, max_mult_Fmsy = max_mult, min_mult_Fmsy = min_mult)
  om_inputs[[i]]$par$catch_paa_pars[,1] = log(L_N_sigma[df.oms$obs_error[i]])
  om_inputs[[i]]$data$agg_catch_sigma[] = 0.1
  om_inputs[[i]]$data$catch_Neff[] = 1
  om_inputs[[i]]$par$index_paa_pars[,1] = log(L_N_sigma[df.oms$obs_error[i]])
  om_inputs[[i]]$data$agg_index_sigma[] = index_sigma[df.oms$obs_error[i]]
  om_inputs[[i]]$data$index_Neff[] = 1
}
om_msy = make_om(Fhist = "Fmsy", N1_state = "Fmsy", selectivity = gf_selectivity, 
  M = gf_M, NAA_re = gf_NAA_re, brp_year = 1, eq_F_init = 0.3, )

#check equilibrium
input = om$input
input$par$log_NAA_sigma[] = -100 #no process error
temp <- fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE)
sim = temp$simulate(complete=TRUE)
sim$NAA #NAA should stay the same throughout time series

#start out overfished and overfishing the whole time series
om_overfished = make_om(Fhist = "H", N1_state = "overfished", selectivity = gf_selectivity, 
  M = gf_M, NAA_re = gf_NAA_re, brp_year = 1, eq_F_init = 0.3)
om_overfished$rep$F
input = om_overfished$input
input$par$log_NAA_sigma[] = -100 #no process error
temp <- fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE)
sim = temp$simulate(complete=TRUE)
sim$NAA #NAA should stay the same throughout time series
  
om_input = om$input


