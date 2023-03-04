

## if(file.exists("c:/Users/timothy.j.miller")) {
##   library(wham, lib.loc = "c:/work/wham/old_packages/77bbd94")
## } else library(wham) #make sure to use the right version of wham

## requires special growth branch
## devtools::install_github("gmoroncorrea/wham", ref='growth')
## devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="77bbd94")

library(tidyr)
library(dplyr)
library(here)
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

## verify_version()

write.dir <- file.path(here(),"Ecov_study", "growth", "inputs")

#### This breaks b/c the input are missing some growth updates,
#### so just manually specify SRab which apparently all these are
#### used for?
## naa_om_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "NAA_om_inputs.RDS"))
## ##SR parameters are the same for all naa_om models
## temp = fit_wham(naa_om_inputs[[1]], do.fit = FALSE, MakeADFun.silent = TRUE)
## SRab = exp(c(temp$rep$log_SR_a[1], temp$rep$log_SR_b[1]))
SRab <- c(5.955694e-01, 2.404283e-05)


#if(!exists("write.dir")) write.dir = getwd()
if(!dir.exists(write.dir)) dir.create(write.dir, recursive = T)
#setwd(write.dir)


#number of simulations for each scenario
#nsim = 1000
#nsim = 100

#Operating model factors
#NAA sigmas for each scenario
R_sig <- c(0.5)
NAA_sig <- c(0.3)
#F time series
M_sig <- 0.3
M_cor <- 0
Ecov_obs_sig <- c(0.1, 0.5)
Ecov_re_sig <- c(0.1,0.5)
Ecov_re_cor <- c(0, 0.5)
Ecov_effect <- c(0, 0.25, 0.5)
Fhist = c("H-MSY","MSY")
#how much observation error
obs_error = c("L", "H")
NAA_M_re <- c("rec","rec+1", "rec+M")
df.oms <- expand.grid(NAA_M_re = NAA_M_re,
  Ecov_obs_sig=Ecov_obs_sig, Ecov_re_sig=Ecov_re_sig, Ecov_re_cor=Ecov_re_cor, Ecov_effect = Ecov_effect,
  Fhist = Fhist, obs_error = obs_error, stringsAsFactors = FALSE)

df.oms <- df.oms[1,]
#logistic-normal age comp SDs for L/H observation error (both indices AND CATCH!!!????)
L_N_sigma = c(L = 0.3, H = 1.5)
#(log) index SDs for L/H observation error
index_sigma = c(L = 0.1, H = 0.4)

n.mods <- dim(df.oms)[1] #288 operating model scenarios
df.oms$Model <- paste0("om_",1:n.mods)
df.oms <- df.oms %>% select(Model, everything()) # moves Model to first col
# look at model table
df.oms
saveRDS(df.oms, file.path(here(),"Ecov_study", "growth", "inputs", "df.oms.RDS"))


gf_info = make_basic_info()



#selectivity is not changing
gf_selectivity = list(
  ## model = c(rep("logistic", gf_info$n_fleets),rep("logistic", gf_info$n_indices)),
  model = c("logistic", "logistic", "len-logistic"),
  initial_pars = list(c(5,1), c(5,1), c(20,1)))



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

### ------------------------------------------------------------
## Setup the growth structure of the model. Assume survey two
## collects marginal length compositions and uses logistic selex
## at length. Growth (k, Linf, CV1, CV10) is estimated
## internally. WAA is then derived from weight-length
## relationship (not estimated).

## These are are the growth values used to create WAA in the
## other models, so use them here
Linf <- 85
k <- 0.3
t0 <- 0
a_LW <- exp(-12.1)
b_LW <- 3.2
L <- Linf*(1-exp(-k*(1:10 - t0)))
W <- a_LW*L^b_LW
plot(1:10, L, ylim=(c(0,Linf)))
L[1] ## length at reference age (1)
CV <- .1
## growth CVs at age 1 and 10
CV*L[1]
CV*L[10]
gf_growth <- list(model='vB_classic', init_vals=c(k, Linf, L[1]),
                  est_pars=1:2, SD_vals=c(CV*L[1], CV*L[10]),
                  SD_est=1:2)
gf_LW <- list(init_vals=c(a_LW, b_LW))
## end of growth changes
### --------------------------------------------------

beta_vals <- list(rep(list(matrix(0,1,length(gf_info$ages))), 4))
# base_om = make_om(Fhist = "Fmsy", N1_state = "overfished", selectivity = gf_selectivity,
#     M = gf_M, NAA_re = gf_NAA_re, age_comp = "logistic-normal-miss0", brp_year = 1, eq_F_init = 0.3,
#     om_input = TRUE, max_mult_Fmsy = 1, min_mult_Fmsy = 1)

#make inputs for operating model (smaller objects to save, can recreate simulated data sets)
om_inputs = list()
for(i in 1:NROW(df.oms)){
  print(paste0("row ", i))
  NAA_re = gf_NAA_re
  #M
  M_i = gf_M
  if(df.oms$NAA_M_re[i] == "rec"){
    NAA_re$sigma = "rec"
    NAA_re$sigma_vals = R_sig
  }
  if(df.oms$NAA_M_re[i] == "rec+1"){
    NAA_re$sigma = "rec+1"
    NAA_re$sigma_vals = c(R_sig, NAA_sig)
  }
  if(df.oms$NAA_M_re[i] == "rec+M"){
    NAA_re$sigma = "rec"
    NAA_re$sigma_vals = R_sig
    M_i$re = "ar1_y"
    M_i$sigma_vals = M_sig * sqrt(1-M_cor^2) #defining marginal variance, but wham estimates conditional var.
    M_i$cor_vals = M_cor
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
    ecov_i$where = "M"
    beta_vals_i = beta_vals
    beta_vals_i[[1]][[2]][] <- df.oms$Ecov_effect[i]
    ecov_i$beta_vals = beta_vals_i
  }
  Fhist. = "Fmsy"
  if(df.oms$Fhist[i] == "H-MSY") Fhist. = "H-L"
  max_mult = 2.5 # fishing at 2.5 x Fmsy
  if(Fhist. == "Fmsy") max_mult = 1
  min_mult = 1 # fishing at Fmsy
  om_inputs[[i]] <-
    make_om(Fhist = Fhist., N1_state = "overfished", selectivity = gf_selectivity,
            M = M_i, NAA_re = NAA_re, ecov = ecov_i,
            age_comp = "logistic-normal-miss0", brp_year = 1, eq_F_init = 0.3,
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

### test the new growth stuff
x <- om_inputs[[1]]
x$data$lengths <- seq(1,100, by=2)
nlbins <- length(x$data$lengths)
x$data$n_lengths <- nlbins
x$data$growth_est
x$data$growth_model
x$data$index_pal <- array(data=1/nlbins, dim=c(2,40, nlbins))
x$data$index_NeffL[,2] <- 1 # dummy value to be updated later??

## left off here, crashes due to index error
## test <- fit_wham(input=x, do.fit=FALSE)

plot(lbins, test$rep$selAL[[3]][1,])


#start out at MSY and continue
om_msy = make_om(Fhist = "Fmsy", N1_state = "overfished", selectivity = gf_selectivity,
  M = gf_M, NAA_re = NAA_re, brp_year = 1, eq_F_init = 0.3, max_mult_Fmsy = 1)
input = om_msy#$input
input$par$log_NAA_sigma[] = -100 #no process error
temp <- fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE)
sim = temp$simulate(complete=TRUE)
sim$NAA #NAA should stay the same throughout time series

#start out overfished and continue for 20 years
om_msy = make_om(Fhist = "H-L", N1_state = "overfished", selectivity = gf_selectivity,
  M = gf_M, NAA_re = NAA_re, brp_year = 1, eq_F_init = 0.3, max_mult_Fmsy = 2.5)
input = om_msy#$input
input$par$log_NAA_sigma[] = -100 #no process error
temp <- fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE)
sim = temp$simulate(complete=TRUE)
sim$NAA #NAA should stay the same throughout time series

saveRDS(om_inputs, file.path(here(), "Ecov_study", "growth", "inputs", "om_inputs.RDS"))


#I don't think we want to use the same (e.g. 1000) seeds for everything.
set.seed(8675309)
seeds = sample(x = (-1e9):(1e9), size = NROW(df.oms)*1000, replace = FALSE)
seeds <- lapply(1:NROW(df.oms), function(x) seeds[(1:1000) + 1000*(x-1)])
saveRDS(seeds, file.path(here(), "Ecov_study", "growth", "inputs","seeds.RDS"))
seeds = readRDS(file.path(here::here(), "Ecov_study", "growth", "inputs","seeds.RDS"))
