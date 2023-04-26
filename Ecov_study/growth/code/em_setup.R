## # devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="devel")
## if(file.exists("c:/Users/timothy.j.miller")) {
##   library(wham, lib.loc = "c:/work/wham/old_packages/77bbd94")
## } else library(wham) #make sure to use the right version of wham
library(wham)
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
source(file.path(here(), "Ecov_study", "growth", "code", "make_om.R"))
source(file.path(here(), "Ecov_study", "growth", "code", "sim_management.R"))
## verify_version()


### this seems unnecessary -cole
## naa_om_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "NAA_om_inputs.RDS"))
## # temp = sapply(naa_om_inputs, function(x) {
## #   temp = fit_wham(x, do.fit = FALSE, MakeADFun.silent = TRUE)
## #   return(temp$rep$log_SR_a)
## # })
## #SR parameters are the same for all naa_om models
## temp = fit_wham(naa_om_inputs[[1]], do.fit = FALSE, MakeADFun.silent = TRUE)
## SRab = exp(c(temp$rep$log_SR_a[1], temp$rep$log_SR_b[1]))

###############Estimating model inputs
#Estimating model factors
SR_model = c(2)
growth_est = c(TRUE, FALSE)[1]
re_config = c("rec","rec+1", "rec+M")[1]
Ecov_est = c(TRUE)
growth_method = c('growth', 'LAA')[1]
growth_re_config = c(NA)
# create df em part 1 (only when Ecov on). TODO: Add Ecov on LAA? 
df.ems.1 <- expand.grid(growth_est = growth_est, re_config = re_config, Ecov_est = Ecov_est, growth_method= growth_method, 
                        growth_re_config = growth_re_config, stringsAsFactors = FALSE)
# Config for growth or LAA random effects (only when Ecov = FALSE):
# Caution when creating this data.frame, some combinations are not possible:
growth_method_1 = c('growth', 'growth')
growth_re_config_1 = c(NA, 'growth_re')
growth_method_2 = c('LAA', 'LAA')
growth_re_config_2 = c(NA, 'LAA_re')
df1 = data.frame(growth_method = c(growth_method_1, growth_method_2),
                 growth_re_config = c(growth_re_config_1, growth_re_config_2))
df.ems.2 <- expand.grid(growth_est = growth_est, re_config = re_config, Ecov_est = FALSE, stringsAsFactors = FALSE)
df.ems.2 = df.ems.2 %>% dplyr::cross_join(df1)
# Merge both df:
# NOT SURE IF THIS IS THE MOST EFFICIENT WAY TO CREATE THE EM DF:
df.ems <- rbind(df.ems.1, df.ems.2)
saveRDS(df.ems, file.path(here(),"Ecov_study", "growth", "inputs", "df.ems.RDS"))

#same as naa_om_setup.R
gf_info = make_basic_info()

#same as naa_om_setup.R
#selectivity is not changing
## gf_selectivity = list(
##   model = c(rep("logistic", gf_info$n_fleets),rep("logistic", gf_info$n_indices)),
##   initial_pars = rep(list(c(5,1)), gf_info$n_fleets + gf_info$n_indices)) #fleet, index

gf_selectivity = list(
  ## model = c(rep("logistic", gf_info$n_fleets),rep("logistic", gf_info$n_indices)),
  model = c("logistic", "logistic", "len-logistic"),
  initial_pars = list(c(5,1), c(5,1), c(65,4)))


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
  use_steepness = 0,
  recruit_model = 2 #random effects with a constant mean
  #recruit_model = 3#, #B-H
  #recruit_pars = SRab
)
## this is the same one from the OM, is this right??
## SRab <- c(5.955694e-01, 2.404283e-05)
## gf_NAA_re = list(
##   N1_pars = exp(10)*exp(-(0:(length(gf_info$ages)-1))*gf_M$initial_means[1]),
##   sigma = "rec", #random about mean
##   cor="iid", #random effects are independent
##   use_steepness = 0,
##   #recruit_model = 2, #random effects with a constant mean
##   recruit_model = 3, #B-H
##   recruit_pars = SRab) #defined above from naa_om_inputs


## gf_ecov <- list(
##   label = "Ecov",
##   process_model = "ar1",
##   logsigma = cbind(rep(log(0.1), length(gf_info$year))),
##   lag = 0,
##   mean = cbind(rep(0, length(gf_info$years))),
##   year = gf_info$years,
##   use_obs = cbind(rep(1, length(gf_info$years))),
##   where = "none",
##   how = 0
## )

## updated for growth -cole
gf_ecov <- list(
  label = "Ecov",
  process_model = "ar1",
  logsigma = cbind(rep(log(0.1), length(gf_info$year))),
  lag = 0,
  mean = cbind(rep(0, length(gf_info$years))),
  year = gf_info$years,
  use_obs = cbind(rep(1, length(gf_info$years))),
  ## process_mean_vals = 0,
  where = "none", ## updated to growth below
  where_subindex = 3, # on L1, growth specific input
  how=0
)

Linf <- 85
k <- 0.3
t0 <- 0
a_LW <- exp(-12.1)
b_LW <- 3.2
L <- Linf*(1-exp(-k*(1:10 - t0)))
W <- a_LW*L^b_LW
L[1] ## length at reference age (1)
CV <- .1
## growth CVs at age 1 and 10
gf_growth <- list(model='vB_classic', re = rep('none', times = 3), 
                  init_vals=c(k, Linf, L[1]),
                  est_pars=1:3, SD_vals=c(CV*L[1], CV*L[10]),
                  SD_est=1:2) 
gf_LW <- list(init_vals=c(a_LW, b_LW))
gf_LAA = list(LAA_vals = L, est_pars = 1:10, re = 'none', SD_vals = c(CV*L[1], CV*L[10]), SD_est=1:2) # fixing SD1

#make inputs for estimating model (smaller objects to save, can overwrinte data elements with simulated data)
em_inputs = list()
## df.ems <- df.ems[8,]
for(i in 1:NROW(df.ems)){
  print(paste0("row ", i))
  NAA_re_i = gf_NAA_re
  M_i = gf_M # mapped to single value below and not estimated for now
  growth_i <- gf_growth
  LAA_i = gf_LAA
  M_i$est_ages = 1:length(gf_info$ages)
  ecov_i = gf_ecov
  selectivity = gf_selectivity

  #NAA_re$sigma_vals = df.oms$R_sig[i] #don't start at true values?
  if(df.ems$re_config[i] == "rec+1") { #full random effects for NAA
    NAA_re_i$sigma = "rec+1"
    #NAA_re$sigma_vals[2] = df.oms$NAA_sig[i] #don't start at true values?
  }
  if(df.ems$re_config[i] == "rec+M") { #full random effects for NAA
    M_i$re = "ar1_y"
    #NAA_re$sigma_vals[2] = df.oms$NAA_sig[i] #don't start at true values?
  }
  ## estimate fixed effect growth pars?
  if(!df.ems$growth_est[i]) {
    growth_i$est_pars <- NULL
    growth_i$SD_est <- NULL
    LAA_i$est_pars <- NULL
    LAA_i$SD_est <- NULL
  }
  # Change Ecov information:
  if(df.ems$Ecov_est[i] & df.ems$growth_method[i] == 'growth'){
    ecov_i$how = 1
    ecov_i$where = df.ems$growth_method[i]
  }
  if(df.ems$Ecov_est[i] & df.ems$growth_method[i] == 'LAA'){
    ecov_i$how = 1
    ecov_i$where = df.ems$growth_method[i]
  }  
  # Change growth information:
  if(df.ems$growth_method[i] == 'growth') { 
    LAA_i = NULL
    if(!is.na(df.ems$growth_re_config[i])) { 
      growth_i$re[3] = "ar1_y" # always on L1
    }
  }
  # Change LAA information
  if(df.ems$growth_method[i] == 'LAA') { 
    growth_i = NULL
    if(!is.na(df.ems$growth_re_config[i])) { 
      LAA_i$re = '2dar1'
    }
  }

  basic_info <- make_basic_info()
  basic_info$fracyr_indices[,1] = 0.25
  basic_info$fracyr_indices[,2] = 0.75
  ## cole added for growth, matches OM
  ny <- length(basic_info$years)
  basic_info$lengths <- seq(1,100, by=2)
  nlbins <- length(basic_info$lengths)
  basic_info$n_lengths <- nlbins
  ## Turn on marginal lengths comps for survey two and turn off the
  ## age comps
  basic_info$use_index_paa <- cbind(rep(1,ny), rep(0,ny))
  basic_info$index_pal <- array(data=1/nlbins, dim=c(2,40, nlbins))
  basic_info$use_index_pal <- cbind(rep(0,ny), rep(1,ny))
  ## for now assuming multinomial N=100
  basic_info$index_NeffL <- cbind(rep(1,ny), rep(100,ny))

  em_inputs[[i]] <-
    prepare_wham_input(basic_info = basic_info,
                       growth=growth_i, LW=gf_LW, LAA= LAA_i,
                       selectivity = selectivity, NAA_re = NAA_re_i, M= M_i,
                       ecov = ecov_i, age_comp = "logistic-normal-miss0",
                       len_comp='multinomial')
  #have to do this after the input is made
  #source("c:/work/wham_master/wham/R/set_ecov.R")
  ## em_inputs[[i]] <- set_ecov(em_inputs[[i]], ecov = ecov_i)
  em_inputs[[i]] = set_M(em_inputs[[i]], M = M_i) #this set_M will change parameter values. needed for iid M_re
  #let's estimate AR1 for all these models
  # if(df.ems$re_config[i] == "rec+M") {
  #   em_inputs[[i]]$map$M_repars = factor(c(1,NA,NA)) #still need to fix rho which is set to 0 above
                                        # }
  ## estimate single M or fixed at truth? turned off for now
  em_inputs[[i]]$map$M_a = factor(NA*rep(1,length(em_inputs[[i]]$par$M_a)))
  ## ## turn on growht estimates?
  ## if(df.ems$growth_est[i]) {
  ##    ## em_inputs[[i]]$map$growth_a = factor(rep(1,length(em_inputs[[i]]$par$growth_a)))
  ## }
  ## seems like leaving where="none" shoudl map these off but doesn't..??
  if(!df.ems$Ecov_est[i]){
    em_inputs[[i]]$random = NULL # will turn off Ecov as well
    # is this still required?:
    em_inputs[[i]]$map$Ecov_re <- factor(NA*em_inputs[[i]]$par$Ecov_re)
    em_inputs[[i]]$par$Ecov_re <- 0*em_inputs[[i]]$par$Ecov_re
    em_inputs[[i]]$map$Ecov_process_pars <- factor(NA*em_inputs[[i]]$par$Ecov_process_pars)
    if(!is.na(df.ems$growth_re_config[i])) { # RE on growth or LAA
      em_inputs[[i]]$random <- df.ems$growth_re_config[i]
      ##  em_inputs[[i]]$random <- 'log_NAA'
    }
  } else {
    ## estimate Ecov effect on L1 growth? for now using penalized
    ## ML for speed
    em_inputs[[i]]$random <- 'Ecov_re'
    ## pars are basically: mean, rho, sigma; need to map off
    ## sigma for now
    ## em_inputs[[i]]$map$Ecov_process_pars <- factor(c(1,2,NA))
  }
  ## can't really esitmate growth SD1 so map it off
  #if(df.ems$growth_est[i]) {
  #  em_inputs[[i]]$map$SDgrowth_par <- factor(c(NA,1))
  #}

  #turn off bias correction
  em_inputs[[i]] = set_simulation_options(em_inputs[[i]], simulate_data = TRUE, simulate_process = TRUE, simulate_projection = TRUE,
    bias_correct_pe = FALSE, bias_correct_oe = FALSE)
  #change Neff so scalar doesn't affect L-N SD
  em_inputs[[i]]$data$catch_Neff[] = 1
  ## not sure why this has to be done
  em_inputs[[i]]$data$index_Neff[] <- 1
  #em_inputs[[i]]$data$index_NeffL <-  basic_info$index_NeffL
  em_inputs[[i]]$data$FXSPR_init[] = 0.3
  em_inputs[[i]]$data$FMSY_init[] = 0.3

  em_inputs[[i]]$map$log_NAA_sigma <- factor(NA*em_inputs[[i]]$par$log_NAA_sigma)
}


df.ems
# what is this for?:
#test1 <- fit_wham(em_inputs[[1]], do.fit=FALSE)
#test2 <- fit_wham(em_inputs[[2]], do.fit=FALSE)
#p1 <- test1$par %>% names %>% unique
#p2 <- test2$par %>% names %>% unique
## extra pars estimated when Ecov turned on
#p1[which(! p1 %in% p2)]

## beta set to zero for initial values
# what is this for?:
#em_inputs[[2]]$par$Ecov_beta %>% max
#em_inputs[[1]]$par$Ecov_beta %>% max
#em_inputs[[1]]$map$Ecov_beta ## all ages mapped to same value
#em_inputs[[2]]$map$Ecov_beta ## all ages mapped to same value
#em_inputs[[1]]$map$Ecov_process_pars
#em_inputs[[2]]$map$Ecov_process_pars
#em_inputs[[1]]$random
#em_inputs[[2]]$random

saveRDS(em_inputs, file.path(here(),"Ecov_study", "growth", "inputs", "em_inputs.RDS"))
