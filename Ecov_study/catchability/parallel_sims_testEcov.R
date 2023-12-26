# Run a handful of tests that:
# 1) Change OM ecov_beta values to cover larger range - help select larger values for broader testing
# 2) Change OM ecov observation error levels to include 0.0001 (OM trend follow observations)
# Run 20 sims for each test
# Pick other OM settings to have low sigma, H-L fishing history, and high correlation
# Environmental covariate has an increasing mean

# Load packages & source functions used in simulation testing
## Packages
library(tidyverse)
library(wham)
library(TAF)
library(varhandle)
library(doParallel)
library(here)

## OM/EM setup functions
source(here::here("Ecov_study", "catchability", "make_om.R")) # Use revised copy that has option not to include a S-R relationship
source(file.path(here::here(),"common_code", "make_basic_info.R")) 
source(here::here("common_code", "set_ecov.R")) 
source(here::here("common_code", "set_NAA.R")) 
source(here::here("common_code", "set_q.R"))
source(here::here("common_code", "set_selectivity.R")) 
source(here::here("common_code", "set_M.R")) 
source(here::here("common_code", "get_FMSY.R"))

## Simulation testing functions
source(here::here("Ecov_study","catchability","simOM.R"))
source(here::here("Ecov_study","catchability","simTestWHAM.R"))

## Processing functions
source(here::here("Ecov_study", "catchability", "postprocess_simTestWHAM.R"))
source(here::here("Ecov_study", "catchability", "combinePerfMet.R"))
source(here::here("Ecov_study", "catchability", "plotResults.R"))


### Factorial OM settings
##### OM
# Ecov process & effect magnitude set up
Ecov_process_sig <- c(0.1, 0.5)
Ecov_process_cor <- c(0, 0.5)
Ecov_process_obs_sig <- c(0.0001, 0.1, 0.5) # Add option so obs followed perfectly
Ecov_effect <- c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0) # beta

# Fishing history
F_hist <- c("H-L", "Fmsy") # High-then FMSY vs. FMSY for entire history

# Observation error
ageComp_sig <- c(0.3, 1.5) # Fleet and index age comp treated with same error
log_index_sig <- c(0.1, 0.4) # agg_index_sigma
log_catch_sig <- 0.1 # agg_catch_sigma

# Factorial combinations
OMsetup <- expand.grid(Ecov_process_sig = Ecov_process_sig, Ecov_process_cor = Ecov_process_cor, Ecov_process_obs_sig = Ecov_process_obs_sig, Ecov_effect = Ecov_effect, F_hist = F_hist, ageComp_sig = ageComp_sig, log_index_sig = log_index_sig, log_catch_sig = log_catch_sig)

OMname <- 1:nrow(OMsetup)

OMsetup <- cbind(OMname, OMsetup)

OMsetup[,"F_hist"] <- as.character(OMsetup[,"F_hist"])


### Factorial EM settings
##### EM
# EM factorial settings
miss_season <- c("BOTH", "ONE", "NONE") # Number of seasons with correctly specified effect, assumes 2 survey indices (spring & fall)
miss_q <- c("NoEcov", "Ecov", "qRand", "qRandEcov") # Catchability setup for EM
EMsetup <- expand.grid(miss_season = miss_season, miss_q = miss_q)


# ### OM/EM setup - environmental covariate updated to increase from 0 to 5 over 40 years
# dir.create(here::here("Ecov_study", "catchability", "Results"))
# 
# n_indices <- 2
# n_ages <- 10
# n_years <- 40
# Ecov_mean_trend <- seq(0, 5, by = 5/39)
# n_fleets = 1
# 
# # Assume logistic selectivity for the fleet and all indices
# sel_list <- list(model = c(rep("logistic", n_fleets),rep("logistic", n_indices)),
#                  initial_pars = lapply(1:(n_fleets + n_indices), function(x) c(5,1)),
#                  n_selblocks = (n_fleets + n_indices))
# 
# # Starting mean M estimate at 0.2, estimate constant M for all ages/years
# M_list <-  list(initial_means = rep(0.2, n_ages))
# 
# # NAA_re = list(recruit_model = 3) # This required for Beverton-Holt S-R
# 
# # Loop over testing OMs only
# for(iom in 1:nrow(OMsetup)){
#   print(iom)
#   ##### Set up generic input #####
#   # Set up input without ecov
#   input <- make_om(Fhist = OMsetup[iom,"F_hist"],
#                    N1_state = "Fmsy", # Default, could also pick "overfished" or "unfished"
#                    selectivity = sel_list,
#                    M = M_list,
#                    catchability = NULL,
#                    NAA_re = NULL)
#   # age_comp, # Default = "logistic-normal-miss0"
#   # brp_year, # Default = 1
#   # eq_F_init, # Default = 0.3
#   # om_input = TRUE, # Don't fit model, only structure generic input
#   # max_mult_Fmsy, # Default = 2.5
#   # min_mult_Fmsy, # Default = 1
#   # F_change_time # Default = 0.5, sets when F changes (in H-L or L-H F_hist scenarios)
# 
#   # Observation error - set L-N SD parameters for catch and index age comps # pulled from project 0 q_om_setup.R
#   input$par$catch_paa_pars[,1] <- log(OMsetup[iom, "ageComp_sig"])
#   input$par$index_paa_pars[,1] <- log(OMsetup[iom, "ageComp_sig"])
#   input$data$agg_catch_sigma[] <- OMsetup[iom,"log_catch_sig"]
#   input$data$agg_index_sigma[] <- OMsetup[iom,"log_index_sig"]
#   input$data$catch_Neff[] <- 1 # Change Neff so scalar doesn't affect L-N SD
#   input$data$index_Neff[] <- 1
# 
#   ##### Set up OM ecov #####
#   # Ecov setup
#   Ecov <- list(label = "Ecov",
#                mean = matrix(Ecov_mean_trend, ncol =1), # Mean = 0
#                logsigma = matrix(log(OMsetup$Ecov_process_obs_sig[iom]), n_years, ncol=1),
#                year = input$years,
#                use_obs = matrix(rep(TRUE, n_years), ncol = 1),
#                lag = 0,
#                where = "q", # Where/how/indices settings need to change if we do sensitivity runs
#                how = 1,
#                indices = list(2), # Ecov only impact fall index catchability
#                process_model = "ar1", # End generic inputs for ecov
#                process_mean_vals = Ecov_mean_trend, # Mean = 0
#                process_sig_vals = rep(OMsetup$Ecov_process_sig[iom]),
#                process_cor_vals = rep(OMsetup$Ecov_process_cor[iom]),
#                beta_vals = list(append(rep(list(matrix(rep(0, n_ages), nrow = 1)), 2), # No impact on R,M
#                                        rep(list(matrix(rep(OMsetup$Ecov_effect[iom], n_ages), nrow = 1)), 2)))) # impact on q by index (here 2 indices)
# 
#   inputOM <- set_ecov(input = input, ecov = Ecov) #!!! something wrong with ecov$beta_vals[[j]][[n]] indexing
# 
# 
#   # Finish setup & save OM
#   inputOM$model_name <- paste0("OM_", paste(OMsetup[iom, ], collapse = "_")) # Set up OM name based on OMsetup
#   OM <- NULL # Clear prior OM object just in case
#   OM <- fit_wham(input = inputOM, do.fit = FALSE, MakeADFun.silent = TRUE)
#   mkdir(here::here("Ecov_study", "catchability", "Results", paste0("OM_", OMsetup$OMname[iom]))) # set up storage
#   saveRDS(OM, file = here::here("Ecov_study", "catchability", "Results", paste0("OM_", OMsetup$OMname[iom]), paste0("OM_", OMsetup$OMname[iom], ".Rds"))) # save OM
#   saveRDS(OMsetup[iom,], file = here::here("Ecov_study", "catchability", "Results", paste0("OM_", OMsetup$OMname[iom]), "OMsettings.Rds")) # save OM settings in a separate file since naming is based on arbitraty order in OMsetup table
# 
# 
#   ##### Iterate over EM misspecifications for this OM
#   for(iem in 1:nrow(EMsetup)){
#     inputEM <- NULL # Reset at start of loop just in case
# 
#     # Seasonal misspecification (q random effect implemented with same misspecification as ecov impact)
#     if(EMsetup[iem, "miss_season"] == "BOTH"){
#       index_list <- list(1) # Both seasons misspecified
#       qRand <- c("ar1", "none")
#     } else if(EMsetup[iem, "miss_season"] == "ONE"){
#       index_list <- list(1,2) # Spring misspecified, fall correct
#       qRand <- c("ar1", "ar1")
#     } else{
#       index_list <- list(2) # Both seasons correctly specified
#       qRand <- c("none", "ar1")
#     }
# 
#     # Catchability model assumptions
#     if(EMsetup[iem, "miss_q"] == "NoEcov"){ # No ecov impact on q or q random effect
#       print("NoEcov")
#       Ecov <- list(label = "NoEcov",
#                    mean = matrix(Ecov_mean_trend, ncol =1), # Mean = 0
#                    logsigma = matrix(log(OMsetup$Ecov_process_obs_sig[iom]), n_years, ncol=1),
#                    year = input$years,
#                    use_obs = matrix(rep(TRUE, n_years), ncol = 1),
#                    lag = 0,
#                    where = "none", # Where/how/indices settings need to change if we do sensitivity runs
#                    how = 1,
#                    indices = index_list, # Ecov impact based on sesonal misspecification setting
#                    process_model = "ar1", # End generic inputs for ecov
#                    process_mean_vals = Ecov_mean_trend, # Mean = 0
#                    process_sig_vals = rep(OMsetup$Ecov_process_sig[iom]),
#                    process_cor_vals = rep(OMsetup$Ecov_process_cor[iom]),
#                    beta_vals = list(append(rep(list(matrix(rep(0, n_ages), nrow = 1)), 2), # No impact on R,M
#                                            rep(list(matrix(rep(OMsetup$Ecov_effect[iom], n_ages), nrow = 1)), 2)))) # impact on q by index (here 2 indices)
# 
#       inputEM <- set_ecov(input = input, ecov = Ecov)
# 
#       # Set up EM name based on EMsetup
#       inputEM$model_name <- paste0("EM_", paste(unfactor(EMsetup[iem, ]), collapse = "_"))
# 
#     } else if(EMsetup[iem, "miss_q"] == "Ecov"){ # Ecov impact on q, correctly specified option when seasonal impact correctly specified
#       print("Ecov")
#       Ecov <- list(label = "Ecov",
#                    mean = matrix(Ecov_mean_trend, ncol =1), # Mean = 0
#                    logsigma = matrix(log(OMsetup$Ecov_process_obs_sig[iom]), n_years, ncol=1),
#                    year = input$years,
#                    use_obs = matrix(rep(TRUE, n_years), ncol = 1),
#                    lag = 0,
#                    where = "q", # Where/how/indices settings need to change if we do sensitivity runs
#                    how = 1,
#                    indices = index_list, # Ecov impact based on sesonal misspecification setting
#                    process_model = "ar1", # End generic inputs for ecov
#                    process_mean_vals = Ecov_mean_trend, # Mean = 0
#                    process_sig_vals = rep(OMsetup$Ecov_process_sig[iom]),
#                    process_cor_vals = rep(OMsetup$Ecov_process_cor[iom]),
#                    beta_vals = list(append(rep(list(matrix(rep(0, n_ages), nrow = 1)), 2), # No impact on R,M
#                                            rep(list(matrix(rep(OMsetup$Ecov_effect[iom], n_ages), nrow = 1)), 2)))) # impact on q by index (here 2 indices)
# 
#       inputEM <- set_ecov(input = input, ecov = Ecov)
# 
#       # Set up EM name based on EMsetup
#       inputEM$model_name <- paste0("EM_", paste(unfactor(EMsetup[iem, ]), collapse = "_"))
# 
#     } else if(EMsetup[iem, "miss_q"] == "qRand"){ # q random effect
#       print("qRand")
#       Ecov <- list(label = "qRand",
#                    mean = matrix(Ecov_mean_trend, ncol =1), # Mean = 0
#                    logsigma = matrix(log(OMsetup$Ecov_process_obs_sig[iom]), n_years, ncol=1),
#                    year = input$years,
#                    use_obs = matrix(rep(TRUE, n_years), ncol = 1),
#                    lag = 0,
#                    where = "none", # Where/how/indices settings need to change if we do sensitivity runs
#                    how = 1,
#                    indices = index_list, # Ecov impact based on sesonal misspecification setting
#                    process_model = "ar1", # End generic inputs for ecov
#                    process_mean_vals = Ecov_mean_trend, # Mean = 0
#                    process_sig_vals = rep(OMsetup$Ecov_process_sig[iom]),
#                    process_cor_vals = rep(OMsetup$Ecov_process_cor[iom]),
#                    beta_vals = list(append(rep(list(matrix(rep(0, n_ages), nrow = 1)), 2), # No impact on R,M
#                                            rep(list(matrix(rep(OMsetup$Ecov_effect[iom], n_ages), nrow = 1)), 2)))) # impact on q by index (here 2 indices)
# 
#       inputEM <- set_ecov(input = input, ecov = Ecov)
# 
#       # Catchability random effects
#       inputEM <- set_q(input = inputEM, catchability = list(re = qRand))
# 
#       # Set up EM name based on EMsetup
#       inputEM$model_name <- paste0("EM_", paste(unfactor(EMsetup[iem, ]), collapse = "_"))
# 
#     } else{ # qRandEcov both q random effect and ecov impact on q
#       print("qRandEcov")
#       Ecov <- list(label = "qRandEcov",
#                    mean = matrix(Ecov_mean_trend, ncol =1), # Mean = 0
#                    logsigma = matrix(log(OMsetup$Ecov_process_obs_sig[iom]), n_years, ncol=1),
#                    year = input$years,
#                    use_obs = matrix(rep(TRUE, n_years), ncol = 1),
#                    lag = 0,
#                    where = "q", # Where/how/indices settings need to change if we do sensitivity runs
#                    how = 1,
#                    indices = index_list, # Ecov impact based on sesonal misspecification setting
#                    process_model = "ar1", # End generic inputs for ecov
#                    process_mean_vals = Ecov_mean_trend, # Mean = 0
#                    process_sig_vals = rep(OMsetup$Ecov_process_sig[iom]),
#                    process_cor_vals = rep(OMsetup$Ecov_process_cor[iom]),
#                    beta_vals = list(append(rep(list(matrix(rep(0, n_ages), nrow = 1)), 2), # No impact on R,M
#                                            rep(list(matrix(rep(OMsetup$Ecov_effect[iom], n_ages), nrow = 1)), 2)))) # impact on q by index (here 2 indices)
# 
#       inputEM <- set_ecov(input = input, ecov = Ecov)
# 
#       # Catchability random effects
#       inputEM <- set_q(input = inputEM, catchability = list(re = qRand))
# 
#       # Set up EM name based on EMsetup
#       inputEM$model_name <- paste0("EM_", paste(unfactor(EMsetup[iem, ]), collapse = "_"))
#     }
# 
#     # Save EM specification
#     mkdir(here::here("Ecov_study", "catchability", "Results", paste0("OM_", iom), paste0("EM_missSeason_", EMsetup[iem, "miss_season"], "_missQ_", EMsetup[iem, "miss_q"]))) # set up storage
# 
#     saveRDS(inputEM, file = here::here("Ecov_study", "catchability", "Results", paste0("OM_", iom), paste0("EM_missSeason_", EMsetup[iem, "miss_season"], "_missQ_", EMsetup[iem, "miss_q"]), "EMinput.Rds")) # save EM input
# 
#   } # End loop over EM specifications
# } # End loop over OMs





#### Run parallel simulations
# Set up parallelization
numCore <- detectCores()

registerDoParallel(numCore-10) # Don't use 2 of the cores

# ##### Simulations with range of ecov beta values
# subsetOM <- OMsetup %>% filter(Ecov_process_sig == 0.1, Ecov_process_cor == 0.5, Ecov_process_obs_sig == 0.1, F_hist == "H-L", ageComp_sig == 0.3, log_index_sig == 0.1, log_catch_sig == 0.1)
# subsetEM <- EMsetup %>% filter(miss_season == "NONE") # No seasonal misspecification

##### Simulations with a range of Ecov observation levels
subsetOM <- OMsetup %>% filter(Ecov_process_sig == 0.1, Ecov_process_cor == 0.5, Ecov_effect == 0.5, F_hist == "H-L", ageComp_sig == 0.3, log_index_sig == 0.1, log_catch_sig == 0.1)
subsetEM <- EMsetup %>% filter(miss_season == "NONE") # No seasonal misspecification


# Run simulation tests
for(iom in 1:nrow(subsetOM)){  # Loop over OMs
  
  # Pull OM
  testOM <- readRDS(here::here("Ecov_study", "catchability", "Results", paste0("OM_", subsetOM[iom, "OMname"]), paste0("OM_", subsetOM[iom, "OMname"], ".Rds")))
  
  # Loop over EMs
  for(iem in 1:nrow(subsetEM)){
    dir <- here::here("Ecov_study", "catchability", "Results", paste0("OM_", subsetOM[iom, "OMname"]), paste0("EM_missSeason_", subsetEM[iem,"miss_season"], "_missQ_", subsetEM[iem,"miss_q"]))
    
    # Pull EM
    testEM <- readRDS(here::here("Ecov_study", "catchability", "Results", paste0("OM_", subsetOM[iom, "OMname"]), paste0("EM_missSeason_", subsetEM[iem, "miss_season"], "_missQ_", subsetEM[iem, "miss_q"]), "EMinput.Rds"))
    
    # Run simulation test in parallelized 2 sim intervals to minimize number of resulting files
    foreach(isim = 1:15) %dopar% { # Run 25 times*2 sims each = 50 sims total in parallel
      simTestWHAM(nsim = 2,
                  OM = testOM,
                  inputEMlist = list(testEM), # Run one EM at a time
                  outdir = dir) # Save in OM directory
    } # End foreach loop over sims
  } # End loop over EMs
} # End loop over OMs




##### Check performance of above OMs with a range of ecov beta parameters
# Find all result files
filenames <- NULL
for(iom in c(7,  19,  31,  43,  55,  67,  79,  91, 103, 115, 127, 139)){
  filenames <- c(filenames, list.files(path = here::here(paste0("Ecov_study/catchability/Results_testEcov/OM_", iom)), pattern = "simWHAM_", recursive = TRUE, full.names = TRUE))
}

# Set storage directory
outdir = here::here("Ecov_study/catchability")

# Post-process results
postprocess_simTestWHAM(filenames = c(filenames), outdir = outdir)

# Plot
# perfMet <- readRDS(here::here("Ecov_study", "catchability", "perfMet_2023-12-14_20-31-52.803322.RDS")) # Only 20 simulations
library(TAF)
mkdir(here::here("Ecov_study", "catchability", "plots_testEcovBeta"))
library(DataExplorer)
plotResults(results = perfMet, convergedONLY = TRUE, outfile = here::here("Ecov_study", "catchability", "plots_testEcovBeta"))


##### Check performance of above OMs with a range of ecov obs errors
# Find all result files
filenames <- NULL
for(iom in c(27, 31, 35)){
  filenames <- c(filenames, list.files(path = here::here(paste0("Ecov_study/catchability/Results_testEcov/OM_", iom)), pattern = "simWHAM_", recursive = TRUE, full.names = TRUE))
}

# Set storage directory
outdir = here::here("Ecov_study/catchability")

# Post-process results
postprocess_simTestWHAM(filenames = c(filenames), outdir = outdir)

# Plot
# perfMet <- readRDS(here::here("Ecov_study", "catchability", "perfMet_2023-12-15_15-47-36.085397.RDS")) # Only 20 simulations
library(TAF)
mkdir(here::here("Ecov_study", "catchability", "plots_testEcovObs"))
library(DataExplorer)
plotResults(results = perfMet, convergedONLY = TRUE, outfile = here::here("Ecov_study", "catchability", "plots_testEcovObs"))



# Alex devel branch SHA: b841eb1d530e94f0c8f4b5de68af891e372ad2ee

