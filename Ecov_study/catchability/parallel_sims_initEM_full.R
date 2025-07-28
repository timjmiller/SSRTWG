# Run full set of OM/EM combinations:
# Include a wider range of environmental effect sizes
# Also fit all EMs to the same OM dataset so that AICs are comparable

# Load packages & source functions used in simulation testing
## Packages
library(tidyverse)
library(wham)
library(TAF)
library(varhandle)
library(doParallel)
library(here)

## OM/EM setup functions
source(here::here("Ecov_study", "catchability", "make_om.R")) # Use revised copy that has option not to include a S-R relationship & sets ecov = NULL so only formatted by prepare_wham_input if provided to make_om
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
Ecov_process_obs_sig <- c(0.0001, 0.1, 0.5) # Add option so obs followed perfectly in OM (i.e. almost no observation error)
Ecov_effect <- c(0, 0.5,  1.5, 3.0) # beta

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


# ## OM/EM setup - environmental covariate updated to increase from 0 to 5 over 40 years
# dir.create(here::here("Ecov_study", "catchability", "Results"))
# 
# n_indices <- 2
# n_ages <- 10
# n_years <- 40
# Ecov_mean_trend <- seq(0, 5, by = 5/39)
# n_fleets <- 1
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
# for(iom in 1:nrow(OMsetup)){ # nrow(OMsetup)){
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
#                years = input$years,
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
#                    years = input$years,
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
#       inputEM <- make_om(Fhist = OMsetup[iom,"F_hist"],
#                        N1_state = "Fmsy", # Default, could also pick "overfished" or "unfished"
#                        selectivity = sel_list,
#                        M = M_list,
#                        catchability = NULL,
#                        NAA_re = NULL,
#                        ecov = Ecov)
# 
#       # Observation error - set L-N SD parameters for catch and index age comps # pulled from project 0 q_om_setup.R
#       inputEM$par$catch_paa_pars[,1] <- log(OMsetup[iom, "ageComp_sig"])
#       inputEM$par$index_paa_pars[,1] <- log(OMsetup[iom, "ageComp_sig"])
#       inputEM$data$agg_catch_sigma[] <- OMsetup[iom,"log_catch_sig"]
#       inputEM$data$agg_index_sigma[] <- OMsetup[iom,"log_index_sig"]
#       inputEM$data$catch_Neff[] <- 1 # Change Neff so scalar doesn't affect L-N SD
#       inputEM$data$index_Neff[] <- 1
# 
# 
#       # Set up EM name based on EMsetup
#       inputEM$model_name <- paste0("EM_", paste(unfactor(EMsetup[iem, ]), collapse = "_"))
# 
#     } else if(EMsetup[iem, "miss_q"] == "Ecov"){ # Ecov impact on q, correctly specified option when seasonal impact correctly specified
#       print("Ecov")
#       Ecov <- list(label = "Ecov",
#                    mean = matrix(Ecov_mean_trend, ncol =1), # Mean = 0
#                    logsigma = matrix(log(OMsetup$Ecov_process_obs_sig[iom]), n_years, ncol=1),
#                    years = input$years,
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
#       inputEM <- make_om(Fhist = OMsetup[iom,"F_hist"],
#                          N1_state = "Fmsy", # Default, could also pick "overfished" or "unfished"
#                          selectivity = sel_list,
#                          M = M_list,
#                          catchability = NULL,
#                          NAA_re = NULL,
#                          ecov = Ecov)
# 
#       # Observation error - set L-N SD parameters for catch and index age comps # pulled from project 0 q_om_setup.R
#       inputEM$par$catch_paa_pars[,1] <- log(OMsetup[iom, "ageComp_sig"])
#       inputEM$par$index_paa_pars[,1] <- log(OMsetup[iom, "ageComp_sig"])
#       inputEM$data$agg_catch_sigma[] <- OMsetup[iom,"log_catch_sig"]
#       inputEM$data$agg_index_sigma[] <- OMsetup[iom,"log_index_sig"]
#       inputEM$data$catch_Neff[] <- 1 # Change Neff so scalar doesn't affect L-N SD
#       inputEM$data$index_Neff[] <- 1
# 
#       # Set up EM name based on EMsetup
#       inputEM$model_name <- paste0("EM_", paste(unfactor(EMsetup[iem, ]), collapse = "_"))
# 
#     } else if(EMsetup[iem, "miss_q"] == "qRand"){ # q random effect
#       print("qRand")
#       Ecov <- list(label = "qRand",
#                    mean = matrix(Ecov_mean_trend, ncol =1), # Mean = 0
#                    logsigma = matrix(log(OMsetup$Ecov_process_obs_sig[iom]), n_years, ncol=1),
#                    years = input$years,
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
#       inputEM <- make_om(Fhist = OMsetup[iom,"F_hist"],
#                          N1_state = "Fmsy", # Default, could also pick "overfished" or "unfished"
#                          selectivity = sel_list,
#                          M = M_list,
#                          catchability = NULL,
#                          NAA_re = NULL,
#                          ecov = Ecov)
# 
#       # Observation error - set L-N SD parameters for catch and index age comps # pulled from project 0 q_om_setup.R
#       inputEM$par$catch_paa_pars[,1] <- log(OMsetup[iom, "ageComp_sig"])
#       inputEM$par$index_paa_pars[,1] <- log(OMsetup[iom, "ageComp_sig"])
#       inputEM$data$agg_catch_sigma[] <- OMsetup[iom,"log_catch_sig"]
#       inputEM$data$agg_index_sigma[] <- OMsetup[iom,"log_index_sig"]
#       inputEM$data$catch_Neff[] <- 1 # Change Neff so scalar doesn't affect L-N SD
#       inputEM$data$index_Neff[] <- 1
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
#                    years = input$years,
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
#       inputEM <- make_om(Fhist = OMsetup[iom,"F_hist"],
#                          N1_state = "Fmsy", # Default, could also pick "overfished" or "unfished"
#                          selectivity = sel_list,
#                          M = M_list,
#                          catchability = NULL,
#                          NAA_re = NULL,
#                          ecov = Ecov)
# 
#       # Observation error - set L-N SD parameters for catch and index age comps # pulled from project 0 q_om_setup.R
#       inputEM$par$catch_paa_pars[,1] <- log(OMsetup[iom, "ageComp_sig"])
#       inputEM$par$index_paa_pars[,1] <- log(OMsetup[iom, "ageComp_sig"])
#       inputEM$data$agg_catch_sigma[] <- OMsetup[iom,"log_catch_sig"]
#       inputEM$data$agg_index_sigma[] <- OMsetup[iom,"log_index_sig"]
#       inputEM$data$catch_Neff[] <- 1 # Change Neff so scalar doesn't affect L-N SD
#       inputEM$data$index_Neff[] <- 1
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

##### Simulations with range of ecov beta values
subsetOM <- OMsetup 
subsetEM <- EMsetup %>% filter(miss_season == "NONE") # No seasonal misspecification



# Run simulation tests
for(iom in 1:nrow(subsetOM)){  # Loop over OMs
  
  # Pull OM
  testOM <- readRDS(here::here("Ecov_study", "catchability", "Results", paste0("OM_", subsetOM[iom, "OMname"]), paste0("OM_", subsetOM[iom, "OMname"], ".Rds")))
  
    omdir <- here::here("Ecov_study", "catchability", "Results", paste0("OM_", subsetOM[iom, "OMname"]))
    
    # Pull EM
    EM_qRand <- readRDS(here::here(omdir, paste0("EM_missSeason_", subsetEM[3, "miss_season"], "_missQ_", subsetEM[3, "miss_q"]), "EMinput.Rds"))
    EM_qRandEcov <- readRDS(here::here(omdir, paste0("EM_missSeason_", subsetEM[4, "miss_season"], "_missQ_", subsetEM[4, "miss_q"]), "EMinput.Rds"))
    EM_Ecov <- readRDS(here::here(omdir, paste0("EM_missSeason_", subsetEM[2, "miss_season"], "_missQ_", subsetEM[2, "miss_q"]), "EMinput.Rds"))
    EM_NoEcov <- readRDS(here::here(omdir, paste0("EM_missSeason_", subsetEM[1, "miss_season"], "_missQ_", subsetEM[1, "miss_q"]), "EMinput.Rds"))
    
    # Run simulation test in parallelized 2 sim intervals to minimize number of resulting files
    foreach(isim = 1:1) %dopar% { # Run 25 times*2 sims each = 50 sims total in parallel
      simTestWHAM(nsim = 2,
                  OM = testOM,
                  inputEMlist = list(EM_qRand, EM_qRandEcov, EM_Ecov, EM_NoEcov), # Run all EMs fit to same OM data
                  outdir = omdir) # Save in OM directory
    } # End foreach loop over sims
} # End loop over OMs




##### Check performance of above OMs with a range of ecov beta parameters
# Find all result files 


filenames <- list.files(path = here::here(paste0("Ecov_study/catchability/Results")), pattern = "simWHAM_", recursive = TRUE, full.names = TRUE)

# "perfmet_2023-12-29_20-13-34.505775.RDS",
# "perfmet_2023-12-29_20-27-36.640791.RDS",
# "perfmet_2023-12-29_20-41-22.442021.RDS",
# "perfmet_2023-12-29_20-55-09.203809.RDS",
# "perfmet_2023-12-29_21-09-09.780349.RDS",
# "perfmet_2023-12-29_21-23-20.431991.RDS",
# "perfmet_2023-12-29_21-37-16.985615.RDS",
# "perfmet_2023-12-29_21-51-28.116545.RDS",
# "perfmet_2023-12-29_22-05-38.41417.RDS",
# "perfmet_2023-12-29_22-19-57.653278.RDS",
# "perfmet_2023-12-29_22-34-20.532968.RDS",
# "perfmet_2023-12-29_22-48-30.489408.RDS",
# "perfmet_2023-12-29_23-02-46.834601.RDS",
# "perfmet_2023-12-29_23-17-07.542494.RDS",
# "perfmet_2023-12-29_23-31-14.37506.RDS",
# "perfmet_2023-12-29_23-45-33.36481.RDS",
# "perfmet_2023-12-29_23-59-44.775598.RDS",
# "perfmet_2023-12-30-00-13-50.700073.RDS",
# "perfmet_2023-01-03_17-30-47.349568.RDS",
# "perfmet_2023-01-03_17-45-25.488529.RDS",
# "perfmet_2023-01-03_18-55-18.280357.RDS",
# "perfmet_2023-01-03_19-09-45.568978.RDS",
# "perfmet_2023-01-03_19-23-45.055622.RDS",
# "perfmet_2023-01-03_19-38-28.88513.RDS",
# "perfmet_2023-01-03_19-53-01.396128.RDS",
# "perfmet_2023-01-03_20-07-36.667071.RDS",
# "perfmet_2023-01-03_20-22-23.60639.RDS",
# "perfmet_2023-01-03_20-37-03.490185.RDS",
# "perfmet_2023-01-04_17-43-37.659487.RDS",
# "perfmet_2023-01-04_17-57-53.103325.RDS",
# "perfmet_2023-01-04_18-12-25.628939.RDS",
# "perfmet_2023-01-04_18-27-07.647875.RDS",
# "perfmet_2023-01-04_18-43-34.467153.RDS",
# "perfmet_2023-01-04_19-01-25.640493.RDS",
# "perfmet_2023-01-04_19-18-03.303149.RDS",
# "perfmet_2023-01-04_19-33-21.512589.RDS"


#filenames <- list.files(path = here::here(paste0("Ecov_study/catchability/Results_initEM_full")), pattern = "simWHAM_", recursive = TRUE, full.names = TRUE)


# Set storage directory
outdir = here::here("Ecov_study/catchability")

# Post-process results
# postprocess_simTestWHAM(filenames = c(filenames[2601:2800]), outdir = outdir)
# postprocess_simTestWHAM(filenames = c(filenames[2801:3000]), outdir = outdir)
# postprocess_simTestWHAM(filenames = c(filenames[3001:3200]), outdir = outdir)
# postprocess_simTestWHAM(filenames = c(filenames[3201:3400]), outdir = outdir)
# postprocess_simTestWHAM(filenames = c(filenames[3401:3600]), outdir = outdir)
# postprocess_simTestWHAM(filenames = c(filenames[3601:3800]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[3801:4000]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[4001:4200]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[4201:4400]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[4401:4600]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[4601:4800]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[4801:5000]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[5001:5200]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[5201:5400]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[5401:5600]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[5601:5800]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[5801:6000]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[6001:6200]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[6202:6400], filenames[6201]), outdir = outdir) # first file didn't converge so formatting is bad
postprocess_simTestWHAM(filenames = c(filenames[6401:6600]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[6601:6800]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[6802:7000],filenames[6801]), outdir = outdir) 
postprocess_simTestWHAM(filenames = c(filenames[7001:7200]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[7201:7400]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[7401:7600]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[7601:7800]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[7801:8000]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[8001:8200]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[8201:8400]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[8402:8600],filenames[8401]), outdir = outdir) 
postprocess_simTestWHAM(filenames = c(filenames[8601:8800]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[8801:9000]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[9001:9200]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[9201:9400]), outdir = outdir)
postprocess_simTestWHAM(filenames = c(filenames[9401:length(filenames)]), outdir = outdir)


combinePerfMet(filenames = c("perfMet_2023-12-29_14-51-37.886169.RDS", #1-100
                             "perfMet_2023-12-29_14-58-44.526053.RDS", #101-200
                             "perfMet_2023-12-29_15-29-39.466686.RDS", #201-300
                             "perfMet_2023-12-29_15-42-49.276661.RDS", #301-400
                             "perfMet_2023-12-29_15-49-50.635176.RDS", #401-500
                             "perfMet_2023-12-29_15-55-16.080999.RDS", #501-600
                             "perfMet_2023-12-29_16-08-09.706646.RDS", #601-700
                             "perfMet_2023-12-29_16-16-22.931958.RDS", #701-800
                             "perfMet_2023-12-29_16-25-37.404685.RDS", #801-900
                             "perfMet_2023-12-29_16-32-04.423158.RDS", #901-1000
                             "perfMet_2023-12-29_16-42-09.57716.RDS", #1001-1100
                             "perfMet_2023-12-29_17-20-47.592459.RDS", #1101-1200
                             "perfMet_2023-12-29_17-28-34.745905.RDS", #1201-1300
                             "perfMet_2023-12-29_17-37-12.177755.RDS", #1301-1400
                             "perfMet_2023-12-29_18-10-45.534984.RDS", #1401-1500
                             "perfMet_2023-12-29_18-19-23.192474.RDS", #1501-1600
                             "perfMet_2023-12-29_18-34-59.439541.RDS", #1601-1800
                             "perfMet_2023-12-29_18-55-49.345083.RDS", #1801-2000
                             "perfMet_2023-12-29_19-23-04.991313.RDS", #2001-2200
                             "perfMet_2023-12-29_19-36-19.603546.RDS", #2201-2400
                             "perfMet_2023-12-29_19-54-20.785649.RDS", #2401-2600
                             ), outdir = outdir)

# # Plot
# perfMet <- readRDS(here::here("Ecov_study", "catchability", "perfMet_2023-12-21_16-12-53.647033.RDS")) 
# library(TAF)
# mkdir(here::here("Ecov_study", "catchability", "plots_testEcovBeta_initEM"))
# library(DataExplorer)
# plotResults(results = perfMet, convergedONLY = TRUE, outfile = here::here("Ecov_study", "catchability", "plots_testEcovBeta_initEM"))
# 
# 
# ##### Check performance of above OMs with a range of ecov obs errors
# # Find all result files
# filenames <- NULL
# for(iom in c(27, 31, 35)){
#   filenames <- c(filenames, list.files(path = here::here(paste0("Ecov_study/catchability/Results_initEM/OM_", iom)), pattern = "simWHAM_", recursive = TRUE, full.names = TRUE))
# }
# 
# # Set storage directory
# outdir = here::here("Ecov_study/catchability")
# 
# # Post-process results
# postprocess_simTestWHAM(filenames = c(filenames), outdir = outdir)
# 
# # Plot
# perfMet <- readRDS(here::here("Ecov_study", "catchability", "perfMet_2023-12-19_19-05-25.884854.RDS")) # Only 20 simulations
# library(TAF)
# mkdir(here::here("Ecov_study", "catchability", "plots_testEcovObs_initEM"))
# library(DataExplorer)
# plotResults(results = perfMet, convergedONLY = TRUE, outfile = here::here("Ecov_study", "catchability", "plots_testEcovObs_initEM"))
