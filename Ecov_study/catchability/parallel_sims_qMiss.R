# Run full set of OM/EM combinations:
# Don't provide EM ecov process arguments, correct spelling mistake for year ecov argument (shouldn't matter since all years have ecov data)
# Fit EMs with both seasons misspecified

# Load packages & source functions used in simulation testing
## Packages
library(tidyverse)
library(wham)
library(TAF)
library(varhandle)
library(doParallel)
library(here)
library(DataExplorer)

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
#                    process_model = "ar1") # End generic inputs for ecov
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
#                    year = input$years,
#                    use_obs = matrix(rep(TRUE, n_years), ncol = 1),
#                    lag = 0,
#                    where = "q", # Where/how/indices settings need to change if we do sensitivity runs
#                    how = 1,
#                    indices = index_list, # Ecov impact based on sesonal misspecification setting
#                    process_model = "ar1") # End generic inputs for ecov
# 
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
#                    year = input$years,
#                    use_obs = matrix(rep(TRUE, n_years), ncol = 1),
#                    lag = 0,
#                    where = "none", # Where/how/indices settings need to change if we do sensitivity runs
#                    how = 1,
#                    indices = index_list, # Ecov impact based on sesonal misspecification setting
#                    process_model = "ar1") # End generic inputs for ecov
# 
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
#                    year = input$years,
#                    use_obs = matrix(rep(TRUE, n_years), ncol = 1),
#                    lag = 0,
#                    where = "q", # Where/how/indices settings need to change if we do sensitivity runs
#                    how = 1,
#                    indices = index_list, # Ecov impact based on sesonal misspecification setting
#                    process_model = "ar1") # End generic inputs for ecov
# 
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

registerDoParallel(25) # Don't use 2 of the cores

##### Simulations with range of ecov beta values
subsetOM <- OMsetup %>% filter(OMname== 249) # Incomplete runs: 247  # 35, 51, 59, 65, 227, 335, 215, 205, 195, 229, 161, 117, 69, 263, 192, 186, 189, 136, 96, 257, 165, 252, 254 check 36/37 # Start 165 and 253 in different runs
# subsetOM <- OMsetup %>% filter(OMname > 247 & OMname < 250)
subsetEM <- EMsetup #%>% filter(miss_season == "NONE" | miss_season == "BOTH") # none or both seasons q misspecified



# Run simulation tests
for(iom in 1:nrow(subsetOM)){  # Loop over OMs
  
  # Pull OM
  testOM <- readRDS(here::here("Ecov_study", "catchability", "Results", paste0("OM_", subsetOM[iom, "OMname"]), paste0("OM_", subsetOM[iom, "OMname"], ".Rds")))
  
  omdir <- here::here("Ecov_study", "catchability", "Results", paste0("OM_", subsetOM[iom, "OMname"]))
  
  # Pull EM
  EM_BOTH_NoEcov <- readRDS(here::here(omdir, paste0("EM_missSeason_", subsetEM[1, "miss_season"], "_missQ_", subsetEM[1, "miss_q"]), "EMinput.Rds"))
  EM_ONE_NoEcov <- readRDS(here::here(omdir, paste0("EM_missSeason_", subsetEM[2, "miss_season"], "_missQ_", subsetEM[2, "miss_q"]), "EMinput.Rds"))
  EM_NONE_NoEcov <- readRDS(here::here(omdir, paste0("EM_missSeason_", subsetEM[3, "miss_season"], "_missQ_", subsetEM[3, "miss_q"]), "EMinput.Rds"))
  EM_BOTH_Ecov <- readRDS(here::here(omdir, paste0("EM_missSeason_", subsetEM[4, "miss_season"], "_missQ_", subsetEM[4, "miss_q"]), "EMinput.Rds"))
  EM_ONE_Ecov <- readRDS(here::here(omdir, paste0("EM_missSeason_", subsetEM[5, "miss_season"], "_missQ_", subsetEM[5, "miss_q"]), "EMinput.Rds"))
  EM_NONE_Ecov <- readRDS(here::here(omdir, paste0("EM_missSeason_", subsetEM[6, "miss_season"], "_missQ_", subsetEM[6, "miss_q"]), "EMinput.Rds"))
  EM_BOTH_qRand <- readRDS(here::here(omdir, paste0("EM_missSeason_", subsetEM[7, "miss_season"], "_missQ_", subsetEM[7, "miss_q"]), "EMinput.Rds"))
  EM_ONE_qRand <- readRDS(here::here(omdir, paste0("EM_missSeason_", subsetEM[8, "miss_season"], "_missQ_", subsetEM[8, "miss_q"]), "EMinput.Rds"))
  EM_NONE_qRand <- readRDS(here::here(omdir, paste0("EM_missSeason_", subsetEM[9, "miss_season"], "_missQ_", subsetEM[9, "miss_q"]), "EMinput.Rds"))
  EM_BOTH_qRandEcov <- readRDS(here::here(omdir, paste0("EM_missSeason_", subsetEM[10, "miss_season"], "_missQ_", subsetEM[10, "miss_q"]), "EMinput.Rds"))
  EM_ONE_qRandEcov <- readRDS(here::here(omdir, paste0("EM_missSeason_", subsetEM[11, "miss_season"], "_missQ_", subsetEM[11, "miss_q"]), "EMinput.Rds"))
  EM_NONE_qRandEcov <- readRDS(here::here(omdir, paste0("EM_missSeason_", subsetEM[12, "miss_season"], "_missQ_", subsetEM[12, "miss_q"]), "EMinput.Rds"))
  
  # Run simulation test in parallelized 2 sim intervals to minimize number of resulting files
  foreach(isim = 1:10) %dopar% { # Run 25 times*2 sims each = 50 sims total in parallel
    simTestWHAM(nsim = 2,
                OM = testOM,
                inputEMlist = list(EM_BOTH_NoEcov, EM_ONE_NoEcov, EM_NONE_NoEcov,
                                   EM_BOTH_Ecov, EM_ONE_Ecov, EM_NONE_Ecov,
                                   EM_BOTH_qRand, EM_ONE_qRand, EM_NONE_qRand,
                                   EM_BOTH_qRandEcov, EM_ONE_qRandEcov, EM_NONE_qRandEcov), # Run all EMs fit to same OM data
                outdir = omdir) # Save in OM directory
  } # End foreach loop over sims
} # End loop over OMs




##### Post-process results
# Find all result files 
filenames <- list.files(path = here::here(paste0("Ecov_study/catchability/Results")), pattern = "simWHAM_", recursive = TRUE, full.names = TRUE)

# Set storage directory
outdir = here::here("Ecov_study/catchability/Results/parallel_process")

# Post-process results
registerDoParallel(40)
# for(idir in 1:length(seq(1150, 4100, by = 50))){ # Store results in sequence named folders to allow debugging if one parallel run doesn't complete
#   dirname <- paste("Ecov_study", "catchability", "Results", "parallel_process", paste0("process_", seq(1150, 4100, by = 50)), sep="/")[idir]
#   dir.create(here::here(dirname))
# }
# for(idir in 1:length(c(seq(5800, 6800, by = 50), head(seq(7850, length(filenames), by = 50),-1)))){
#   dirname <- paste("Ecov_study", "catchability", "Results", "parallel_process", paste0("process_", c(seq(5800, 6800, by = 50), head(seq(7850, length(filenames), by = 50),-1))), sep="/")[idir]
#   dir.create(here::here(dirname))
# }
# for(idir in 1:length(c(seq(1, 1100, by = 50), seq(4150, 5750, by = 50), seq(6850, 7800, by = 50)))){
#   dirname <- paste("Ecov_study", "catchability", "Results", "parallel_process", paste0("process_", c(seq(1, 1100, by = 50), seq(4150, 5750, by = 50), seq(6850, 7800, by = 50))), sep="/")[idir]
#   dir.create(here::here(dirname))
# }
# foreach(ifile = seq(1150, 4100, by = 50)) %dopar% { 
# foreach(ifile = c(seq(5800, 6800, by = 50), head(seq(7850, length(filenames), by = 50),-1))) %dopar% {
#foreach(ifile = c(seq(1, 1100, by = 50), seq(4150, 5750, by = 50), seq(6850, 7800, by = 50))) %dopar% {

# for(idir in 1:length(head(seq(1, length(filenames), by=50), -1))){
#   dirname <- paste("Ecov_study", "catchability", "Results", "parallel_process", paste0("process_", head(seq(1, length(filenames), by=50), -1)), sep="/")[idir]
#   dir.create(here::here(dirname))
# }
#foreach(ifile = head(seq(1, length(filenames), by=50), -1)) %dopar% {
foreach(ifile = c(9201, 9251, 9301, 9351, 9401, 9451, 9501, 9551)) %dopar% {
  print(paste0("ifile ", ifile))
  postprocess_simTestWHAM(filenames =  c(filenames[ifile:(ifile+49)]), outdir = paste0(outdir,"/process_", ifile))
}
postprocess_simTestWHAM(filenames =  c(filenames[tail(seq(1, length(filenames), by = 50), n=1):length(filenames)]), outdir = paste0(outdir,"/process_", 9601)) # Last files may be less than 50 interval so processed independently
postprocess_simTestWHAM(filenames = c(filenames[9152:9200], filenames[9151]), outdir = paste0(outdir,"/process_", 9151))
postprocess_simTestWHAM(filenames = c(filenames[9552:9600], filenames[9551]), outdir = paste0(outdir,"/process_", 9551))

parallel_processed <-  list.files(path = here::here(paste0("Ecov_study/catchability/Results/parallel_process")), pattern = "perfMet_", recursive = TRUE, full.names = TRUE)

# Check that all blocks of files have been processed
check <- parallel_processed %>% strsplit(., "/", fixed = TRUE) %>% unlist() %>% matrix(ncol = 11, byrow = T) %>% as.data.frame()
check <- check[,10] %>% strsplit("process_", fixed = TRUE) %>% unlist() %>% matrix(ncol = 2, byrow = TRUE)
check <- check[,2] %>% as.numeric() 
seq(1, length(filenames), by=50)[-which(check %in% seq(1, length(filenames), by=50))] # if numeric(0) then finished processing

# Aggregate processed files
combinePerfMet(filenames = parallel_processed, outdir = outdir)

##### Aggregate processed files by seasonal miss-specification (required because full results file too large)
batchCombine <- function(filenames=NULL, manualLast = NULL){
  aggPerfMet <- NULL
  
  for(ifile in 1:length(filenames)){
    print(ifile)
    
    # Read in RDS file
    perfMet <- readRDS(file = here::here(filenames[ifile]))
    
    if(is.null(manualLast) == TRUE){
      if(ifile == 1){ # If this is the first file, append to the aggPerfMet storage
        aggPerfMet <- rbind(aggPerfMet, perfMet)
      } else{
        # ID number of sims to append
        nsim <- perfMet$sim %>% unique() %>% length()
        # ID last sim number in aggPerfMet 
          lastsim <- aggPerfMet$sim %>% unique() %>% as.numeric() %>% max()
        # New sim numbers
        newSim <- lastsim + 1:nsim
        # Replace old sim numbers with new sim numbers
        for(isim in 1:nsim){
          perfMet[which(perfMet$sim == isim), "sim"] <- newSim[isim]
        }
        # Append newly renumbered simulations to aggPerfMet storage
        aggPerfMet <- rbind(aggPerfMet, perfMet)
      }
    } else{ 
      if(ifile == 1){ # Use manualLast number to order sims within the batch combined files
        # ID number of sims to append
        nsim <- perfMet$sim %>% unique() %>% length()
        # ID last sim number in aggPerfMet 
        lastsim <- manualLast
        
        # New sim numbers
        newSim <- lastsim + 1:nsim
        # Replace old sim numbers with new sim numbers
        for(isim in 1:nsim){
          perfMet[which(perfMet$sim == isim), "sim"] <- newSim[isim]
        }
        # Append newly renumbered simulations to aggPerfMet storage
        aggPerfMet <- rbind(aggPerfMet, perfMet)
      } else{
        # ID number of sims to append
        nsim <- perfMet$sim %>% unique() %>% length()
        # ID last sim number in aggPerfMet 
        lastsim <- aggPerfMet$sim %>% unique() %>% as.numeric() %>% max()
        # New sim numbers
        newSim <- lastsim + 1:nsim
        # Replace old sim numbers with new sim numbers
        for(isim in 1:nsim){
          perfMet[which(perfMet$sim == isim), "sim"] <- newSim[isim]
        }
        # Append newly renumbered simulations to aggPerfMet storage
        aggPerfMet <- rbind(aggPerfMet, perfMet)
      }
        
    } # End if statement for manualLast
  } # End loop over files
  
  # Save results by seasonal misspecification 
  timeStamp <- Sys.time() %>% gsub(" ", "_",.) %>% gsub(":", "-", .) # Change millisecond half of Sys.time() output to avoid having spaces/weird characters in filenames
  
  aggPerfMet %>% filter(EM_miss_season == "NONE") %>% filter(F_hist == "H-L") %>% saveRDS(., file =  here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE_HL", paste0("perfMet_missSeason_NONE_HL_",timeStamp, ".Rds")))
  aggPerfMet %>% filter(EM_miss_season == "NONE") %>% filter(F_hist == "Fmsy") %>% saveRDS(., file =  here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE_Fmsy", paste0("perfMet_missSeason_NONE_Fmsy_",timeStamp, ".Rds")))
  aggPerfMet %>% filter(EM_miss_season == "ONE") %>% filter(F_hist == "H-L") %>% saveRDS(., file =  here::here("Ecov_study", "catchability", "Results", "plots_missSeason_ONE_HL", paste0("perfMet_missSeason_ONE_HL_",timeStamp, ".Rds")))
  aggPerfMet %>% filter(EM_miss_season == "ONE") %>% filter(F_hist == "Fmsy") %>% saveRDS(., file =  here::here("Ecov_study", "catchability", "Results", "plots_missSeason_ONE_Fmsy", paste0("perfMet_missSeason_ONE_Fmsy_",timeStamp, ".Rds")))
  aggPerfMet %>% filter(EM_miss_season == "BOTH") %>% filter(F_hist == "H-L") %>% saveRDS(., file =  here::here("Ecov_study", "catchability", "Results", "plots_missSeason_BOTH_HL", paste0("perfMet_missSeason_BOTH_HL_",timeStamp, ".Rds")))
  aggPerfMet %>% filter(EM_miss_season == "BOTH") %>% filter(F_hist == "Fmsy") %>% saveRDS(., file =  here::here("Ecov_study", "catchability", "Results", "plots_missSeason_BOTH_Fmsy", paste0("perfMet_missSeason_BOTH_Fmsy_",timeStamp, ".Rds")))
}

batchCombine(filenames = parallel_processed[1:10], manualLast = NULL)
readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE_Fmsy", "perfMet_missSeason_NONE_Fmsy_2024-01-18_14-48-45.404633.RDS"))$sim %>% unique() %>% as.numeric() %>% max() # NONE option is last in each file so should be max sim from the file
batchCombine(filenames = parallel_processed[11:20], manualLast = 12000)
readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE_Fmsy", "perfMet_missSeason_NONE_Fmsy_2024-01-18_14-55-07.363647.RDS"))$sim %>% unique() %>% as.numeric() %>% max() # NONE option is last in each file so should be max sim from the file
batchCombine(filenames = parallel_processed[21:30], manualLast = 24000)
# readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "perfMet_missSeason_NONE_2024-01-18_05-24-04.481042.RDS"))$sim %>% unique() %>% as.numeric() %>% max()
batchCombine(filenames = parallel_processed[31:40], manualLast = 36000)
# readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "perfMet_missSeason_NONE_2024-01-18_05-29-31.42681.RDS"))$sim %>% unique() %>% as.numeric() %>% max()
batchCombine(filenames = parallel_processed[41:50], manualLast = 48000)
# readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "perfMet_missSeason_NONE_2024-01-18_05-35-51.096573.RDS"))$sim %>% unique() %>% as.numeric() %>% max()
batchCombine(filenames = parallel_processed[51:60], manualLast = 60000)
# readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "perfMet_missSeason_NONE_2024-01-18_05-40-11.302437.RDS"))$sim %>% unique() %>% as.numeric() %>% max()
batchCombine(filenames = parallel_processed[61:70], manualLast = 72000)
# readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "perfMet_missSeason_NONE_2024-01-18_05-43-54.779406.RDS"))$sim %>% unique() %>% as.numeric() %>% max()
batchCombine(filenames = parallel_processed[71:80], manualLast = 84000)
# readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "perfMet_missSeason_NONE_2024-01-18_05-45-13.199043.RDS"))$sim %>% unique() %>% as.numeric() %>% max()
batchCombine(filenames = parallel_processed[81:90], manualLast = 96000)
# readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "perfMet_missSeason_NONE_2024-01-18_05-48-00.055611.RDS"))$sim %>% unique() %>% as.numeric() %>% max()
batchCombine(filenames = parallel_processed[91:100], manualLast = 108000)
# readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "perfMet_missSeason_NONE_2024-01-18_05-48-47.913159.RDS"))$sim %>% unique() %>% as.numeric() %>% max()
batchCombine(filenames = parallel_processed[101:110], manualLast = 120000)
# readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "perfMet_missSeason_NONE_2024-01-18_05-52-07.416731.RDS"))$sim %>% unique() %>% as.numeric() %>% max()
batchCombine(filenames = parallel_processed[111:120], manualLast = 132000)
# readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "perfMet_missSeason_NONE_2024-01-18_05-52-41.688762.RDS"))$sim %>% unique() %>% as.numeric() %>% max()
batchCombine(filenames = parallel_processed[121:130], manualLast = 144000)
# readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "perfMet_missSeason_NONE_2024-01-18_05-56-28.449862.RDS"))$sim %>% unique() %>% as.numeric() %>% max()
batchCombine(filenames = parallel_processed[131:140], manualLast = 156000)
# readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "perfMet_missSeason_NONE_2024-01-18_05-56-52.446157.RDS"))$sim %>% unique() %>% as.numeric() %>% max()
batchCombine(filenames = parallel_processed[141:150], manualLast = 168000)
# readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "perfMet_missSeason_NONE_2024-01-18_06-00-14.716907.RDS"))$sim %>% unique() %>% as.numeric() %>% max()
batchCombine(filenames = parallel_processed[151:160], manualLast = 180000)
# readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "perfMet_missSeason_NONE_2024-01-18_06-01-05.860818.RDS"))$sim %>% unique() %>% as.numeric() %>% max()
batchCombine(filenames = parallel_processed[161:170], manualLast = 192000)
# readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "perfMet_missSeason_NONE_2024-01-18_06-04-15.81736.RDS"))$sim %>% unique() %>% as.numeric() %>% max()
batchCombine(filenames = parallel_processed[171:180], manualLast = 204000)
# readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "perfMet_missSeason_NONE_2024-01-18_06-04-36.916695.RDS")) %>% select(sim) %>% unique() %>% as.numeric() %>% max()
batchCombine(filenames = parallel_processed[181:length(parallel_processed)], manualLast = 216000)
# readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "perfMet_missSeason_NONE_2024-01-18_06-09-16.129359.RDS")) %>% select(sim) %>% unique() %>% as.numeric() %>% max()

# # Append batches for no misspecification
# missSeason_NONE <- list.files(path = here::here(paste0("Ecov_study/catchability/Results/plots_missSeason_NONE")), pattern = "perfMet_missSeason_NONE", recursive = FALSE, full.names = TRUE)
# for(ifile in 1:8){
#   print(ifile)
#   if(ifile == 1){
#     file1 <- readRDS(missSeason_NONE[1])
#     file2 <- readRDS(missSeason_NONE[2])
#     rbind(file1, file2) %>% saveRDS(., file =  here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "aggPerfMet_missSeason_NONE.Rds"))
#   } else if(ifile > 2){
#     aggFile <- readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "aggPerfMet_missSeason_NONE.Rds"))
#     tempfile <- readRDS(missSeason_NONE[ifile])
#     # append to file
#     rbind(aggFile, tempfile) %>% saveRDS(., file =  here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "aggPerfMet_missSeason_NONE.Rds"))
#   }
# }
# #!!! start here
# for(ifile in 9:length(missSeason_NONE)){
#   print(ifile)
#   if(ifile == 9){
#     file1 <- readRDS(missSeason_NONE[9])
#     file2 <- readRDS(missSeason_NONE[10])
#     rbind(file1, file2) %>% saveRDS(., file =  here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "aggPerfMet_missSeason_NONE2.Rds"))
#   } else if(ifile > 10){
#     aggFile <- readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "aggPerfMet_missSeason_NONE2.Rds"))
#     tempfile <- readRDS(missSeason_NONE[ifile])
#     # append to file
#     rbind(aggFile, tempfile) %>% saveRDS(., file =  here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "aggPerfMet_missSeason_NONE2.Rds"))
#   }
# }
# agg1 <- readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "aggPerfMet_missSeason_NONE.Rds")) 
# agg2 <- readRDS(here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "aggPerfMet_missSeason_NONE2.Rds"))
# rbind(agg1, agg2) %>% saveRDS(., file =  here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE", "aggPerfMet_missSeason_NONE.Rds"))
# 
# # Function to expedite combining files (I hope)
# rbindCombine <- function(filenames = NULL, outfile = NULL){
#   for(ifile in 1:length(filenames)){
#     print(ifile)
#     if(ifile == 1){
#       file1 <- readRDS(filenames[1])
#       file2 <- readRDS(filenames[2])
#       rbind(file1, file2) %>% saveRDS(., file =  outfile)
#     } else if(ifile > 2){
#       aggFile <- readRDS(outfile)
#       tempfile <- readRDS(missSeason_ONE[ifile])
#       # append to file
#       rbind(aggFile, tempfile) %>% saveRDS(., file =  outfile)
#     }
#   }
# }

testrbindCombine <- function(filenames = NULL, outfile = NULL){
      file1 <- readRDS(filenames[1])
      file2 <- readRDS(filenames[2])
      file3 <- readRDS(filenames[3])
      file4 <- readRDS(filenames[4])
      if(length(filenames) == 5){
        file5 <- readRDS(filenames[5])
        rbind(file1, file2, file3, file4, file5) %>% saveRDS(., file =  outfile)
      } else{
        rbind(file1, file2, file3, file4) %>% saveRDS(., file =  outfile)
      }
}

# Append batches for one season misspecified Fhist H-L
missSeason_ONE <- list.files(path = here::here(paste0("Ecov_study/catchability/Results/plots_missSeason_ONE_HL")), pattern = "perfMet_missSeason_ONE", recursive = FALSE, full.names = TRUE)
testrbindCombine(filenames = missSeason_ONE[1:5], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_ONE_HL", "aggPerfMet_missSeason_ONE1-5.Rds"))
testrbindCombine(filenames = missSeason_ONE[6:10], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_ONE_HL", "aggPerfMet_missSeason_ONE6-10.Rds"))
testrbindCombine(filenames = missSeason_ONE[11:15], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_ONE_HL", "aggPerfMet_missSeason_ONE11-15.Rds"))
testrbindCombine(filenames = missSeason_ONE[16:length(missSeason_ONE)], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_ONE_HL", "aggPerfMet_missSeason_ONE16-end.Rds"))
final_ONE <- missSeason_ONE <- list.files(path = here::here(paste0("Ecov_study/catchability/Results/plots_missSeason_ONE_HL")), pattern = "aggPerfMet_missSeason_ONE", recursive = FALSE, full.names = TRUE)
testrbindCombine(filenames = final_ONE, outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_ONE_HL", "aggPerfMet_missSeason_ONE_HL.Rds"))
# Check that all OM/EM pairs have been processed
path = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_ONE_HL")
results <- readRDS(paste0(path,"/aggPerfMet_missSeason_ONE_HL.Rds"))
check <- cbind(rep(seq(1,nrow(OMsetup)),each=4), rep(paste0("EM_ONE_", miss_q), nrow(OMsetup))) %>% as.data.frame()
colnames(check) <- c("OMshortName",	"EMshortName") # Set up full combinations
check2 <- results %>% count(OMshortName, EMshortName) %>% mutate(nsim = n/40) # Check if completed
check3 <- full_join(check, check2) 
check3$OMshortName <- check3$OMshortName %>% as.numeric()
colnames(OMsetup)[1] <- "OMshortName"
full_join(check3, OMsetup, by = "OMshortName") %>% write.csv(., paste0(path,"/simSummary.csv")) # /40 years so nsim = number of full simulations for each OM/EM
# Plot
plotResults(results = results, convergedONLY = TRUE, outfile = path)

# Append batches for one season misspecified Fhist Fmsy
missSeason_ONE <- list.files(path = here::here(paste0("Ecov_study/catchability/Results/plots_missSeason_ONE_Fmsy")), pattern = "perfMet_missSeason_ONE", recursive = FALSE, full.names = TRUE)
testrbindCombine(filenames = missSeason_ONE[1:5], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_ONE_Fmsy", "aggPerfMet_missSeason_ONE1-5.Rds"))
testrbindCombine(filenames = missSeason_ONE[6:10], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_ONE_Fmsy", "aggPerfMet_missSeason_ONE6-10.Rds"))
testrbindCombine(filenames = missSeason_ONE[11:15], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_ONE_Fmsy", "aggPerfMet_missSeason_ONE11-15.Rds"))
testrbindCombine(filenames = missSeason_ONE[16:length(missSeason_ONE)], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_ONE_Fmsy", "aggPerfMet_missSeason_ONE16-end.Rds"))
final_ONE <- missSeason_ONE <- list.files(path = here::here(paste0("Ecov_study/catchability/Results/plots_missSeason_ONE_Fmsy")), pattern = "aggPerfMet_missSeason_ONE", recursive = FALSE, full.names = TRUE)
testrbindCombine(filenames = final_ONE, outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_ONE_Fmsy", "aggPerfMet_missSeason_ONE_Fmsy.Rds"))
# Check that all OM/EM pairs have been processed
path = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_ONE_Fmsy")
results <- readRDS(paste0(path,"/aggPerfMet_missSeason_ONE_Fmsy.Rds"))
check <- cbind(rep(seq(1,nrow(OMsetup)),each=4), rep(paste0("EM_ONE_", miss_q), nrow(OMsetup))) %>% as.data.frame()
colnames(check) <- c("OMshortName",	"EMshortName") # Set up full combinations
check2 <- results %>% count(OMshortName, EMshortName) %>% mutate(nsim = n/40) # Check if completed
check3 <- full_join(check, check2) 
check3$OMshortName <- check3$OMshortName %>% as.numeric()
colnames(OMsetup)[1] <- "OMshortName"
full_join(check3, OMsetup, by = "OMshortName") %>% write.csv(., paste0(path,"/simSummary.csv")) # /40 years so nsim = number of full simulations for each OM/EM
# Plot
plotResults(results = results, convergedONLY = TRUE, outfile = path)


# Append batches for both seasons misspecified Fhist H-L
missSeason_BOTH <- list.files(path = here::here(paste0("Ecov_study/catchability/Results/plots_missSeason_BOTH_HL")), pattern = "perfMet_missSeason_BOTH", recursive = FALSE, full.names = TRUE)
testrbindCombine(filenames = missSeason_BOTH[1:5], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_BOTH_HL", "aggPerfMet_missSeason_BOTH1-5.Rds"))
testrbindCombine(filenames = missSeason_BOTH[6:10], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_BOTH_HL", "aggPerfMet_missSeason_BOTH6-10.Rds"))
testrbindCombine(filenames = missSeason_BOTH[11:15], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_BOTH_HL", "aggPerfMet_missSeason_BOTH11-15.Rds"))
testrbindCombine(filenames = missSeason_BOTH[16:length(missSeason_BOTH)], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_BOTH_HL", "aggPerfMet_missSeason_BOTH16-end.Rds"))
final_BOTH <- missSeason_BOTH <- list.files(path = here::here(paste0("Ecov_study/catchability/Results/plots_missSeason_BOTH_HL")), pattern = "aggPerfMet_missSeason_BOTH", recursive = FALSE, full.names = TRUE)
testrbindCombine(filenames = final_BOTH, outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_BOTH_HL", "aggPerfMet_missSeason_BOTH_HL.Rds"))
# Check that all OM/EM pairs have been processed
path = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_BOTH_HL")
results <- readRDS(paste0(path,"/aggPerfMet_missSeason_BOTH_HL.Rds"))
check <- cbind(rep(seq(1,nrow(OMsetup)),each=4), rep(paste0("EM_BOTH_", miss_q), nrow(OMsetup))) %>% as.data.frame()
colnames(check) <- c("OMshortName",	"EMshortName") # Set up full combinations
check2 <- results %>% count(OMshortName, EMshortName) %>% mutate(nsim = n/40) # Check if completed
check3 <- full_join(check, check2) 
check3$OMshortName <- check3$OMshortName %>% as.numeric()
colnames(OMsetup)[1] <- "OMshortName"
full_join(check3, OMsetup, by = "OMshortName") %>% write.csv(., paste0(path,"/simSummary.csv")) # /40 years so nsim = number of full simulations for each OM/EM
# Plot
plotResults(results = results, convergedONLY = TRUE, outfile = path)

# Append batches for both seasons misspecified Fhist Fmsy
missSeason_BOTH <- list.files(path = here::here(paste0("Ecov_study/catchability/Results/plots_missSeason_BOTH_Fmsy")), pattern = "perfMet_missSeason_BOTH", recursive = FALSE, full.names = TRUE)
testrbindCombine(filenames = missSeason_BOTH[1:5], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_BOTH_Fmsy", "aggPerfMet_missSeason_BOTH1-5.Rds"))
testrbindCombine(filenames = missSeason_BOTH[6:10], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_BOTH_Fmsy", "aggPerfMet_missSeason_BOTH6-10.Rds"))
testrbindCombine(filenames = missSeason_BOTH[11:15], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_BOTH_Fmsy", "aggPerfMet_missSeason_BOTH11-15.Rds"))
testrbindCombine(filenames = missSeason_BOTH[16:length(missSeason_BOTH)], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_BOTH_Fmsy", "aggPerfMet_missSeason_BOTH16-end.Rds"))
final_BOTH <- missSeason_BOTH <- list.files(path = here::here(paste0("Ecov_study/catchability/Results/plots_missSeason_BOTH_Fmsy")), pattern = "aggPerfMet_missSeason_BOTH", recursive = FALSE, full.names = TRUE)
testrbindCombine(filenames = final_BOTH, outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_BOTH_Fmsy", "aggPerfMet_missSeason_BOTH_Fmsy.Rds"))
# Check that all OM/EM pairs have been processed
path = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_BOTH_Fmsy")
results <- readRDS(paste0(path,"/aggPerfMet_missSeason_BOTH_Fmsy.Rds"))
check <- cbind(rep(seq(1,nrow(OMsetup)),each=4), rep(paste0("EM_BOTH_", miss_q), nrow(OMsetup))) %>% as.data.frame()
colnames(check) <- c("OMshortName",	"EMshortName") # Set up full combinations
check2 <- results %>% count(OMshortName, EMshortName) %>% mutate(nsim = n/40) # Check if completed
check3 <- full_join(check, check2) 
check3$OMshortName <- check3$OMshortName %>% as.numeric()
colnames(OMsetup)[1] <- "OMshortName"
full_join(check3, OMsetup, by = "OMshortName") %>% write.csv(., paste0(path,"/simSummary.csv")) # /40 years so nsim = number of full simulations for each OM/EM
# Plot
plotResults(results = results, convergedONLY = TRUE, outfile = path)

# Append batches for NONE seasons misspecified Fhist H-L
missSeason_NONE <- list.files(path = here::here(paste0("Ecov_study/catchability/Results/plots_missSeason_NONE_HL")), pattern = "perfMet_missSeason_NONE", recursive = FALSE, full.names = TRUE)
testrbindCombine(filenames = missSeason_NONE[1:5], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE_HL", "aggPerfMet_missSeason_NONE1-5.Rds"))
testrbindCombine(filenames = missSeason_NONE[6:10], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE_HL", "aggPerfMet_missSeason_NONE6-10.Rds"))
testrbindCombine(filenames = missSeason_NONE[11:15], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE_HL", "aggPerfMet_missSeason_NONE11-15.Rds"))
testrbindCombine(filenames = missSeason_NONE[16:length(missSeason_NONE)], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE_HL", "aggPerfMet_missSeason_NONE16-end.Rds"))
final_NONE <- missSeason_NONE <- list.files(path = here::here(paste0("Ecov_study/catchability/Results/plots_missSeason_NONE_HL")), pattern = "aggPerfMet_missSeason_NONE", recursive = FALSE, full.names = TRUE)
testrbindCombine(filenames = final_NONE, outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE_HL", "aggPerfMet_missSeason_NONE_HL.Rds"))
# Check that all OM/EM pairs have been processed
path = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE_HL")
results <- readRDS(paste0(path,"/aggPerfMet_missSeason_NONE_HL.Rds"))
check <- cbind(rep(seq(1,nrow(OMsetup)),each=4), rep(paste0("EM_NONE_", miss_q), nrow(OMsetup))) %>% as.data.frame()
colnames(check) <- c("OMshortName",	"EMshortName") # Set up full combinations
check2 <- results %>% count(OMshortName, EMshortName) %>% mutate(nsim = n/40) # Check if completed
check3 <- full_join(check, check2) 
check3$OMshortName <- check3$OMshortName %>% as.numeric()
colnames(OMsetup)[1] <- "OMshortName"
full_join(check3, OMsetup, by = "OMshortName") %>% write.csv(., paste0(path,"/simSummary.csv")) # /40 years so nsim = number of full simulations for each OM/EM
# Plot
plotResults(results = results, convergedONLY = TRUE, outfile = path)

# Append batches for NONE seasons misspecified Fhist Fmsy
missSeason_NONE <- list.files(path = here::here(paste0("Ecov_study/catchability/Results/plots_missSeason_NONE_Fmsy")), pattern = "perfMet_missSeason_NONE", recursive = FALSE, full.names = TRUE)
testrbindCombine(filenames = missSeason_NONE[1:5], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE_Fmsy", "aggPerfMet_missSeason_NONE1-5.Rds"))
testrbindCombine(filenames = missSeason_NONE[6:10], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE_Fmsy", "aggPerfMet_missSeason_NONE6-10.Rds"))
testrbindCombine(filenames = missSeason_NONE[11:15], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE_Fmsy", "aggPerfMet_missSeason_NONE11-15.Rds"))
testrbindCombine(filenames = missSeason_NONE[16:length(missSeason_NONE)], outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE_Fmsy", "aggPerfMet_missSeason_NONE16-end.Rds"))
final_NONE <- missSeason_NONE <- list.files(path = here::here(paste0("Ecov_study/catchability/Results/plots_missSeason_NONE_Fmsy")), pattern = "aggPerfMet_missSeason_NONE", recursive = FALSE, full.names = TRUE)
testrbindCombine(filenames = final_NONE, outfile = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE_Fmsy", "aggPerfMet_missSeason_NONE_Fmsy.Rds"))
# Check that all OM/EM pairs have been processed
path = here::here("Ecov_study", "catchability", "Results", "plots_missSeason_NONE_Fmsy")
results <- readRDS(paste0(path,"/aggPerfMet_missSeason_NONE_Fmsy.Rds"))
check <- cbind(rep(seq(1,nrow(OMsetup)),each=4), rep(paste0("EM_NONE_", miss_q), nrow(OMsetup))) %>% as.data.frame()
colnames(check) <- c("OMshortName",	"EMshortName") # Set up full combinations
check2 <- results %>% count(OMshortName, EMshortName) %>% mutate(nsim = n/40) # Check if completed
check3 <- full_join(check, check2) 
check3$OMshortName <- check3$OMshortName %>% as.numeric()
colnames(OMsetup)[1] <- "OMshortName"
full_join(check3, OMsetup, by = "OMshortName") %>% write.csv(., paste0(path,"/simSummary.csv")) # /40 years so nsim = number of full simulations for each OM/EM
# Plot
plotResults(results = results, convergedONLY = TRUE, outfile = path)



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



