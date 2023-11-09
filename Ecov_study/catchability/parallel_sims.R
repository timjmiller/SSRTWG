# Try parallelizing code

# Load packages & source functions used in simulation testing
## Packages
library(tidyverse)
library(wham)
library(TAF)
library(varhandle)
library(doParallel)
library(here)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Don't forget to create a Results directory before running!

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
Ecov_process_obs_sig <- c(0.1, 0.5)
Ecov_effect <- c(0, 0.25, 0.5) # beta

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


### Could have OM/EM setup here if needed to run remotely but since I am testing this locally all I need is the 




#### Run parallel simulations
# Set up parallelization
numCore <- detectCores()

registerDoParallel(numCore-20) # Don't use 2 of the cores

# Alex - to run your Fmsy simulations adjusts the number of simulations in line 76 and change F_hist in line 79 to = "Fmsy"
# Set number of simulations to run
nsim <- 50

# ##### Run settings for separate Fmsy and H-L simulations with environmental effect implemented
# # Subset of OMs to run
# # subsetOM <- OMsetup %>% filter(F_hist == "H-L") %>% filter(Ecov_effect != 0)# %>% filter(OMname > 153) # Only run 3 OMs for 2 simulations as a test
# subsetOM <- OMsetup %>% filter(F_hist == "Fmsy") %>% filter(Ecov_effect != 0) %>% filter(OMname > 176) # Only run 3 OMs for 2 simulations as a test
# 
# # Subset of EMs to run
# subsetEM <- EMsetup %>% as.data.frame() %>% filter(miss_season == "NONE")

##### Run settings for OMs with no environmental effects - smaller subset of EMs run for these OMs
# Subset of OMs to run
subsetOM <- OMsetup %>% filter(Ecov_effect == 0) #%>% filter(OMname > 124)
# EMs 
subsetEM <- EMsetup %>% filter(miss_season == "NONE") %>% filter(miss_q == "qRand" | miss_q == "NoEcov")



# Run simulation tests
foreach(iom = 1:nrow(subsetOM)) %dopar% { # Run foreach loop in parallel
  omdir <- here::here("Ecov_study", "catchability", "Results", paste0("OM_", subsetOM[iom, "OMname"]))
  
  # Pull OM
  testOM <- readRDS(here::here(omdir, paste0("OM_", subsetOM[iom, "OMname"], ".Rds")))
  
  for(iem in 1:nrow(subsetEM)){
    # Pull EM
    testEM <- readRDS(here::here(omdir, paste0("EM_missSeason_", subsetEM[iem, "miss_season"], "_missQ_", subsetEM[iem, "miss_q"]), "EMinput.Rds"))
    
    # Run simulation test
    simTestWHAM(nsim = nsim,
                OM = testOM,
                inputEMlist = list(testEM), # Run one EM at a time
                outdir = here::here(omdir, paste0("EM_missSeason_", subsetEM[iem, "miss_season"], "_missQ_", subsetEM[iem, "miss_q"])))
  }
}





#Error in unserialize(socklist[[n]]) : error reading from connection
# Gives up on OM 22 missSeason_NONE_missQ_qRandEcov