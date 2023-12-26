#' @title Post-process simulation results
#' 
#' @param filenames A vector of strings (including .RData extensions) indicating what results files to read in and post process
#' @param outdir A string for the directory where a plot folder will be generated
#' @param earlySim Boolean, TRUE = assumes EM results are stored in corresponding folder and EMshortName pulled from folder name (early sims used this option due to error with EM naming that has now been resoleved), if FALSE use EM_miss_season and EM_miss_q from simulation results. Default = FALSE.
#' 
#' @return A dataframe containing the following columns describing model performance
#' \itemize{
#'   \item{Year - Year associated with performance metric}
#'   \item{sim - Simulation number (based on total number of simulations in all files processed)}
#'   \item{EM - Name of estimation model used to calculate performance metric}
#'   \item{converged - Boolean indicating whether model converged with invertible hessian ("Yes") or not ("No")}
#'   \item{relSSB - Ratio of SSB from EM:OM}
#'   \item{relF - Ratio of F from EM:OM}
#'   \item{relR - Ratio of R from EM:OM}
#'   \item{relEcov_x - Ratio of Ecov_x from EM:OM}
#' }

# # Find all result files
# filenames1 <- list.files(path = here::here("Ecov_study/catchability/remote1_Results"), pattern = "simWHAM_", recursive = TRUE, full.names = TRUE)
# filenames2 <- list.files(path = here::here("Ecov_study/catchability/remote2_Results"), pattern = "simWHAM_", recursive = TRUE, full.names = TRUE)
# filenames3 <- list.files(path = here::here("Ecov_study/catchability/remote3_Results_supplementalRuns"), pattern = "simWHAM_", recursive = TRUE, full.names = TRUE) # simulations to fill gaps so 100 simulations for each OM/EM pair
# outdir = here::here("Ecov_study/catchability")
# postprocess_simTestWHAM(filenames = c(filenames1, filenames2), outdir = outdir)
# # file.remove(filenames) # will delete incorrect files for debugging purposes

postprocess_simTestWHAM <- function(filenames = NULL, outdir = here::here(), earlySim = FALSE){
  # Set up storage for processed results
  perfMet <- NULL
  
  # Start sim counter at 0, will increment once at the innermost loop
  simNum <- 0
  
  # Establish total number of simulations across all files (based on nsim in filename)
  nsim <- rep(NA, length(filenames))
  for(ifile in 1:length(filenames)){
    fileSpecs <- strsplit(filenames[ifile], "_")
    nsim[ifile] <- fileSpecs[[1]][which(grepl("simWHAM", fileSpecs[[1]]))+1] %>% as.numeric()
  }
  
  #nsim <- c(0, nsim) # Add 0 so isim indexing in file processed output goes from 1:total nsim
  
  # Loop over files to calculate performance metrics
  for(ifile in 1:length(filenames)){
    print(ifile)
    
    # Read in RData object
    load(file=filenames[ifile]) # loads object named "results"
    
    # Loop over simulations in each file
    for(isim in 1:nsim[ifile]){ # start at ifile+1 since nsim[ifile] = 0 for sim numbering reasons
      # print(paste0("sim ", isim))
      
      # Pull out EM names 
      EMs <- names(results)[2:length(names(results))]
      OMname <- names(results)[which(grepl("OM", names(results))==TRUE | grepl("WHAM", names(results))==TRUE)] # model with "OM" or "WHAM" in name (default OM label is "WHAM for unnamed stock")
      
      # # Set up sim numbering for ifile
      # lastSim <- ifile*nsim[ifile]*length(EMs)
      # firstSim <- lastSim+1-(nsim[ifile]*length(EMs))
      # simStorage <- firstSim:lastSim
      # simIncrement <- length(EMs)*(isim-1)
      
      # Loop over EM in each ifile
      for(iEM in 1:length(EMs)){
        # print(paste0("EM ", iEM))
        
        # Increase simNum by 1 #!!! May need to confirm that this still works if multiple EMs fit to same OM
        simNum <- simNum + 1
        
        # Dimensions
        nyear <- results[OMname][[1]][isim][[1]]$dataOM$n_years_model
        
        ##### Simulation settings #####
        # OM + EM
        seed <- rep(results[OMname][[1]][isim][[1]]$seed, nyear) 
        F_hist <- rep(results[OMname][[1]][isim][[1]]$F_hist, nyear)
        ageComp_sig <- rep(results[OMname][[1]][isim][[1]]$ageComp_sig, nyear)
        log_catch_sig <- rep(results[OMname][[1]][isim][[1]]$log_catch_sig, nyear)
        log_index_sig <- rep(results[OMname][[1]][isim][[1]]$log_index_sig, nyear)
        Year <- 1:nyear
        sim <- rep(simNum, nyear)  
        # OM
        OM_ecov_effect <- rep(results[OMname][[1]][isim][[1]]$Ecov_effect, nyear)
        OM_ecov_process_cor <- rep(results[OMname][[1]][isim][[1]]$Ecov_process_cor, nyear)
        OM_ecov_process_obs_sig <- rep(results[OMname][[1]][isim][[1]]$Ecov_process_obs_sig, nyear)
        OM_ecov_process_sig <- rep(results[OMname][[1]][isim][[1]]$Ecov_process_sig, nyear)
        OMshortName <- rep(results[OMname][[1]][isim][[1]]$OMname)
        # EM
        # EM_ecov_effect <- rep(results[EMs[iEM]][[1]][isim][[1]]$Ecov_effect, nyear) # These are initial settings only so may not need them
        # EM_ecov_process_cor <- rep(results[EMs[iEM]][[1]][isim][[1]]$Ecov_process_cor, nyear)
        # EM_ecov_process_obs_sig <- rep(results[EMs[iEM]][[1]][isim][[1]]$Ecov_process_obs_sig, nyear)
        # EM_ecov_process_sig <- rep(results[EMs[iEM]][[1]][isim][[1]]$Ecov_process_sig, nyear)
        
        if(earlySim == TRUE){
          # !!! Several model names were misslabeled so instead rely on file path name which is correct
          filepath <- filenames[ifile] %>% str_split(., "/") %>% unlist()
          filepath <- filepath[length(filepath)-1] %>% str_split(., "_")
          EM_miss_season <- filepath[[1]][3]
          EM_miss_q <- filepath[[1]][5]
          EMshortName <- paste("EM", EM_miss_season, EM_miss_q, sep="_")
        } else{ # Use model_name
          EM_miss_season <- results[EMs[iEM]][[1]][isim][[1]]$EM_miss_season
          EM_miss_q <- results[EMs[iEM]][[1]][isim][[1]]$EM_miss_q
          EMshortName <- paste("EM", EM_miss_season, EM_miss_q, sep="_")
        }
        
        
        # EM_miss_season <- rep(results[EMs[iEM]][[1]][isim][[1]]$EM_miss_season, nyear)
        # EM_miss_q <- rep(results[EMs[iEM]][[1]][isim][[1]]$EM_miss_q, nyear)
        # EMshortName <-  rep(EMs[iEM], nyear)
        # # For debugging misslabeling look at the following:
        # EMinput$data$Ecov_where #shows seasonal misspecification
        # EMinput$data$Ecov_how # shows if ecov implemented
        # EMinput$data$use_q_re # shows if qrand implemented
        
        ##### Store raw results #####
        # OM
        OM_SSB <- results[OMname][[1]][isim][[1]]$dataOM$SSB
        OM_F <-  results[OMname][[1]][isim][[1]]$dataOM$F
        colnames(OM_F) <- "OM_F"
        OM_FAA <-  results[OMname][[1]][isim][[1]]$dataOM$FAA_tot
        colnames(OM_FAA) <- paste("OM_FAA", 1:10, sep="_")
        OM_R <-  results[OMname][[1]][isim][[1]]$dataOM$NAA[,1]
        OM_NAA <-  results[OMname][[1]][isim][[1]]$dataOM$NAA
        colnames(OM_NAA) <- paste("OM_NAA", 1:10, sep="_")
        OM_Catch <-  results[OMname][[1]][isim][[1]]$dataOM$pred_catch
        colnames(OM_Catch) <- "OM_Catch"
        OM_CAA <- results[OMname][[1]][isim][[1]]$dataOM$pred_CAA[,1,]
        colnames(OM_CAA) <- paste("OM_CAA", 1:10, sep="_")
        OM_FMSY <-  rep(exp(results[OMname][[1]][isim][[1]]$dataOM$log_FXSPR_static), nyear) # F40
        OM_SSBMSY <- rep(exp(results[OMname][[1]][isim][[1]]$dataOM$log_SSB_FXSPR_static), nyear) # at F40
        OM_MSY <- rep(exp(results[OMname][[1]][isim][[1]]$dataOM$log_Y_FXSPR_static), nyear)
        OM_selAA_cat <- results[OMname][[1]][isim][[1]]$dataOM$selAA[[1]]
        colnames(OM_selAA_cat) <- paste("OM_selCat", 1:10, sep="_")
        OM_selAA_ind1 <- results[OMname][[1]][isim][[1]]$dataOM$selAA[[2]]
        colnames(OM_selAA_ind1) <- paste("OM_selInd1", 1:10, sep="_")
        OM_selAA_ind2 <- results[OMname][[1]][isim][[1]]$dataOM$selAA[[3]]
        colnames(OM_selAA_ind2) <- paste("OM_selInd2", 1:10, sep="_")
        OM_q <- results[OMname][[1]][isim][[1]]$dataOM$q # For spring and fall surveys respectively
        colnames(OM_q) <- c("OM_q_index1", "OM_q_index2")
        OM_Ecov_obs <- results[OMname][[1]][isim][[1]]$dataOM$Ecov_obs # Observations
        colnames(OM_Ecov_obs) <- "OM_Ecov_obs"
        
        
        # EM
        if(results[EMs[iEM]][[1]][isim][[1]]$whamConverge == TRUE){
          Converged <- results[EMs[iEM]][[1]][isim][[1]]$whamConverge
          AIC <- results[EMs[iEM]][[1]][isim][[1]]$AIC
          EM_SSB <- results[EMs[iEM]][[1]][isim][[1]]$SSB
          EM_F <-  results[EMs[iEM]][[1]][isim][[1]]$F
          colnames(EM_F) <- "EM_F"
          EM_FAA <-  results[EMs[iEM]][[1]][isim][[1]]$FAA
          colnames(EM_FAA) <- paste("EM_FAA", 1:10, sep="_")
          EM_R <-  results[EMs[iEM]][[1]][isim][[1]]$R
          EM_NAA <-  results[EMs[iEM]][[1]][isim][[1]]$NAA
          colnames(EM_NAA) <- paste("EM_NAA", 1:10, sep="_")
          EM_Catch <-  results[EMs[iEM]][[1]][isim][[1]]$Catch
          colnames(EM_Catch) <- "EM_Catch"
          EM_CAA <- results[EMs[iEM]][[1]][isim][[1]]$CAA #!!! Double check that this gets saved correctly
          colnames(EM_CAA) <- paste("EM_CAA", 1:10, sep="_")
          EM_FMSY <-  rep(results[EMs[iEM]][[1]][isim][[1]]$FMSY, nyear)
          EM_SSBMSY <- rep( results[EMs[iEM]][[1]][isim][[1]]$SSBMSY, nyear)
          EM_MSY <- rep(results[EMs[iEM]][[1]][isim][[1]]$MSY, nyear)
          EM_selAA_cat <- results[EMs[iEM]][[1]][isim][[1]]$SelAA[[1]]
          colnames(EM_selAA_cat) <- paste("EM_selCat", 1:10, sep="_")
          EM_selAA_ind1 <- results[EMs[iEM]][[1]][isim][[1]]$SelAA[[2]]
          colnames(EM_selAA_ind1) <- paste("EM_selInd1", 1:10, sep="_")
          EM_selAA_ind2 <- results[EMs[iEM]][[1]][isim][[1]]$SelAA[[3]]
          colnames(EM_selAA_ind2) <- paste("EM_selInd2", 1:10, sep="_")
          EM_q <- results[EMs[iEM]][[1]][isim][[1]]$pars_q # For spring and fall surveys respectively
          colnames(EM_q) <- c("EM_q_index1", "EM_q_index2")
          EM_ecovBeta_ind1 <- results[EMs[iEM]][[1]][isim][[1]]$pars_ecovBeta_ind1
          EM_ecovBeta_ind2 <- results[EMs[iEM]][[1]][isim][[1]]$pars_ecovBeta_ind2
          EM_q_re <- results[EMs[iEM]][[1]][isim][[1]]$q_re # q random effect
          colnames(EM_q_re) <- c("EM_qre_Ind1", "EM_qre_Ind2")
          EM_Ecov_re <- results[EMs[iEM]][[1]][isim][[1]]$Ecov_re # Ecov random effect
          colnames(EM_Ecov_re) <- "EM_Ecov_re"
          #!!! pars_Ecov_process = 3 parameters, what is order for labeling purposes here?
          EM_Ecov_pred <- results[EMs[iEM]][[1]][isim][[1]]$Ecov_x # EM predictions
          colnames(EM_Ecov_pred) <- "EM_Ecov_pred"
          # EM_Ecov_pred <- NA #!!! placeholder for EM predictions for now since models prior to 12/20/23 at 4pm do not have the data saved to support this
          
          # Combine settings and raw results into a single storage data.frame 
          storage <- cbind(seed, F_hist, ageComp_sig, log_catch_sig, log_index_sig, Year, sim, 
                           OM_ecov_effect, OM_ecov_process_cor, OM_ecov_process_obs_sig, OM_ecov_process_sig, OMshortName,
                           EM_miss_season, EM_miss_q, EMshortName,
                           OM_SSB, OM_F, OM_FAA, OM_R, OM_NAA, OM_Catch, OM_CAA, OM_FMSY, OM_SSBMSY, OM_MSY, OM_selAA_cat, OM_selAA_ind1, OM_selAA_ind2, OM_q, OM_Ecov_obs,
                           Converged, AIC, EM_SSB, EM_F, EM_FAA, EM_R, EM_NAA, EM_Catch, EM_CAA, EM_FMSY, EM_SSBMSY, EM_MSY, EM_selAA_cat, EM_selAA_ind1, EM_selAA_ind2, EM_q, EM_ecovBeta_ind1, EM_ecovBeta_ind2, EM_q_re, EM_Ecov_re, EM_Ecov_pred) %>% 
            as.data.frame() 
          numericIndex <- which(colnames(storage) %in% c("OMshortName", 'EMshortName', "EM_miss_season", "EM_miss_q", "F_hist", "Converged") == FALSE)
          storage[,numericIndex] <- sapply(storage[,numericIndex], as.numeric)
          
          
          
          ##### Calculate performance metrics #####
          
          # Mohn's rho
          EM_MohnsRho_SSB <- rep(results[EMs[iEM]][[1]][isim][[1]]$MohnsRho_SSB, nyear)
          EM_MohnsRho_F <- rep(results[EMs[iEM]][[1]][isim][[1]]$MohnsRho_F, nyear)
          EM_MohnsRho_R <- rep(results[EMs[iEM]][[1]][isim][[1]]$MohnsRho_R, nyear)
          storage <- cbind(storage, EM_MohnsRho_SSB, EM_MohnsRho_F, EM_MohnsRho_R) # Append to storage
          
          # EM convergence check
          EM_converged <- rep(results[EMs[iEM]][[1]][isim][[1]]$whamConverge, nyear)
          storage <- cbind(storage, EM_converged) # Append to storage
          
          # Calculate all other performance metrics
          simPerfMet <- storage %>%      # Output settings
            mutate(Converged = results[EMs[iEM]][[1]][isim][[1]]$whamConverge,
                   
                   ##### EM results #####
                   EM_MohnsRho_SSB = EM_MohnsRho_SSB,
                   EM_MohnsRho_F = EM_MohnsRho_F,
                   EM_MohnsRho_R = EM_MohnsRho_R,
                   EM_converged = EM_converged,
                   status_SSB = EM_SSB/EM_SSBMSY,
                   status_F = EM_F/EM_FMSY,
                   status_Y = EM_Catch/EM_MSY,
                   EM_SSB = EM_SSB,
                   EM_R = EM_R,
                   EM_F = EM_F,
                   EM_ecovBeta_ind1 = EM_ecovBeta_ind1,
                   EM_ecovBeta_ind2 = EM_ecovBeta_ind2,
                   OM_ecov_effect_beta = OM_ecov_effect,
                   
                   ##### EM/OM relative statistics #####
                   relSSB = EM_SSB/OM_SSB,
                   relF = EM_F/OM_F,
                   relR = EM_R/OM_R,
                   relSSBMSY = EM_SSBMSY/OM_SSBMSY,
                   relFMSY = EM_FMSY/OM_FMSY,
                   relMSY = EM_MSY/OM_MSY,
                   relq_index1 = EM_q_index1/OM_q_index1,
                   relq_index2 = EM_q_index2/OM_q_index2,
                   relEcovBeta_ind1 = EM_ecovBeta_ind1/OM_ecov_effect,
                   relEcovBeta_ind2 = EM_ecovBeta_ind2/OM_ecov_effect,
                   # Relative FAA
                   relFAA_1 = EM_FAA_1/OM_FAA_1, 
                   relFAA_2 = EM_FAA_2/OM_FAA_2,
                   relFAA_3 = EM_FAA_3/OM_FAA_3,
                   relFAA_4 = EM_FAA_4/OM_FAA_4,
                   relFAA_5 = EM_FAA_5/OM_FAA_5,
                   relFAA_6 = EM_FAA_6/OM_FAA_6,
                   relFAA_7 = EM_FAA_7/OM_FAA_7,
                   relFAA_8 = EM_FAA_8/OM_FAA_8, 
                   relFAA_9 = EM_FAA_9/OM_FAA_9, 
                   relFAA_10 = EM_FAA_10/OM_FAA_10,
                   # Relative NAA
                   relNAA_1 = EM_NAA_1/OM_NAA_1,
                   relNAA_2 = EM_NAA_2/OM_NAA_2,
                   relNAA_3 = EM_NAA_3/OM_NAA_3,
                   relNAA_4 = EM_NAA_4/OM_NAA_4,
                   relNAA_5 = EM_NAA_5/OM_NAA_5,
                   relNAA_6 = EM_NAA_6/OM_NAA_6,
                   relNAA_7 = EM_NAA_7/OM_NAA_7,
                   relNAA_8 = EM_NAA_8/OM_NAA_8,
                   relNAA_9 = EM_NAA_9/OM_NAA_9,
                   relNAA_10 = EM_NAA_10/OM_NAA_10,
                   # Relative CAA
                   relCAA_1 = EM_CAA_1/OM_CAA_1, 
                   relCAA_2 = EM_CAA_2/OM_CAA_2,
                   relCAA_3 = EM_CAA_3/OM_CAA_3,
                   relCAA_4 = EM_CAA_4/OM_CAA_4,
                   relCAA_5 = EM_CAA_5/OM_CAA_5,
                   relCAA_6 = EM_CAA_6/OM_CAA_6,
                   relCAA_7 = EM_CAA_7/OM_CAA_7,
                   relCAA_8 = EM_CAA_8/OM_CAA_8,
                   relCAA_9 = EM_CAA_9/OM_CAA_9,
                   relCAA_10 = EM_CAA_10/OM_CAA_10,
                   # Relative Catch selAA
                   relselCat_1 = EM_selCat_1/OM_selCat_1,
                   relselCat_2 = EM_selCat_2/OM_selCat_2,
                   relselCat_3 = EM_selCat_3/OM_selCat_3,
                   relselCat_4 = EM_selCat_4/OM_selCat_4,
                   relselCat_5 = EM_selCat_5/OM_selCat_5,
                   relselCat_6 = EM_selCat_6/OM_selCat_6,
                   relselCat_7 = EM_selCat_7/OM_selCat_7,
                   relselCat_8 = EM_selCat_8/OM_selCat_8,
                   relselCat_9 = EM_selCat_9/OM_selCat_9,
                   relselCat_10 = EM_selCat_10/OM_selCat_10,
                   # Relative index 1 selAA
                   relselInd1_1 = EM_selInd1_1/OM_selInd1_1,
                   relselInd1_2 = EM_selInd1_2/OM_selInd1_2,
                   relselInd1_3 = EM_selInd1_3/OM_selInd1_3,
                   relselInd1_4 = EM_selInd1_4/OM_selInd1_4,
                   relselInd1_5 = EM_selInd1_5/OM_selInd1_5,
                   relselInd1_6 = EM_selInd1_6/OM_selInd1_6,
                   relselInd1_7 = EM_selInd1_7/OM_selInd1_7,
                   relselInd1_8 = EM_selInd1_8/OM_selInd1_8,
                   relselInd1_9 = EM_selInd1_9/OM_selInd1_9,
                   relselInd1_10 = EM_selInd1_10/OM_selInd1_10,
                   # Relative index 2 sel AA
                   relselInd2_1 = EM_selInd2_1/OM_selInd2_1,
                   relselInd2_2 = EM_selInd2_2/OM_selInd2_2,
                   relselInd2_3 = EM_selInd2_3/OM_selInd2_3,
                   relselInd2_4 = EM_selInd2_4/OM_selInd2_4,
                   relselInd2_5 = EM_selInd2_5/OM_selInd2_5,
                   relselInd2_6 = EM_selInd2_6/OM_selInd2_6,
                   relselInd2_7 = EM_selInd2_7/OM_selInd2_7,
                   relselInd2_8 = EM_selInd2_8/OM_selInd2_8,
                   relselInd2_9 = EM_selInd2_9/OM_selInd2_9,
                   relselInd2_10 = EM_selInd2_10/OM_selInd2_10) # Currently don't use EM_qre_Ind1, EM_qre_Ind2, EM_Ecov_re, or pars_Ecov_process
          
          # Append perfMets 
          perfMet <- rbind(perfMet, simPerfMet)
        } else{ # If EM did not converge for the simulation save only OM results and EM convergence status
          Converged <- results[EMs[iEM]][[1]][isim][[1]]$whamConverge
          AIC <- NA
          EM_SSB = NA
          EM_F = NA
          EM_FAA = NA
          EM_R = NA
          EM_NAA = NA 
          EM_Catch = NA
          EM_CAA = NA
          EM_FMSY = NA
          EM_SSBMSY = NA
          EM_MSY = NA
          EM_selAA_cat = NA
          EM_selAA_ind1 = NA
          EM_selAA_ind2 = NA
          EM_q = NA
          EM_ecovBeta_ind1 = NA
          EM_ecovBeta_ind2 = NA
          EM_q_re = NA
          EM_Ecov_re = NA
          EM_Ecov_pred <- NA
          
          # storage <- cbind(seed, F_hist, ageComp_sig, log_catch_sig, log_index_sig, Year, sim, 
          #                  OM_ecov_effect, OM_ecov_process_cor, OM_ecov_process_obs_sig, OM_ecov_process_sig, OMshortName,
          #                  EM_miss_season, EM_miss_q, EMshortName,
          #                  Converged, AIC,
          #                  matrix(rep(NA, (ncol(perfMet)-17)*length(F_hist)), ncol = (ncol(perfMet)-17))) %>% # Fill remaining performance metrics with NAs
          #   as.data.frame() 
          storage <- cbind(seed, F_hist, ageComp_sig, log_catch_sig, log_index_sig, Year, sim, 
                           OM_ecov_effect, OM_ecov_process_cor, OM_ecov_process_obs_sig, OM_ecov_process_sig, OMshortName,
                           EM_miss_season, EM_miss_q, EMshortName,
                           OM_SSB, OM_F, OM_FAA, OM_R, OM_NAA, OM_Catch, OM_CAA, OM_FMSY, OM_SSBMSY, OM_MSY, OM_selAA_cat, OM_selAA_ind1, OM_selAA_ind2, OM_q, OM_Ecov_obs,
                           Converged, AIC,  #EM_SSB, EM_F, EM_FAA, EM_R, EM_NAA, EM_Catch, EM_CAA, EM_FMSY, EM_SSBMSY, EM_MSY, EM_selAA_cat, EM_selAA_ind1, EM_selAA_ind2, EM_q, EM_ecovBeta_ind1, EM_ecovBeta_ind2, EM_q_re, EM_Ecov_re, EM_Ecov_pred, 
                           matrix(rep(NA, (ncol(perfMet)-87)*length(F_hist)), ncol = (ncol(perfMet)-87))) %>% # Fill remaining performance metrics with NAs
            as.data.frame() 
          
          names(storage) <- names(perfMet)
          # numericIndex <- which(colnames(storage) %in% c("OMshortName", 'EMshortName', "EM_miss_season", "EM_miss_q", "F_hist", "Converged") == FALSE)
          # storage[,numericIndex] <- sapply(storage[,numericIndex], as.numeric) # This introduces NAs by coercion since using NAs as placeholders
          
          
          # Append perfMets
          perfMet <- rbind(perfMet, storage)
          
        } # End handling for unconverged EMs    
      } # End loop over EM
      
    } # End loop over simulations in ifile
    
  } # End loop over filenames
  
  # Save returned perfMet in outdir
  timeStamp <- Sys.time() %>% gsub(" ", "_",.) %>% gsub(":", "-", .) # Change millisecond half of Sys.time() output to avoid having spaces/weird characters in filenames
  saveRDS(perfMet, file = paste0(outdir, paste0("/perfMet_",timeStamp, ".RDS")))
  
  # Return
  return(perfMet)
}

