#' @title  Fit specified EMs to simulated OM data and store results
#' 
#' @param nsim A number specifying how many simulations to run, default = 1.
#' @param OM A fitted WHAM model to use as an operating model, no default.
#' @param inputEMlist A list of WHAM input data objects for all estimation models, no default.
#' @param outdir A filename specifying the directory where time-stamped results should be stored, if NULL (default) then the results are not saved as RData,
#' @param seeds A vector of length nsim containing seed numbers for OM simulated data, no default. If not provided a random number used to set the seed and is returned for reproducibility.
#' 
#' @return A nested list with names and results corresponding to the operating models and all estimation models. For each model results contains a sub-list with items for each simulation's results stored as sub-sub lists (because WHAM results are stored as a list). If outdir != NULL then results are saved as RData and returned as an object from the function.

# Function based on simulation testing vignette:https://github.com/timjmiller/wham/blob/devel/vignettes/ex10_simulation.Rmd

# simTestWHAM(nsim = 1, OM = OM1, inputEMlist = list(genericInput), outdir = here::here("Ecov_study", "catchability"))

simTestWHAM <- function(nsim = 1,
                        OM = NULL,
                        inputEMlist = NULL,
                        outdir = NULL,
                        seeds = NULL){
        # Simulation set-up
        obs_names = c("agg_indices","agg_catch","catch_paa","index_paa", "Ecov_obs", "obsvec") # Data overwritten by OM
        namesEM <- lapply(inputEMlist, FUN = function(inputEMlist){inputEMlist$model_name}) %>% unlist() # EM model names
        
        # Storage set-up
        results <- NULL #!!!!! Still need to work on setting this up
        storeEM <- lapply(namesEM, FUN = assign, value = list())
        names(storeEM) <- namesEM
        sim_inputs <- vector(mode='list', length = nsim)
        
        # Set seeds for all simulation 
        if(is.null(seeds)==TRUE){ # If vector of seeds (of length nsim) is not provided, set up random seed for OM data generation
          seeds <- sample(1:1000000, size = nsim, replace = FALSE)
        }
        
        
        # Loop over simulations
        for(isim in 1:nsim){ 
                # Generate OM data for each simulation
                sim_inputs[[isim]] <- simOM(OM, seed = seeds, isim = isim) #!!! this now returns is a list of dataOM and seed for each replicate
                
                #!!! store specification details here !!!!!!!!!!!!
                
                
                print(paste0("Simulation_", isim))
                for(iEM in 1:length(inputEMlist)){ # Loop over estimation models
                  # Set up temporary storage for this EM/sim
                  tempStore <- NULL
                  
                  # Store EM misspecification and simulation seed (also in OM results)
                  tempStore$seed <- sim_inputs[[isim]]$seed
                  # !!! store misspecification here !!!!!!!!
                  
                        # Pull EM initial input from list
                        inputEM <- inputEMlist[[iEM]] 
                        
                        # Update EM input with OM simulated data
                        inputEM$data[obs_names] = sim_inputs[[isim]]$dataOM[obs_names]
                        
                        # Fit WHAM model to updated EM input and save results
                        fitEM <- fit_wham(inputEM, do.osa = FALSE, retro.silent = TRUE, MakeADFun.silent = TRUE, save.sdrep = FALSE) #!!!! Figure out how to store
                        
                        # Calculate Mohn's rho
                        MohnsRho <- mohns_rho(fitEM)
                        
                        # # Check convergence !!! get this working, look back at cod code
                        # if(check_convergence(fitEM, ret=TRUE)$is_sdrep == FALSE){ # If sdrep had error
                        #         converged <- "No"
                        # }
                        # if(check_convergence(fitEM, ret=TRUE)$is_sdrep == TRUE & !is.null(check_convergence(fitEM, ret=TRUE)$na_sdrep)){
                        #         if(check_convergence(fitEM, ret=TRUE)$na_sdrep == TRUE){
                        #                 converged <- "No" # Not converged due to NA on diagonal of hessian
                        #         } else{
                        #                 converged <- "Yes"
                        #         }
                        # }
                        # fitEM$converged <- converged
                        
                        # Pull together EM results for this EM/sim
                        tempStore$SSB <- fitEM$rep$SSB
                        tempStore$F <- fitEM$rep$F
                        tempStore$FAA <- fitEM$rep$FAA_tot
                        tempStore$R <- fitEM$rep$NAA[,1]
                        tempStore$NAA <- fitEM$rep$NAA
                        tempStore$Catch <- fitEM$rep$pred_catch
                        tempStore$CAA <- fitEM$rep$rep$pred_CAA[,1,]
                        tempStore$FMSY <- exp(fitEM$rep$log_FXSPR_static) # F40
                        tempStore$SSBMSY <- exp(fitEM$rep$log_SSB_FXSPR_static) # at F40
                        tempStore$MSY <- exp(fitEM$rep$log_Y_FXSPR_static) # at F40
                        tempStore$SelAA <- fitEM$rep$selAA
                        # check_convergence !!! get this working and store here
                        tempStore$MohnsRho_SSB <- MohnsRho["SSB"]
                        tempStore$MohnsRho_F <- MohnsRho["Fbar"]
                        tempStore$MohnsRho_R <- MohnsRho["R"]
                        tempStore$pars_Ecov_beta <- fitEM$rep$Ecov_beta[3,,1,]
                        tempStore$pars_Ecov_process <- fitEM$rep$Ecov_process_pars
                        tempStore$pars_q <- fitEM$rep$q
                        
                        # Store model-specific results for given isim here (append tempStore EM results to storage list for that EM)
                        storeEM[which(namesEM == namesEM[iEM])][[1]] <- append(storeEM[which(namesEM == namesEM[iEM])][[1]], list(tempStore))
                        
                } # End loop over EM
        } # End loop over nsims
        
        # Store all results
        results$OM <- sim_inputs # Includes seed numbers for each sim
        names(results) <- OM$model_name # Name OM column based on model_name
        results <- append(results, storeEM) # Should add named EM results to the final results object
        
        # Save time-stamped results as RData if output directory provided as argument
        timeStamp <- Sys.time() %>% gsub(" ", "_",.) %>% gsub(":", "-", .) # Change millisecond half of Sys.time() output to avoid having spaces/weird characters in filenames
        if(is.null(outdir)==FALSE){
                save(results, file = paste(outdir, paste0("simWHAM_", nsim, "_nsim_", OM$model_name, "_OM_", timeStamp, ".RData"), sep="/"))
        }
        
        return(results) # Always return results as an object
}





