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
                #!!! Save the seed somewhere here
                
                
                print(paste0("Simulation_", isim))
                for(iEM in 1:length(inputEMlist)){ # Loop over estimation models
                        # Pull EM initial input from list
                        inputEM <- inputEMlist[[iEM]] 
                        
                        # Update EM input with OM simulated data
                        inputEM$data[obs_names] = sim_inputs[[isim]]$dataOM[obs_names]
                        
                        # Fit WHAM model to updated EM input and save results
                        resultEM <- fit_wham(inputEM, do.osa = FALSE, retro.silent = TRUE, MakeADFun.silent = TRUE, save.sdrep = FALSE) #!!!! Figure out how to store
                        
                        # # Check convergence
                        # if(check_convergence(resultEM, ret=TRUE)$is_sdrep == FALSE){ # If sdrep had error
                        #         converged <- "No"
                        # }
                        # if(check_convergence(resultEM, ret=TRUE)$is_sdrep == TRUE & !is.null(check_convergence(resultEM, ret=TRUE)$na_sdrep)){
                        #         if(check_convergence(resultEM, ret=TRUE)$na_sdrep == TRUE){
                        #                 converged <- "No" # Not converged due to NA on diagonal of hessian
                        #         } else{
                        #                 converged <- "Yes"
                        #         }
                        # }
                        # resultEM$converged <- converged
                        
                        # Store model-specific results for given isim here (append EM results to storage list for that EM)
                        storeEM[which(namesEM == namesEM[iEM])][[1]] <- append(storeEM[which(namesEM == namesEM[iEM])][[1]], list(resultEM))
                        
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



