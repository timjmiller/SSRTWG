#' @param filenames A vector of strings (including .RData extensions) indicating what results files to read in and post process
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


postprocess_simTestWHAM <- function(filenames = NULL){
        # Set up storage for processed results
        perfMet <- NULL
        
        # Establish total number of simulations across all files (based on nsim in filename)
        nsim <- rep(NA, length(filenames))
        for(ifile in 1:length(filenames)){
                fileSpecs <- strsplit(filenames[ifile], "_")
                nsim[ifile] <- fileSpecs[[1]][which(grepl("simWHAM", fileSpecs[[1]]))+1] %>% as.numeric()
        }
        
        nsim <- c(0, nsim) # Add 0 so isim indexing in file processed output goes from 1:total nsim
        
        # Loop over files to calculate performance metrics
        for(ifile in 1:length(filenames)){
                # Read in RData object
                load(file=filenames[ifile]) # loads object named "results"
                
                # Loop over simulations in each file
                for(isim in 1:nsim[ifile+1]){ # start at ifile+1 since nsim[ifile] = 0 for sim numbering reasons
                        
                        # Pull out EM names 
                        EMs <- names(results)[2:length(names(results))]
                        OMname <- names(results)[which(grepl("OM", names(results))==TRUE | grepl("WHAM", names(results))==TRUE)] # model with "OM" or "WHAM" in name (default OM label is "WHAM for unnamed stock")
                        
                        
                        # Loop over EM in each ifile
                        for(iEM in 1:length(EMs)){
                                
                                
                                # Store Years, isim, iEM, convergence (mod$opt$convergence == 0 and hessian invertible $na_sdrep == FALSE)
                                Year <- results[EMs[iEM]][[1]][isim][[1]]$years
                                sim <- rep((isim+nsim[ifile]), length(Year)) # + 0 for first file, starts sim numbering at nsim[ifile=1]+1 for second file
                                EM <- rep(EMs[iEM], length(Year))
                                if(is.na(results[EMs[iEM]][[1]][[isim]]$na_sdrep) | results[EMs[iEM]][[1]][[isim]]$na_sdrep == TRUE) { # If NA in hessian diagonal or na_sdrep not evaluated (NA) model didn't converge
                                        converged <- rep("No", length(Year))
                                } else{
                                        converged <- rep("Yes", length(Year))
                                }
                                # converged <- rep(results[EMs[iEM]][[1]][[isim]]$converged, length(Year))
                                
                                
                                # Proportion of EM:OM SSB
                                # EM SSB calculated based on plotting script https://github.com/timjmiller/wham/blob/6c30301d5b68c2bb0c0330cc4f112432c7950dfd/R/wham_plots_tables.R#L725-L738
                                # tempSSB_EM <- results[EMs[iEM]][[1]][isim][[1]]$sdrep %>% summary() #EM SSB must be calculated, this line triggers an NaNs produced warning for EMs that didn't converge with invertible hessians
                                # tempSSB_EM <- tempSSB_EM[rownames(tempSSB_EM) == "log_SSB",]
                                # tempSSB_EM <- exp(cbind(tempSSB_EM, tempSSB_EM[,1]+qnorm(0.975)*cbind(-tempSSB_EM[,2], tempSSB_EM[,2])))#/1000 Don't divide by 1000 so units in mt rather than kmt
                                # tempSSB_EM <- tempSSB_EM[,1] 
                                tempSSB_EM <- results[EMs[iEM]][[1]][isim][[1]]$rep$SSB
                                # OM SSB
                                tempSSB_OM <- results[OMname][[1]][isim][[1]]$dataOM$SSB # Pull SSB directly from simulated data
                                # relative SSB
                                relSSB <- tempSSB_EM/tempSSB_OM
                                
                                
                                # Proportion of EM:OM F
                                # EM calculated F based on plotting scripts
                                # mod <- results[EMs[iEM]][[1]][isim][[1]]
                                # std = results[EMs[iEM]][[1]][isim][[1]]$sdrep %>% summary()
                                # years_full <- mod$years_full
                                # n_ages = mod$env$data$n_ages
                                # faa.ind <- which(rownames(std) == "log_FAA_tot")
                                # log.faa <- matrix(std[faa.ind,1], length(years_full), n_ages)
                                # faa.cv <- matrix(std[faa.ind,2], length(years_full), n_ages)
                                # age.full.f <- apply(log.faa,1, function(x) max(which(x == max(x))))
                                # full.f.ind = cbind(1:length(years_full), age.full.f)
                                # log.full.f <- log.faa[full.f.ind]
                                # full.f.cv <- faa.cv[full.f.ind]
                                # full.f <- exp(log.full.f)
                                full.f <- results[EMs[iEM]][[1]][isim][[1]]$rep$F
                                # OM F
                                tempF_OM <- c(results[OMname][[1]][isim][[1]]$dataOM$F)
                                # relative F
                                relF <- full.f/tempF_OM
                                
                                
                                # Proportion of EM:OM Recruitment 
                                #   # EM R calculated based on plotting script: https://github.com/timjmiller/wham/blob/6c30301d5b68c2bb0c0330cc4f112432c7950dfd/R/wham_plots_tables.R#L765-L780
                                # tempR_EM <-  results[EMs[iEM]][[1]][isim][[1]]$sdrep %>% summary()
                                # ind <- rownames(tempR_EM) == "log_NAA_rep"
                                # tempR_EM <-  exp(array(tempR_EM[ind,1], dim = c(length(results[EMs[iEM]][[1]][isim][[1]]$years), 8)))
                                # tempR_EM <- tempR_EM[,1]
                                tempR_EM <- results[EMs[iEM]][[1]][isim][[1]]$rep$NAA[,1]
                                # OM recruitment
                                tempR_OM <- results[OMname][[1]][isim][[1]]$dataOM$NAA[,1]
                                # relative R
                                relR <- tempR_EM/tempR_OM
                                
                                
                                # # Relative Ecov_x values (EM:OM ratio) #!!!! revisit
                                # if("TRUE" %in% (results[EMs[iEM]][[1]][isim][[1]]$input$data$Ecov_where > 0)){ # If link to ecov specified in EM calculate EM:OM ratio for Ecov_x values
                                #     # EM
                                # tempEcovx_EM <- results[EMs[iEM]][[1]][isim][[1]]$sdrep %>% summary() 
                                # tempEcovx_EM <- tempEcovx_EM[rownames(tempEcovx_EM)=="Ecov_x",]
                                #   # OM
                                # tempEcovx_OM <- results$OM[isim][[1]]$Ecov_x
                                #   # Relative Ecov_x
                                # relEcov_x <- tempEcovx_EM[,"Estimate"]/tempEcovx_OM
                                # } else{ # If no ecov link specified then save an NA for that model in all years 
                                #   relEcov_x <- rep(NA, length(Year))
                                # }
                                
                                
                                # # Relative Ecov_beta (EM:OM ratio) for all four indices
                                #   # EM
                                # tempEcovbeta_EM <- results[EMs[iEM]][[1]][isim][[1]]$sdrep %>% summary()
                                # tempEcovbeta_EM <- tempEcovbeta_EM[rownames(tempEcovbeta_EM)=="Ecov_beta"]
                                #   # OM
                                # tempEcovbeta_OM <- results$OM[isim][[1]]$Ecov_beta # by age rather than by index?
                                # 
                                
                                
                                # Store results
                                # EM perfMets
                                EMperfMets <- cbind(Year, sim, EM, converged, relSSB, relF, relR)#, relEcov_x) #!!! revisit
                                
                                # Append perfMets to storage
                                perfMet <- rbind(perfMet, EMperfMets)
                                
                        } # End loop over EM
                        
                } # End loop over simulations in ifile
                
        } # End loop over filenames
        
        # Clean table rownames and ensure that numeric results are not treated as strings
        rownames(perfMet) <- NULL
        colnames(perfMet) <- c("Year", "sim", "EM", "converged", "relSSB", "relF", "relR")# , "relEcov_x") #!!! revisit
        perfMet <- perfMet %>% as.data.frame() %>% 
                dplyr::summarise(Year = as.numeric(Year),
                                 sim = as.numeric(sim),
                                 EM = EM,
                                 converged = converged,
                                 relSSB = as.numeric(relSSB),
                                 relF = as.numeric(relF),
                                 relR = as.numeric(relR))# , # Revisit !!!
                                 # relEcov_x = as.numeric(relEcov_x)) # Revisit !!!
        
        
        return(perfMet)
}

