#' @title Combine processed performance metrics
#' @description
#' Combine multiple files with post-processed performance metrics, the major action this function performs is renumbering sims so each sim has a unique number
#' 
#' @param filenames A vector of strings (including full file path and .RData extensions) for performance metrics files to 
#' @param outdir A string for the directory where a plot folder will be generated



combinePerfMet <- function(filenames=NULL, outdir = here::here()){
  
  aggPerfMet <- NULL
  
  for(ifile in 1:length(filenames)){
    print(ifile)
    
    # Read in RDS file
    perfMet <- readRDS(file = here::here(filenames[ifile]))
    
    
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
  }
 
  # Save returned aggPerfMet in outdir
  timeStamp <- Sys.time() %>% gsub(" ", "_",.) %>% gsub(":", "-", .) # Change millisecond half of Sys.time() output to avoid having spaces/weird characters in filenames
  saveRDS(aggPerfMet, file = paste0(outdir, paste0("/aggPerformanceMetric_",timeStamp, ".Rds")))
  
}