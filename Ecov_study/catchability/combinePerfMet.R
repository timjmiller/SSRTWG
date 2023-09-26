#' @title Combine processed performance metrics
#' @description
#' Combine multiple files with post-processed performance metrics, the major action this function performs is renumbering sims so each sim has a unique number
#' 
#' @param filenames A vector of strings (including full file path and .RData extensions) for performance metrics files to 



combinePerfMet <- function(filenames=NULL){
  
  # Read in file
  
  # if first file do nothing
  
  # For subsequent files:
    # Find unique sim
    # generate unique sims starting 1 after largest sim in first file
    # replace original unique sim with regenerated values that continue numbering from prior file
  
  # save results
  
}