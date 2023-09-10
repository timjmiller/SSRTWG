#' @title  Generate simulated OM data
#' 
#' @param OM A fitted WHAM model to use as the operating model, no default.
#' @param seed A number to set the seed for OM simulated data, no default. 
#' @param isim A number indicating what simulation data is generated for, no default.
#' 
#' @return A list containing:
#' \itemize{
#'   \item{dataOM - A simulated OM dataset, may be used to overwrite input$data either entirely or in part in simulations}
#'   \item{seed - A number used to set the seed for OM simulated data}
#' }

# Function based on simulation testing vignette:https://github.com/timjmiller/wham/blob/devel/vignettes/ex10_simulation.Rmd
        
simOM <- function(OM=NULL,seed=NULL, isim = NULL){
        # input <- OM$input
        set.seed(seed[isim])
        dataOM <- OM$simulate(complete=TRUE) # Overwrite input with simulated data & return
        
        # Return
        returnList <- NULL
        returnList$dataOM <- dataOM
        returnList$seed <- seed[isim]
        return(returnList)
}
