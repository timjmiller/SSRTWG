#' Set simulation options
#'
#' @param input model input list that would be provided to \code{\link{fit_wham}}.
#' @param simulate_data TRUE/FALSE, 1 value or a vector of 3 for Catch, Indices, Ecov_obs in that order. Default = \code{TRUE}.
#' @param simulate_process TRUE/FALSE, 1 value or a vector of 5 for NAA, MAA, sel, Ecov, q in that order. Default = \code{TRUE}.
#' @param simulate_projection TRUE/FALSE, Whether to simulate observations and processes during the projection period (if any). Default = \code{TRUE}.
#' @param bias_correct_oe TRUE/FALSE, Whether to bias correct simulated aggregate catch and index to be mean-unbiased on natural scale. Default = \code{FALSE}.
#' @param bias_correct_pe TRUE/FALSE, Whether to bias correct simulated log-normal process (random effects) errors for numbers at age and M. Default = \code{FALSE}.
#' 
#' @return the modified input that would be provided ot \code{\link{fit_wham}}.
#'
#' @export
#' 
#' @seealso \code{\link{fit_wham}}, \code{\link{prepare_wham_input}}
#'
set_simulation_options <- function(input, simulate_data = TRUE, simulate_process = TRUE, simulate_projection = TRUE, 
  bias_correct_pe = FALSE, bias_correct_oe = FALSE){
  
  data = input$data
  if(!all(is.logical(simulate_data))) stop("simulate_data must be TRUE or FALSE values")
  if(!(length(simulate_data) %in% c(1, length(data$simulate_data)))){
    stop(paste0("simulate_data must be a single value or length = ", length(data$simulate_data)))
  }
  if(!all(is.logical(simulate_process))) stop("simulate_process must be TRUE or FALSE values")
  if(!(length(simulate_process) %in% c(1, length(data$simulate_process)))){
    stop(paste0("simulate_process must be a single value or length = ", length(data$simulate_process)))
  }
  if(!all(is.logical(simulate_projection))) stop("simulate_projection must be TRUE or FALSE")
  if(!(length(simulate_projection) == 1)){
    stop(paste0("simulate_process must be a single value"))
  }
  if(!all(is.logical(bias_correct_oe))) stop("bias_correct_oe must be TRUE or FALSE")
  if(!(length(bias_correct_oe) == 1)){
    stop(paste0("bias_correct_oe must be a single value"))
  }
  if(!all(is.logical(bias_correct_pe))) stop("bias_correct_pe must be TRUE or FALSE")
  if(!(length(bias_correct_pe) == 1)){
    stop(paste0("bias_correct_pe must be a single value"))
  }
  simulate_data = as.integer(simulate_data)
  simulate_process = as.integer(simulate_process)
  simulate_projection = as.integer(simulate_projection)
  bias_correct_pe = as.integer(bias_correct_pe)
  bias_correct_oe = as.integer(bias_correct_oe)
  
  data$simulate_data[] = simulate_data
  data$simulate_process[] = simulate_process
  data$simulate_period[2] = simulate_projection
  input$data = data
  return(input)

}