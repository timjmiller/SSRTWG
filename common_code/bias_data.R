# bias_data.R

# function to bias data sent to WHAM
# agg_catch_multiplier can be either a constant (applied to all years and fleets) or else a matrix with dimensions nyears by nfleets
bias_data <- function(data, multiply_agg_catch_flag=FALSE, agg_catch_multiplier=1.0){
  if(multiply_agg_catch_flag == TRUE){
    data$agg_catch <-  data$agg_catch * agg_catch_multiplier
    data$obsvec[data$keep_C+1] <-  data$agg_catch
  }
  return(data)
}

# helper function to define agg_catch_multiplier
create_agg_catch_multiplier <- function(input, multiplier=0.5, first_year=floor(input$data$n_years_model/2), last_year=length(input$years),  all_fleets_flag=TRUE){
  agg_catch_multiplier <- input$data$agg_catch-input$data$agg_catch+1.0 # to ensure correct dimensions
  if(all_fleets_flag == TRUE){
    agg_catch_multiplier[first_year:last_year, ] <- multiplier
  }else{
    # not defined yet
  }
  return(agg_catch_multiplier)
}
