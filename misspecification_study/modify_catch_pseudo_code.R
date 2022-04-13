# modify_catch_pseudo_code.R
# pseudo code for now to help understand what is needed to explore 
# missing (or extra) catch in SSRTWG explorations

###
# control code used in modify_catch directory
# all the usual libraries
# start with an input object 
# can be created by reading an ASAP file or generating de novo
input <- 

# create operating model
om <- fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE)

# modify_catch function parameters
# determine which years will have catch modified
# should this be in 1:nyears form or actual years?
years_catch_changed <- 
  
# supply multiplier for changed aggregate catch
multiplier_aggcatch <- 
  
# do we want the ability to deal with age proportions?
# for now leave age proportions alone
  
  
# for now assume no projections
# will need to adjust catch advice in MSE loops to account for change
  
#simulate data from operating model
set.seed(8675309) #use same seed for all operating models?
#all RE and data are simulated
sim_input[[m]] = lapply(1:nsim, function(x) {
  input_i = input
  sim = om$simulate(complete=TRUE)
  # here is where the modify catch function would be used
  sim = modify_catch(sim, other parameters defined below)
  input_i$data = sim
  return(input_i)
})


###
# modify_catch function
# input is a data object and parameters to control how the catch 
# is modified
# output is a data object with the catch modifications applied
# apply the adjustment to the catch
# note agg_catch has dimensions n_model_years by nfleets
data$agg_catch[years_catch_changed, ] <-  data$agg_catch[years_catch_changed, ] * multiplier_aggcatch
return(data)
