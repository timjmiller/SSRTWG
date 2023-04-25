#' Specify model and parameter configuration for natural mortality
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{wham::prepare_wham_input}})
#' @param M (optional) list specifying natural mortality options: model, random effects, initial values, and parameters to fix (see details)
#' 
#' \code{M} specifies estimation options for natural mortality and can overwrite M-at-age values specified in the ASAP data file.
#' If \code{NULL}, the M-at-age matrix from the ASAP data file is used (M fixed, not estimated). To estimate M-at-age
#' shared/mirrored among some but not all ages, modify \code{input$map$M_a} after calling \code{prepare_wham_input}
#' (see vignette for more details). \code{M} is a list with the following entries:
#'   \describe{
#'     \item{$model}{Natural mortality model options are:
#'                    \describe{
#'                      \item{"constant"}{(default) estimate a single mean M shared across all ages}
#'                      \item{"age-specific"}{estimate M_a independent for each age}
#'                      \item{"weight-at-age"}{specifies M as a function of weight-at-age, \eqn{M_y,a = exp(b0 + b1*log(W_y,a))}, as in
#'                        \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1095-8649.1996.tb00060.x}{Lorenzen (1996)} and
#'                        \href{https://www.nrcresearchpress.com/doi/10.1139/cjfas-2017-0035}{Miller & Hyun (2018)}.}
#'                    }
#'                  }
#'     \item{$re}{Time- and age-varying (random effects) on M. Note that random effects can only be estimated if
#'                mean M-at-age parameters are (\code{$est_ages} is not \code{NULL}).
#'                 \describe{
#'                   \item{"none"}{(default) M constant in time and across ages.}
#'                   \item{"iid"}{M varies by year and age, but uncorrelated.}
#'                   \item{"ar1_a"}{M correlated by age (AR1), constant in time.}
#'                   \item{"ar1_y"}{M correlated by year (AR1), constant all ages.}
#'                   \item{"2dar1"}{M correlated by year and age (2D AR1), as in \href{https://www.nrcresearchpress.com/doi/10.1139/cjfas-2015-0047}{Cadigan (2016)}.}
#'                 }
#'               }
#'     \item{$initial_means}{Initial/mean M-at-age vector, with length = number of ages (if \code{$model = "age-specific"})
#'                          or 1 (if \code{$model = "weight-at-age" or "constant"}). If \code{NULL}, initial mean M-at-age values are taken
#'                          from the first row of the MAA matrix in the ASAP data file.}
#'     \item{$est_ages}{Vector of ages to estimate M (others will be fixed at initial values). E.g. in a model with 6 ages,
#'                      \code{$est_ages = 5:6} will estimate M for the 5th and 6th ages, and fix M for ages 1-4. If \code{NULL},
#'                      M at all ages is fixed at \code{M$initial_means} (if not \code{NULL}) or row 1 of the MAA matrix from the ASAP file (if \code{M$initial_means = NULL}).}
#'     \item{$logb_prior}{(Only if \code{$model = "weight-at-age"}) TRUE or FALSE (default), should a N(0.305, 0.08) prior be
#'                        used on log_b? Based on Fig. 1 and Table 1 (marine fish) in \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1095-8649.1996.tb00060.x}{Lorenzen (1996)}.}
#'     \item{$b_size}{if \code{$model = "weight-at-age"}, initial value for random effect with prior defined by logb_prior (intended for simulating data).
#'     \item{$sigma_val}{Initial standard deviation value to use for the M deviations. Values are not used if \code{M$re} = "none". Otherwise, a single value.
#'     \item{$cor_vals}{Initial correlation values to use for the M deviations. If unspecified all initial values are 0. When \code{M$re} = 
#'                  \describe{
#'                    \item{"iid" or "none"}{values are not used.}
#'                    \item{"ar1_a" or "ar1_y"}{cor_vals must be a single value.}
#'                    \item{"2dar1"}{2 values must be specified. First is for "age", second is for "year".}
#'                  }
#'                }
#'   }

set_M = function(input, M)
{
  data = input$data
  par = input$par
  map = input$map

  #clear any map definitions that may exist. necessary because some configurations may not define map elements.
  map <- map[(!names(map) %in% c("log_b", "M_repars", "M_a","M_re"))]
  
  data$n_M_a = data$n_ages
  data$M_model = 2
  data$use_b_prior = 0
  data$M_re_model = 1 # default = no RE / 'none'
  data$M_est <- rep(0, data$n_M_a) # default = don't estimate M
  M_first_est = NA  
  M_re_ini = matrix(NA, data$n_years_model, data$n_M_a)
  if(is.null(input$asap3)) asap3 = NULL
  else {
    asap3 = input$asap3
    M_a_ini <- log(asap3$M[1,])
    M_re_ini[] <- matrix(log(asap3$M), data$n_years_model, data$n_M_a) - matrix(M_a_ini, data$n_years_model, data$n_M_a, byrow=T)
  }

  # natural mortality options, default = use values from ASAP file, no estimation
  if(is.null(M))
  {
    data$M_model = 2
    data$use_b_prior = 0
    data$M_re_model = 1 # default = no RE / 'none'
    data$M_est <- rep(0, data$n_M_a) # default = don't estimate M
    M_first_est = NA
    # Use ASAP file M values
    if(!is.null(asap3)){
    }
    else {
      M_a_ini <- log(rep(0.2, data$n_ages))
      M_re_ini[] <- 0
    }
  }
  if(!is.null(M)){
    if(!is.null(M$model)){ # M model options
      if(!(M$model %in% c("constant","age-specific","weight-at-age"))) stop("M$model must be either 'constant', 'age-specific', or 'weight-at-age'")
      if(!is.null(M$re)) if(M$model == "age-specific" & M$re == "ar1_a") stop("Cannot estimate age-specific mean M and AR1 deviations M_a.
If you want an AR1 process on M-at-age, set M$model = 'constant' and M$re = 'ar1_a'.")
      data$M_model <- match(M$model, c("constant","age-specific","weight-at-age"))
      if(M$model %in% c("constant","weight-at-age")){
        data$n_M_a = 1
        data$M_est = 0
        if(!is.null(asap3)) M_a_ini = log(asap3$M[1,1])
        else M_a_ini = log(0.2)
        if(!is.null(asap3)) {
          if(is.null(M$initial_means) & length(unique(asap3$M[1,])) > 1) warning("Constant or weight-at-age M specified (so only 1 mean M parameter),
but first row of MAA matrix has > 1 unique value.
Initializing M at age-1 MAA values. To avoid this warning
without changing ASAP file, specify M$initial_means.")
        }
        if(!is.null(M$logb_prior)){
          if(!is.logical(M$logb_prior)) stop("M$logb_prior must be either TRUE or FALSE")
          if(M$logb_prior) data$use_b_prior = 1
        }
      }
    }
    if(!is.null(M$re)){
      if(length(M$re) != 1) stop("Length(M$re) must be 1")
      if(!(M$re %in% c("none","iid","ar1_a","ar1_y","2dar1"))) stop("M$re must be one of the following: 'none','iid','ar1_a','ar1_y','2dar1'")
      data$M_re_model <- match(M$re, c("none","iid","ar1_a","ar1_y","2dar1"))
    }
    if(!is.null(M$initial_means)){
      if(length(M$initial_means) != data$n_M_a) stop("Length(M$initial_means) must be # ages (if age-specific M) or 1 (if constant or weight-at-age M)")
      M_a_ini <- log(M$initial_means)
      # overwrite ASAP file M values
      M_re_ini[] <- 0# if estimating mean M for any ages, initialize yearly deviations at 0
    }
    if(!is.null(M$est_ages)){
      if(!all(M$est_ages %in% 1:data$n_M_a)) stop("All M$est_ages must be in 1:n.ages (if age-specific M) or 1 (if constant or weight-at-age M)")
      data$M_est[M$est_ages] = 1 # turn on estimation for specified M-at-age
      M_first_est <- M$est_ages[1]
      M_re_ini[] <- 0# if estimating mean M for any ages, initialize yearly deviations at 0
    }
  }
  data$n_M_est <- sum(data$M_est)

  # natural mortality pars
  par$M_a <- M_a_ini # deviations by age
  par$M_re <- M_re_ini # deviations from mean M_a on log-scale, PARAMETER_ARRAY
  par$M_repars <- rep(0, 3)
  if(is.null(M$sigma_val)){
    par$M_repars[1] <- log(0.1) # start sigma at 0.1, rho at 0
  } else {
    if(length(M$sigma_val)==1) par$M_repar[1] = log(M$sigma_val)
    else stop("M$sigma_val can be only 1 value.")
  }
  if(is.null(M$cor_vals)) {
    #already set to 0 above
    #if(data$M_re_model == 3) par$M_repars[3] <- 0 # if ar1 over ages only, fix rho_y = 0
    #if(data$M_re_model == 4) par$M_repars[2] <- 0 # if ar1 over years only, fix rho_a = 0
    } else{
    #inv_trans_rho <- function(rho) 0.5 * log(rho+1) - log(1-rho) # 0.5 because needed transformation on cpp side is unusual.
    if(data$M_re_model == 3) par$M_repars[2] <- wham:::gen.logit(M$cor_vals, -1, 1, 2) #inv_trans_rho(M$cor_vals) # if ar1 over ages only, fix rho_y = 0
    if(data$M_re_model == 4) par$M_repars[3] <- wham:::gen.logit(M$cor_vals, -1, 1, 2) #inv_trans_rho(M$cor_vals) # if ar1 over years only, fix rho_a = 0
    if(data$M_re_model == 5) par$M_repars[2:3] <- wham:::gen.logit(M$cor_vals, -1, 1, 2) #inv_trans_rho(M$cor_vals) # if 2dar1 years and ages
  }
  # M_repars: sigma_M, rho_M_a, rho_M_y
  if(data$M_re_model == 1) tmp <- rep(NA,3) # no RE pars to estimate
  if(data$M_re_model == 2) tmp <- c(1,NA,NA) # estimate sigma
  if(data$M_re_model == 3) tmp <- c(1,2,NA) # ar1_a: estimate sigma, rho_a
  if(data$M_re_model == 4) tmp <- c(1,NA,2) # ar1_y: estimate sigma, rho_y
  if(data$M_re_model == 5) tmp <- 1:3 # 2dar1: estimate all
  map$M_repars = factor(tmp)

  # check if only 1 estimated mean M (e.g. because weight-at-age M or if all but 1 age is fixed), can't estimate rho_a
  # if(data$n_M_est < 2) par$M_repars[2] <- 0
  if(is.null(M$b_size)) par$log_b = log(0.305)
  else par$log_b = log(M$b_size)
  

  tmp <- par$M_a
  tmp[data$M_est==0] = NA
  ind.notNA <- which(!is.na(tmp))
  tmp[ind.notNA] <- 1:length(ind.notNA)
  map$M_a <- factor(tmp)
  if(data$M_model != 3) map$log_b = factor(rep(NA,length(par$log_b)))

  # M_re: "none","iid","ar1_a","ar1_y","2dar1"
  tmp <- par$M_re
  if(data$M_re_model == 1) tmp[] = NA # no RE (either estimate RE for all ages or none at all)
  if(data$M_re_model %in% c(2,5)){ # 2d ar1
    tmp[] = 1:(dim(tmp)[1]*dim(tmp)[2]) # all y,a estimated
  }
  if(data$M_re_model == 3){ # ar1_a (devs by age, constant by year)
    for(i in 1:dim(tmp)[2]) tmp[,i] = i
  }
  if(data$M_re_model == 4){ # ar1_y (devs by year, constant by age)
    for(i in 1:dim(tmp)[1]) tmp[i,] = i
  }
  map$M_re <- factor(tmp)


  input$data = data
  input$par = par
  input$map = map
  
  #may need to update these 
	# projection data will always be modified by 'prepare_projection'
	input = wham:::set_proj(input, proj.opts = NULL) #proj options are used later after model fit, right?

	#set any parameters as random effects
	input$random = NULL
	input = wham:::set_random(input)
  return(input)

}