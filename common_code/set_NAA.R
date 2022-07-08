#' Specify model and parameter configuration for numbers at age
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{prepare_wham_input}})
#' @param NAA_re (optional) list specifying options for numbers-at-age random effects, initial parameter values, and recruitment model (see details)
#' 
#' If \code{NAA_re = NULL}, a traditional statistical catch-at-age model is fit (NAA = pred_NAA for all ages, deterministic). Otherwise,
#' \code{NAA_re} specifies numbers-at-age configuration. It is a list with the following possible entries:
#'   \describe{
#'     \item{$sigma}{Which ages allow deviations from pred_NAA? Common options are specified with the strings:
#'                    \describe{
#'                      \item{"rec"}{Random effects on recruitment (deviations), all other ages deterministic}
#'                      \item{"rec+1"}{"Full state space" model with 2 estimated \code{sigma_a}, one for recruitment and one shared among other ages}
#'                    }
#'                   Alternatively, you can specify a more complex structure by entering a vector with length = n.ages, where each entry points to the
#'                   NAA_sigma to use for that age. E.g. c(1,2,2,3,3,3) will estimate 3 \code{sigma_a}, with recruitment (age-1) deviations having their
#'                   own \code{sigma_R}, ages 2-3 sharing \code{sigma_2}, and ages 4-6 sharing \code{sigma_3}.
#'                  }
#'     \item{$sigma_vals}{Initial standard deviation values to use for the NAA deviations. Values are not used if recruit_model = 1 and NAA_re$sigma is
#'                  not specifed. Otherwise when \code{NAA_re$sigma} =
#'                  \describe{
#'                    \item{"rec"}{must be a single value.}
#'                    \item{"rec+1"}{2 values must be specified. First is for the first age class (recruits), second is for all other ages.}
#'                    \item{vector of values (length = number of age classes)}{either 1 value or the number of values is equal to the number of unique values provided to \code{NAA_re$sigma}.}
#'                  }
#'                }
#'     \item{$cor}{Correlation structure for the NAA deviations. Options are:
#'                  \describe{
#'                    \item{"iid"}{NAA deviations vary by year and age, but uncorrelated.}
#'                    \item{"ar1_a"}{NAA deviations correlated by age (AR1).}
#'                    \item{"ar1_y"}{NAA deviations correlated by year (AR1).}
#'                    \item{"2dar1"}{NAA deviations correlated by year and age (2D AR1).}
#'                  }
#'                }
#'     \item{$cor_vals}{Initial correlation values to use for the NAA deviations. If unspecified all initial values are 0. When \code{NAA_re$cor} = 
#'                  \describe{
#'                    \item{"iid"}{values are not used.}
#'                    \item{"ar1_a" or "ar1_y"}{cor_vals must be a single value.}
#'                    \item{"2dar1"}{2 values must be specified. First is for "age", second is for "year".}
#'                  }
#'                }
#'     \item{$N1_model}{Integer determining which way to model the initial numbers at age:
#'       \describe{
#'          \item{0}{(default) age-specific fixed effects parameters}
#'          \item{1}{2 fixed effects parameters: an initial recruitment and an instantaneous fishing mortality rate to generate an equilibruim abundance at age.}
#'       }
#'     }
#'     \item{$N1_pars}{if N1_model = 0, then these would be the initial values to use for abundance at age in the first year. If N1_model = 1, This would be the
#'        initial numbers in the first age class and the equilibrium fishing mortality rate generating the rest of the numbers at age in the first year.
#'     }
#'     \item{$recruit_model}{Integer determining how to model recruitment. Overrides \code{recruit_model} argument to \code{prepare_wham_input}. Must make sure \code{NAA_re$sigma}, \code{NAA_re$cor}
#'        and \code{ecov} are properly specified.
#'       \describe{
#'           \item{1}{SCAA, estimating all recruitements as fixed effects or a random walk if NAA_re$sigma specified}
#'           \item{2}{estimating a mean recruitment with yearly recruitements as random effects}
#'           \item{3}{Beverton-Holt stock-recruitment with yearly recruitements as random effects}
#'           \item{4}{Ricker stock-recruitment with yearly recruitements as random effects}
#'       }
#'     }
#'     \item{$use_steepness}{T/F determining whether to use a steepness parameterization for a stock-recruit relationship. Only used if recruit_model>2}.
#'     \item{$recruit_pars}{vector of initial parameters for recruitment model. If use_steepness=F, parameters are "alpha" and "beta"
#'        otherwise they are steepness and R0.
#'     }
#'   }

set_NAA = function(input, NAA_re=NULL)
{

  data = input$data
  par = input$par
  map = input$map
  
  #clear any map definitions that may exist. necessary because some configurations may not define map elements.
  map <- map[(!names(map) %in% c("mean_rec_pars", "log_NAA_sigma", "trans_NAA_rho","logR_proj", "log_NAA"))]
  
  if(is.null(input$asap3)) asap3 = NULL
  else asap3 = input$asap3

  #set up initial NAA
  if(is.null(NAA_re$N1_model)) {
    data$N1_model = 0 #0: just age-specific numbers at age
  } else {
    data$N1_model = NAA_re$N1_model 
  }
  if(is.null(NAA_re$N1_pars)){
    if(data$N1_model == 0){
      if(!is.null(asap3)) par$log_N1_pars = log(asap3$N1_ini) # use N1_ini values from asap3 file
      else par$log_N1_pars = log(exp(10)*exp(-(0:(data$n_ages-1))*0.2))
    }
    if(data$N1_model == 1) par$log_N1_pars = c(10,log(0.1)) # allowed in wham.cpp but no option to set here (must be hard-coded after calling prepare_wham_input)
  }
  else par$log_N1_pars = log(NAA_re$N1_pars)

  if(!(data$N1_model %in% 0:1)) stop("NAA_re$N1_model can only be 0 or 1 currently")

  # NAA_re options
  if(is.null(NAA_re$sigma)){ # default = SCAA
    data$n_NAA_sigma <- 0
    data$NAA_sigma_pointers <- rep(1,data$n_ages)
  } else {
    if(NAA_re$sigma == "rec"){
      data$n_NAA_sigma <- 1
      data$NAA_sigma_pointers <- rep(1,data$n_ages)
    } else {
      if(NAA_re$sigma == "rec+1"){ # default state-space model with two NAA_sigma (one for recruits, one for ages > 1)
        data$n_NAA_sigma <- 2
        data$NAA_sigma_pointers <- c(1,rep(2,data$n_ages-1))
      } else {
        if(length(NAA_re$sigma) != data$n_ages) stop("NAA_re$sigma must either be 'rec' (random effects on recruitment only), 
'rec+1' (random effects on all NAA with ages > 1 sharing sigma_a,
or a vector with length == n.ages specifying which sigma_a to use for each age.")
        if(length(NAA_re$sigma) == data$n_ages){
          if(is.na(unique(NAA_re$sigma))){
            data$n_NAA_sigma <- 0
            data$NAA_sigma_pointers <- rep(1,data$n_ages)            
          } else {
            data$n_NAA_sigma <- max(unique(NAA_re$sigma), na.rm=T)
            data$NAA_sigma_pointers <- NAA_re$sigma
          }
        }
      }
    }
  }
  if(data$recruit_model > 2 & data$n_NAA_sigma == 0) warning("SCAA model specified, yearly recruitment deviations estimated as fixed effects. Stock-recruit function also specified. WHAM will fit the SCAA model but without estimating a stock-recruit function.
This message will not appear if you set recruit_model = 2 (random about mean).")

  # NAA_re pars
  par$log_NAA_sigma = numeric(length = ifelse(data$n_NAA_sigma==0,1,data$n_NAA_sigma))
  if(is.null(NAA_re$sigma_vals)) par$log_NAA_sigma[] <- 0
  else {
    if(length(NAA_re$sigma_vals) %in% c(1,data$n_NAA_sigma)) par$log_NAA_sigma[] <- log(NAA_re$sigma_vals)
    else stop(paste0("length of NAA_re$sigma_vals is not 1 or ", data$n_NAA_sigma, ". It must be consistent with other elements of NAA_re."))
  }
  if(data$n_NAA_sigma == 0) map$log_NAA_sigma <- factor(NA)

  if(!is.null(NAA_re$cor)){
    if(!NAA_re$cor %in% c("iid","ar1_a","ar1_y","2dar1")) stop("NAA_re$cor must be one of 'iid','ar1_a','ar1_y','2dar1'")
  } else {
    NAA_re$cor <- 'iid'
  }
  inv_trans_rho <- function(rho, s = 2) (log(rho+1) - log(1-rho))/s # 0.5 because needed transformation on cpp side is unusual.
  if(is.null(NAA_re$cor_vals)) par$trans_NAA_rho <- c(0,0)
  else {
    if(length(NAA_re$cor_vals) == 2) par$trans_NAA_rho <- inv_trans_rho(NAA_re$cor_vals)
    if(length(NAA_re$cor_vals) == 1) {
      if(NAA_re$cor == "ar1_a") {
        par$trans_NAA_rho <- c(inv_trans_rho(NAA_re$cor_vals), 0)
      }
      if(NAA_re$cor == "ar1_y") {
        par$trans_NAA_rho <- c(0, inv_trans_rho(NAA_re$cor_vals))
      }
    }
    if(!length(NAA_re$cor_vals) %in% 1:2) stop(paste0("length of NAA_re$cor_vals is not consistent with other elements of NAA_re$cor."))
  }

  #NAA_rho map
  if(!is.null(NAA_re$cor)) {
    if(NAA_re$cor == "iid") map$trans_NAA_rho <- factor(c(NA,NA))
    if(NAA_re$cor == "ar1_a") map$trans_NAA_rho <- factor(c(1,NA))
    if(NAA_re$cor == "ar1_y") map$trans_NAA_rho <- factor(c(NA,1))
  } else map$trans_NAA_rho <- factor(c(NA,NA))

  par$log_NAA = matrix(10, data$n_years_model-1, data$n_ages)
  par$logR_proj <- 0 # will be set by prepare_projection if SCAA
  map$logR_proj <- factor(NA)
  
  #NAA_re map
  tmp <- par$log_NAA
  if(data$n_NAA_sigma < 2) tmp[,-1] <- NA # always estimate Rec devs (col 1), whether random effect or not
  ind.notNA <- which(!is.na(tmp))
  tmp[ind.notNA] <- 1:length(ind.notNA)
  map$log_NAA = factor(tmp)

  #set up recruitment
  if(!is.null(NAA_re$recruit_model)) {
    data$recruit_model = NAA_re$recruit_model #overrides recruit_model argument to wham::prepare_wham_input
    if(data$recruit_model > 1 & length(NAA_re$sigma) == 0) stop("NAA_re$recruit_model > 1 has been specified, but NAA_re$sigma must either be 'rec' (random effects on recruitment only), 
      'rec+1' (random effects on all NAA with ages > 1 sharing sigma_a,
      or a vector with length == n.ages specifying which sigma_a to use for each age.")
  }
  par$mean_rec_pars = numeric(c(0,1,2,2)[data$recruit_model])

  if(!is.null(NAA_re$use_steepness)) data$use_steepness = sum(NAA_re$use_steepness)
  else data$use_steepness = 0 #use regular SR parameterization by default, steepness still can be estimated as derived par.
  
  if(!is.null(NAA_re$recruit_pars)){
    if(data$recruit_model == 2) par$mean_rec_pars[] = log(NAA_re$recruit_pars[1])
    if(data$recruit_model %in% 3:4){
      if(data$use_steepness == 1){
        if(data$recruit_model == 3) par$mean_rec_pars[1] = log(NAA_re$recruit_pars[1] - 0.2) - log(1-NAA_re$recruit_pars[1])
        if(data$recruit_model == 4) par$mean_rec_pars[1] = log(NAA_re$recruit_pars[1] - 0.2)
        par$mean_rec_pars[2] = log(NAA_re$recruit_pars[2])
      }
      else par$mean_rec_pars[] = log(NAA_re$recruit_pars)
    }
  }
  else{
    par$mean_rec_pars[] = 0	
    if(data$recruit_model==2) {
  		if(!is.null(asap3)) par$mean_rec_pars = log(asap3$N1_ini[1]) # initialize R0 at initial age-1
  		else par$mean_rec_pars = 10
  	}  
    if(data$recruit_model==4) par$mean_rec_pars[2] = -10
    if(data$recruit_model == 1) map$mean_rec_pars = factor(rep(NA, length(par$mean_rec_pars)))
    if(data$n_NAA_sigma == 0) map$mean_rec_pars = factor(rep(NA, length(par$mean_rec_pars))) #SCAA with recruitments as fixed effects
  }

  input$data = data
  input$par = par
  input$map = map

  #may need to update these 
	# projection data will always be modified by 'prepare_projection'
	#input = wham:::set_proj(input, proj.opts = NULL) #proj options are used later after model fit, right?

	#set any parameters as random effects
	input$random = NULL
	input = wham:::set_random(input)

  return(input)
}