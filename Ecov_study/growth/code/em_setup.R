# remotes::install_github(repo = 'GiancarloMCorrea/wham', ref='growth', INSTALL_opts = c("--no-docs", "--no-multiarch", "--no-demo"))
library(wham)
library(tidyr)
library(dplyr)
library(here)

source(file.path(here(), "common_code", "make_basic_info.R"))
source(file.path(here(), "common_code", "set_NAA.R"))
source(file.path(here(), "common_code", "set_M.R"))
source(file.path(here(), "common_code", "set_q.R"))
source(file.path(here(), "common_code", "set_ecov.R"))
source(file.path(here(), "common_code", "set_selectivity.R"))
source(file.path(here(), "common_code", "set_simulation_options.R"))
source(file.path(here(), "common_code", "get_FMSY.R"))
source(file.path(here(), "Ecov_study", "growth", "code", "make_om.R"))
source(file.path(here(), "Ecov_study", "growth", "code", "sim_management.R"))
## verify_version()

# Create data.frame with EM configurations:
em_method = c('EWAA', rep(c('WAA', 'WAA', 'WAA'), times = 2))
em_re_method = c(NA, rep(c('iid', '2dar1', '3dgmrf'), times = 2))
em_est_fixed = c(NA, rep(FALSE, times = 3), rep(TRUE, times = 3))
df.ems = data.frame(em_method = em_method, em_re_method = em_re_method, em_est_fixed = em_est_fixed)
saveRDS(df.ems, file.path(here(),"Ecov_study", "growth", "inputs", "df.ems.RDS"))

# Make basic inputs (age comps)
gf_info = make_basic_info()

# Selectivity configuration (age-based for now)
gf_selectivity = list(
  model = c("double-normal", "logistic"),
  initial_pars = list(c(4,-1,1.5,1.4,0,0.4), c(1,0.5)))

# Natural mortality
gf_M = list(initial_means = 0.2, model = "constant")

# NAA configuration
gf_NAA_re = list(
  N1_pars = c(1e+05, 0),
  sigma = "rec", #random about mean
  cor = "iid", #random effects are independent
  recruit_model = 2,
  N1_model = 1
)

# Ecov configuration (for OM):
gf_ecov <- list(
  label = "Ecov",
  process_model = "ar1",
  logsigma = cbind(rep(log(0.1), length(gf_info$year))),
  lag = 0,
  mean = cbind(rep(0, length(gf_info$years))),
  year = gf_info$years,
  use_obs = cbind(rep(1, length(gf_info$years))),
  where = "none", ## updated to growth below
  where_subindex = 3, # on L1, growth specific input
  how=0)

# OM parameters:
Linf <- 85
k <- 0.3
t0 <- 0
a_LW <- exp(-12.1)
b_LW <- 3.2
L_a <- Linf*(1-exp(-k*(1:10 - t0)))
W_a <- a_LW*L_a^b_LW
CV_a <- .1

# gf_growth <- list(model='vB_classic', re = rep('none', times = 3), 
#                   init_vals=c(k, Linf, L[1]),
#                   est_pars=1:3, SD_vals=c(CV_a*L[1], CV_a*L[10]),
#                   SD_est=1:2) 
# gf_LW <- list(init_vals=c(a_LW, b_LW))
# gf_LAA = list(LAA_vals = L, est_pars = 1:10, re = 'none', 
#               SD_vals = c(CV_a*L[1], CV_a*L[10]), SD_est=1:2) # fixing SD1
# WAA configuration:
gf_WAA = list(WAA_vals = W_a,
              re = c('none'))


#make inputs for estimating model (smaller objects to save, can overwrinte data elements with simulated data)
em_inputs = list()
## df.ems <- df.ems[8,]
for(i in 1:NROW(df.ems)){
  print(paste0("EM config ", i))
  NAA_re_i = gf_NAA_re
  M_i = gf_M
  #growth_i <- gf_growth
  #LAA_i = gf_LAA
  WAA_i = gf_WAA
  ecov_i = gf_ecov
  selectivity = gf_selectivity

  ## estimate fixed effect growth pars?
  # if(!df.ems$growth_est[i]) {
  #   growth_i$est_pars <- NULL
  #   growth_i$SD_est <- NULL
  #   LAA_i$est_pars <- NULL
  #   LAA_i$SD_est <- NULL
  # }
  # Change Ecov information:
  # if(df.ems$Ecov_est[i] & df.ems$growth_method[i] == 'growth'){
  #   ecov_i$how = 1
  #   ecov_i$where = df.ems$growth_method[i]
  # }
  # if(df.ems$Ecov_est[i] & df.ems$growth_method[i] == 'LAA'){
  #   ecov_i$how = 1
  #   ecov_i$where = df.ems$growth_method[i]
  # }

  # Change growth information:
  # if(df.ems$growth_method[i] == 'growth') { 
  #   LAA_i = NULL
  #   if(!is.na(df.ems$growth_re_config[i])) { 
  #     growth_i$re[3] = "ar1_y" # always on L1
  #   }
  # }
  # Change LAA information
  # if(df.ems$growth_method[i] == 'LAA') { 
  #   growth_i = NULL
  #   if(!is.na(df.ems$growth_re_config[i])) { 
  #     LAA_i$re = '2dar1'
  #   }
  # }
  # Change WAA information -------------------------------
  if(df.ems$em_method[i] == 'WAA') { # nonparametric approach
    if(!is.na(df.ems$em_re_method[i])) { # random effects structure
      WAA_i$re = df.ems$em_re_method[i]
    }
    if(df.ems$em_est_fixed[i]) { # estimate fixed effects?
      WAA_i$est_pars = 1:10
    }
  }
  if(df.ems$em_method[i] == 'EWAA') { # EWAA approach
    WAA_i = NULL
  }
  # Change growth and LAA information (semiparametric)
  # if(df.ems$growth_method[i] == 'semiparametric') { 
  #   if(!is.na(df.ems$growth_re_config[i])) { 
  #     LAA_i$re = '2dar1'
  #   }
  # }

  basic_info <- make_basic_info()
  # Add length information:
  ny <- length(basic_info$years)
  basic_info$lengths <- seq(1, 120, by=2) # length bins
  nlbins <- length(basic_info$lengths)
  basic_info$n_lengths <- nlbins
  # Age comps:
  basic_info$use_catch_paa <- matrix(1, ncol = basic_info$n_fleets, nrow = ny)
  basic_info$use_index_paa <- matrix(1, ncol = basic_info$n_indices, nrow = ny)

  em_inputs[[i]] <-
    prepare_wham_input(basic_info = basic_info,
                       selectivity = selectivity, NAA_re = NAA_re_i, M= M_i,
                       WAA = WAA_i, 
                       age_comp = "logistic-normal-miss0",
                       len_comp = 'multinomial')

  # if(!df.ems$Ecov_est[i]){ # if we do not estimate Ecov:
  #   em_inputs[[i]]$random = NULL # will turn off Ecov as well
  #   # is this still required?:
  #   em_inputs[[i]]$map$Ecov_re <- factor(NA*em_inputs[[i]]$par$Ecov_re)
  #   em_inputs[[i]]$par$Ecov_re <- 0*em_inputs[[i]]$par$Ecov_re
  #   em_inputs[[i]]$map$Ecov_process_pars <- factor(NA*em_inputs[[i]]$par$Ecov_process_pars)
  #   if(!is.na(df.ems$growth_re_config[i])) { # RE on growth or LAA
  #     em_inputs[[i]]$random <- df.ems$growth_re_config[i]
  #   }
  # } else { # if we estimate Ecov
  #   em_inputs[[i]]$random <- 'Ecov_re'
  # }

  #turn off bias correction
  em_inputs[[i]] = set_simulation_options(em_inputs[[i]], simulate_data = TRUE, simulate_process = TRUE, simulate_projection = FALSE,
    bias_correct_pe = FALSE, bias_correct_oe = FALSE)
  #change Neff so scalar doesn't affect L-N SD
  em_inputs[[i]]$data$catch_Neff[] = 1
  ## not sure why this has to be done (me either, Giancarlo)
  em_inputs[[i]]$data$index_Neff[] <- 1
  em_inputs[[i]]$data$FXSPR_init[] = 0.3
  em_inputs[[i]]$data$FMSY_init[] = 0.3
  em_inputs[[i]]$map$log_NAA_sigma <- factor(NA*em_inputs[[i]]$par$log_NAA_sigma)
}

saveRDS(em_inputs, file.path(here(),"Ecov_study", "growth", "inputs", "em_inputs.RDS"))
