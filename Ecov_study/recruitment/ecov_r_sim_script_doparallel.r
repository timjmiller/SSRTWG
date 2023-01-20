# devtools::load_all("~/work/wham/wham")
# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="devel")
library(wham)
library(tidyr)
library(dplyr)
library(here)
library(doParallel)
x <- detectCores()      
registerDoParallel(x-1) #leave one core for other tasks
#writeLines(paste(x), "cores_detected.txt") #print how many cores were used   

source(file.path(here(),"common_code", "make_basic_info.R"))
source(file.path(here(),"common_code", "set_ecov.R")) #load set_ecov.r function
source(file.path(here(),"common_code", "get_FMSY.R")) #load set_ecov.r function
source(file.path(here(), "Project_0", "code", "make_om.R"))

write.dir <- file.path(here(),"Ecov_study", "recruitment", "results") # create directory for analysis

if(!exists("write.dir")) write.dir = getwd()  #if we don't specify above, set as current wd
if(!dir.exists(write.dir)) dir.create(write.dir, recursive = T)  #if the write.dir directory doesn't exist, create it
setwd(write.dir)

nsim = 25 #number of simulations for each scenario

################################################################
##--FUNCTIONS--#################################################
################################################################
#Function to modify inputs
mod_input <- function(input,NAA_re,ecov,df.mods,m){
  NAA_re$sigma_vals[2] = df.mods$NAA_sig[m]  #only makes sense with 'rec+1'?
  
  #modify ecov  
  ecov$logsigma          <- matrix(log(df.mods$obs_sig[m]), length(gf_info$years), 1)
  ecov$process_mean_vals <- df.mods$Ecov_mean[m]
  ecov$process_sig_vals  <- df.mods$Ecov_sig[m]
  ecov$process_cor_vals  <- df.mods$Ecov_phi[m]
  ecov$beta_vals         <- list(lapply(1:ecov$n_effects, function(x) matrix(df.mods$beta[m],1,input0$data$n_ages)))
  return(set_ecov(input,ecov))  
}

get_FXSPR <- function(input,NAA_re,brp_year=1){
  temp <- fit_wham(input0,do.fit=FALSE,MakeADFun.silent=FALSE)
  return(exp(temp$rep$log_FXSPR[brp_year]))
}
################################################################
##--EXPERIMENTAL FACTORS--######################################
################################################################
Ecov_sig   <- c(0.1,0.5)     
ar1_y      <- c(0,0.95)
beta       <- c(0.3,1.0) 
obs_sig    <- c(1e-5,0.25) 
NAA_sig    <- c(1e-5,0.25)
#R_sig      <- 
#F_hist     <- c("Fmsy","H")
#NEED: 1) obs error on NAA; 2) process error on NAA; F history

df.mods       <- expand.grid(Ecov_sig=Ecov_sig, 
                             Ecov_phi = ar1_y, 
                             beta = beta, 
                             obs_sig = obs_sig,
                             NAA_sig = NAA_sig)
n.mods        <- dim(df.mods)[1]
df.mods$Model <- paste0("m_",1:n.mods)
df.mods       <- df.mods %>% select(Model, everything()) # moves Model to first col
df.mods       <- cbind(Recruitment=3, NAA_re="rec+1", df.mods) #recruit model codes: 1=RW, 2=random about mean, 3=BH, 4=Ricker; 'rec+1' is full state space model
df.mods$nsim  <- rep(nsim,nrow(df.mods)) #number of simulations per 

saveRDS(df.mods,file.path(write.dir, "om_sim_inputs_GLB_recruitment_doparallel.RDS"))

##################################################################
##--BASIC SETTINGS--#####################################F#########
##################################################################
gf_info <- make_basic_info()

#Define selectivity; not changing
selectivity <- list(model       =c(rep("logistic", gf_info$n_fleets),rep("logistic", gf_info$n_indices)),
                    initial_pars=rep(list(c(5,1)), gf_info$n_fleets + gf_info$n_indices)) #fleet, index

#Define M; not changing
M <- list(initial_means=rep(0.2, length(gf_info$ages)))

#Define NAA re to be modified
NAA_re <- list(
  N1_pars = exp(10)*exp(-(0:(length(gf_info$ages)-1))*M$initial_means[1]),
  sigma = "rec", #random about mean
  cor="iid", #random effects are independent
  #use_steepness = 1, #GLB: also don't need for BH?
  recruit_model = 2, #random effects with a constant mean
  recruit_pars <- exp(10)
  #recruit_model = 3, #B-H
  #recruit_pars = c(0.75,exp(10)) #GLB: this is for BH?
)

#Define ecov to be modified within loop
ecov <- list(label = "AR1_ecov",
  years = gf_info$years,
  lag = 0,
  use_obs = matrix(TRUE, length(gf_info$years), 1),
  link_model = "linear", 
  process_model = "ar1",
  mean = matrix(0, length(gf_info$years), 1),
  process_mean_values = 0,
  #n_effects = 2+input$data$n_indices, #GLB: what does this do?
  n_effects = 3,
  where = "recruit",
  how = 1) 

input0 <- prepare_wham_input(basic_info = gf_info, selectivity = selectivity, M = M, NAA_re = NAA_re, age_comp="logistic-normal-miss0")

################################################################
##--SIMULATIONS--###############################################
################################################################
sim_input = em_input <-list()
for(m in 1:n.mods){
  input = mod_input(input0,NAA_re,ecov,df.mods,m)
  
  ##--SIMULATE WITH WHAM--######################
  print(paste0("m: ", m))
  om = fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE) 
  #set.seed(8675309) #use same seed for all operating models?
  sim_input[[m]] = lapply(1:nsim, function(x) {      #all RE and data are simulated
    input_i = input
    sim = om$simulate(complete=TRUE)
    input_i$data = sim
    return(input_i)
  })
  em_input[[m]] = lapply(1:nsim, function(x) {     #put in data simulated from operating model
    input_i = input
    input_i$data = sim_input[[m]][[x]]$data #put in simulated operating model data
    return(input_i)
  })
}
saveRDS(sim_input, file.path(write.dir, "om_sim_data_GLB_recruitment_doparallel.RDS"))
saveRDS(em_input,  file.path(write.dir, "em_input_GLB_recruitment_doparallel.RDS"))

########################################################
##--FIT MODELS--########################################
########################################################
em_fits <- foreach(m=1:n.mods) %dopar% {
  lapply(1:nsim, function(x){
    cat(paste("model:",m, "fit:", x, "start \n"))
    out = fit_wham(em_input[[m]][[x]], do.osa = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE, save.sdrep = FALSE)
    out$rep <- out$report()
    cat(paste("model:",m, "fit:", x, "done \n"))
    return(out)
  })
}
saveRDS(em_fits, file.path(write.dir, "em_fits_GLB_recruitment_doparallel.RDS"))


