# devtools::load_all("~/work/wham/wham")
# devtools::install_github("timjmiller/wham", dependencies=TRUE, ref="devel")
library(wham)
library(tidyr)
library(dplyr)
library(here)
source(file.path(here(),"common_code", "set_Ecov.R"))
# create directory for analysis, e.g.


write.dir <- file.path(here(),"Ecov_study", "results")

if(!exists("write.dir")) write.dir = getwd()  #if we don't specify above, set as current wd
if(!dir.exists(write.dir)) dir.create(write.dir, recursive = T)  #if the write.dir directory doesn't exist, create it
setwd(write.dir)

#number of simulations for each scenario
#nsim = 1000
nsim = 10

#M sigmas and correlation parameters for each scenario
Ecov_sig   <- c(0.1, 0.3, 0.5) #units?
ar1_y      <- c(0,0.25,0.5,0.75)
Ecov_mean  <- 0
beta       <- c(0,0.3,0.7) #units?
Ecov_where <- c("recruit")
obs_sig    <- c(1e-5,0.1,0.25) #units?

df.mods <- expand.grid(Ecov_sig = Ecov_sig, Ecov_phi = ar1_y, Ecov_mean = Ecov_mean, beta = beta, Ecov_where = Ecov_where, obs_sig = obs_sig)

n.mods <- dim(df.mods)[1] #108 scenarios
df.mods$Model <- paste0("m_",1:n.mods)
df.mods <- df.mods %>% select(Model, everything()) # moves Model to first col
df.mods <- cbind(Recruitment = 2, NAA_re = "rec+1", df.mods) #recruit model not yet discussed by WG.
# look at model table: 28 scenarios
df.mods

#make a list of input components that prepare_wham_input can use to generate an input for fit_wham
#this will be the generic flatfish/groundfish life histor information
source(file.path(here(),"common_code", "make_basic_info.R"))
groundfish_info = make_basic_info()

selectivity = list(model = c(rep("logistic", groundfish_info$n_fleets),rep("logistic", groundfish_info$n_indices)),
    initial_pars = rep(list(c(5,1)), groundfish_info$n_fleets + groundfish_info$n_indices)) #fleet, index

M = list(initial_means = rep(0.2, length(groundfish_info$ages)))


sim_input = list()
# sim the 324 models!
for(m in 1:n.mods){
    #initial numbers at age
    NAA_re = list(N1_pars = exp(10)*exp(-(0:(length(groundfish_info$ages)-1))*M$initial_means[1]))
    NAA_re$recruit_pars = exp(10) #mean recruitment, if recruit model changes, need to change this
    NAA_re$recruit_model = df.mods$Recruitment[m] #random effects with a constant mean
    NAA_re$sigma = "rec+1"
    NAA_re$cor = "iid"
    input = prepare_wham_input(basic_info = groundfish_info, selectivity = selectivity, M = M, NAA_re = NAA_re)
    
    ecov = list(label = "AR1_ecov")
    ecov$mean = matrix(df.mods$Ecov_mean[m], input$data$n_years_model, 1)
    ecov$logsigma = matrix(log(df.mods$obs_sig[m]), input$data$n_years_model, 1)
    ecov$years = input$years
    ecov$lag = 0
    ecov$use_obs = matrix(TRUE, input$data$n_years_model, 1)
    ecov$link_model = "linear"  #GLB: this is always true, can it come outside the loop?
    ecov$process_model = "ar1"  
    ecov$process_mean_vals = df.mods$Ecov_mean[m]
    ecov$process_sig_vals = df.mods$Ecov_sig[m]
    ecov$process_cor_vals = df.mods$Ecov_phi[m]
    n_effects = 2+input$data$n_indices
    ecov$beta_vals = list(lapply(1:n_effects, function(x) matrix(df.mods$beta[m],1,input$data$n_ages)))
    ecov$where = df.mods$Ecov_where[m]
    ecov$how = 1 #not yet discussed by WG.

    input = set_ecov(input, ecov)

    print(paste0("m: ", m))
    om = fit_wham(input, do.fit = FALSE, MakeADFun.silent = TRUE)

    #simulate data from operating model
    set.seed(8675309) #use same seed for all operating models?
    #all RE and data are simulated
    sim_input[[m]] = lapply(1:nsim, function(x) {
        input_i = input
        sim = om$simulate(complete=TRUE)
        input_i$data = sim
        return(input_i)
    })

}
#save simulated data sets to google drive?
saveRDS(sim_input, file.path(write.dir, "om_sim_data_recruitment.RDS"))
sim_input <- readRDS(file.path(write.dir, "om_sim_data_recruitment.RDS"))


#This will now generate a list of inputs to estimate models that match the operating model.
#initial values are commented out to start at "generic" starting values for estimation.
#this could be modified to estimate models that do not match the operating model.
em_input = list()
# sim the 90 models!
for(m in 1:n.mods){

    #initial numbers at age
    NAA_re = list(N1_pars = exp(10)*exp(-(0:(length(groundfish_info$ages)-1))*M$initial_means[1]))
    NAA_re$recruit_pars = exp(10) #mean recruitment, if recruit model changes, need to change this
    NAA_re$recruit_model = df.mods$Recruitment[m] #random effects with a constant mean
    NAA_re$sigma = "rec+1"
    NAA_re$cor = "iid"
    input = prepare_wham_input(basic_info = groundfish_info, selectivity = selectivity, M = M, NAA_re = NAA_re)
    
    ecov = list(label = "AR1_ecov")
    ecov$mean = matrix(df.mods$Ecov_mean[m], input$data$n_years_model, 1)
    ecov$logsigma = matrix(log(df.mods$obs_sig[m]), input$data$n_years_model, 1)
    ecov$years = input$years
    ecov$lag = 0
    ecov$use_obs = matrix(TRUE, input$data$n_years_model, 1)
    ecov$link_model = "linear"
    ecov$process_model = "ar1"
    #ecov$process_mean_vals = df.mods$Ecov_mean[m]   #GLB: why are these commented out? For estimation?
    #ecov$process_sig_vals = df.mods$Ecov_sig[m]
    #ecov$process_cor_vals = df.mods$Ecov_phi[m]
    n_effects = 2+input$data$n_indices
    #ecov$beta_vals = list(lapply(1:n_effects, function(x) matrix(df.mods$beta[m],1,input$data$n_ages)))
    ecov$where = "recruit"
    ecov$how = 1 #not yet discussed by WG.
    
    input = set_ecov(input, ecov)
    
    #put in data simulated from operating model
    em_input[[m]] = lapply(1:nsim, function(x) {
        input_i = input
        input_i$data = sim_input[[m]][[x]]$data #put in simulated operating model data
        return(input_i)
    })
}
#save 
saveRDS(em_input, file.path(write.dir, "em_input_GLB_recruitment.RDS"))
em_input <- readRDS(file.path(write.dir, "em_input_GLB_recruitment.RDS"))



#NOT RUN YET
#em_input <- readRDS(file.path(write.dir, "em_input.RDS"))
em_fits = list()
#for(m in 1:n.mods){ #just do the first scenario until we are ready to do all of them.
ms <- c(8,20,32)
for(m in ms){
    em_fits[[m]] = lapply(1:nsim, function(x){
        cat(paste("model:",m, "fit:", x, "start \n"))
        out = fit_wham(em_input[[m]][[x]], do.osa = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE, save.sdrep = FALSE)
        cat(paste("model:",m, "fit:", x, "done \n"))
        return(out)
        })
}


