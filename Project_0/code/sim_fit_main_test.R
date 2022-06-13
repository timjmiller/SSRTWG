library(here)
library(wham)
#sim_fit_main.R


#naa_oms

om_inputs = readRDS(file.path(here(),"Project_0","inputs", "NAA_om_inputs.RDS"))
em_inputs = readRDS(file.path(here(),"Project_0","inputs", "em_inputs.RDS"))
script.full.path = file.path(here(),"Project_0", "code", "naa_om_sim_fit_script.R")

#this file must be created on the server
fname = "naa_om_sim_fit_test_commands.txt"
write("#commands to run on server.", file = fname, append = FALSE)
write("#create simulated data from an operating model and fit with an estimating model.", file = fname, append = TRUE)
for(i in 1:length(om_inputs)) for(j in 1:length(em_inputs)){ 
  write(paste0("Rscript --vanilla " script.full.path, " " , i, " ",  j, " 1 &"), file = fname, append = TRUE)
}

df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.oms.RDS"))
#######################################################
#need to have matching assumptions about CVs for catch and indices, too
obs_names = c("agg_catch","agg_catch_sigma", "agg_indices", "agg_index_sigma", "catch_paa", "index_paa", 
  "Ecov_obs", "obs", "obsvec")
#######################################################

#######################################################
#Do we want to use the same (e.g. 1000) seeds for everything?
#set.seed(1234)
#seeds = sample(x = (-1e9):(1e9), size = 1000, replace = FALSE)
#saveRDS(seeds, file.path(here(), "Project_0", "inputs","seeds.RDS"))
seeds = readRDS(file.path(here(), "Project_0", "inputs","seeds.RDS"))
#can change sims to 101:200, etc. to do subsets of simulations and stop < 1000 if time is an issue.
sims = 1:100
sims = 1
#######################################################

naa_om_results = list()
naa_om_truth = list()
n_oms = length(om_inputs)
#n_oms = 1 #testing
#for(this_om in 1:n_oms){
for(this_om in 2:n_oms){
  naa_om_results[[this_om]] = list()
  naa_om_truth[[this_om]] = list()
  for(this_sim in sims){
    print(paste0("OM: ", this_om, " Sim: ", this_sim))
    # Set seed
    om = fit_wham(om_inputs[[this_om]], do.fit = FALSE, MakeADFun.silent = TRUE)
    set.seed(seeds[this_sim])
    sim_data = om$simulate(complete=TRUE)
    naa_om_truth[[this_om]][[this_sim]] = sim_data
    #save the version for reproducibility
    naa_om_truth[[this_om]][[this_sim]]$wham_version = om$wham_version
    naa_om_results[[this_om]][[this_sim]] = list()
    for(this_em in 1:length(em_inputs)){
      library(wham)
      print(paste0("OM: ", this_om, " Sim: ", this_sim, " EM: ", this_em))
      res = list()      
      EM_input <- em_inputs[[this_em]] # Read in the EM 
      #put simulated data into the em input
      EM_input$data[obs_names] = sim_data[obs_names]

      #do fit withouth sdreport first
      fit <- tryCatch(fit_wham(EM_input, do.sdrep=F, do.osa=F, do.retro=T, do.proj=F, MakeADFun.silent=FALSE),
        error = function(e) conditionMessage(e))
				
      # Deal with issues fitting EM to non-matching OM data
      # empty elements below can be used to summarize convergence information
      if(!'err' %in% names(fit) & class(fit) != "character"){
        res$wham_version <- fit$wham_version
        res$TMB_version <- fit$TMB_version
        res$opt <- fit$opt
        res$reps <- fit$rep
        res$mohns_rho <- tryCatch(mohns_rho(fit),
          error = function(e) conditionMessage(e))
        fit$sdrep <- tryCatch(TMB::sdreport(fit), # no bc
                error = function(e) conditionMessage(e))
        if(class(fit$sdrep) == "sdreport"){ 
          res$sdreps <- list(
            "Estimate_par" = as.list(fit$sdrep, what = "Est"),
            "SE_par" = as.list(fit$sdrep, what = "Std"),
            "Estimate_rep" = as.list(fit$sdrep, what = "Est", report = TRUE),
            "SE_rep" = as.list(fit$sdrep, what = "Std", report = TRUE))
        }
      }
      naa_om_results[[this_om]][[this_sim]][[this_em]] = res
      #saveRDS(naa_om_results[[this_om]][[this_sim]][[this_em]], file = file.path(here(),"results", 
      #  paste0("naa_om_results_om", this_om, "_em", this_em, "_sim", this_sim, "_test.RDS"))) 
    }
  }
}



saveRDS(naa_om_results, file = file.path(here(),"results", "naa_om_results_test.RDS"))
