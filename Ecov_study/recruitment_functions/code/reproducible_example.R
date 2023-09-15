library(here)

source(file.path(here(),'Ecov_study','recruitment_functions','code','om_setup.r'))
source(file.path(here(),'Ecov_study','recruitment_functions','code','em_setup.r'))

obs_names <- c("agg_catch","agg_catch_sigma", "agg_indices", "agg_index_sigma", "catch_paa", "index_paa", 
               "Ecov_obs", "obs", "obsvec")

#######################################################
# Example 1; works
#######################################################
omj <- 1; emk <- 2
df.oms[omj,] #ecov$how=1
df.ems[emk,] #ecov$how=1

om                       <- fit_wham(om_inputs[[omj]], do.fit = FALSE, MakeADFun.silent = TRUE)
sim_data                 <- om$simulate(complete=TRUE)
EM_input                 <- em_inputs[[emk]] 
EM_input$data[obs_names] <- sim_data[obs_names]

fit <- fit_wham(EM_input, do.sdrep=F, do.osa=F, do.retro=F, do.proj=F, MakeADFun.silent=TRUE)

#######################################################
# Example 2; this causes R to crash due to segfault error
#######################################################
omj <- 1; emk <- 1  #both om and em have ecov$how=1
df.oms[omj,] #ecov_how=1
df.ems[emk,] #ecov_how=0

om                       <- fit_wham(om_inputs[[omj]], do.fit = FALSE, MakeADFun.silent = TRUE)
sim_data                 <- om$simulate(complete=TRUE)
EM_input                 <- em_inputs[[emk]] 
EM_input$data[obs_names] <- sim_data[obs_names]

fit <- fit_wham(EM_input, do.sdrep=F, do.osa=F, do.retro=F, do.proj=F, MakeADFun.silent=TRUE)

#######################################################
# Example 3; this errors out due to optimizer returning NANs
#######################################################
omj <- 1; emk <- 3  #both om and em have ecov$how=1
df.oms[omj,] #ecov_how=1
df.ems[emk,] #ecov_how=0

om                       <- fit_wham(om_inputs[[omj]], do.fit = FALSE, MakeADFun.silent = TRUE)
sim_data                 <- om$simulate(complete=TRUE)
EM_input                 <- em_inputs[[emk]] 
EM_input$data[obs_names] <- sim_data[obs_names]

fit <- fit_wham(EM_input, do.sdrep=F, do.osa=F, do.retro=F, do.proj=F, MakeADFun.silent=TRUE)
