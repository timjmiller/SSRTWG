# script to generate and plot the simulated data
# liz brooks
# feb 21, 2023
# example using the Ecov_study for Mortality on branch gregdev



# load required libraries ====

if(file.exists("C:/Users/liz.brooks")) {
  library(wham, lib.loc = "C:/Users/liz.brooks/Documents/R/win-library/4.1")
} else library(wham) #make sure to use the right version of wham
library(tidyr)
library(dplyr)
library(here)
library(ggplot2)

# source code dependencies ====
source(here("common_code", "fn_sim_plot_data.R"))  #functions to simulate data and make some plots
source(file.path(here(), "Ecov_study", "mortality", "code", "sim_management.R")) #this has verify_version()

# check wham version ====
verify_version()
# The right commit: 77bbd94 of wham is loaded! #this should be the output of verify_version()


# specify output directory for simulated data (rds file) and plots (png) ====
sim.out=here("Ecov_study", "mortality", "simdata_RperS")


# specify location of inputs for function calls ====
seed.rds=file.path(here(), "Ecov_study", "mortality", "inputs","seeds.RDS")           #random seeds
om.inputs.rds=file.path(here(), "Ecov_study", "mortality", "inputs", "M_om_inputs.RDS") #inputs for oms
df.oms.rds=file.path(here(),"Ecov_study", "mortality", "inputs", "df.oms.RDS")        #matrix of om specs


# call function to simulate_data  ====
## some oms that i've compared: 
##    om.mod.nums=c(10, 11, 12, 34, 35, 36)   
##    om.mod.nums=c(10, 11, 12, 58, 59, 60)  
##    om.mod.nums=c(22, 23, 24, 70, 71, 72)

sim_data <- simulate_data(seed.rds=seed.rds,   #pass list of seeds
                          same.seeds=TRUE,  # if TRUE, uses the same 1000 seeds for all oms
                          om.inputs.rds=om.inputs.rds, #pass name of rds file with om.inputs
                          df.oms.rds=df.oms.rds,    #pass name of rds file with df.om specs
                          om.mod.nums=c(22, 23, 24, 70, 71, 72) , #pass vector of om model numbers to be simulated (1-288)
                          sim.nums=c(1,22,116),     #pass simulation iteration numbers (corresponds to seeds)
                          simdata_dir=sim.out,      #pass the directory for simulated output
                          write.sim.data=FALSE)     #if TRUE, writes rds file with simulated data



# call function to plot simulate_data  ====
#vec.ts.obj.names <- c( "SSB" )               #can be just one time series object
vec.ts.obj.names <- c( "agg_catch", "agg_indices", 
                       "Ecov_x" , "Ecov_obs", "Ecov_obs_sigma",  "Ecov_re",
                        "SSB", "F", "MAA", "NAA", "q" )  #or can be all time series objects
                    #note: MAA just plots M[1], since re are shared over all ages

plot_sim_data(ts.obj.names=vec.ts.obj.names,   # vector with any of the above obj names (i haven't tested any other output)
              sim_data=sim_data,      #list of simulated data (from function sim_data)
              simdata_dir=sim.out,    #directory for saving output
              df.oms.rds=df.oms.rds,  #pass name of rds file with df.om specs
              pheight=12, pwidth=6, text.size=3)    #plot height and plot width of saved png plots



###  END

# recruitment inputs  ====
wd <- file.path(here(),"Ecov_study", "recruitment", "inputs_20230410")
seed.rds=file.path(here(), "Ecov_study", "recruitment", "inputs","seeds.RDS")           #random seeds
om.inputs.rds=file.path(wd, "om_inputs.RDS") #inputs for oms
df.oms.rds=file.path(wd, "df.oms.RDS")        #matrix of om specs
sim.out=here(wd, "simdata1")

om.nums <-  seq(41,44)  #seq(5,8)   #seq(197,200)  #c(9,33,393,417)  #seq(9,12)
sim.nums <- c(22,47,116)  #c(25, 44, 95)
sim_data <- simulate_data(seed.rds=seed.rds,   #pass list of seeds
                          same.seeds=TRUE,  # if TRUE, uses the same 1000 seeds for all oms
                          om.inputs.rds=om.inputs.rds, #pass name of rds file with om.inputs
                          df.oms.rds=df.oms.rds,    #pass name of rds file with df.om specs
                          om.mod.nums=om.nums , #pass vector of om model numbers to be simulated (1-288)
                          sim.nums=sim.nums,     #pass simulation iteration numbers (corresponds to seeds)
                          simdata_dir=sim.out,      #pass the directory for simulated output
                          write.sim.data=TRUE)     #if TRUE, writes rds file with simulated data


vec.ts.obj.names <- c( "agg_catch", "agg_indices", 
                       "Ecov_x" , "Ecov_obs", "Ecov_obs_sigma",  "Ecov_re",
                       "SSB", "F", "MAA", "NAA", "q" )  #or can be all time series objects
#note: MAA just plots M[1], since re are shared over all ages

plot_sim_data(ts.obj.names=vec.ts.obj.names,   # vector with any of the above obj names (i haven't tested any other output)
              sim_data=sim_data,      #list of simulated data (from function sim_data)
              simdata_dir=sim.out,    #directory for saving output
              df.oms.rds=df.oms.rds,  #pass name of rds file with df.om specs
              pheight=12, pwidth=6, text.size=3)    #plot height and plot width of saved png plots

