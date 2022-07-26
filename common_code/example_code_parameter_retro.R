# simple code to run Ex1 vignette and test functions 
#     to calculate mohn's rho for wham parameter estimates and make plots
# liz brooks



library(wham)
library(here)
library(tidyverse)


source(here("common_code", "fn_retro_par_estimates.R") )


wham.dir <- find.package("wham")
file.path(wham.dir, "example_scripts")
write.dir <- here("common_code", "test_retro_par")
if(dir.exists(write.dir)==F) dir.create(write.dir)
file.copy(from=file.path(wham.dir,"extdata","ex1_SNEMAYT.dat"), to=write.dir, overwrite=FALSE)

asap3 <- read_asap3_dat(here(write.dir,"ex1_SNEMAYT.dat"  )  )

input1 <- prepare_wham_input(asap3, recruit_model=2, model_name="Ex 1: SNEMA Yellowtail Flounder",
                             selectivity=list(model=rep("age-specific",3), 
                                              re=rep("none",3), 
                                              initial_pars=list(c(0.5,0.5,0.5,1,1,0.5),c(0.5,0.5,0.5,1,0.5,0.5),c(0.5,1,1,1,0.5,0.5)),  
                                              fix_pars=list(4:5,4,2:4)),
                             NAA_re = list(sigma="rec", cor="iid"))




m1 <- fit_wham(input1, do.osa = F) # turn off OSA residuals to save time in ex



plot.retro.pars(m=m1, yr.plot.start=2005, od=here(write.dir, "Ex1"), save.name="Ex1", plot.f='png' )
# m is a fitted wham model
# this function produces two plots for each parameter: one of the full time series, 
#          and one that only plots years yr.plot.start - final model year
# od is the output directory where plots and csv files are saved
# save.name is prepended in front of each saved file (enter "" if none desired)
# plot.f is the file type used to save plots 
#         (a combined pdf of all plots is automatically created)
# the value of rho is printed in the facet for all plots except index selectivity, 
#         for which rho is printed on the plot above the selectivity at age curve
