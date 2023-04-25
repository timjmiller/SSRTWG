#.libPaths("~/Rlib/")
library(here)
df.ems = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.ems.RDS"))
df.oms = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.oms.RDS"))
om_inputs = readRDS(file.path(here::here(),"Ecov_study","mortality","inputs", "om_inputs.RDS"))
em_inputs = readRDS(file.path(here::here(),"Ecov_study","mortality","inputs", "em_inputs.RDS"))
oms = 1:length(om_inputs)
ems = 1:length(em_inputs)
source(file.path(here::here(),"Ecov_study","mortality","code","sim_management.R"))


aggregate_hpcc_results(sims = 1:100, oms = 1:288, ems = 1:12)
aggregate_hpcc_results(sims = 1:100, oms = 101:288, ems = 1:12)

aggregate_hpcc_results(sims = 1:100, oms = 8, ems = 3)
aggregate_hpcc_results(sims = 1:100, oms = 53, ems = 4)
aggregate_hpcc_results(sims = 1:100, oms = 26, ems = 1)
aggregate_hpcc_results(sims = 1:100, oms = 37, ems = 2)
