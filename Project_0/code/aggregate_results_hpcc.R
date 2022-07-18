library(here)
library(wham)
df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.oms.RDS"))
om_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "NAA_om_inputs.RDS"))
em_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "em_inputs.RDS"))
oms = 1:length(om_inputs)
ems = 1:length(em_inputs)
source(file.path(here::here(),"Project_0","code","sim_management.R")

aggregate_hpcc_results(sim = 7, oms = 1, ems = ems, res_dir = file.path(here::here(),"Project_0", "results", "naa_om"))

x = readRDS(file.path(here::here(),"Project_0", "results", "naa_om", "om_1", "sim_7.RDS"))
x = readRDS(file.path(here::here(),"Project_0", "results", "naa_om", "om_1", "sim_1.RDS"))

aic = 2*sapply(x,function(y) y$fit$opt$obj + length(y$fit$opt$par))
aic - min(aic)
