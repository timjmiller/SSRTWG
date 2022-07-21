.libPaths("~/Rlib/")
library(here)
library(wham)
df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.oms.RDS"))
om_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "NAA_om_inputs.RDS"))
em_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "em_inputs.RDS"))
oms = 1:length(om_inputs)
ems = 1:length(em_inputs)
source(file.path(here::here(),"Project_0","code","sim_management.R"))

sapply(1:24, function(x) aggregate_hpcc_results(sim = 1, oms = x, ems = ems, res_dir = file.path(here::here(),"Project_0", "results", "naa_om")))
#aggregate_hpcc_results(sim = 1, oms = 1, ems = ems, res_dir = file.path(here::here(),"Project_0", "results", "naa_om"))
#aggregate_hpcc_results(sim = 1, oms = 3, ems = ems, res_dir = file.path(here::here(),"Project_0", "results", "naa_om"))

sapply(1:100, function(x) aggregate_hpcc_results(sim = x, oms = 1, ems = ems, res_dir = file.path(here::here(),"Project_0", "results", "naa_om"))

x = readRDS(file.path(here::here(),"Project_0", "results", "naa_om", "om_1", "sim_7.RDS"))
x = readRDS(file.path(here::here(),"Project_0", "results", "naa_om", "om_1", "sim_1.RDS"))
x = readRDS(file.path(here::here(),"Project_0", "results", "naa_om", "om_3", "sim_1.RDS"))

aggregate_hpcc_results(sim = 1, oms = 10, ems = ems, res_dir = file.path(here::here(),"Project_0", "results", "naa_om"))
x = readRDS(file.path(here::here(),"Project_0", "results", "naa_om", "om_10", "sim_1.RDS"))

aic = 2*sapply(x,function(y) {
  out = NA
  if(!is.null(y)) out = y$fit$opt$obj + length(y$fit$opt$par)
  return(out)
})
aic - min(aic, na.rm = TRUE)
