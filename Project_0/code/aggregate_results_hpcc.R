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
sapply(1:100, function(x) aggregate_hpcc_results(sim = x, oms = 1, ems = 1:20, res_dir = file.path(here::here(),"Project_0", "results", "naa_om")))
sapply(1:100, function(x) aggregate_hpcc_results(sim = x, oms = 2, ems = 1:20, res_dir = file.path(here::here(),"Project_0", "results", "naa_om")))
sapply(1:100, function(x) aggregate_hpcc_results(sim = x, oms = 3, ems = 1:20, res_dir = file.path(here::here(),"Project_0", "results", "naa_om")))
sapply(4:24, function(y) sapply(1:100, function(x) aggregate_hpcc_results(sim = x, oms = y, ems = 1:20, res_dir = file.path(here::here(),"Project_0", "results", "naa_om"))))
#aggregate_hpcc_results(sim = 1, oms = 1, ems = ems, res_dir = file.path(here::here(),"Project_0", "results", "naa_om"))
#aggregate_hpcc_results(sim = 1, oms = 3, ems = ems, res_dir = file.path(here::here(),"Project_0", "results", "naa_om"))

#M operating models
df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.M.oms.RDS"))
om_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "M_om_inputs.RDS"))
em_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "em_inputs.RDS"))
oms = 1:length(om_inputs)
ems = 1:length(em_inputs)
ems = 5:24
aggregate_hpcc_results(sim = 1, oms = 1, ems = 5:24, res_dir = file.path(here::here(),"Project_0", "results", "M_om"))
sapply(1:100, function(x) aggregate_hpcc_results(sim = x, oms = 1, ems = 5:24, res_dir = file.path(here::here(),"Project_0", "results", "M_om")))
sapply(2:16, function(y) sapply(1:100, function(x) aggregate_hpcc_results(sim = x, oms = y, ems = 5:24, res_dir = file.path(here::here(),"Project_0", "results", "M_om"))))


#Sel operating models
df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.Sel.oms.RDS"))
om_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "Sel_om_inputs.RDS"))
em_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "em_inputs.RDS"))
oms = 1:length(om_inputs)
ems = 1:length(em_inputs)
ems = c(5:20,25:28)
aggregate_hpcc_results(sim = 1, oms = 1, ems = ems, res_dir = file.path(here::here(),"Project_0", "results", "Sel_om"))
sapply(1:100, function(x) aggregate_hpcc_results(sim = x, oms = 1, ems = ems, res_dir = file.path(here::here(),"Project_0", "results", "Sel_om")))
sapply(2:16, function(y) sapply(1:100, function(x) aggregate_hpcc_results(sim = x, oms = y, ems = ems, res_dir = file.path(here::here(),"Project_0", "results", "Sel_om"))))



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
