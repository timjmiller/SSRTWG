#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
sim = as.integer(args[1])
om = as.integer(args[2])
em = as.integer(args[3])
.libPaths("~/Rlib/")
library(wham)
source(file.path(here::here(), "Project_0","code", "sim_management.R"))
verify_version()
#cat(paste0("number of cores available: ", parallel::detectCores(), "/n"))
#stop()
run_hpcc_jobs_rev(this_sim = sim, this_om = om, this_em = em,
  script.full.path = file.path(here::here(),"Project_0", "code", "Sel_om_sim_fit_script_hpcc.R"), 
  df.ems = readRDS(file.path(here::here(),"Project_0","inputs", "df.ems.RDS")), 
  df.oms = readRDS(file.path(here::here(),"Project_0","inputs", "df.Sel.oms.RDS")), 
  om_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "Sel_om_inputs.RDS")),
  em_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "em_inputs.RDS")),
  write.dir = file.path(here::here(),"Project_0", "results", "Sel_om")
)
