
x = readRDS(file.path(here::here(),"Project_0", "results", "naa_om", paste0("om_", 1), paste0("sim_",1,".RDS")))
x = readRDS(file.path(here::here(),"Project_0", "results", "naa_om", paste0("om_", 1), paste0("sim_",27,".RDS")))


all_naa_relR = lapply(1:NROW(df.oms), function(y){
  res = lapply(1:100, function(x){
    print(paste0("sim_",x,".RDS"))
    sim = readRDS(file.path(here::here(),"Project_0", "results", "naa_om", paste0("om_", y), paste0("sim_",x,".RDS")))
    relSSB = sapply(sim,function(y) {
      out <- rep(NA, 40)
      if(length(y)) if(length(y$fit)) out = y$fit$rep$NAA[,1]/y$truth$NAA[,1]
      return(out)
    })
    return(relSSB)
  })
  return(res)
})
saveRDS(all_naa_relR, file = here("Project_0","results", "all_naa_relR_results.RDS"))

all_naa_relF = lapply(1:NROW(df.oms), function(y){
  res = lapply(1:100, function(x){
    print(paste0("sim_",x,".RDS"))
    sim = readRDS(file.path(here::here(),"Project_0", "results", "naa_om", paste0("om_", y), paste0("sim_",x,".RDS")))
    relSSB = sapply(sim,function(y) {
      out <- rep(NA, 40)
      if(length(y)) if(length(y$fit)) out = y$fit$rep$F[,1]/y$truth$F[,1]
      return(out)
    })
    return(relSSB)
  })
  return(res)
})
saveRDS(all_naa_relF, file = here("Project_0","results", "all_naa_relF_results.RDS"))

use.df.ems[SR_Mest,]




df.ems = readRDS(file.path(here::here(),"Project_0","inputs", "df.ems.RDS"))
df.oms = readRDS(file.path(here::here(),"Project_0","inputs", "df.oms.RDS"))
use.df.ems = df.ems[1:20,]
ynames = paste0(
  c("MR","SR")[match(use.df.ems$SR_model, c(2,3))],
  "_", c("ME","MF")[match(use.df.ems$M_est,c(TRUE,FALSE))],
  "_", use.df.ems$re_config)
xnames = paste0(
  c("HM","MSY")[match(df.oms$Fhist, c("H-MSY","MSY"))],
  "_", c("RsigL","RsigH")[match(df.oms$R_sig,c(0.5,1.5))],
  "_", c("Nsig0", "NsigL", "NsigH")[match(df.oms$NAA_sig, c(NA,0.25,0.5))], 
  "_", c("OEL", "OEH")[match(df.oms$obs_error, c("L","H"))])



naa_ems_aic = readRDS(file = file.path(here::here(),"Project_0","results", "NAA_om_aic_results.RDS"))
source(file.path(here::here(),"Project_0","code", "heat.plot.fn.R"))
png(file.path(here::here(), "Project_0", "results", "naa_om_aic.png"), 
  width = 10*144, height = 10*144, res = 144, pointsize = 12)
heat.plot.fn(naa_ems_aic, xlabs = xnames, ylabs = ynames, main.title = "AIC", palette = "viridis")
dev.off()

df.oms = readRDS(file.path(here::here(),"Project_0","inputs", "df.M.oms.RDS"))
use.df.ems = df.ems[5:24,]
use.df.ems$re_config[c(5:8,17:20)] = paste0("M_",c("iid","ar1")[match(use.df.ems$M_re_cor[c(5:8,17:20)], c("iid","ar1_y"))])
use.df.ems = use.df.ems[c(5:8,17:20,1:4,9:16),]
ynames = paste0(
  c("MR","SR")[match(use.df.ems$SR_model, c(2,3))],
  "_", c("ME","MF")[match(use.df.ems$M_est,c(TRUE,FALSE))],
  "_", use.df.ems$re_config)
xnames = paste0(
  c("HM","MSY")[match(df.oms$Fhist, c("H-MSY","MSY"))],
  "_", c("MsigL","MsigH")[match(df.oms$M_sig,c(0.1,0.5))],
  "_", c("Mcor0", "McorH")[match(df.oms$M_cor, c(0,0.9))], 
  "_", c("OEL", "OEH")[match(df.oms$obs_error, c("L","H"))])

M_ems_aic = readRDS(file = file.path(here::here(),"Project_0","results", "M_om_aic_results.RDS"))[c(5:8,17:20,1:4,9:16),]
source(file.path(here::here(),"Project_0","code", "heat.plot.fn.R"))
png(file.path(here::here(), "Project_0", "results", "M_om_aic.png"), 
  width = 10*144, height = 10*144, res = 144, pointsize = 12)
heat.plot.fn(M_ems_aic, xlabs = xnames, ylabs = ynames, main.title = "AIC", palette = "viridis")
dev.off()


df.oms = readRDS(file.path(here::here(),"Project_0","inputs", "df.Sel.oms.RDS"))
use.df.ems = df.ems[c(5:20,25:28),]
use.df.ems$re_config[c(9:12,17:20)] = paste0("Sel_",c("iid","ar1")[match(use.df.ems$sel_re_cor[c(9:12,17:20)], c("iid","ar1_y"))])
use.df.ems = use.df.ems[c(9:12,17:20,1:8,13:16),]
ynames = paste0(
  c("MR","SR")[match(use.df.ems$SR_model, c(2,3))],
  "_", c("ME","MF")[match(use.df.ems$M_est,c(TRUE,FALSE))],
  "_", use.df.ems$re_config)
xnames = paste0(
  c("HM","MSY")[match(df.oms$Fhist, c("H-MSY","MSY"))],
  "_", c("SelsigL","SelsigH")[match(df.oms$Sel_sig,c(0.1,0.5))],
  "_", c("Selcor0", "SelcorH")[match(df.oms$Sel_cor, c(0,0.9))], 
  "_", c("OEL", "OEH")[match(df.oms$obs_error, c("L","H"))])

Sel_ems_aic = readRDS(file = file.path(here::here(),"Project_0","results", "Sel_om_aic_results.RDS"))[c(9:12,17:20,1:8,13:16),]
source(file.path(here::here(),"Project_0","code", "heat.plot.fn.R"))
png(file.path(here::here(), "Project_0", "results", "Sel_om_aic.png"), 
  width = 10*144, height = 10*144, res = 144, pointsize = 12)
heat.plot.fn(Sel_ems_aic, xlabs = xnames, ylabs = ynames, main.title = "AIC", palette = "viridis")
dev.off()

df.oms = readRDS(file.path(here::here(),"Project_0","inputs", "df.q.oms.RDS"))
use.df.ems = df.ems[c(5:20,29:32),]
use.df.ems$re_config[c(13:20)] = paste0("q_",c("iid","ar1")[match(use.df.ems$q_re_cor[c(13:20)], c("iid","ar1"))])
use.df.ems = use.df.ems[c(13:20,1:12),]
ynames = paste0(
  c("MR","SR")[match(use.df.ems$SR_model, c(2,3))],
  "_", c("ME","MF")[match(use.df.ems$M_est,c(TRUE,FALSE))],
  "_", use.df.ems$re_config)
xnames = paste0(
  c("HM","MSY")[match(df.oms$Fhist, c("H-MSY","MSY"))],
  "_", c("qsigL","qsigH")[match(df.oms$q_sig,c(0.1,0.5))],
  "_", c("qcor0", "qcorH")[match(df.oms$q_cor, c(0,0.9))], 
  "_", c("OEL", "OEH")[match(df.oms$obs_error, c("L","H"))])

q_ems_aic = readRDS(file = file.path(here::here(),"Project_0","results", "q_om_aic_results.RDS"))[c(13:20,1:12),]
source(file.path(here::here(),"Project_0","code", "heat.plot.fn.R"))
png(file.path(here::here(), "Project_0", "results", "q_om_aic.png"), 
  width = 10*144, height = 10*144, res = 144, pointsize = 12)
heat.plot.fn(q_ems_aic, xlabs = xnames, ylabs = ynames, main.title = "AIC", palette = "viridis")
dev.off()
