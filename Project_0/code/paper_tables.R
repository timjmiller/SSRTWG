library(here)
df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))


#NAA oms: ems = 1-20
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.oms.RDS"))
naa_outer_res = sapply(1:NROW(df.oms), function(y){
  res = sapply(1:100, function(x){
    print(paste0("sim_",x,".RDS"))
    sim = readRDS(file.path(here::here(),"Project_0", "results", "naa_om", paste0("om_", y), paste0("sim_",x,".RDS")))
    aic = 2*sapply(sim,function(y) {
      out = NA
      if(length(y$fit)) out = y$fit$opt$obj + length(y$fit$opt$par)
      return(out)
    })
    return(aic)
  })
  tmp = apply(res,2, function(x) {
    if(any(!is.na(x))) {
      return(x == min(x,na.rm=T))
    } else return(rep(NA, length(x)))
  })
  return(apply(tmp,1,sum,na.rm=T))
})
saveRDS(naa_outer_res, file = file.path(here(),"Project_0","results", "NAA_om_aic_results.RDS"))



#M oms: ems = 5-24
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.M.oms.RDS"))
M_outer_res = sapply(1:NROW(df.oms), function(y){
  res = sapply(1:100, function(x){
    print(paste0("sim_",x,".RDS"))
    sim = readRDS(file.path(here::here(),"Project_0", "results", "M_om", paste0("om_", y), paste0("sim_",x,".RDS")))
    aic = 2*sapply(sim,function(y) {
      out = NA
      if(length(y$fit)) out = y$fit$opt$obj + length(y$fit$opt$par)
      return(out)
    })
    return(aic)
  })
  tmp = apply(res,2, function(x) {
    if(any(!is.na(x))) {
      return(x == min(x,na.rm=T))
    } else return(rep(NA, length(x)))
  })
  return(apply(tmp,1,sum,na.rm=T))
})
saveRDS(M_outer_res, file = file.path(here(),"Project_0","results", "M_om_aic_results.RDS"))

#Sel oms: ems = 5-20, 25-28
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.Sel.oms.RDS"))
Sel_outer_res = sapply(1:NROW(df.oms), function(y){
  res = sapply(1:100, function(x){
    print(paste0("sim_",x,".RDS"))
    sim = readRDS(file.path(here::here(),"Project_0", "results", "Sel_om", paste0("om_", y), paste0("sim_",x,".RDS")))
    aic = 2*sapply(sim,function(y) {
      out = NA
      if(length(y$fit)) out = y$fit$opt$obj + length(y$fit$opt$par)
      return(out)
    })
    return(aic)
  })
  tmp = apply(res,2, function(x) {
    if(any(!is.na(x))) {
      return(x == min(x,na.rm=T))
    } else return(rep(NA, length(x)))
  })
  return(apply(tmp,1,sum,na.rm=T))
})
saveRDS(Sel_outer_res, file = file.path(here(),"Project_0","results", "Sel_om_aic_results.RDS"))


#q oms: ems = 5-20, 29-32
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.q.oms.RDS"))
q_outer_res = sapply(1:NROW(df.oms), function(y){
  res = sapply(1:100, function(x){
    print(paste0("sim_",x,".RDS"))
    sim = readRDS(file.path(here::here(),"Project_0", "results", "q_om", paste0("om_", y), paste0("sim_",x,".RDS")))
    aic = 2*sapply(sim,function(y) {
      out = NA
      if(length(y$fit)) out = y$fit$opt$obj + length(y$fit$opt$par)
      return(out)
    })
    return(aic)
  })
  tmp = apply(res,2, function(x) {
    if(any(!is.na(x))) {
      return(x == min(x,na.rm=T))
    } else return(rep(NA, length(x)))
  })
  return(apply(tmp,1,sum,na.rm=T))
})
saveRDS(q_outer_res, file = file.path(here(),"Project_0","results", "q_om_aic_results.RDS"))

df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.oms.RDS"))
which(is.na(df.oms$NAA_sig))
em_cats = c("rec","rec+1","M_re","sel_re","q_re")
rec_aic_row = rep(NA, length(em_cats))
names(rec_aic_row) = em_cats
temp = df.ems[1:20,]
for(i in c("rec","rec+1","M_re","sel_re","q_re")){
  rec_aic_row[i] = sum(naa_outer_res[which(temp$re_config==i),which(is.na(df.oms$NAA_sig))])
}

recplus1_aic_row = rep(NA, length(em_cats))
names(recplus1_aic_row) = em_cats
temp = df.ems[1:20,]
for(i in c("rec","rec+1","M_re","sel_re","q_re")){
  recplus1_aic_row[i] = sum(naa_outer_res[which(temp$re_config==i),which(!is.na(df.oms$NAA_sig))])
}

df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.M.oms.RDS"))
em_cats = c("rec+1","M_re","sel_re","q_re")
M_aic_row = rep(NA, length(em_cats))
names(M_aic_row) = em_cats
temp = df.ems[5:24,]
for(i in c("rec+1","M_re","sel_re","q_re")){
  M_aic_row[i] = sum(M_outer_res[which(temp$re_config==i),])
}
M_aic_row = c(NA,M_aic_row)
names(M_aic_row)[1] = "rec"

df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.Sel.oms.RDS"))
em_cats = c("rec+1","M_re","sel_re","q_re")
sel_aic_row = rep(NA, length(em_cats))
names(sel_aic_row) = em_cats
temp = df.ems[c(5:20,25:28),]
for(i in c("rec+1","M_re","sel_re","q_re")){
  sel_aic_row[i] = sum(Sel_outer_res[which(temp$re_config==i),])
}
sel_aic_row = c(NA,sel_aic_row)
names(sel_aic_row)[1] = "rec"

df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.q.oms.RDS"))
em_cats = c("rec+1","M_re","sel_re","q_re")
q_aic_row = rep(NA, length(em_cats))
names(q_aic_row) = em_cats
temp = df.ems[c(5:20,29:32),]
for(i in c("rec+1","M_re","sel_re","q_re")){
  q_aic_row[i] = sum(q_outer_res[which(temp$re_config==i),])
}
q_aic_row = c(NA,q_aic_row)
names(q_aic_row)[1] = "rec"

aic_table = rbind(rec_aic_row, recplus1_aic_row, M_aic_row, sel_aic_row, q_aic_row)
saveRDS(aic_table, file = file.path(here(),"Project_0","results", "aic_table_coarse.RDS"))
aic_table = aic_table/apply(aic_table,1,sum, na.rm = T)
x = latex(out, file = paste0(parentdir, "/paper/base_model_description.tex"), 
  table.env = FALSE, col.just = c("l","l","l","r","p{0.5\\textwidth}"), rowname = NULL)


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
