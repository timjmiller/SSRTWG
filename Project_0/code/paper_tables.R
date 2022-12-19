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


all_naa_aic = lapply(1:NROW(df.oms), function(y){
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
  return(res)
})
saveRDS(all_naa_aic, file = file.path(here(),"Project_0","results", "all_naa_aic_results.RDS"))

x = readRDS(file.path(here::here(),"Project_0", "results", "naa_om", paste0("om_", 1), paste0("sim_",1,".RDS")))
x = readRDS(file.path(here::here(),"Project_0", "results", "naa_om", paste0("om_", 1), paste0("sim_",27,".RDS")))

all_naa_relssb = lapply(1:NROW(df.oms), function(y){
  res = lapply(1:100, function(x){
    print(paste0("sim_",x,".RDS"))
    sim = readRDS(file.path(here::here(),"Project_0", "results", "naa_om", paste0("om_", y), paste0("sim_",x,".RDS")))
    relSSB = sapply(sim,function(y) {
      out <- rep(NA, 40)
      if(length(y)) if(length(y$fit)) out = y$fit$rep$SSB/y$truth$SSB
      return(out)
    })
    return(relSSB)
  })
  return(res)
})
saveRDS(all_naa_relssb, file = here("Project_0","results", "all_naa_relssb_results.RDS"))

df.M.oms = readRDS(file.path(here(),"Project_0","inputs", "df.M.oms.RDS"))

all_M_om_relssb = lapply(1:NROW(df.M.oms), function(y){
  res = lapply(1:100, function(x){
    print(paste0("sim_",x,".RDS"))
    sim = readRDS(file.path(here::here(),"Project_0", "results", "M_om", paste0("om_", y), paste0("sim_",x,".RDS")))
    relSSB = sapply(sim,function(y) {
      out <- rep(NA, 40)
      if(length(y)) if(length(y$fit)) out = y$fit$rep$SSB/y$truth$SSB
      return(out)
    })
    return(relSSB)
  })
  return(res)
})
saveRDS(all_M_om_relssb, file = here("Project_0","results", "all_M_om_relssb_results.RDS"))


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

all_naa_aic <- readRDS(file = file.path(here(),"Project_0","results", "all_naa_aic_results.RDS"))
aic_fn <- function(all, est_ind, om_ind = NULL){
  if(!is.null(om_ind)) all <- all[om_ind]
  out <- sapply(all, function(res){
    tmp = apply(res[est_ind,],2, function(x) {
      if(any(!is.na(x))) {
        return(x == min(x,na.rm=T))
      } else return(rep(NA, length(x)))
    })
    return(apply(tmp,1,sum,na.rm=T))
  })  
  return(out)
}

temp <- df.ems[1:20,]
use.df.ems <- df.ems[1:20,]
ynames = paste0(
  c("MR","SR")[match(use.df.ems$SR_model, c(2,3))],
  "-", c("ME","MF")[match(use.df.ems$M_est,c(TRUE,FALSE))],
  "-", use.df.ems$re_config)
library(Hmisc)

#EM: estimate SR, M fixed
SR_Mfixed = which(use.df.ems$SR_model ==3 & use.df.ems$M_est == FALSE)
naa_em_SR_Mfixed <- t(aic_fn(all_naa_aic, SR_Mfixed))
out <- cbind(df.oms[-1], naa_em_SR_Mfixed)
colnames(out) <- c("$\\sigma_R$", "$\\sigma_N$", "F-history", "Obs Error", "R only", "NAA", "M", "Sel", "q") 
x = latex(out, file = here("Project_0","paper","naa_om_em_SR_MF_aic_table.tex"), 
  table.env = FALSE, col.just = rep("r", dim(out)[2]), rowname = NULL)

#EM: estimate SR, M estimated
SR_Mest = which(use.df.ems$SR_model ==3 & use.df.ems$M_est == TRUE)
naa_em_SR_Mest <- t(aic_fn(all_naa_aic, SR_Mest))
out <- cbind(df.oms[-1], naa_em_SR_Mest)
colnames(out) <- c("$\\sigma_R$", "$\\sigma_N$", "F-history", "Obs Error", "R only", "NAA", "M", "Sel", "q") 
x = latex(out, file = here("Project_0","paper","naa_om_em_SR_ME_aic_table.tex"), 
  table.env = FALSE, col.just = rep("r", dim(out)[2]), rowname = NULL)

#EM: No SR, M fixed
R_Mfixed = which(use.df.ems$SR_model ==2 & use.df.ems$M_est == FALSE)
naa_em_R_Mfixed <- t(aic_fn(all_naa_aic, R_Mfixed))
out <- cbind(df.oms[-1], naa_em_R_Mfixed)
colnames(out) <- c("$\\sigma_R$", "$\\sigma_N$", "F-history", "Obs Error", "R only", "NAA", "M", "Sel", "q") 
x = latex(out, file = here("Project_0","paper","naa_om_em_R_MF_aic_table.tex"), 
  table.env = FALSE, col.just = rep("r", dim(out)[2]), rowname = NULL)

#EM: No SR, M estimated
R_Mest = which(use.df.ems$SR_model ==2 & use.df.ems$M_est == TRUE)
naa_em_R_Mest <- t(aic_fn(all_naa_aic, R_Mest))
out <- cbind(df.oms[-1], naa_em_R_Mest)
colnames(out) <- c("$\\sigma_R$", "$\\sigma_N$", "F-history", "Obs Error", "R only", "NAA", "M", "Sel", "q") 
x = latex(out, file = here("Project_0","paper","naa_om_em_R_ME_aic_table.tex"), 
  table.env = FALSE, col.just = rep("r", dim(out)[2]), rowname = NULL)


use.df.ems[SR_Mest,]


#OM: just rec re, EM: Assume R only, Estimating R and SR, M fixed
om_ind <- which(is.na(df.oms$NAA_sig))
SR_rec_Mfixed = which(use.df.ems$re_config == "rec" & use.df.ems$M_est == FALSE)
temp <- t(aic_fn(all_naa_aic, SR_rec_Mfixed, om_ind = om_ind))
out <- cbind(df.oms[om_ind,-1], temp)
colnames(out) <- c("$\\sigma_R$", "$\\sigma_N$", "F-history", "Obs Error", "R only", "BH") 
x = latex(out, file = here("Project_0","paper","rec_om_em_R_SR_MF_aic_table.tex"), 
  table.env = FALSE, col.just = rep("r", dim(out)[2]), rowname = NULL)

#OM: just rec re, EM: Assume R only, Estimating R and SR, M estimated
om_ind <- which(is.na(df.oms$NAA_sig))
SR_rec_Mest = which(use.df.ems$re_config == "rec" & use.df.ems$M_est == TRUE)
temp <- t(aic_fn(all_naa_aic, SR_rec_Mest, om_ind = om_ind))
out <- cbind(df.oms[om_ind,-1], temp)
colnames(out) <- c("$\\sigma_R$", "$\\sigma_N$", "F-history", "Obs Error", "R only", "BH") 
x = latex(out, file = here("Project_0","paper","rec_om_em_R_SR_ME_aic_table.tex"), 
  table.env = FALSE, col.just = rep("r", dim(out)[2]), rowname = NULL)

#OM: NAA re, EM: Assume R only, Estimating R and SR, M fixed
om_ind <- which(!is.na(df.oms$NAA_sig))
SR_rec_Mfixed = which(use.df.ems$re_config == "rec+1" & use.df.ems$M_est == FALSE)
temp <- t(aic_fn(all_naa_aic, SR_rec_Mfixed, om_ind = om_ind))
out <- cbind(df.oms[om_ind,-1], temp)
colnames(out) <- c("$\\sigma_R$", "$\\sigma_N$", "F-history", "Obs Error", "R only", "BH") 
x = latex(out, file = here("Project_0","paper","naa_om_em_R_SR_MF_aic_table.tex"), 
  table.env = FALSE, col.just = rep("r", dim(out)[2]), rowname = NULL)

#OM: NAA re, EM: Assume R only, Estimating R and SR, M estimated
om_ind <- which(!is.na(df.oms$NAA_sig))
SR_rec_Mest = which(use.df.ems$re_config == "rec+1" & use.df.ems$M_est == TRUE)
temp <- t(aic_fn(all_naa_aic, SR_rec_Mest, om_ind = om_ind))
out <- cbind(df.oms[om_ind,-1], temp)
colnames(out) <- c("$\\sigma_R$", "$\\sigma_N$", "F-history", "Obs Error", "R only", "BH") 
x = latex(out, file = here("Project_0","paper","naa_om_em_R_SR_ME_aic_table.tex"), 
  table.env = FALSE, col.just = rep("r", dim(out)[2]), rowname = NULL)



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
#aic_table <- readRDS(file.path(here(),"Project_0","results", "aic_table_coarse.RDS"))
aic_table = aic_table/apply(aic_table,1,sum, na.rm = T)


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


all_naa_om_relssb <- readRDS(file = here("Project_0","results", "all_naa_relssb_results.RDS"))
median_relbias_ssb <- lapply(all_naa_om_relssb, function(x){
  sapply(1:20, function(y) {
    relssb_omx_emy <- sapply(1:100, function(z) {
      #print(length(x))
      #print(dim(x[[z]]))
      #print(c(y,z))
      return(x[[z]][,y])
    })
    return(apply(relssb_omx_emy,1,median, na.rm = TRUE))
  })
})

all_naa_om_relF <- readRDS(file = here("Project_0","results", "all_naa_relF_results.RDS"))
median_relbias_F <- lapply(all_naa_om_relF, function(x){
  sapply(1:20, function(y) {
    rel_omx_emy <- sapply(1:100, function(z) {
      #print(length(x))
      #print(dim(x[[z]]))
      #print(c(y,z))
      return(x[[z]][,y])
    })
    return(apply(rel_omx_emy,1,median, na.rm = TRUE))
  })
})

all_naa_om_relR <- readRDS(file = here("Project_0","results", "all_naa_relR_results.RDS"))
median_relbias_R <- lapply(all_naa_om_relR, function(x){
  sapply(1:20, function(y) {
    rel_omx_emy <- sapply(1:100, function(z) {
      #print(length(x))
      #print(dim(x[[z]]))
      #print(c(y,z))
      return(x[[z]][,y])
    })
    return(apply(rel_omx_emy,1,median, na.rm = TRUE))
  })
})

all_M_om_relssb <- readRDS(file = here("Project_0","results", "all_M_om_relssb_results.RDS"))
M_om_median_relbias_ssb <- lapply(all_M_om_relssb, function(x){
  sapply(1:20, function(y) {
    relssb_omx_emy <- sapply(1:100, function(z) {
      #print(length(x))
      #print(dim(x[[z]]))
      #print(c(y,z))
      return(x[[z]][,y])
    })
    return(apply(relssb_omx_emy,1,median, na.rm = TRUE))
  })
})

relbias_fn <- function(median_relbias, em_ind, yrange = c(-1,2)){

  line.labs = c("rec", "rec+1", "M", "sel", "q")
  xvals <- 1:40
  par(mfrow = c(4,6), oma = c(4,4,3,3), mar = c(1,1,1,1))
  for(om_ind in 1:24){
    yvals <- median_relbias[[om_ind]][,em_ind]-1
    ncolors = length(em_ind)
    colors.unique = hcl.colors(ncolors, palette = "viridis")
    plot(range(xvals), range(yvals), xlim = range(xvals), ylim = yrange, 
      xlab = "", ylab = "", type = "n", axes = FALSE)
    if(om_ind < 19) axis(1, labels = FALSE)
    else axis(1)
    if(om_ind %in% c(1,7,13,19)) axis(2)#, cex.axis = 0.75, las = 1)
    else axis(2, labels = FALSE)
    if(om_ind == 1) legend("topright", bty = "none", col = colors.unique, lty = 1, legend = line.labs)
    ulabs <- c("rec0.5", "rec1.5","rec0.5N0.25","rec1.5N0.25","rec0.5N0.5","rec1.5N0.5")
    rlabs <- c("F: H-MSY, OE Low","F: MSY, OE Low", "F: H-MSY, OE Hi","F: MSY, OE Hi")
    if(om_ind %in% 1:6) mtext(text = ulabs[om_ind], side = 3, line = 1)
    if(om_ind %in% c(6,12,18,24)) mtext(text = rlabs[om_ind/6], side = 4, line = 1)
    box()
    ind <- integer()
    for (j in 1:NCOL(yvals)) {
      #vals = as.integer(z[j,])+1
      lines(xvals, yvals[,j], col = colors.unique[j], bg = colors.unique[j])
      if(all(yvals[,j] < yrange[1] || yvals[,j] > yrange[2]))
      ind <- c(ind, j)
    }
    if(length(ind)) {
      print(ind)
      legend("bottomright", legend = paste0(line.labs[ind], " out"))
    }
  }
}

em_ind <- which(use.df.ems$SR_model == 2 & use.df.ems$M_est == FALSE)
use.df.ems[em_ind,]
png(file.path(here::here(), "Project_0", "paper", "naa_om_R_MF_relbias_ssb.png"), 
  width = 10*144, height = 8*144, res = 144, pointsize = 12)
relbias_fn(median_relbias_ssb, em_ind)
dev.off()

relbias_fn(median_relbias_F, em_ind)
relbias_fn(median_relbias_R, em_ind)

em_ind <- which(use.df.ems$SR_model == 2 & use.df.ems$M_est == TRUE)
use.df.ems[em_ind,]
png(file.path(here::here(), "Project_0", "paper", "naa_om_R_ME_relbias_ssb.png"), 
  width = 10*144, height = 8*144, res = 144, pointsize = 12)
relbias_fn(median_relbias_ssb, em_ind)
dev.off()

em_ind <- which(use.df.ems$SR_model == 3 & use.df.ems$M_est == FALSE)
use.df.ems[em_ind,]
png(file.path(here::here(), "Project_0", "paper", "naa_om_SR_MF_relbias_ssb.png"), 
  width = 10*144, height = 8*144, res = 144, pointsize = 12)
relbias_fn(median_relbias_ssb, em_ind)
dev.off()

em_ind <- which(use.df.ems$SR_model == 3 & use.df.ems$M_est == TRUE)
use.df.ems[em_ind,]
png(file.path(here::here(), "Project_0", "paper", "naa_om_SR_ME_relbias_ssb.png"), 
  width = 10*144, height = 8*144, res = 144, pointsize = 12)
relbias_fn(median_relbias_ssb, em_ind)
dev.off()

terminal_year_median_relbias_ssb <- sapply(median_relbias_ssb, function(x) x[40,])

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
