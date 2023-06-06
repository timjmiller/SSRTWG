library("here")
df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))
all_naa_aic <- readRDS(file = file.path(here(),"Project_0","results", "all_naa_aic_results.RDS"))

#first NAA operating models
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.oms.RDS"))
use.df.ems <- df.ems[1:20,]
library(Hmisc)
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

#OM: just rec re, EM: Assume R only, Estimating R and SR, M fixed
om_ind <- which(is.na(df.oms$NAA_sig))
SR_rec_Mfixed = which(use.df.ems$re_config == "rec" & use.df.ems$M_est == FALSE)
temp <- t(aic_fn(all_naa_aic, SR_rec_Mfixed, om_ind = om_ind))
# out <- cbind(df.oms[om_ind,-1], temp)
# colnames(out) <- c("$\\sigma_R$", "$\\sigma_N$", "F-history", "Obs Error", "R only", "BH") 
# x = latex(out, file = here("Project_0","paper","rec_om_em_R_SR_MF_aic_table.tex"), 
#   table.env = FALSE, col.just = rep("r", dim(out)[2]), rowname = NULL)

#OM: just rec re, EM: Assume R only, Estimating R and SR, M estimated
om_ind <- which(is.na(df.oms$NAA_sig))
SR_rec_Mest = which(use.df.ems$re_config == "rec" & use.df.ems$M_est == TRUE)
temp <- cbind(temp, t(aic_fn(all_naa_aic, SR_rec_Mest, om_ind = om_ind)))
colnames(temp) <- c("R (M fix)", "BH (M fix)", "R (M est)", "BH (M est)")
#out <- cbind(out, temp)
out <- cbind(df.oms[om_ind,-1], temp)
#colnames(out) <- c("$\\sigma_R$", "$\\sigma_N$", "F-history", "Obs Error", "R only", "BH") 
# x = latex(out, file = here("Project_0","paper","rec_om_em_R_SR_ME_aic_table.tex"), 
#   table.env = FALSE, col.just = rep("r", dim(out)[2]), rowname = NULL)

#OM: NAA re, EM: Assume rec+1, Estimating R and SR, M fixed
om_ind <- which(!is.na(df.oms$NAA_sig))
SR_rec_Mfixed = which(use.df.ems$re_config == "rec+1" & use.df.ems$M_est == FALSE)
temp <- t(aic_fn(all_naa_aic, SR_rec_Mfixed, om_ind = om_ind))
#out <- cbind(df.oms[om_ind,-1], temp)
# colnames(out) <- c("$\\sigma_R$", "$\\sigma_N$", "F-history", "Obs Error", "R only", "BH") 
# x = latex(out, file = here("Project_0","paper","naa_om_em_R_SR_MF_aic_table.tex"), 
#   table.env = FALSE, col.just = rep("r", dim(out)[2]), rowname = NULL)

#OM: NAA re, EM: Assume rec+1, Estimating R and SR, M estimated
om_ind <- which(!is.na(df.oms$NAA_sig))
SR_rec_Mest = which(use.df.ems$re_config == "rec+1" & use.df.ems$M_est == TRUE)
temp <- cbind(temp, t(aic_fn(all_naa_aic, SR_rec_Mest, om_ind = om_ind)))
colnames(temp) <- c("R (M fix)", "BH (M fix)", "R (M est)", "BH (M est)")
temp <- cbind(df.oms[om_ind,-1], temp)

out <- rbind(out, temp)
colnames(out)[1:4] <- c("$\\sigma_R$", "$\\sigma_N$", "F-history", "Obs Error")
x = latex(out, file = here("Project_0","paper","naa_om_em_R_BH_aic_table.tex"), 
  table.env = FALSE, col.just = rep("r", dim(out)[2]), rowname = NULL)


#OM: M re, EM: Estimating R and SR, M fixed
ems <- c(5:24)
use.df.ems = df.ems[ems,]
ind <- which(use.df.ems$re_config == "M_re")
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.M.oms.RDS"))
all_M_aic <- readRDS(file = file.path(here(),"Project_0","results", "all_M_aic_results.RDS"))
use.df.ems$re_config[ind] = paste0("M_",c("iid","ar1")[match(use.df.ems$M_re_cor[ind], c("iid","ar1_y"))])

#OM: M re are iid and EM assumption matches
om_ind <- which(df.oms$M_cor == 0)
SR_rec_Mfixed = which(use.df.ems$re_config %in% c("M_iid") & use.df.ems$M_est == FALSE)
temp <- t(aic_fn(all_M_aic, SR_rec_Mfixed, om_ind = om_ind))
# M estimated in EMs
om_ind <- which(df.oms$M_cor == 0)
SR_rec_Mfixed = which(use.df.ems$re_config %in% c("M_iid") & use.df.ems$M_est == TRUE)
temp <- cbind(temp, t(aic_fn(all_M_aic, SR_rec_Mfixed, om_ind = om_ind)))
colnames(temp) <- c("R (M fix)", "BH (M fix)", "R (M est)", "BH (M est)")
out <- cbind(df.oms[om_ind,c("M_sig", "M_cor", "Fhist", "obs_error")], temp)

#OM: M re are cor and EM assumption matches, M fixed
om_ind <- which(df.oms$M_cor != 0)
SR_rec_Mfixed = which(use.df.ems$re_config %in% c("M_ar1") & use.df.ems$M_est == FALSE)
temp <- t(aic_fn(all_M_aic, SR_rec_Mfixed, om_ind = om_ind))
#OM: M re are cor and EM assumption matches, M est
SR_rec_Mfixed = which(use.df.ems$re_config %in% c("M_ar1") & use.df.ems$M_est == TRUE)
temp <- cbind(temp, t(aic_fn(all_M_aic, SR_rec_Mfixed, om_ind = om_ind)))
colnames(temp) <- c("R (M fix)", "BH (M fix)", "R (M est)", "BH (M est)")

out <- rbind(out, cbind(df.oms[om_ind,c("M_sig", "M_cor", "Fhist", "obs_error")], temp))
colnames(out)[1:4] <- c("$\\sigma_M$", "$\\rho_M$", "F-history", "Obs Error")
x = latex(out, file = here("Project_0","paper","M_om_em_R_BH_aic_table.tex"), 
  table.env = FALSE, col.just = rep("r", dim(out)[2]), rowname = NULL)

#OM: Sel re, EM: Estimating R and SR, M fixed
ems <- c(5:20,25:28)
use.df.ems = df.ems[ems,]
ind <- which(use.df.ems$re_config == "sel_re")
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.Sel.oms.RDS"))
all_Sel_aic <- readRDS(file = file.path(here(),"Project_0","results", "all_Sel_aic_results.RDS"))
use.df.ems$re_config[ind] = paste0("Sel_",c("iid","ar1")[match(use.df.ems$sel_re_cor[ind], c("iid","ar1_y"))])

om_ind <- which(df.oms$Sel_cor == 0)
SR_rec_Mfixed = which(use.df.ems$re_config %in% c("Sel_iid") & use.df.ems$M_est == FALSE)
temp <- t(aic_fn(all_Sel_aic, SR_rec_Mfixed, om_ind = om_ind))
# M estimated in EMs
om_ind <- which(df.oms$Sel_cor == 0)
SR_rec_Mfixed = which(use.df.ems$re_config %in% c("Sel_iid") & use.df.ems$M_est == TRUE)
temp <- cbind(temp, t(aic_fn(all_Sel_aic, SR_rec_Mfixed, om_ind = om_ind)))
colnames(temp) <- c("R (M fix)", "BH (M fix)", "R (M est)", "BH (M est)")
out <- cbind(df.oms[om_ind,c("Sel_sig", "Sel_cor", "Fhist", "obs_error")], temp)

#OM: Sel re are cor and EM assumption matches, M fixed
om_ind <- which(df.oms$Sel_cor != 0)
SR_rec_Mfixed = which(use.df.ems$re_config %in% c("Sel_ar1") & use.df.ems$M_est == FALSE)
temp <- t(aic_fn(all_Sel_aic, SR_rec_Mfixed, om_ind = om_ind))
#OM: Sel re are cor and EM assumption matches, M est
SR_rec_Mfixed = which(use.df.ems$re_config %in% c("Sel_ar1") & use.df.ems$M_est == TRUE)
temp <- cbind(temp, t(aic_fn(all_Sel_aic, SR_rec_Mfixed, om_ind = om_ind)))
colnames(temp) <- c("R (M fix)", "BH (M fix)", "R (M est)", "BH (M est)")

out <- rbind(out, cbind(df.oms[om_ind,c("Sel_sig", "Sel_cor", "Fhist", "obs_error")], temp))
colnames(out)[1:4] <- c("$\\sigma_{Sel}$", "$\\rho_{Sel}$", "F-history", "Obs Error")
x = latex(out, file = here("Project_0","paper","Sel_om_em_R_BH_aic_table.tex"), 
  table.env = FALSE, col.just = rep("r", dim(out)[2]), rowname = NULL)

#OM: q re, EM: Estimating R and SR, M fixed
ems <- c(5:20,29:32)
use.df.ems = df.ems[ems,]
ind <- which(use.df.ems$re_config == "q_re")
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.q.oms.RDS"))
all_q_aic <- readRDS(file = file.path(here(),"Project_0","results", "all_q_aic_results.RDS"))
use.df.ems$re_config[ind] = paste0("q_",c("iid","ar1")[match(use.df.ems$q_re_cor[ind], c("iid","ar1"))])

om_ind <- which(df.oms$q_cor == 0)
SR_rec_Mfixed = which(use.df.ems$re_config %in% c("q_iid") & use.df.ems$M_est == FALSE)
temp <- t(aic_fn(all_q_aic, SR_rec_Mfixed, om_ind = om_ind))
# M estimated in EMs
om_ind <- which(df.oms$q_cor == 0)
SR_rec_Mfixed = which(use.df.ems$re_config %in% c("q_iid") & use.df.ems$M_est == TRUE)
temp <- cbind(temp, t(aic_fn(all_q_aic, SR_rec_Mfixed, om_ind = om_ind)))
colnames(temp) <- c("R (M fix)", "BH (M fix)", "R (M est)", "BH (M est)")
out <- cbind(df.oms[om_ind,c("q_sig", "q_cor", "Fhist", "obs_error")], temp)

#OM: q re are cor and EM assumption matches, M fixed
om_ind <- which(df.oms$q_cor != 0)
SR_rec_Mfixed = which(use.df.ems$re_config %in% c("q_ar1") & use.df.ems$M_est == FALSE)
temp <- t(aic_fn(all_q_aic, SR_rec_Mfixed, om_ind = om_ind))
#OM: q re are cor and EM assumption matches, M est
SR_rec_Mfixed = which(use.df.ems$re_config %in% c("q_ar1") & use.df.ems$M_est == TRUE)
temp <- cbind(temp, t(aic_fn(all_q_aic, SR_rec_Mfixed, om_ind = om_ind)))
colnames(temp) <- c("R (M fix)", "BH (M fix)", "R (M est)", "BH (M est)")

out <- rbind(out, cbind(df.oms[om_ind,c("q_sig", "q_cor", "Fhist", "obs_error")], temp))
colnames(out)[1:4] <- c("$\\sigma_{q}$", "$\\rho_{q}$", "F-history", "Obs Error")
x = latex(out, file = here("Project_0","paper","q_om_em_R_BH_aic_table.tex"), 
  table.env = FALSE, col.just = rep("r", dim(out)[2]), rowname = NULL)
