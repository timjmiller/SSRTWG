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

#OM: NAA re, EM: Assume rec+1, Estimating R and SR, M fixed
om_ind <- which(!is.na(df.oms$NAA_sig))
SR_rec_Mfixed = which(use.df.ems$re_config == "rec+1" & use.df.ems$M_est == FALSE)
temp <- t(aic_fn(all_naa_aic, SR_rec_Mfixed, om_ind = om_ind))
out <- cbind(df.oms[om_ind,-1], temp)
colnames(out) <- c("$\\sigma_R$", "$\\sigma_N$", "F-history", "Obs Error", "R only", "BH") 
x = latex(out, file = here("Project_0","paper","naa_om_em_R_SR_MF_aic_table.tex"), 
  table.env = FALSE, col.just = rep("r", dim(out)[2]), rowname = NULL)

#OM: NAA re, EM: Assume rec+1, Estimating R and SR, M estimated
om_ind <- which(!is.na(df.oms$NAA_sig))
SR_rec_Mest = which(use.df.ems$re_config == "rec+1" & use.df.ems$M_est == TRUE)
temp <- t(aic_fn(all_naa_aic, SR_rec_Mest, om_ind = om_ind))
out <- cbind(df.oms[om_ind,-1], temp)
colnames(out) <- c("$\\sigma_R$", "$\\sigma_N$", "F-history", "Obs Error", "R only", "BH") 
x = latex(out, file = here("Project_0","paper","naa_om_em_R_SR_ME_aic_table.tex"), 
  table.env = FALSE, col.just = rep("r", dim(out)[2]), rowname = NULL)

#OM: M re, EM: Estimating R and SR, M fixed
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.M.oms.RDS"))
use.df.ems = df.ems[5:24,]
use.df.ems$re_config[c(5:8,17:20)] = paste0("M_",c("iid","ar1")[match(use.df.ems$M_re_cor[c(5:8,17:20)], c("iid","ar1_y"))])
use.df.ems = use.df.ems[c(5:8,17:20,1:4,9:16),]

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
