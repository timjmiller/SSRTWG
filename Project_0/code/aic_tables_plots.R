library(here)
df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))

make_results <- FALSE
if(make_results) {
  #NAA oms: ems = 1-20
  df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.oms.RDS"))
  all_naa_aic = lapply(1:NROW(df.oms), function(y){
    res = sapply(1:100, function(x){
      print(paste0("om_", y, ", sim_",x))
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

  naa_outer_res = sapply(all_naa_aic, function(y){
    res <- y
    tmp = apply(res,2, function(x) {
      if(any(!is.na(x))) {
        return(x == min(x,na.rm=T))
      } else return(rep(NA, length(x)))
    })
    return(apply(tmp,1,sum,na.rm=T))
  })
  saveRDS(naa_outer_res, file = file.path(here(),"Project_0","results", "naa_om_aic_choice_results.RDS"))

  #M oms: ems = 5-20, 21-24
  df.M.oms = readRDS(file.path(here(),"Project_0","inputs", "df.M.oms.RDS"))
  all_M_aic = lapply(1:NROW(df.M.oms), function(y){
    res = sapply(1:100, function(x){
      print(paste0("om_", y, ", sim_",x))
      sim = readRDS(file.path(here::here(),"Project_0", "results", "M_om", paste0("om_", y), paste0("sim_",x,".RDS")))
      aic = 2*sapply(sim,function(y) {
        out = NA
        if(length(y$fit)) out = y$fit$opt$obj + length(y$fit$opt$par)
        return(out)
      })
      return(aic)
    })
    return(res)
  })
  saveRDS(all_M_aic, file = file.path(here(),"Project_0","results", "all_M_aic_results.RDS"))

  M_outer_res = sapply(all_M_aic, function(y){
    res <- y
    tmp = apply(res,2, function(x) {
      if(any(!is.na(x))) {
        return(x == min(x,na.rm=T))
      } else return(rep(NA, length(x)))
    })
    return(apply(tmp,1,sum,na.rm=T))
  })
  saveRDS(M_outer_res, file = file.path(here(),"Project_0","results", "M_om_aic_choice_results.RDS"))

  #Sel oms: ems = 5-20, 25-28
  df.Sel.oms = readRDS(file.path(here(),"Project_0","inputs", "df.Sel.oms.RDS"))
  all_Sel_aic = lapply(1:NROW(df.Sel.oms), function(y){
    res = sapply(1:100, function(x){
      print(paste0("om_", y, ", sim_",x))
      sim = readRDS(file.path(here::here(),"Project_0", "results", "Sel_om", paste0("om_", y), paste0("sim_",x,".RDS")))
      aic = 2*sapply(sim,function(y) {
        out = NA
        if(length(y$fit)) out = y$fit$opt$obj + length(y$fit$opt$par)
        return(out)
      })
      return(aic)
    })
    return(res)
  })
  saveRDS(all_Sel_aic, file = file.path(here(),"Project_0","results", "all_Sel_aic_results.RDS"))

  Sel_outer_res = sapply(all_Sel_aic, function(y){
    res <- y
    tmp = apply(res,2, function(x) {
      if(any(!is.na(x))) {
        return(x == min(x,na.rm=T))
      } else return(rep(NA, length(x)))
    })
    return(apply(tmp,1,sum,na.rm=T))
  })
  saveRDS(Sel_outer_res, file = file.path(here(),"Project_0","results", "Sel_om_aic_choice_results.RDS"))

  #q oms: ems = 5-20, 29-32
  df.q.oms = readRDS(file.path(here(),"Project_0","inputs", "df.q.oms.RDS"))
  all_q_aic = lapply(1:NROW(df.q.oms), function(y){
    res = sapply(1:100, function(x){
      print(paste0("om_", y, ", sim_",x))
      sim = readRDS(file.path(here::here(),"Project_0", "results", "q_om", paste0("om_", y), paste0("sim_",x,".RDS")))
      aic = 2*sapply(sim,function(y) {
        out = NA
        if(length(y$fit)) out = y$fit$opt$obj + length(y$fit$opt$par)
        return(out)
      })
      return(aic)
    })
    return(res)
  })
  saveRDS(all_q_aic, file = file.path(here(),"Project_0","results", "all_q_aic_results.RDS"))
  q_outer_res = sapply(all_q_aic, function(y){
    res <- y
    tmp = apply(res,2, function(x) {
      if(any(!is.na(x))) {
        return(x == min(x,na.rm=T))
      } else return(rep(NA, length(x)))
    })
    return(apply(tmp,1,sum,na.rm=T))
  })
  saveRDS(q_outer_res, file = file.path(here(),"Project_0","results", "q_om_aic_choice_results.RDS"))
}

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


df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.oms.RDS"))
which(is.na(df.oms$NAA_sig))
em_cats = c("rec","rec+1","M_re","sel_re","q_re")
rec_aic_row = rep(NA, length(em_cats))
names(rec_aic_row) = em_cats
temp = df.ems[1:20,]
#aic_fn(all_naa_aic, which(temp$re_config=="rec"), which(is.na(df.oms$NAA_sig)))
#aic_fn(all_naa_aic, which(temp$re_config=="rec+1"), which(is.na(df.oms$NAA_sig)))
aic_choice <- aic_fn(all_naa_aic, 1:20)
for(i in em_cats){
  rec_aic_row[i] <- sum(aic_choice[which(temp$re_config==i),which(is.na(df.oms$NAA_sig))])
}
recplus1_aic_row = rep(NA, length(em_cats))
names(recplus1_aic_row) = em_cats
temp = df.ems[1:20,]
for(i in em_cats){
  recplus1_aic_row[i] = sum(aic_choice[which(temp$re_config==i),which(!is.na(df.oms$NAA_sig))])
}

df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.M.oms.RDS"))
em_cats = c("rec+1","M_re","sel_re","q_re")
M_aic_row = rep(NA, length(em_cats))
names(M_aic_row) = em_cats
temp = df.ems[5:24,]
aic_choice <- aic_fn(all_M_aic, 1:20)
for(i in c("rec+1","M_re","sel_re","q_re")){
  M_aic_row[i] = sum(aic_choice[which(temp$re_config==i),])
}
M_aic_row = c(NA,M_aic_row)
names(M_aic_row)[1] = "rec"

df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.Sel.oms.RDS"))
em_cats = c("rec+1","M_re","sel_re","q_re")
sel_aic_row = rep(NA, length(em_cats))
names(sel_aic_row) = em_cats
temp = df.ems[c(5:20,25:28),]
aic_choice <- aic_fn(all_Sel_aic, 1:20)
for(i in c("rec+1","M_re","sel_re","q_re")){
  sel_aic_row[i] = sum(aic_choice[which(temp$re_config==i),])
}
sel_aic_row = c(NA,sel_aic_row)
names(sel_aic_row)[1] = "rec"

df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.q.oms.RDS"))
em_cats = c("rec+1","M_re","sel_re","q_re")
q_aic_row = rep(NA, length(em_cats))
names(q_aic_row) = em_cats
temp = df.ems[c(5:20,29:32),]
aic_choice <- aic_fn(all_q_aic, 1:20)
for(i in c("rec+1","M_re","sel_re","q_re")){
  q_aic_row[i] = sum(aic_choice[which(temp$re_config==i),])
}
q_aic_row = c(NA,q_aic_row)
names(q_aic_row)[1] = "rec"

aic_table = rbind(rec_aic_row, recplus1_aic_row, M_aic_row, sel_aic_row, q_aic_row)
aic_table = aic_table/apply(aic_table,1,sum, na.rm = T)
saveRDS(aic_table, file = file.path(here(),"Project_0","results", "aic_table_coarse.RDS"))
#aic_table <- readRDS(file.path(here(),"Project_0","results", "aic_table_coarse.RDS"))

make_plot_df <- function(aic_res, df.ems, df.oms, em_ind = 1:20, SR_est =FALSE, M_est = FALSE){
  use.df.ems <- df.ems[em_ind,]
  SR_M_ind = which(use.df.ems$SR_model == ifelse(SR_est,3,2) & use.df.ems$M_est == M_est)
  aic_res <- t(aic_fn(aic_res, SR_M_ind))
  print(aic_res)
  out <- cbind(df.oms[-1], aic_res)
  print(out)
  out <- out %>% pivot_longer(cols = as.character(1:5), names_to = "EM", values_to = "n")
  print(head(out))
  if(max(em_ind)==20) {
    out$NAA_sig[which(is.na(out$NAA_sig))] <- 0
    out <- out %>% mutate(NAA_sig = recode(NAA_sig,
        "0" = "sigma['2+'] == 0",
        "0.25" = "sigma['2+'] == 0.25",
        "0.5" = "sigma['2+'] == 0.5")) %>% as.data.frame
    out <- out %>% mutate(R_sig = recode(R_sig,
        "0.5" = "sigma[R] == 0.5",
        "1.5" = "sigma[R] == 1.5")) %>% as.data.frame
    out <- out %>% mutate(EM = recode(EM,
        "1" = "R",
        "2" = "R+S",
        "3" = "R+M",
        "4" = "R+Sel",
        "5" = "R+q"
      )) %>% as.data.frame
  }
  if(max(em_ind)==24) {
    out <- out %>% 
      mutate(M_sig = recode(M_sig,
        "0.1" = "sigma['M'] == 0.1",
        "0.5" = "sigma['M'] == 0.5")) %>% 
      mutate(M_cor = recode(M_cor,
        "0" = "rho['M'] == 0",
        "0.9" = "rho['M'] == 0.9")) %>%
      mutate(EM = recode(EM,
        "1" = "R+S",
        "2" = "R+M(iid)",
        "3" = "R+Sel",
        "4" = "R+q",
        "5" = "R+M(AR1)"
      )) %>% as.data.frame
  }
  if(max(em_ind)==28) {
    out <- out %>% 
      mutate(Sel_sig = recode(Sel_sig,
        "0.1" = "sigma['Sel'] == 0.1",
        "0.5" = "sigma['Sel'] == 0.5")) %>% 
      mutate(Sel_cor = recode(Sel_cor,
        "0" = "rho['Sel'] == 0",
        "0.9" = "rho['Sel'] == 0.9")) %>%
      mutate(EM = recode(EM,
        "1" = "R+S",
        "2" = "R+M",
        "3" = "R+Sel(iid)",
        "4" = "R+q",
        "5" = "R+Sel(AR1)"
      )) %>% as.data.frame
  }
  if(max(em_ind)==32) {
    out <- out %>% 
      mutate(q_sig = recode(q_sig,
        "0.1" = "sigma['q'] == 0.1",
        "0.5" = "sigma['q'] == 0.5")) %>%
      mutate(q_cor = recode(q_cor,
        "0" = "rho['q'] == 0",
        "0.9" = "rho['q'] == 0.9")) %>%
      mutate(EM = recode(EM,
        "1" = "R+S",
        "2" = "R+M",
        "3" = "R+Sel",
        "4" = "R+q(iid)",
        "5" = "R+q(AR1)"
      )) %>% as.data.frame
  }
  return(out)
}


df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))
em_ind <- c(5:20,25:28)
use.df.ems <- df.ems[em_ind,]
types <- c("naa", "M", "Sel", "q")
types.plt <- c("R, R+S","R+M", "R+Sel", "R+q")
est_conds <- cbind(c(FALSE,FALSE,TRUE,TRUE), c(FALSE,TRUE,FALSE,TRUE))
em_inds <- cbind(1:20, 5:24, c(5:20,25:28), c(5:20, 29:32))
for(i in 1:4) {
  all_aic <- readRDS(file = here("Project_0","results", paste0("all_", types[i], "_aic_results.RDS")))
  df.oms = readRDS(here("Project_0","inputs", paste0("df.", ifelse(i==1, "", paste0(types[i],".")),"oms.RDS")))
  if(i>1) re_nms <- paste0(types[i],c("_sig","_cor"))
  if(i ==1) re_nms <- c("R_sig","NAA_sig")
  for(j in 1:4){
    print(c(i,j))
    out <- make_plot_df(all_aic, df.ems, df.oms, em_ind = em_inds[,i], SR_est = est_conds[j,1], M_est = est_conds[j,2])
    totals <- out %>% group_by(obs_error, Fhist, get(re_nms[1]), get(re_nms[2])) %>% summarize(n = sum(n))
    plt <- ggplot(out, aes(x = obs_error, y = n, fill = EM)) + scale_fill_viridis_d() + 
        #facet_grid(get(re_nms[1]) + get(re_nms[2]) ~ Fhist) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
        facet_grid(get(re_nms[1]) + get(re_nms[2]) ~ Fhist, labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
        geom_col(position = "fill") + theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("Best EM proportion") + xlab("Level of observation error") +
        ggtitle(paste0(types.plt[i], " OMs: M ", ifelse(est_conds[j,2], "estimated", "= 0.2"), " and ", 
          ifelse(est_conds[j,1], "B-H ", "no SRR "), "assumed")) + theme(plot.title = element_text(hjust = 0.5)) + labs(fill = "EM Process Error")
        #geom_text(stat = "count", aes(label = ..count..), vjust = 1.5, colour = "black")
    plt
    ggsave(here("Project_0", "paper", paste0(types[i],"_om_proportion_best_aic_", ifelse(est_conds[j,1], "SR", "R"), "_", 
        ifelse(est_conds[j,2], "ME", "MF"),".png")), plt, width = 20, height = 12, units = "in")
    # png(here("Project_0", "paper", paste0(types[i],"_om_proportion_best_aic_", ifelse(est_conds[j,1], "SR", "R"), "_", 
    #     ifelse(est_conds[j,2], "ME", "MF"),".png")), width = 20, height = 12, units = "in", res = 300)
    # print(plt)
    # dev.off()

  }
}

out_all <- cbind.data.frame(Fhist = character(), obs_error = character(), n = integer(), OM = character(), EM = character(), M_est = logical(), SR_est = logical())
for(i in 1:4) {
  all_aic <- readRDS(file = here("Project_0","results", paste0("all_", types[i], "_aic_results.RDS")))
  df.oms = readRDS(here("Project_0","inputs", paste0("df.", ifelse(i==1, "", paste0(types[i],".")),"oms.RDS")))
  if(i>1) re_nms <- paste0(types[i],c("_sig","_cor"))
  if(i ==1) re_nms <- c("R_sig","NAA_sig")
  for(j in 1:4){
    #j <- 1 #SR = F, M = F
    out <- make_plot_df(all_aic, df.ems, df.oms, em_ind = em_inds[,i], SR_est = est_conds[j,1], M_est = est_conds[j,2])
    out <- cbind(out, M_est = est_conds[j,2], SR_est = est_conds[j,1])
    if(i ==1) {
      out$OM <- "R"
      out$OM[out$NAA_sig != unique(out$NAA_sig)[1]] <- "R+S"
    }
    if(i !=1) {
      if(i ==2) {
        out$OM <- "R+M"
        out$EM[which(out$EM %in% c("R+M(iid)", "R+M(AR1)"))] <- "R+M"
      }
      if(i ==3) {
        out$OM <- "R+Sel"
        out$EM[which(out$EM %in% c("R+Sel(iid)", "R+Sel(AR1)"))] <- "R+Sel"
      }
      if(i ==4) {
        out$OM <- "R+q"
        out$EM[which(out$EM %in% c("R+q(iid)", "R+q(AR1)"))] <- "R+q"
      }
    }
    out_all <- rbind(out_all, out[c("Fhist","obs_error","n","OM","EM", "M_est", "SR_est")])
  }
}
out_all$obs_error <- factor(out_all$obs_error, levels = c("L","H"))
out_all$EM <- factor(out_all$EM, levels = c("R","R+S", "R+M", "R+Sel", "R+q"))
out_all$OM <- factor(out_all$OM, levels = c("R","R+S", "R+M", "R+Sel", "R+q"))

for(j in 1:4){
  temp <- out_all %>% filter(M_est == est_conds[j,2] & SR_est == est_conds[j,1]) %>% group_by(Fhist,obs_error,OM, EM) %>% summarize(N = sum(n)) %>% as.data.frame
  #for(i in c("Fhist","obs_error","OM","EM")) temp[i] <- factor(temp[i])
  plt <- ggplot(temp, aes(x = obs_error, y = N, fill = EM)) + scale_fill_viridis_d() + 
      facet_grid(Fhist ~ OM) + #, labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
      geom_col(position = "fill") + theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("Best EM proportion") + xlab("Level of observation error") +
      ggtitle(paste0("M ", ifelse(est_conds[j,2], "estimated", "= 0.2"), " and ", 
        ifelse(est_conds[j,1], "B-H ", "no SRR "), "assumed")) + 
      theme(plot.title = element_text(hjust = 0.5)) + labs(fill = "EM Process Error")
  plt
  ggsave(here("Project_0", "paper", paste0("om_em_proportion_best_aic_coarse_", ifelse(est_conds[j,1], "SR", "R"), "_", 
      ifelse(est_conds[j,2], "ME", "MF"),".png")), plt, width = 20, height = 12, units = "in")
}

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






############################################
#AIC weights
############################################

make_plot_df <- function(aic_res, df.ems, df.oms, em_ind = 1:20, SR_est =FALSE, M_est = FALSE){
  use.df.ems <- df.ems[em_ind,]
  SR_M_ind = which(use.df.ems$SR_model == ifelse(SR_est,3,2) & use.df.ems$M_est == M_est)
  aic_res <- t(aic_wt_fn(aic_res, SR_M_ind))
  print(aic_res)
  cnames <- paste0(rep(1:length(SR_M_ind), each = 3), "_", c("lo", "med", "hi"))
  colnames(aic_res) <- cnames
  out <- cbind(df.oms[-1], aic_res)
  out <- out %>% pivot_longer(
    cols = all_of(cnames),
    names_to = c("EM", "Type"),
    names_pattern = "(.)_(.*)" #magic
  )
  out <- out %>% pivot_wider(names_from = Type, values_from = value) %>% as.data.frame

  print(head(out))
  if(max(em_ind)==20) {
    out$NAA_sig[which(is.na(out$NAA_sig))] <- 0
    out <- out %>% mutate(NAA_sig = recode(NAA_sig,
        "0" = "sigma['2+'] == 0",
        "0.25" = "sigma['2+'] == 0.25",
        "0.5" = "sigma['2+'] == 0.5")) %>% as.data.frame
    out <- out %>% mutate(R_sig = recode(R_sig,
        "0.5" = "sigma[R] == 0.5",
        "1.5" = "sigma[R] == 1.5")) %>% as.data.frame
    out <- out %>% mutate(EM = recode(EM,
        "1" = "R",
        "2" = "R+S",
        "3" = "R+M",
        "4" = "R+Sel",
        "5" = "R+q"
      )) %>% as.data.frame
  }
  if(max(em_ind)==24) {
    out <- out %>% 
      mutate(M_sig = recode(M_sig,
        "0.1" = "sigma['M'] == 0.1",
        "0.5" = "sigma['M'] == 0.5")) %>% 
      mutate(M_cor = recode(M_cor,
        "0" = "rho['M'] == 0",
        "0.9" = "rho['M'] == 0.9")) %>%
      mutate(EM = recode(EM,
        "1" = "R+S",
        "2" = "R+M(iid)",
        "3" = "R+Sel",
        "4" = "R+q",
        "5" = "R+M(AR1)"
      )) %>% as.data.frame
  }
  if(max(em_ind)==28) {
    out <- out %>% 
      mutate(Sel_sig = recode(Sel_sig,
        "0.1" = "sigma['Sel'] == 0.1",
        "0.5" = "sigma['Sel'] == 0.5")) %>% 
      mutate(Sel_cor = recode(Sel_cor,
        "0" = "rho['Sel'] == 0",
        "0.9" = "rho['Sel'] == 0.9")) %>%
      mutate(EM = recode(EM,
        "1" = "R+S",
        "2" = "R+M",
        "3" = "R+Sel(iid)",
        "4" = "R+q",
        "5" = "R+Sel(AR1)"
      )) %>% as.data.frame
  }
  if(max(em_ind)==32) {
    out <- out %>% 
      mutate(q_sig = recode(q_sig,
        "0.1" = "sigma['q'] == 0.1",
        "0.5" = "sigma['q'] == 0.5")) %>%
      mutate(q_cor = recode(q_cor,
        "0" = "rho['q'] == 0",
        "0.9" = "rho['q'] == 0.9")) %>%
      mutate(EM = recode(EM,
        "1" = "R+S",
        "2" = "R+M",
        "3" = "R+Sel",
        "4" = "R+q(iid)",
        "5" = "R+q(AR1)"
      )) %>% as.data.frame
  }
  return(out)
}

custom_boxplot_stat <- function(x){#, n.yrs = 1, n.sim = 100) {
  x <- x[which(!is.na(x))]
  n <- length(x)
  bnds95 <- qbinom(c(0.025,0.975), n, 0.5)/n # 95% CI bounds for median (sims x years)
  #bnds80 <- qbinom(c(0.1,0.9), n, 0.5)/n # 80% CI bounds for median (sims x years)
  #r <- quantile(x, probs = c(bnds95[1], bnds80[1], 0.5, bnds80[2], bnds95[2]))
  r <- quantile(x, probs = c(bnds95[1], 0.5, bnds95[2]))
  names(r) <- c("lo", "med", "hi")
  r
}

aic_wt_fn <- function(all, est_ind, om_ind = NULL){
  if(!is.null(om_ind)) all <- all[om_ind]
  out <- sapply(all, function(res){
    tmp <- apply(res[est_ind,],2, function(x) return(x))
    tmp <- apply(tmp,2, function(x) {
      if(any(!is.na(x))) {
        out <- rep(0, length(x))
        out <- exp(-0.5*(x-min(x,na.rm =TRUE)))
        out[which(is.na(out))] <- 0
        out <- out/sum(out)
        return(out)
      } else return(rep(NA, length(x)))
    })
    return(apply(tmp,1,custom_boxplot_stat))
  })  
  return(out)
}
types <- c("naa", "M", "Sel", "q")
types.plt <- c("R, R+S","R+M", "R+Sel", "R+q")
est_conds <- cbind(c(FALSE,FALSE,TRUE,TRUE), c(FALSE,TRUE,FALSE,TRUE))
em_inds <- cbind(1:20, 5:24, c(5:20,25:28), c(5:20, 29:32))
SR_M_ind = which(use.df.ems$SR_model == ifelse(SR_est,3,2) & use.df.ems$M_est == M_est)
df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))
for(i in 1:4) {
  all_aic <- readRDS(file = here("Project_0","results", paste0("all_", types[i], "_aic_results.RDS")))
  df.oms = readRDS(here("Project_0","inputs", paste0("df.", ifelse(i==1, "", paste0(types[i],".")),"oms.RDS")))
  if(i>1) re_nms <- paste0(types[i],c("_sig","_cor"))
  if(i ==1) re_nms <- c("R_sig","NAA_sig")
  for(j in 1:4){
    print(c(i,j))
    out <- make_plot_df(all_aic, df.ems, df.oms, em_ind = em_inds[,i], SR_est = est_conds[j,1], M_est = est_conds[j,2])
    #totals <- out %>% group_by(obs_error, Fhist, get(re_nms[1]), get(re_nms[2])) %>% summarize(n = sum(n))
plt <- ggplot(out, aes(x = obs_error, y = med, colour = EM))  + scale_colour_viridis_d() + 
    geom_point(position = position_dodge(0.1), size = 4) + #geom_line(position = position_dodge(0.1), linewidth = 1) + 
    facet_grid(get(re_nms[1]) + get(re_nms[2]) ~ Fhist, labeller = label_parsed) + #label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(0, 1)) + ylab("AIC weight") + xlab("Level of observation error") +
    ggtitle(paste0(types.plt[i], " OMs: M ", ifelse(est_conds[j,2], "estimated", "= 0.2"), " and ", 
      ifelse(est_conds[j,1], "B-H ", "no SRR "), "assumed")) + theme(plot.title = element_text(hjust = 0.5)) + labs(fill = "EM Process Error") +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, position = position_dodge(0.1), linewidth = 1) 
  plt
    ggsave(here("Project_0", "paper", paste0(types[i],"_om_aic_weights_", ifelse(est_conds[j,1], "SR", "R"), "_", 
        ifelse(est_conds[j,2], "ME", "MF"),".png")), plt, width = 8, height = 12, units = "in")

  }
}

