mohns_rho_fn <- function(df.oms, om_type = "naa"){
  all_out <- lapply(1:NROW(df.oms), function(y){
    res <- lapply(1:100, function(x){
      print(paste0("om_", y, ", sim_",x))
      sim = readRDS(file.path(here::here(),"Project_0", "results", paste0(om_type,"_om"), paste0("om_", y), paste0("sim_",x,".RDS")))
      simres <- lapply(sim,function(z) {
        out <-rep(NA,3)
        if(length(z)) if(length(z$fit)) {
          if(!is.null(z$fit$mohns_rho)) if(!is.character(z$fit$mohns_rho)) out <- z$fit$mohns_rho[1:3]
        }
        return(out)
      })
      return(simres)
    })
    return(res)
  })
  return(all_out)
}
make_results <- FALSE
if(make_results) {
  library(here)
  df.oms <- readRDS(here("Project_0","inputs", "df.oms.RDS"))
  all_naa_mohns_rho <- mohns_rho_fn(df.oms,"naa")
  saveRDS(all_naa_mohns_rho, file = here("Project_0","results", "all_naa_mohns_rho_results.RDS"))
  df.M.oms <- readRDS(here("Project_0","inputs", "df.M.oms.RDS"))
  all_M_mohns_rho <- mohns_rho_fn(df.M.oms,"M")
  saveRDS(all_M_mohns_rho, file = here("Project_0","results", "all_M_mohns_rho_results.RDS"))
  df.Sel.oms <- readRDS(here("Project_0","inputs", "df.Sel.oms.RDS"))
  all_Sel_mohns_rho <- mohns_rho_fn(df.Sel.oms,"Sel")
  saveRDS(all_Sel_mohns_rho, file = here("Project_0","results", "all_Sel_mohns_rho_results.RDS"))
  df.q.oms <- readRDS(here("Project_0","inputs", "df.q.oms.RDS"))
  all_q_mohns_rho <- mohns_rho_fn(df.q.oms,"q")
  saveRDS(all_q_mohns_rho, file = here("Project_0","results", "all_q_mohns_rho_results.RDS"))
}
custom_boxplot_stat <- function(x){#, n.yrs = 1, n.sim = 100) {
  x <- x[which(!is.na(x))]
  n <- length(x)
  bnds95 <- qbinom(c(0.025,0.975), n, 0.5)/n # 95% CI bounds for median (sims x years)
  bnds80 <- qbinom(c(0.1,0.9), n, 0.5)/n # 80% CI bounds for median (sims x years)
  r <- quantile(x, probs = c(bnds95[1], bnds80[1], 0.5, bnds80[2], bnds95[2]))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
get_df_ems <- function(type = "naa"){
  df.ems = readRDS(here::here("Project_0","inputs", "df.ems.RDS"))
  df.ems$pe <- c("R","R+S","R+M","R+Sel","R+q")[match(df.ems$re_config, c("rec","rec+1","M_re", "sel_re","q_re"))]
  em_ind <- 1:20
  if(type != "naa"){
    if(type == "M") {
      em_ind <- 5:24
      ind <- which(df.ems$re_config =="M_re")
      df.ems$pe[ind] <- paste0(df.ems$pe[ind], "(", df.ems$M_re_cor[ind], ")") 
    }
    if(type == "Sel") {
      em_ind <- c(5:20,25:28)
      ind <- which(df.ems$re_config =="sel_re")
      df.ems$pe[ind] <- paste0(df.ems$pe[ind], "(", df.ems$sel_re_cor[ind], ")") 
    }
    if(type == "q") {
      em_ind <- c(5:20,29:32)
      ind <- which(df.ems$re_config =="q_re")
      df.ems$pe[ind] <- paste0(df.ems$pe[ind], "(", df.ems$q_re_cor[ind], ")") 
    }
  }
  df.ems <- df.ems[em_ind,]
  df.ems$meanR <- "Mean R"
  df.ems$meanR[df.ems$SR_model ==3] <- "B-H"
  return(df.ems)
}
naa_om_df.ems <- get_df_ems()
M_om_df.ems <- get_df_ems("M")
Sel_om_df.ems <- get_df_ems("Sel")
q_om_df.ems <- get_df_ems("q")

temp <- melt()
make_plot_df <- function(type = "naa", res = all_naa_om_relssb) {
  library(reshape2)
  res <- melt(res)
  names(res) <- c("rho", "em", "sim","om")
  res$Type = c("SSB","F","R")
  # res <- res %>% mutate(column = recode(column,
  #     "1" = "relerror",
  #     "2" = "cv",
  #     "3" = "in_ci",
  #   )) %>% as.data.frame   
  # df <- res %>% pivot_wider(names_from = column, values_from = value) %>% as.data.frame
  # print(dim(df))
  # if(is_SE) df <- filter(df, !is.na(cv))
  # print(dim(df))
  # df$relerror = df$relerror - 1

  df <- res %>% group_by(om, em, Type) %>%
    reframe(stats = custom_boxplot_stat(rho)) %>% as.data.frame
  df$type <- c("ymin", "lower", "middle", "upper", "ymax")
  library(tidyr)
  df.ems <- get_df_ems(type)
  if(type == "naa") {
    df.oms <- readRDS(here::here("Project_0", "inputs", "df.oms.RDS"))
  } else {
    df.oms <- readRDS(here::here("Project_0", "inputs", paste0("df.", type, ".oms.RDS")))
  }
  print(head(df.oms))
  df <- df %>% pivot_wider(names_from = type, values_from = stats) %>% as.data.frame
  df$em_pe = df.ems$pe#[df$em]
  df <- cbind(df, df.ems[df$em,])
  df <- cbind(df, df.oms[df$om,])
  if(type == "naa") {
    df$NAA_sig[which(is.na(df$NAA_sig))] <- 0
    df <- df %>% mutate(NAA_sig = recode(NAA_sig,
        "0" = "sigma['2+'] == 0",
        "0.25" = "sigma['2+'] == 0.25",
        "0.5" = "sigma['2+'] == 0.5")) %>%
      mutate(R_sig = recode(R_sig,
        "0.5" = "sigma[R] == 0.5",
        "1.5" = "sigma[R] == 1.5"))
  }
  if(type == "M") {
    df <- df %>% mutate(M_sig = recode(M_sig,
        "0.1" = "sigma['M'] == 0.1",
        "0.5" = "sigma['M'] == 0.5")) %>% 
      mutate(M_cor = recode(M_cor,
        "0" = "rho['M'] == 0",
        "0.9" = "rho['M'] == 0.9"))
  }
  if(type == "Sel") {
    df <- df %>% mutate(Sel_sig = recode(Sel_sig,
        "0.1" = "sigma['Sel'] == 0.1",
        "0.5" = "sigma['Sel'] == 0.5")) %>% 
      mutate(Sel_cor = recode(Sel_cor,
        "0" = "rho['Sel'] == 0",
        "0.9" = "rho['Sel'] == 0.9"))
  }
  if(type == "q") {
    df <- df %>% mutate(q_sig = recode(q_sig,
        "0.1" = "sigma['q'] == 0.1",
        "0.5" = "sigma['q'] == 0.5")) %>%
      mutate(q_cor = recode(q_cor,
        "0" = "rho['q'] == 0",
        "0.9" = "rho['q'] == 0.9"))  
  }
  df <- df %>% mutate(obs_error = recode(obs_error,
      "L" = "Low obs error (indices, age comp)",
      "H" = "High obs error (indices, age comp)"))
  df <- df %>% mutate(Fhist = recode(Fhist,
      "H-MSY" = "F history: High->FMSY",
      "MSY" = "F history: FMSY"))
  df <- df %>% as.data.frame
  facs <- c("Fhist","NAA_sig", "R_sig", "M_sig", "M_cor", "Sel_sig", "Sel_cor", "q_sig", "q_cor", "obs_error")
  df[names(df) %in% facs] <- lapply(df[names(df) %in% facs], as.factor)
  return(df)
}

all_naa_mohns_rho <-  readRDS(file = here("Project_0","results", "all_naa_mohns_rho_results.RDS"))
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
df <- make_plot_df(type = "naa", res = all_naa_mohns_rho)
# length(all_naa_mohns_rho) #om
# length(all_naa_mohns_rho[[1]]) #sim
# length(all_naa_mohns_rho[[1]][[1]]) #em
# length(all_naa_mohns_rho[[1]][[1]][[1]]) #rho ssb, f, r
# sim = readRDS(here("Project_0", "results", paste0("naa","_om"), paste0("om_", 1), paste0("sim_",1,".RDS")))

types <- c("naa", "M", "Sel", "q")
types.plt <- c("R, R+S","R+M", "R+Sel", "R+q")
#est_conds <- cbind(c(FALSE,FALSE,TRUE,TRUE), c(FALSE,TRUE,FALSE,TRUE))
estM <- c(FALSE,FALSE,TRUE,TRUE)
estSR<- c("Mean R", "B-H", "Mean R", "B-H")
em_inds <- cbind(1:20, 5:24, c(5:20,25:28), c(5:20, 29:32))
df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))
for(i in 1:4) {
  all_mohns_rho <- readRDS(file = here("Project_0","results", paste0("all_", types[i], "_mohns_rho_results.RDS")))
  df.oms = readRDS(here("Project_0","inputs", paste0("df.", ifelse(i==1, "", paste0(types[i],".")),"oms.RDS")))
  if(i>1) re_nms <- paste0(types[i],c("_sig","_cor"))
  if(i ==1) re_nms <- c("R_sig","NAA_sig")
  out <- make_plot_df(type = types[i], res = all_mohns_rho)
  names(out)[which(names(out) %in% re_nms)] <- c("re_name_1", "re_name_2")
  for(j in 1:4){
    print(c(i,j))
  plt <- ggplot(filter(out, M_est== estM[j] & meanR == estSR[j]), aes(x = Type, y = middle, colour = pe))  + scale_colour_viridis_d() + 
      geom_hline(aes(yintercept=0), linewidth = 1, linetype = "dashed", colour = "black") +
      geom_point(position = position_dodge(0.5), size = 4) + 
      facet_grid(re_name_1 + re_name_2 ~ Fhist + obs_error, labeller = labeller(re_name_1 = label_parsed, re_name_2 = label_parsed)) + 
      theme_bw() + coord_cartesian(ylim = c(-0.5, 0.5)) + ylab("Mohn's rho") + xlab("Population attribute") +
      ggtitle(paste0(types.plt[i], " OMs: M ", ifelse(estM[j], "estimated", "= 0.2"), " and ", 
        ifelse(estSR[j]=="B-H", "BH assumed", "no SRR"))) + theme(plot.title = element_text(hjust = 0.5)) + labs(color = "EM Process Error") +
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.25, position = position_dodge(0.5), linewidth = 1) 
  plt
    ggsave(here("Project_0", "paper", paste0(types[i],"_om_mohns_rho_", ifelse(estSR[j]=="B-H","BH","R"), "_", 
        ifelse(estM[j], "ME", "MF"),".png")), plt, width = 12, height = 8, units = "in")

  }
}

