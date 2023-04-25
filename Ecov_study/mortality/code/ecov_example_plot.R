library(here)
library(ggplot2)
library(dplyr)
df.ems = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.ems.RDS"))
df.oms = readRDS(file.path(here(),"Ecov_study","mortality","inputs", "df.oms.RDS"))
om_inputs <- readRDS(file.path(here(),"Ecov_study","mortality","inputs", "om_inputs.RDS"))

Ecov_obs_sigs <- unique(df.oms$Ecov_obs_sig)
Ecov_re_sigs <- unique(df.oms$Ecov_re_sig)
Ecov_re_cors <- unique(df.oms$Ecov_re_cor)

Ecov_types <- expand.grid(Ecov_obs_sig = Ecov_obs_sigs,Ecov_re_sig = Ecov_re_sigs,Ecov_re_cor =Ecov_re_cors)
#par(mfcol = c(4,2), oma = c(5,5,0,0), mar = c(1,1,1,1))
ecov_df <- data.frame(Ecov_obs_sig = numeric(), Ecov_re_sig = numeric(), Ecov_re_cor = numeric(), Ecov = numeric(), type = character())
cn <- colnames(ecov_df)
for(i in 1:NROW(Ecov_types)){
   om <- which(df.oms$Ecov_obs_sig == Ecov_types[i,1] & df.oms$Ecov_re_sig == Ecov_types[i,2] & df.oms$Ecov_re_cor == Ecov_types[i,3])[1]
   print(om)
   res <- readRDS(file.path(here::here(),"Ecov_study","mortality", "results", "aggregated_results", paste0("om_", om, "_em_", 1, ".RDS")))
   good_sims <- which(sapply(res, length) == 2)
   tmp <- cbind.data.frame(Ecov_types[i,1], Ecov_types[i,2], Ecov_types[i,3], c(res[[good_sims[1]]]$truth$Ecov_x,res[[good_sims[1]]]$truth$Ecov_obs),
    rep(c("True", "Observation"), each = 40))
   colnames(tmp) <- cn
   ecov_df <- rbind(ecov_df, tmp)
}
ecov_df = cbind(ecov_df, Year = om_inputs[[1]]$years)
facs <- c("Ecov_obs_sig", "Ecov_re_sig","Ecov_re_cor", "type")
ecov_df[facs] <- lapply(all_res[facs], factor)
ecov_df_mod <- ecov_df %>%
  mutate(Ecov_obs_sig = recode(Ecov_obs_sig,
    "0.1" = "Ecov obs SD = 0.1",
    "0.5" = "Ecov obs SD = 0.5"
  ))
ecov_df_mod <- ecov_df_mod %>%
  mutate(Ecov_re_sig = recode(Ecov_re_sig,
    "0.1" = "Ecov AR1 SD = 0.1",
    "0.5" = "Ecov AR1 SD = 0.5"
  ))
ecov_df_mod <- ecov_df_mod %>%
  mutate(Ecov_re_cor = recode(Ecov_re_cor,
    "0" = "Ecov is iid",
    "0.5" = "Ecov AR1 Cor = 0.5"
  ))

plt <- ggplot(ecov_df_mod, aes(x = Year, y = Ecov, colour = type)) + 
    geom_line(position = position_dodge(0.1), linewidth = 1) + geom_point(position = position_dodge(0.1), size = 4) + 
    facet_grid(Ecov_obs_sig ~ Ecov_re_sig + Ecov_re_cor, labeller = label_wrap_gen(width = 30)) + #, labeller = label_parsed) + 
    theme_bw() + coord_cartesian(ylim = c(-1, 1)) + ylab("Environmental Covariate") + xlab("Year") #+
plt
ggsave(here("Ecov_study","mortality", "paper", "Ecov_true_obs_example.png"), plt)
