library(here)
library(ggplot2)
library(ggh4x)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
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
ecov_df[facs] <- lapply(ecov_df[facs], factor)
ecov_df_mod <- ecov_df %>%
    mutate(Ecov_obs_sig = recode(Ecov_obs_sig,
      "0.1" = "sigma[italic(e)] == 0.1",
      "0.5" = "sigma[italic(e)] == 0.5"
    ))
ecov_df_mod <- ecov_df_mod %>%
  mutate(Ecov_re_sig = recode(Ecov_re_sig,
                              "0.1" = "sigma[italic(E)] == 0.1",
                              "0.5" = "sigma[italic(E)] == 0.5"
  ))
ecov_df_mod <- ecov_df_mod %>%
  mutate(Ecov_re_cor = recode(Ecov_re_cor,
                              "0" = "rho[italic(E)] == 0",
                              "0.5" = "rho[italic(E)] == 0.5"
  ))

theme_set(theme_bw())
theme_update(strip.text.x = element_text(size = rel(2)), strip.text.y = element_text(size = rel(2)), strip.placement = "outside", strip.background = element_rect(), #fill = "transparent"), 
             axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5)), legend.text = element_text(size = rel(2)), 
             legend.title = element_text(size = rel(2)))

plt <- ggplot(ecov_df_mod, aes(x = Year, y = Ecov, colour = type)) + 
    scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
    geom_line(position = position_dodge(0.1), linewidth = 1, alpha = 0.7) + 
    geom_point(position = position_dodge(0.1), size = 4, alpha = 0.7) + 
    facet_nested(Ecov_obs_sig ~ Ecov_re_sig + Ecov_re_cor, labeller = label_parsed) +
    labs(colour = element_blank()) + 
    scale_x_continuous(breaks = seq(1985,2015,10)) +
    coord_cartesian(ylim = c(-1, 1)) + ylab("Covariate") + xlab("Year") +
    guides(col = guide_legend(override.aes = list(shape = 15, size = 10, linetype = 0)))
plt
cairo_pdf(here::here("Ecov_study","mortality", "manuscript", "Ecov_true_obs_example.pdf"), width = 30*2/3, height = 20*2/3)
plt
dev.off()

ggsave(here("Ecov_study","mortality", "manuscript", "Ecov_true_obs_example.png"), plt)
