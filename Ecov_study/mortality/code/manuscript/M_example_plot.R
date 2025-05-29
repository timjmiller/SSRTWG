library(here)
library(ggplot2)
library(ggh4x)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
df.ems = readRDS(here::here("Ecov_study","mortality","inputs", "df.ems.RDS"))
df.oms = readRDS(here::here("Ecov_study","mortality","inputs", "df.oms.RDS"))
om_inputs <- readRDS(here::here("Ecov_study","mortality","inputs", "om_inputs.RDS"))

Ecov_obs_sigs <- unique(df.oms$Ecov_obs_sig)
Ecov_re_sigs <- unique(df.oms$Ecov_re_sig)
Ecov_re_cors <- unique(df.oms$Ecov_re_cor)
Ecov_effect <- unique(df.oms$Ecov_effect)
PE <- c("rec", "rec+M")
temp <- readRDS(here::here("Ecov_study","mortality", "results", "aggregated_results", paste0("om_", 1, "_em_", 1, ".RDS")))

Ecov_types <- expand.grid(Ecov_re_sig = Ecov_re_sigs,Ecov_re_cor =Ecov_re_cors, Ecov_effect = Ecov_effect, process_error = PE)
#par(mfcol = c(4,2), oma = c(5,5,0,0), mar = c(1,1,1,1))
M_df <- data.frame(Ecov_re_sig = numeric(), Ecov_re_cor = numeric(), Ecov_effect = numeric(), process_error = character(), M = numeric())
cn <- colnames(M_df)
for(i in 1:NROW(Ecov_types)){
   om <- which(df.oms$Ecov_re_sig == Ecov_types[i,1] & df.oms$Ecov_re_cor == Ecov_types[i,2] & df.oms$Ecov_effect == Ecov_types[i,3] & df.oms$NAA_M_re == Ecov_types[i,4])[1]
   print(om)
   res <- readRDS(file.path(here::here(),"Ecov_study","mortality", "results", "aggregated_results", paste0("om_", om, "_em_", 1, ".RDS")))
   good_sims <- which(sapply(res, length) == 2)
   tmp <- cbind.data.frame(Ecov_types[i,1], Ecov_types[i,2], Ecov_types[i,3], Ecov_types[i,4], res[[good_sims[1]]]$truth$MAA[,1])
   colnames(tmp) <- cn
   M_df <- rbind(M_df, tmp)
}
M_df = cbind(M_df, Year = om_inputs[[1]]$years)
facs <- c("Ecov_re_sig","Ecov_re_cor", "Ecov_effect", "process_error")
M_df[facs] <- lapply(M_df[facs], factor)
M_df_mod <- M_df %>%
  mutate(M_re_present = recode(process_error,
    "rec" = "without M process error",
    "rec+M" = "with M process error"
  ))
M_df_mod <- M_df_mod %>%
  mutate(Ecov_effect = recode(Ecov_effect,
    "0" = "beta[E] == 0",
    "0.25" = "beta[E] == 0.25",
    "0.5" = "beta[E] == 0.5"
  ))
M_df_mod <- M_df_mod %>%
  mutate(Ecov_re_sig = recode(Ecov_re_sig,
                              "0.1" = "sigma[italic(E)] == 0.1",
                              "0.5" = "sigma[italic(E)] == 0.5"
  ))
M_df_mod <- M_df_mod %>%
  mutate(Ecov_re_cor = recode(Ecov_re_cor,
                              "0" = "rho[italic(E)] == 0",
                              "0.5" = "rho[italic(E)] == 0.5"
  ))

theme_update(strip.text.x = element_text(size = rel(2)), strip.text.y = element_text(size = rel(2)), strip.placement = "outside", strip.background = element_rect(), #fill = "transparent"), 
             axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(1.5)), legend.text = element_text(size = rel(2)), 
             legend.title = element_text(size = rel(2)))

plt <- ggplot(M_df_mod, aes(x = Year, y = M)) + 
  scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "turbo", drop = FALSE) +
  geom_line(position = position_dodge(0.1), linewidth = 1, alpha = 0.7) + geom_point(position = position_dodge(0.1), size = 4, alpha = 0.7) + 
  facet_nested(Ecov_re_sig + Ecov_re_cor ~ M_re_present + Ecov_effect, labeller = labeller(Ecov_effect = label_parsed, Ecov_re_sig = label_parsed, Ecov_re_cor = label_parsed)) + 
  scale_x_continuous(breaks = seq(1985,2015,10)) +
  ylab(expression(italic(M))) + xlab("Year")
plt
cairo_pdf(here::here("Ecov_study","mortality","manuscript", "M_example.pdf"), width = 30*2/3, height = 20*2/3)
plt
dev.off()
