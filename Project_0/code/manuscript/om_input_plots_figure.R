library(here)
library(dplyr)
library(wham, lib.loc = "c:/work/wham/old_packages/77bbd94")
df.oms <- readRDS(file.path(here(),"Project_0","inputs", "df.oms.RDS"))
om_inputs <- readRDS(file.path(here(),"Project_0","inputs", "NAA_om_inputs.RDS"))
x <- fit_wham(om_inputs[[1]], do.fit =F)
library(ggplot2)
library(patchwork)
theme_set(theme_bw())
theme_update(strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size = rel(2)), 
      axis.title = element_text(size = rel(2)), axis.text = element_text(size = rel(2)), legend.text = element_text(size = rel(2)), #text = element_text(size = rel(2)), 
      legend.title = element_text(size = rel(2)))

df <- cbind.data.frame(age_lab = x$ages, age = 1:om_inputs[[1]]$data$n_ages, sel = x$rep$selAA[[1]][1,], mat = x$input$data$mature[1,], wgt = x$input$data$waa[1,1,])

sel_plt <- ggplot(df, aes(x = age, y = sel)) + geom_line() + ylab("Selectivity") + xlab("Age") + scale_x_continuous(breaks=seq(1,10,1),labels=df$age_lab)
sel_plt

mat_plt <- ggplot(df, aes(x = age, y = mat)) + geom_line() + ylab("Proportion Mature") + xlab("Age") + scale_x_continuous(breaks=seq(1,10,1),labels=df$age_lab)
mat_plt

wgt_plt <- ggplot(df, aes(x = age, y = wgt)) + geom_line() + ylab("Mass (kg)") + xlab("Age") + scale_x_continuous(breaks=seq(1,10,1),labels=df$age_lab)
wgt_plt

png(here("Project_0","paper","om_sr.png"), width=7, height=7, units='in', res = 300)
sr.fn <- function(ssb) exp(x$rep$log_SR_a[1])*ssb/(1+exp(x$rep$log_SR_b[1])*ssb)

ssb <- seq(0,10e5,10)
df <- cbind.data.frame(ssb = ssb/1000, R = sr.fn(ssb)/1000)
sr_plt <- ggplot(df, aes(x = ssb, y = R)) + geom_line() + ylab(expression(Recruitment~(10^6))) + xlab("Spawning biomass (kmt)")
sr_plt

cairo_pdf(here("manuscript", "om_input_plots_figure.pdf"), width = 20, height = 20)
mat_plt+wgt_plt+sel_plt+sr_plt
dev.off()
