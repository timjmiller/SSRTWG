library(here)
library(dplyr)
library(Hmisc)
df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))

#NAA oms: ems = 1-20
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.oms.RDS"))
om_inputs = readRDS(file.path(here(),"Project_0","inputs", "NAA_om_inputs.RDS"))

naa_om_tab <- df.oms
naa_om_tab$Model <- paste0("NAA_", 1:24)
naa_om_tab <- naa_om_tab %>%
  mutate(obs_error = recode(obs_error,
    "L" = "Index sigma (log scale) = 0.1, Age composition sigma (logistic normal) = 0.3",
    "H" = "Index sigma (log scale) = 0.4, Age composition sigma (logistic normal) = 1.5"
  ))
naa_om_tab <- naa_om_tab %>%
  mutate(Fhist = recode(Fhist,
    "H-MSY" = "$2.5 F_{\\text{MSY}} \\rightarrow F_{\\text{MSY}}$",
    "MSY" = "F_{\\text{MSY}}"
  ))

names(naa_om_tab)[2] <- "$\\sigma_R$"
names(naa_om_tab)[3] <- "$\\sigma_{2+}$"
names(naa_om_tab)[4] <- "Fishing History"
names(naa_om_tab)[5] <- "Observation Uncertainty"

x = latex(naa_om_tab, file = here("Project_0","paper","naa_om_tab.tex"), 
  table.env = FALSE, col.just = rep("r", dim(naa_om_tab)[2]), rowname = NULL)


#M oms: ems = 5-24
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.M.oms.RDS"))
om_inputs = readRDS(file.path(here(),"Project_0","inputs", "M_om_inputs.RDS"))

om_tab <- df.oms
om_tab$Model <- paste0("M_", 1:16)
om_tab <- om_tab %>%
  mutate(obs_error = recode(obs_error,
    "L" = "Index sigma (log scale) = 0.1, Age composition sigma (logistic normal) = 0.3",
    "H" = "Index sigma (log scale) = 0.4, Age composition sigma (logistic normal) = 1.5"
  ))
om_tab <- om_tab %>%
  mutate(Fhist = recode(Fhist,
    "H-MSY" = "$2.5 F_{\\text{MSY}} \\rightarrow F_{\\text{MSY}}$",
    "MSY" = "F_{\\text{MSY}}"
  ))

names(om_tab)[2] <- "$\\sigma_R$"
names(om_tab)[3] <- "$\\sigma_{M}$"
names(om_tab)[4] <- "$\\rho_{M}$"
names(om_tab)[5] <- "Fishing History"
names(om_tab)[6] <- "Observation Uncertainty"

x = latex(om_tab, file = here("Project_0","paper","M_om_tab.tex"), 
  table.env = FALSE, col.just = rep("r", dim(om_tab)[2]), rowname = NULL)


#Sel oms: ems = 5-20, 25-28
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.Sel.oms.RDS"))
om_inputs = readRDS(file.path(here(),"Project_0","inputs", "Sel_om_inputs.RDS"))

om_tab <- df.oms
om_tab$Model <- paste0("Sel_", 1:16)
om_tab <- om_tab %>%
  mutate(obs_error = recode(obs_error,
    "L" = "Index sigma (log scale) = 0.1, Age composition sigma (logistic normal) = 0.3",
    "H" = "Index sigma (log scale) = 0.4, Age composition sigma (logistic normal) = 1.5"
  ))
om_tab <- om_tab %>%
  mutate(Fhist = recode(Fhist,
    "H-MSY" = "$2.5 F_{\\text{MSY}} \\rightarrow F_{\\text{MSY}}$",
    "MSY" = "F_{\\text{MSY}}"
  ))

names(om_tab)[2] <- "$\\sigma_R$"
names(om_tab)[3] <- "$\\sigma_{M}$"
names(om_tab)[3] <- "$\\rho_{M}$"
names(om_tab)[4] <- "Fishing History"
names(om_tab)[5] <- "Observation Uncertainty"

x = latex(om_tab, file = here("Project_0","paper","M_om_tab.tex"), 
  table.env = FALSE, col.just = rep("r", dim(om_tab)[2]), rowname = NULL)
