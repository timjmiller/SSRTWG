library(here)
library(dplyr)
library(Hmisc)
df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))

#NAA oms: ems = 1-20
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.oms.RDS"))
om_inputs = readRDS(file.path(here(),"Project_0","inputs", "NAA_om_inputs.RDS"))

naa_om_tab <- df.oms
naa_om_tab$Model <- paste0("$\\text{NAA}_{", 1:24, "}$")
naa_om_tab <- naa_om_tab %>%
  mutate(obs_error = recode(obs_error,
    "L" = "Index SD = 0.1, Age composition SD = 0.3",
    "H" = "Index SD = 0.4, Age composition SD = 1.5"
  ))
naa_om_tab <- naa_om_tab %>%
  mutate(Fhist = recode(Fhist,
    "H-MSY" = "$2.5 F_{\\text{MSY}} \\rightarrow F_{\\text{MSY}}$",
    "MSY" = "$F_{\\text{MSY}}$"
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
om_tab$Model <- paste0("$M_{", 1:16, "}$")
om_tab <- om_tab %>%
  mutate(obs_error = recode(obs_error,
    "L" = "Index SD = 0.1, Age composition SD = 0.3",
    "H" = "Index SD = 0.4, Age composition SD = 1.5"
  ))
om_tab <- om_tab %>%
  mutate(Fhist = recode(Fhist,
    "H-MSY" = "$2.5 F_{\\text{MSY}} \\rightarrow F_{\\text{MSY}}$",
    "MSY" = "$F_{\\text{MSY}}$"
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
om_tab$Model <- paste0("$\\text{Sel}_{", 1:16, "}$")
om_tab <- om_tab %>%
  mutate(obs_error = recode(obs_error,
    "L" = "Index SD = 0.1, Age composition SD = 0.3",
    "H" = "Index SD = 0.4, Age composition SD = 1.5"
  ))
om_tab <- om_tab %>%
  mutate(Fhist = recode(Fhist,
    "H-MSY" = "$2.5 F_{\\text{MSY}} \\rightarrow F_{\\text{MSY}}$",
    "MSY" = "$F_{\\text{MSY}}$"
  ))

names(om_tab)[2] <- "$\\sigma_R$"
names(om_tab)[3] <- "$\\sigma_{\\text{Sel}}$"
names(om_tab)[4] <- "$\\rho_{\\text{Sel}}$"
names(om_tab)[5] <- "Fishing History"
names(om_tab)[6] <- "Observation Uncertainty"

x = latex(om_tab, file = here("Project_0","paper","Sel_om_tab.tex"), 
  table.env = FALSE, col.just = rep("r", dim(om_tab)[2]), rowname = NULL)

#q oms: ems = 5-20, 29-32
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.q.oms.RDS"))
om_inputs = readRDS(file.path(here(),"Project_0","inputs", "q_om_inputs.RDS"))

om_tab <- df.oms
om_tab$Model <- paste0("$q_{", 1:16, "}$")
om_tab <- om_tab %>%
  mutate(obs_error = recode(obs_error,
    "L" = "Index SD = 0.1, Age composition SD = 0.3",
    "H" = "Index SD = 0.4, Age composition SD = 1.5"
  ))
om_tab <- om_tab %>%
  mutate(Fhist = recode(Fhist,
    "H-MSY" = "$2.5 F_{\\text{MSY}} \\rightarrow F_{\\text{MSY}}$",
    "MSY" = "$F_{\\text{MSY}}$"
  ))

names(om_tab)[2] <- "$\\sigma_R$"
names(om_tab)[3] <- "$\\sigma_{q}$"
names(om_tab)[4] <- "$\\rho_{q}$"
names(om_tab)[5] <- "Fishing History"
names(om_tab)[6] <- "Observation Uncertainty"

x = latex(om_tab, file = here("Project_0","paper","q_om_tab.tex"), 
  table.env = FALSE, col.just = rep("r", dim(om_tab)[2]), rowname = NULL)


em_tab <- df.ems
em_tab$Model <- paste0("EM$_{", 1:32, "}$")
em_tab$pe <- em_tab$re_config
em_tab$pe[which(em_tab$pe == "M_re")] = paste0(em_tab$pe[which(em_tab$pe == "M_re")], "_", em_tab$M_re_cor[which(em_tab$pe == "M_re")])
em_tab$pe[which(em_tab$pe == "sel_re")] = paste0(em_tab$pe[which(em_tab$pe == "sel_re")], "_", em_tab$sel_re_cor[which(em_tab$pe == "sel_re")])
em_tab$pe[which(em_tab$pe == "q_re")] = paste0(em_tab$pe[which(em_tab$pe == "q_re")], "_", em_tab$q_re_cor[which(em_tab$pe == "q_re")])
em_tab <- em_tab %>%
  mutate(SR_model = recode(SR_model,
    "2" = "Mean recruitment",
    "3" = "Beverton-Holt"
  ))
em_tab[["Mean $M$"]] <- "Estimated"
em_tab[["Mean $M$"]][which(!em_tab$M_est)] <- "0.2"

em_tab <- em_tab %>%
  mutate(pe = recode(pe,
    "rec" = "Recruitment ($\\sigma_R$ estimated)",
    "rec+1" = "Recruitment and survival ($\\sigma_R$, $\\sigma_{2+}$ estimated)",
    "M_re_iid" = "Recruitment and uncorrelated natural mortality ($\\sigma_R$, $\\sigma_{M}$ estimated, $\\rho_{M} = 0$)",
    "M_re_ar1_y" = "Recruitment and AR1 natural mortality ($\\sigma_R$, $\\sigma_{M}$, $\\rho_{M}$ estimated)",
    "sel_re_iid" = "Recruitment and uncorrelated fleet selectivity ($\\sigma_R$, $\\sigma_{\\text{Sel}}$ estimated, $\\rho_{\\text{Sel}} = 0$)",
    "sel_re_ar1_y" = "Recruitment and AR1 selectivity ($\\sigma_R$, $\\sigma_{\\text{Sel}}$, $\\rho_{\\text{Sel}}$ estimated)",
    "q_re_iid" = "Recruitment and uncorrelated catchability (spring index) ($\\sigma_R$, $\\sigma_{q}$ estimated, $\\rho_{q} = 0$)",
    "q_re_ar1" = "Recruitment and AR1 catchability (spring index) ($\\sigma_R$, $\\sigma_{q}$, $\\rho_{q}$ estimated)",
  ))
em_tab <- em_tab[c("Model","SR_model","Mean $M$", "pe")]
names(em_tab)[2] <- "Recruitment model"
names(em_tab)[4] <- "Process error assumption"

x = latex(em_tab, file = here("Project_0","paper","em_tab.tex"), 
  table.env = FALSE, col.just = rep("r", dim(em_tab)[2]), rowname = NULL)

