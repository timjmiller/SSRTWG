library(here)
library(dplyr)
library(Hmisc)
df.ems <- readRDS(file.path(here(),"Ecov_study","mortality", "inputs", "df.ems.RDS"))
df.oms <- readRDS(file.path(here(),"Ecov_study","mortality", "inputs", "df.oms.RDS"))
om_inputs <- readRDS(file.path(here(),"Ecov_study","mortality", "inputs", "om_inputs.RDS"))

factor_tab <- cbind.data.frame(
  Factor = c(
    "OM Process error", 
    "OM Fishing History", 
    "OM Covariate effect size ($\\beta_E$)", 
    "OM Population Observation Error", 
    "OM Covariate Observation Error", 
    "OM Covariate Process Error SD ($\\sigma_E$)", 
    "OM Covariate Process Error Correlation ($\\rho_E$)",
    "EM Process Error",
    "EM Covariate Effect",
    "EM Median Natural Mortality Rate Parameter"),
  Levels = 
    c("R, R+S, R+M", 
    "$F_{\\text{MSY}}$, $2.5F_{\\text{MSY}} \\rightarrow F_{\\text{MSY}}$", 
    "0, 0.25, 0.5", 
    "Low, High", 
    "Low ($\\sigma_e = 0.1$), High ($\\sigma_e = 0.5$)",
    "0.1, 0.5",
    "0, 0.5",
    "R, R+S, R+M",
    "None ($\\beta_E = 0$), $\\beta_E$ estimated",
    "Known, Estimated"  ))


x = latex(factor_tab, file = here("Ecov_study","mortality","manuscript","factor_table.tex"), 
  table.env = FALSE, col.just = c("l","c"), rowname = NULL)

