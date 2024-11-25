library(here)
library(tidyverse)
library(wham)
library(rpart)
library(rpart.plot)
source(file.path(here::here(), "Ecov_study","recruitment_functions", "code", "analysis", "functions_for_analysis.R" ) )

#dir <- file.path(here::here(), "Ecov_study","recruitment_functions", "results_pile2")
res.dir     <- 'results'   # 'results'     'results_beta_fix'   # results folder where AIC dataframes are
plot.dir    <- 'plots'    # 'plots_lizruns'  'plots_beta_fix'  
table.dir   <- 'tables'   # 'results'     'results_beta_fix'
plot.suffix <- ''      # '_beta_fix'   '' 

df.oms    <- readRDS(file.path(here::here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS"))

rho.df <- readRDS(mrho.df,  file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0("mrho.df.RDS") ) )
