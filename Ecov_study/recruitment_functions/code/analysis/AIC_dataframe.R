library(here)
library(wham)
#folder link to id
# dir <- '~/dropbox/working/state_space_assessments/cluster_download/results_noSR//'
#df.oms          <- readRDS(file.path(here(),"Ecov_study","recruitment_functions", "inputs", "df.oms_noSR.RDS"))
source(file.path(here::here(), "Ecov_study","recruitment_functions", "code", "analysis", "functions_for_analysis.R" ) )

#res.path <- 'E:/results_beta_fix'  # directory where simulation runs are (beta unstandardized)
res.path    <- file.path(here::here(),"Ecov_study","recruitment_functions","results")  # directory where simulation runs are (beta unstandardized)
save.path   <- file.path(here::here(),"Ecov_study",'recruitment_functions', 'results')  
save.suffix <- ''

df.oms    <- readRDS(file.path(here::here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS"))
df.ems    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))

n_oms <- nrow(df.oms)
n_ems <- nrow(df.ems)
n_sims <- 100

################################################################################
# make 2 dataframes (AIC_all, AIC_weight)  =====================================
## AIC_all has info on the EM with lowest AIC per OM-Sim
## AIC_weight has info on all EMs per OM-Sim
## both dfs contain convergence info
################################################################################
get_aic_convergence_info(df.oms, df.ems, 
                         nsims=n_sims, 
                         res.path=res.path, 
                         save.rds=TRUE, 
                         save.suffix=save.suffix, 
                         save.path=save.path) 

