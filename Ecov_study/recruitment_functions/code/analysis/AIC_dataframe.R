library(here)
library(wham)
#folder link to id
# dir <- '~/dropbox/working/state_space_assessments/cluster_download/results_noSR//'
#df.oms          <- readRDS(file.path(here(),"Ecov_study","recruitment_functions", "inputs", "df.oms_noSR.RDS"))
source(file.path(here::here(), "Ecov_study","recruitment_functions", "code", "analysis", "functions_for_analysis.R" ) )

#res.path <- 'E:/results_beta_fix'  # directory where simulation runs are (beta unstandardized)
res.path <- file.path(here::here(),"Ecov_study","recruitment_functions","results")  # directory where simulation runs are (beta unstandardized)

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
t4 <- Sys.time()
get_aic_convergence_info(df.oms, df.ems, nsims=n_sims, res.path=res.path, save.rds=TRUE, save.suffix='_beta_fix', save.path=here::here("Ecov_study", 'recruitment_functions', 'results_beta_fix')  ) 
# "This took 20.7794454693794" minutes (liz's runs)!  
cd


t5 <- Sys.time()
t5-t4   # Time difference of 4.075107 hours