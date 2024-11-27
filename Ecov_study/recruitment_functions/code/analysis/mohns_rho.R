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

## specify bad.grad.label and bad.se.value (these are the thresholds set in convergence_summaries.R to determine convergence) 
bad.grad.value <- 1E-6 #tim used 1E-6 (the abs of the exponent will be used for output suffix; ex: filename_grad_6.png)
bad.grad.label <- as.numeric(strsplit(as.character(bad.grad.value), split="-")[[1]][2])
bad.se.value   <- 100 #tim used 100 (this value will be used for output suffix; ex: filename_se_100.png)

df.oms    <- readRDS(file.path(here::here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS"))
df.ems    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))
conv.runs <- readRDS(file.path(here::here(),'Ecov_study','recruitment_functions',res.dir, paste0("conv.runs_grad_", bad.grad.label, "_SE_", bad.se.value, plot.suffix, ".RDS")  ) )

conv.runs <- conv.runs %>%
  rename(Sim=sim) %>%
  mutate(ok.run =(bad.opt+ bad.conv+ bad.sdrep+ bad.grad+ bad.se.big) ) %>%
  replace_na(list(ok.run=1)) %>%
  relocate(OM, EM, Sim, ok.run, bad.opt, bad.conv, bad.sdrep, bad.grad, bad.se.big) %>%
  filter(ok.run==0)  # drop unconverged runs (there shouldn't be any at this point)

n_oms <- nrow(df.oms)
n_ems <- nrow(df.ems)
n_sims <- 100
nyears <- 40

n.conv.runs <- nrow(conv.runs)
ten.pct <- round(n.conv.runs/10, 0)


mrho.df <- matrix(NA,ncol=16,nrow=nrow(conv.runs))
colnames(mrho.df) <- c('OM','EM','Sim','ssb','fbar','n1','n2','n3','n4','n5','n6','n7','n8','n9','n10','rdev')

for(irun in 1:n.conv.runs) {    #   1:n.conv.runs
print(irun/n.conv.runs)
  dat <- try(readRDS(file.path(here::here(),'Ecov_study','recruitment_functions','results', paste0("om", conv.runs$OM[irun], '/','sim',conv.runs$Sim[irun],'_','em', conv.runs$EM[irun],'.RDS') ) ) )
  if(class(dat)!='try-error'){
    mrho        <- mohns_rho_set_peel(model=dat, npeels=7, ny=40, na=10)
    mrho_raneff <- mohns_rho_randeff_peel(model=dat, npeels=7, ny=40)

    mrho.df[irun,] <- c(OM=conv.runs$OM[irun],
                          EM=conv.runs$EM[irun],
                          Sim=conv.runs$Sim[irun],
                          mrho_ssb=mrho[1],
                          mrho_fbar=mrho[2],
                          mrho_n1=mrho[3],
                          mrho_n2=mrho[4],
                          mrho_n3=mrho[5],
                          mrho_n4=mrho[6],
                          mrho_n5=mrho[7],
                          mrho_n6=mrho[8],
                          mrho_n7=mrho[9],
                          mrho_n8=mrho[10],
                          mrho_n9=mrho[11],
                          mrho_n10=mrho[12],
                          mrho_rdev=mrho_raneff)
  }
}

print('saving mrho...')         
saveRDS(mrho.df,  file.path(here::here(),'Ecov_study','recruitment_functions',res.dir , paste0("mrho.df.RDS") ) )
