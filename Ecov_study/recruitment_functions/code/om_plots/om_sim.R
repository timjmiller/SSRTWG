
library(here)

om_sim <- function(om_input, nsim){
  om  <- fit_wham(om_input, do.fit = FALSE, MakeADFun.silent = TRUE)
  lapply(1:nsim, function(x) return(om$simulate(complete=TRUE)))
}

n_sim <- 100
n_om  <- length(om_inputs)

lapply(1:n_om, function(x){
  sim <- om_sim(om_inputs[[x]],nsim=n_sim)
  saveRDS(file=paste0(here(),'/Ecov_study/recruitment_functions/results/om_',x,'_','n_sims_',n_sim),sim)
}) 

