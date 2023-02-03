args <- as.numeric(commandArgs(trailingOnly=TRUE))

library(wham)
library(here)
library(doParallel)
ncores <- detectCores()      
registerDoParallel(ncores-1) #leave one core for other tasks

write.dir <- file.path(here(),"Ecov_study", "recruitment", "results") # create directory for analysis
nsim = 100 #number of simulations for each scenario

em_input <- readRDS(file.path(write.dir, "em_input_GLB_recruitment_doparallel.RDS"))

########################################################
##--FIT MODELS--########################################
########################################################
em_fits <- foreach(x = 1:nsim) %dopar% {
    cat(paste("model:",args[1], "fit:", x, "start \n"))
    out = fit_wham(em_input[[args[1]]][[x]], do.osa = FALSE, MakeADFun.silent = TRUE, retro.silent = TRUE, save.sdrep = FALSE)
    out$rep <- out$report()
    cat(paste("model:",args[1],"fit:", x, "done \n"))
    return(out)
}
saveRDS(em_fits, file.path(write.dir, paste0("em_fits_m_",args[1],".RDS")))
