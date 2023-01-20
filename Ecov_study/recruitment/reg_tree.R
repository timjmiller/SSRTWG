library(rpart)
library(rpart.plot)
library(wham)

##--Read data--###############

dir <- "~/dropbox/working/state_space_assessments/cluster/results/"
dir <- "~/dropbox/working/state_space_assessments/SSRTWG/Ecov_study/recruitment/results/"
em_fits <- readRDS(file.path(dir,'em_fits_GLB_recruitment_doparallel.RDS'))
df.mods <- readRDS(file.path(dir,'om_sim_inputs_GLB_recruitment_doparallel.RDS'))
n.mods <- nrow(df.mods)
nsim <- length(em_fits[[1]])


temp <- list()
for(m in 1:n.mods){
    temp[[m]] <- lapply(1:nsim, function(x){
    ##SSB
    ssb_sim <- em_fits[[m]][[x]]$input$data$SSB
    ssb_est <- em_fits[[m]][[x]]$report()$SSB
    
    r_sim <- em_fits[[m]][[x]]$input$data$NAA[,1]
    #r_est <- em_fits[[m]][[x]]$report()$NAA[,1]
    r_est <- exp(em_fits[[m]][[x]]$parList$log_NAA[,1])
    
    F_sim <- em_fits[[m]][[x]]$input$data$F
    F_est <- em_fits[[m]][[x]]$report()$F
    
    return(list(dssb =mean(ssb_est - ssb_sim),
                d2ssb=sqrt(mean((ssb_est - ssb_sim)^2)),
                dr   =mean(r_est - r_sim),
                d2r  =sqrt(mean((r_est - r_sim)^2)),
                dF   =mean(F_est - F_sim),
                d2F  =sqrt(mean((F_est - F_sim)^2)),
                ))
  })
}
results <- unlist(temp)

################################################
## MODEL BUILDING ##############################
################################################
df <- df.mods[rep(1:nrow(df.mods),each=2),] %>%
  mutate(dssb =results[names(results)=='dssb'],
              d2ssb=results[names(results)=='d2ssb'],
              dr=results[names(results)=='dr'],
              d2r=results[names(results)=='d2r'],
              dF=results[names(results)=='dF'],
              d2F=results[names(results)=='d2F'])

fit <- rpart(dr ~ Ecov_sig + Ecov_phi + beta + obs_sig + NAA_sig, data=df)
rpart.plot(fit)


plot(r_sim)
points(r_est,col='red')
