library(here)
ems <- c('em1','em2','em3')

#folder link to id
dir <- '~/dropbox/working/state_space_assessments/cluster_download/results_google/results/'

folders <- list.files(dir)

ff <- lapply(1:10,function(om){
  lapply(1:3, function(em){
    sapply(1:100, function(sim){
      print(paste0("om ",om," em ",em," sim ",sim))
      
      dat <- try(readRDS(paste0(dir, "om", om, '/','sim',sim,'_','em',em,'.RDS')),silent=TRUE)
      if(class(dat)!='try-error'){
        #dat$truth$SSB - dat$fit$rep$SSB
        #dat$truth$SSB - dat$fit$rep$SSB
        dat$truth$NAA[,1] - dat$fit$rep$NAA
      }
    })
  })
})


om <- 1; em <- 1; sim <- 1


om <- 4; em <- 3; sim <- 26
