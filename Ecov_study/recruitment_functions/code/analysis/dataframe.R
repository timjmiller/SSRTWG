library(here)
ems <- c('em1','em2','em3')
n_ems <- length(ems)
n_oms <- 96
#folder link to id
dir <- '~/dropbox/working/state_space_assessments/cluster_download/results_noSR//'

#df.oms          <- readRDS(file.path(here(),"Ecov_study","recruitment_functions", "inputs", "df.oms.RDS"))
df.oms          <- readRDS(file.path(here(),"Ecov_study","recruitment_functions", "inputs", "df.oms_noSR.RDS"))
n_oms <- nrow(df.oms)

folders <- list.files(dir)

ff <- lapply(1:10,function(om){
  lapply(1:3, function(em){
    sapply(1:100, function(sim){
      print(paste0("om ",om," em ",em," sim ",sim))
      
      dat <- try(readRDS(paste0(dir, "om", om, '/','sim',sim,'_','em',em,'.RDS')),silent=TRUE)
      if(class(dat)!='try-error'){
        #dat$truth$SSB - dat$fit$rep$SSB
        dat$truth$NAA[,1] - dat$fit$rep$NAA[,1]
      }
    })
  })
})

###################################
## AIC ############################
###################################
nsim <- n_oms*100

AIC           <- data.frame(matrix(nrow=nsim, ncol=ncol(df.oms)+4))
colnames(AIC) <- c('sim',colnames(df.oms),'aic_pick','correct_form','correct_SR')

k <- 1
for(om in 1:n_oms){
  print(paste0("OM = ",om))
  for(sim in 1:100){
    DAT <- sapply(1:4, function(em){
      dat <- tryCatch(readRDS(paste0(dir, "om", om, '/','sim',sim,'_','em',em,'.RDS')),
                      error = function(e) conditionMessage(e))
      aic <- NA
      if(class(dat)!='try-error' & class(dat)!='character'){
        if(length(dat$fit)>0) aic = as.numeric(dat$fit$opt$objective) + as.numeric(length(dat$fit$opt$par))
      }
      aic
    })
    
    em_match <- which(df.ems$ecov_how==df.oms$Ecov_how[om] &
          df.ems$r_mod==df.oms$recruit_mod[om])
    
    aic_pick <- which(DAT==min(DAT,na.rm=TRUE))
    AIC[k,]  <- data.frame(sim=sim,df.oms[om,],aic_pick=aic_pick,
                           correct_form=ifelse(aic_pick==em_match,1,0))#,
                           #correct_SR=ifelse(aic_pick%in%c(1,2,3),1,0))
    k <- k + 1
  }  
}

saveRDS(AIC,file.path(here(), 'Ecov_study','recruitment_functions','results','AIC.rds'))
