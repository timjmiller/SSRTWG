library(wham)
library(here)

df.ems    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))
df.oms    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.oms.RDS"))

n.em <- nrow(df.ems)
n.om <- nrow(df.oms)
om.names <- paste0("om", seq(1,n.om))
#iters <- seq(1,50) #liz ran 50
iters <- seq(1,100) #greg ran 100
n.iter <- length(iters)
n.iter.em <- n.iter*n.em

sim.names <- paste0(rep("sim",n.iter.em), rep(iters, each=n.em), "_em", rep(seq(1,n.em), n.iter), ".RDS")

om.fails <- c()
em.fails <- c()
iter.fails <- c()

fail.list <- list()
t1 <- Sys.time()

for (i in 1:n.om) {
  om=om.names[i]
  
  #sims.run <- list.files(file.path(here::here(),"Ecov_study", "recruitment_functions", "results", om), pattern = ".RDS")
  sims.run <- list.files(file.path('E:/results_beta_fix', om), pattern = ".RDS")
  fails <- which((sim.names %in% sims.run)==FALSE)
  
  if(length(fails)==0) fail.list[[i]] <- NA
  
  if(length(fails>0)) {
    om.fails <- as.numeric(substr(om, 3, nchar(om) ) )
    tmp <- strsplit(sim.names[fails], split="_")
    iter.fails <- as.numeric(c(unlist(lapply(tmp, function(x) {substr(x[1], 4, nchar(x[1]) ) }))) )
    em.fails <- as.numeric(c(unlist( lapply(tmp, function(x) {substr(x[2], 3, 3 ) }))  ) )
    
    fail.list[[i]] <- list(om=om.fails, nfails = length(tmp) , iter_em=cbind(iter.fails, em.fails))
    
  }
    om.fails <- c()
    em.fails <- c()
    iter.fails <- c()

}
t2 <- Sys.time()
t2-t1  # 22.01293 secs

om.reruns <- which(lapply(fail.list, function(x) is.na(x[1]) )==FALSE)



#####################################################################################
# make df of failed runs ====

om.fails.df <- c()  #

for(i in 1:nrow(df.oms)) {
  
  if(is.na(fail.list[[i]][1])==FALSE) om.fails.df <- rbind(om.fails.df,
                                                           cbind(OM=rep(i, fail.list[[i]]$nfails), 
                                                                 N.fail=rep(fail.list[[i]]$nfails, fail.list[[i]]$nfails), 
                                                                 fail.list[[i]]$iter_em) )
}

saveRDS(om.fails.df, file.path(here::here("Ecov_study", 'recruitment_functions', 'results_beta_fix'),'om.fails.df_beta_unstandardized.RDS'))


#####################################################################################

# saveRDS(fail.list, file=file.path(here::here(),"Ecov_study", "recruitment_functions", 
#                                   "results", "failed.list_1.RDS") ) 
# saveRDS(fail.list, file=file.path(here::here(),"Ecov_study", "recruitment_functions", 
#                                   "results", "failed.list_2.RDS") ) 
# saveRDS(fail.list, file=file.path(here::here(),"Ecov_study", "recruitment_functions", 
#                                   "results", "failed.list_3.RDS") ) 
# saveRDS(fail.list, file=file.path(here::here(),"Ecov_study", "recruitment_functions", 
#                                   "results", "failed.list_4.RDS") ) 
# saveRDS(fail.list, file=file.path(here::here(),"Ecov_study", "recruitment_functions", 
#                                   "results", "failed.list_5.RDS") ) 
# saveRDS(fail.list, file=file.path(here::here(),"Ecov_study", "recruitment_functions", 
#                                   "results", "failed.list_6.RDS") ) 
# saveRDS(fail.list, file=file.path(here::here(),"Ecov_study", "recruitment_functions", 
#                                   "results", "failed.list_7.RDS") ) 

saveRDS(fail.list, file=file.path(here::here(),"Ecov_study", "recruitment_functions", 
                                  "results", "failed.list_8.RDS") ) 



################################################################
# re-run the failed sims ==============================
################################################################
#om.folder <- om.reruns


library(wham)
library(here)
source(file.path(here::here(), "Ecov_study", "recruitment_functions", "code", "sim_management.R"))
library(doParallel)
#verify_version()
# The right commit: 77bbd94 of wham is loaded! 
shortRversion <- function() {
  rvs <- R.version.string
  if(grepl("devel", (st <- R.version$status)))
    rvs <- sub(paste0(" ",st," "), "-devel_", rvs, fixed=TRUE)
  gsub("[()]", "", gsub(" ", "_", sub(" version ", "-", rvs)))
}

RVersion.short=shortRversion()
TMBversion=packageVersion('TMB')
MatrixVersion=packageVersion('Matrix')

# NOTES
# 1. om_inputs and df.oms and seeds were created in om_setup.R
# 2. em_inputs and df.ems were created in em_setup.R

om_inputs <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "om_inputs.RDS"))
em_inputs <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "em_inputs.RDS"))
df.ems    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))
df.oms    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.oms.RDS"))


#######################################################
seeds <- readRDS(file.path(here::here(), "Ecov_study", "recruitment_functions", "inputs","seeds.RDS"))
#######################################################

numCore <- detectCores()

myCluster <- makeCluster(round(0.6*numCore,0), # number of cores to use
                         type = "PSOCK") # type of cluster
# print(paste0("Number of Cores: ",detectCores()))
# print(paste0("Cores Registered: ",myCluster))
registerDoParallel(myCluster)


# sims <- seq(1,50)  #create a sequence, so i can respecify and pick up later
# nsims = length(sims)  #50 #100
t1 <- Sys.time()
print(t1)
for(om in 1:length(om.reruns))   {
  
  omi = om.reruns[om]

  write.dir <- file.path(here::here(),"Ecov_study", "recruitment_functions", "results", paste0("om", omi))
  dir.create(write.dir, recursive = T, showWarnings = FALSE)
  
  gradient.txt <- c()
  write.table(gradient.txt, file=file.path(here::here(), "Ecov_study", "recruitment_functions",
                                           'results',  paste0("om", omi),
                                           "gradient.test.txt"), append=FALSE)
  
  warns <- c()
  errs <- c()
  
  n.iter.fails <- fail.list[[om.reruns[om]]]$nfails

  foreach(isim = 1:n.iter.fails )  %dopar% {
    library(wham)
    
    iterj = fail.list[[om.reruns[om]]]$iter_em[isim,1]
    emk = fail.list[[om.reruns[om]]]$iter_em[isim,2]
    #for(emk in 1:nrow(df.ems))  {
      # clear this.err and this.warn for tryCatch
      this.err <- NULL
      this.warn <- NULL
      #extract info
      x        <- data.frame(df.ems[emk,])
      names(x) <- paste0('em_',names(x))
      y        <- data.frame(df.oms[omi,])
      names(y) <- paste0('om_',names(y))
      model    <- cbind(im=iterj, om=omi, em=emk, optimized=FALSE, sdreport=FALSE, y,x)
      
      obs_names <- c("agg_catch","agg_catch_sigma", "agg_indices", "agg_index_sigma", "catch_paa", "index_paa", 
                     "Ecov_obs", "obs", "obsvec")
      
      

      
      om                 <- fit_wham(om_inputs[[omi]], do.fit = FALSE, MakeADFun.silent = TRUE)
      set.seed(seeds[[omi]][iterj])
      sim_data           <- om$simulate(complete=TRUE)
      truth              <- sim_data
      truth$wham_version <- om$wham_version #save the version for reproducibility
      EM_input           <- em_inputs[[emk]] # Read in the EM
      
      #put simulated data into the em input
      EM_input$data[obs_names] = sim_data[obs_names]
      
      #not estimating observation error in Ecov
      EM_input$par$Ecov_obs_logsigma[] <- om_inputs[[omi]]$par$Ecov_obs_logsigma
      
      ## the fixed effects used to generate truth
      ompars <- data.frame(par=names(om$par), value=om$par) |> dplyr::filter(par!='F_devs')
      ompars$par2 <- sapply(unique(ompars$par), function(x) {
        y <- which(ompars$par==x)
        if(length(y)==1) return(x)
        x <- paste(x, 1:length(y), sep='_')
        return(x)
      }) %>% unlist
      res <- list(truth = truth, model=model, ompars=ompars, seed=seeds[[omi]][iterj], RVersion=RVersion.short,
                  TMBversion=TMBversion, MatrixVersion=MatrixVersion)
      res$fit <- list()
      
      ## Build test object to test for initial 0 gradients?
      test <- fit_wham(EM_input, do.fit=FALSE, do.sdrep=F, do.osa=F, do.retro=F, do.proj=F, MakeADFun.silent=TRUE)
      if(any(test$gr()==0)) {
        gradient.txt <- ((paste0("Initial gradients 0 in OM: ", omi, " Sim: ", iterj, 
                                 " EM: ", emk)))
        write.table(gradient.txt, file=file.path(here::here(), "Ecov_study", "recruitment_functions",
                                                 'results',  paste0("om", omi),
                                                 "gradient.test.txt"), append=TRUE)
        
      }
      
      
      #do fit without sdreport first
      fit <- tryCatch(fit_wham(EM_input, do.check=TRUE, do.sdrep=F, do.osa=F, do.retro=F, do.proj=F, 
                               MakeADFun.silent=TRUE), error = function(e) conditionMessage(e) )

      
      # Deal with issues fitting EM to non-matching OM data
      # empty elements below can be used to summarize convergence information
      if(!'err' %in% names(fit) & class(fit) != "character"){
        res$model$optimized <- TRUE
        res$fit <- fit[c("wham_version", "TMB_version", "opt", "final_gradient", "rep")]
        ## res$fit$mohns_rho <- tryCatch(mohns_rho(fit),
        ##                               error = function(e) conditionMessage(e))
        fit$sdrep <- tryCatch(TMB::sdreport(fit), # no bc
                              error = function(e) conditionMessage(e))
        if(class(fit$sdrep) == "sdreport"){
          res$model$sdreport <- TRUE
          res$fit$sdrep <- list(
            "Estimate_par" = as.list(fit$sdrep, what = "Est"),
            "SE_par" = as.list(fit$sdrep, what = "Std"),
            "Estimate_rep" = as.list(fit$sdrep, what = "Est", report = TRUE),
            "SE_rep" = as.list(fit$sdrep, what = "Std", report = TRUE))
        }
        empars <- data.frame(par=names(res$fit$opt$par), value=res$fit$opt$par)%>%
          dplyr::filter(!grepl(x=par,'F_devs|log_NAA'))
        empars$par2 <- sapply(unique(empars$par), function(x) {
          y <- which(empars$par==x)
          if(length(y)==1) return(x)
          x <- paste(x, 1:length(y), sep='_')
          return(x)
        }) %>% unlist
        res$empars <- empars
      }
      
      rds.fn = file.path(here::here(), "Ecov_study", "recruitment_functions", "results", paste0("om", omi), paste0("sim", iterj, "_em", emk, ".RDS"))
      saveRDS(res, file = rds.fn)

  #  } # end emk loop over df.ems
    
    
  } # end iterj loop across iterations (nsims) 
  
  
  
} # end omi loop over df.oms

t2 <- Sys.time()
print(t2)
print(t2 - t1)


stopCluster(cl = myCluster)  # this seemed to fix the error below:
# Error in file(con, "w") : all connections are in use
# In addition: There were 23 warnings (use warnings() to see them)
# warnings()
# 23: In .Internal(gzfile(description, open, encoding, compression)) :
#   closing unused connection 105 (<-lizRmachine:11908)
# > 


