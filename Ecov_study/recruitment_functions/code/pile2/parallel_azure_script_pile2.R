rm(list=ls(all.names = FALSE)) # clear session memory
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

#om_inputs <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "om_inputs.RDS"))
#em_inputs <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "em_inputs.RDS"))
#df.ems    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))
#df.oms    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.oms.RDS"))
om_inputs <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "code", "pile2", "om_inputs.RDS"))
em_inputs <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "code", "pile2", "em_inputs.RDS"))
df.ems    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "code", "pile2", "df.ems.RDS"))
df.oms    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "code", "pile2", "df.oms.RDS"))


#######################################################
#seeds <- readRDS(file.path(here::here(), "Ecov_study", "recruitment_functions", "inputs","seeds.RDS"))
seeds <- readRDS(file.path(here::here(), "Ecov_study", "recruitment_functions", "code", "pile2","seeds.RDS"))
#######################################################

numCore <- detectCores()

myCluster <- makeCluster(round(0.8*numCore,0), # number of cores to use
                         type = "PSOCK") # type of cluster
# print(paste0("Number of Cores: ",detectCores()))
# print(paste0("Cores Registered: ",myCluster))
registerDoParallel(myCluster)


sims <- seq(1,10)  #create a sequence, so i can respecify and pick up later
nsims = length(sims)  #50 #100
om_set <- seq(71,72)    #seq(1,nrow(df.oms))
n_om_set <- length(om_set)
em_set <- seq(1,6)    #seq(1,nrow(df.ems))
n_em_set <- length(em_set)

# specs for the retro peels & projections ===
n.yrs.peel <- 10
n.yrs.proj <- 10
yrs.avg <- 5  # will be used to average ecov and bio pars

t1 <- Sys.time()
print(t1)

for(om_index in 1:n_om_set )   {
  omi <- om_set[om_index]
# for(omi in 7:7)   {    #test with OM 67 (Ecov_how=1, r_mod=3, low corr, low obs err, high Rsig)
  
  write.dir <- file.path(here::here(),"Ecov_study", "recruitment_functions", "results_pile2", paste0("om", omi))
  dir.create(write.dir, recursive = T, showWarnings = FALSE)
  
    gradient.txt <- c()
  write.table(gradient.txt, file=file.path(here::here(), "Ecov_study", "recruitment_functions",
                                          'results_pile2',  paste0("om", omi),
                                           "gradient.test.txt"), append=FALSE)
  
  warns <- c()
  errs <- c()
  
  foreach(iterj = sims[1]:sims[nsims] )  %dopar% {
    library(wham)
    
    # for(emk in 1:1)  {
      for(em_index  in 1:n_em_set  )  {
        emk <- em_set[em_index]
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
      
      
      # txt.start <- (paste0("START OM: ", omi, " Sim: ", iterj, " EM: ", emk, "\n"))
      
      
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
      test <- fit_wham(EM_input, do.fit=FALSE, do.sdrep=F, do.osa=F, do.retro=FALSE, do.proj=F, MakeADFun.silent=TRUE)
      if(any(test$gr()==0)) {
        gradient.txt <- ((paste0("Initial gradients 0 in OM: ", omi, " Sim: ", iterj, 
                                                     " EM: ", emk)))
        write.table(gradient.txt, file=file.path(here::here(), "Ecov_study", "recruitment_functions",
                                                 'results_pile2',  paste0("om", omi),
                                                 "gradient.test.txt"), append=TRUE)
        
      }
      
      
      #do fit without sdreport first
      fit <- tryCatch(fit_wham(EM_input, do.check=TRUE, do.sdrep=F, do.osa=F,  do.retro=FALSE, do.proj=F, 
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
        
           # do peels ====
           
          t1 <- Sys.time()
            fit$peels <- retro(fit, ran = unique(names(fit$env$par[fit$env$random])), do.sdrep = TRUE,
                               n.peels=n.yrs.peel, save.input = TRUE, MakeADFun.silent = TRUE,
                               retro.silent = TRUE, n.newton=1) #n.newton = 0 (2.25 min); n.newton=1 (3.1 min; 4.531353 mins with sdrep); n.newton=2 (3.9 min)
            t2 <- Sys.time()
            t2-t1    # 7.9 secs
            
          
        # project from last peel ====
          #t1 <- Sys.time()
          tmp.proj.cont.ecov <- project_wham(
            model=fit$peels[[n.yrs.peel]],
            proj.opts = list(n.yrs = n.yrs.proj,  
                             use.last.F = FALSE, use.avg.F = FALSE, use.FXSPR = FALSE, use.FMSY = FALSE, 
                             proj.F = tail(truth$F, n.yrs.peel),# use Truth so can compare projections with Trut  
                             #proj.catch = NULL, 
                             avg.yrs = tail(fit$peels[[n.yrs.proj]]$input$years, yrs.avg), # this will cause selectivity to be averaged over last 5 years
                             cont.ecov = TRUE, use.last.ecov = FALSE, avg.ecov.yrs = NULL, proj.ecov = NULL, 
                             cont.Mre = NULL, 
                             #avg.rec.yrs = NULL,  #(only relevant if recr are fixed effects)
                             percentFXSPR = 100, percentFMSY = 100 #,
                             #proj_F_opt = rep(3,n.yrs.proj)
            ),
            n.newton = 3,
            do.sdrep = TRUE,
            MakeADFun.silent = TRUE,
            save.sdrep = TRUE
          )  
          
          # t2 <- Sys.time()
          # t2-t1    # 7.9 secs
          
          #t1 <- Sys.time()
          tmp.proj.avg.ecov <- project_wham(
            model=fit$peels[[n.yrs.peel]],
            proj.opts = list(n.yrs = n.yrs.proj,  
                             use.last.F = FALSE, use.avg.F = FALSE, use.FXSPR = FALSE, use.FMSY = FALSE, 
                             proj.F = tail(truth$F, n.yrs.peel),  # use Truth so can compare projections with Truth
                             #proj.catch = NULL, 
                             avg.yrs = tail(fit$peels[[n.yrs.proj]]$input$years, yrs.avg), # this will cause selectivity to be averaged over last 5 years
                             cont.ecov = FALSE, use.last.ecov = FALSE, 
                             avg.ecov.yrs = tail(fit$peels[[n.yrs.proj]]$input$years, yrs.avg), 
                             proj.ecov = NULL, 
                             cont.Mre = NULL, 
                             #avg.rec.yrs = NULL,  #(only relevant if recr are fixed effects)
                             percentFXSPR = 100, percentFMSY = 100 #,
                             #proj_F_opt = rep(3,n.yrs.proj)
            ),
            n.newton = 3,
            do.sdrep = TRUE,
            MakeADFun.silent = TRUE,
            save.sdrep = TRUE
          )  
          
          # t2 <- Sys.time()
          # t2-t1      # 7.9 secs
          
          
          # t1 <- Sys.time()
          tmp.proj.use.ecov <- project_wham(
            model=fit$peels[[n.yrs.peel]],
            proj.opts = list(n.yrs = n.yrs.proj,  
                             use.last.F = FALSE, use.avg.F = FALSE, use.FXSPR = FALSE, use.FMSY = FALSE, 
                             proj.F = tail(truth$F, n.yrs.peel),  # use Truth so can compare projections with Truth
                             #proj.catch = NULL, 
                             avg.yrs = tail(fit$peels[[n.yrs.proj]]$input$years, yrs.avg), # this will cause selectivity to be averaged over same years as ecov
                             cont.ecov = FALSE, use.last.ecov = FALSE, 
                             avg.ecov.yrs = NULL, 
                             proj.ecov = tail(fit$input$data$Ecov_obs, n.yrs.proj), 
                             cont.Mre = NULL, 
                             #avg.rec.yrs = NULL,  #(only relevant if recr are fixed effects)
                             percentFXSPR = 100, percentFMSY = 100 #,
                             #proj_F_opt = rep(3,n.yrs.proj)
            ),
            n.newton = 3,
            do.sdrep = TRUE,
            MakeADFun.silent = TRUE,
            save.sdrep = TRUE
          )  
          
          # t2 <- Sys.time()
          # t2-t1      # 7.9 secs
          
          
          

          
          #res$peels <- fit$peels
          res$peels <- lapply(fit$peels, function(x) list(rep=x$rep,
                                                        sdrep= list(
                                                      "Estimate_par" = as.list(x$sdrep, what = "Est"),
                                                      "SE_par" = as.list(x$sdrep, what = "Std"),
                                                      "Estimate_rep" = as.list(x$sdrep, what = "Est", report = TRUE),
                                                            "SE_rep" = as.list(x$sdrep, what = "Std", report = TRUE) ) )
          )
          
          
          
          res$proj$cont.ecov$rep <- tmp.proj.cont.ecov$rep
          res$proj$cont.ecov$sdrep <- list(
            "Estimate_par" = as.list(tmp.proj.cont.ecov$sdrep, what = "Est"),
            "SE_par" = as.list(tmp.proj.cont.ecov$sdrep, what = "Std"),
            "Estimate_rep" = as.list(tmp.proj.cont.ecov$sdrep, what = "Est", report = TRUE),
            "SE_rep" = as.list(tmp.proj.cont.ecov$sdrep, what = "Std", report = TRUE)
            )
          
          res$proj$avg.ecov$rep <- tmp.proj.avg.ecov$rep
          res$proj$avg.ecov$sdrep <- list(
            "Estimate_par" = as.list(tmp.proj.avg.ecov$sdrep, what = "Est"),
            "SE_par" = as.list(tmp.proj.avg.ecov$sdrep, what = "Std"),
            "Estimate_rep" = as.list(tmp.proj.avg.ecov$sdrep, what = "Est", report = TRUE),
            "SE_rep" = as.list(tmp.proj.avg.ecov$sdrep, what = "Std", report = TRUE)
            )
          
          res$proj$use.ecov$rep <- tmp.proj.use.ecov$rep
          res$proj$use.ecov$sdrep <- list(
            "Estimate_par" = as.list(tmp.proj.use.ecov$sdrep, what = "Est"),
            "SE_par" = as.list(tmp.proj.use.ecov$sdrep, what = "Std"),
            "Estimate_rep" = as.list(tmp.proj.use.ecov$sdrep, what = "Est", report = TRUE),
            "SE_rep" = as.list(tmp.proj.use.ecov$sdrep, what = "Std", report = TRUE)
            )
          
          
          
  
          
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
      
      rds.fn = file.path(here::here(), "Ecov_study", "recruitment_functions", "results_pile2", paste0("om", omi), paste0("sim", iterj, "_em", emk, ".RDS"))
      saveRDS(res, file = rds.fn)
      
      #  txt.end <- (paste0("END OM: ", omi, " Sim: ", iterj, " EM: ", emk, "\n"))
      #  txt.run <- rbind(txt.start, txt.end)
      # write.table(txt.run, file.path(here::here(), "Ecov_study", "recruitment_functions", "results", paste0("om", omi), 'runinfo.txt'))
      
    } # end emk loop over df.ems
    
    
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



