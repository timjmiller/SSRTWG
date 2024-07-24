# functions used for analysis SSRTWG - Project 3 (Ecov & Recruitment)

#=======================================================================================
#   Grab Convergence Info ====
#=======================================================================================

#tim's convergence check code, additions from liz
# convergence IF out[1] == 1; out[2] ==0; out[3]==0; out[4]<max.grad.threshold; out[5]< not too big;
convergence_fn <- function(em=NULL){
  out <- rep(NA,9) #convergence type  1, and 2
  names(out) <- c("opt", "conv", "sdrep", "max_grad", "SE_par_max", "SE_par_max_name", "SE_par_max2", "SE_par_max2_name", "max_grad_name")
  if(!is.character(em)) {
    if(!is.null(em$fit$opt)) {
      out[1] <- 1 
      out[2] <- em$fit$opt$conv
    }
    if(!is.null(em$fit$sdrep)) out[3] <- as.integer(sum(sapply(em$fit$sdrep$SE_par, function(g) any(is.nan(g)))))
    if(!is.null(em$fit$final_gradient)) out[4] <- max(abs(em$fit$final_gradient))
    if(!is.null(em$fit$sdrep)) {
      maxs <- sapply(em$fit$sdrep$SE_par, function(g) ifelse(any(!is.na(g)), max(g,na.rm=TRUE), NA))
      out[5] <- ifelse(any(!is.na(maxs)), max(maxs, na.rm =TRUE), NA)
      if(!is.na(out[5]) )  {
        out[6] <- names(which(maxs==out[5]))
        out[7] <- sort(maxs, decreasing=TRUE) [2]
        out[8] <- names(sort(maxs, decreasing=TRUE)) [2]
        out[9] <- names(em$fit$opt$par)[which.max( abs(em$fit$final_gradient))]
      }
    }
  }
  return(out)
}

#=======================================================================================
#   Calculate mohn's rho for user-specified number of peels   ====
#=======================================================================================
# (just modifying wham's fn mohns_rho)
mohns_rho_set_peel <- function (model=NULL, npeels=NULL, ny=NULL, na=NULL) 
{
  #npeels = length(model$peels)
  # ny = model$env$data$n_years_model
  # na = model$env$data$n_ages
  ages.lab <- seq(1,na)
  if (npeels) {
    rho = c(mean(sapply(1:npeels, function(x) model$peels[[x]]$rep$SSB[ny - 
                                                                         x]/model$fit$rep$SSB[ny - x] - 1)), mean(sapply(1:npeels, 
                                                                                                                     function(x) model$peels[[x]]$rep$Fbar[ny - x]/model$fit$rep$Fbar[ny - 
                                                                                                                                                                                    x] - 1)))
    names(rho) = c("SSB", "Fbar")
    rho.naa = sapply(1:na, function(y) {
      mean(sapply(1:npeels, function(x) model$peels[[x]]$rep$NAA[ny - 
                                                                   x, y]/model$fit$rep$NAA[ny - x, y] - 1))
    })
    names(rho.naa) = c("R", paste0("N", ages.lab[2:na]))
    rho = c(rho, rho.naa)
    return(rho)
  }
  else stop("There are no peels in this model")
}


#=======================================================================================
#   Calculate mohn's rho for Recruitment random effects for user-specified number of peels   ====
#=======================================================================================
# (modifying my fn_retro_par_estimates.R script in SSRTWG/common_code)
mohns_rho_randeff_peel <- function (model=NULL, npeels=NULL, ny=NULL) 
{
  #npeels = length(model$peels)
  # ny = model$env$data$n_years_model
  
  ny.re = ny-1  # there are ny-1 random effects, first year is mean
  
  if (npeels) {
    rho = mean(sapply(1:npeels, function(x) model$peels[[x]]$rep$NAA_devs[(ny.re -x),1]/model$fit$rep$NAA_devs[(ny.re - x),1] - 1))
    names(rho) = c("R_dev")
    
    return(rho)
  }
  else stop("There are no peels in this model")
}


#=======================================================================================
#   Collect AIC and convergence info across OMs, EMs, Sims  ====
#=======================================================================================

get_aic_convergence_info <- function(df.oms=NULL, df.ems=NULL, nsims=NULL, res.path=NULL, 
                                     save.rds=TRUE, save.suffix=NULL, save.path=NULL)  {
  # df.oms is dataframe of om specs
  # df.ems is dataframe of em specs
  # nsims is the number of simulated iterations of each om
  # res.path is the path to the results directory where each om is a subdirectory (ex:here::here("Ecov_study", 'recruitment_functions', 'results') )
  # save.rds (T/F) save the results as RDS file? (default=TRUE)
  # save.suffix is appended at very end of filename (before extension) to distinguish cases
  # save.path is the path to the directory where you want to save results
  
  if(!dir.exists(save.path))  dir.create(save.path)
  
  nsim <-  nrow(df.oms)*n_sims
  nsim.all <- nrow(df.oms)*nrow(df.ems)*n_sims
  
  
  AIC           <- data.frame(matrix(nrow=nsim, ncol=ncol(df.oms)+20))
  colnames(AIC) <- c('sim',colnames(df.oms),'aic_pick', 'daic_next', 'daic_last', 
                     'correct_form','correct_ecov', 'correct_SR', 'EM_ecov', 'EM_SR',
                     "opt", "conv", "sdrep", "max_grad", "SE_par_max",  "SE_par_max_name",  "SE_par_max2",  "SE_par_max2_name", "max_grad_name","ecov_slope","ssb_cv")
  
  AIC_weight  <- data.frame(matrix(nrow=nsim.all, ncol=(ncol(df.oms)+20)) )
  colnames(AIC_weight) <- c('sim',colnames(df.oms),'OM','EM', 'EM_ecov_how', 'EM_r_mod','AIC', 'dAIC', 'AIC_rank', 'Model_prob',
                            "opt", "conv", "sdrep", "max_grad", "SE_par_max", "SE_par_max_name", "SE_par_max2",  "SE_par_max2_name", "max_grad_name","ecov_slope","ssb_cv")
  

k <- 1
for(om in 1:nrow(df.oms)){
  print(paste0("OM = ",om, " out of ", nrow(df.oms)))
    # get aic ====
  for(sim in 1:n_sims){
  print(paste0("SIM = ",sim, " out of 100"))
    
    DAT <- sapply(1:nrow(df.ems), function(em){
      dat <- tryCatch(readRDS(file.path(res.path, paste0("om", om, '/','sim',sim,'_','em',em,'.RDS')) ),
                      error = function(e) conditionMessage(e))
      aic <- NA
        if(class(dat)!='try-error' & class(dat)!='character'){
          if(length(dat$fit)>0) aic = 2*(as.numeric(dat$fit$opt$objective) + as.numeric(length(dat$fit$opt$par)) )
        }
        aic
      })

    dat <- tryCatch(readRDS(file.path(res.path, paste0("om", om, '/','sim',sim,'_','em',1,'.RDS')) ),
                    error = function(e) conditionMessage(e))
    ecov_slope <- NA
    ssb_cv <- NA
    if(class(dat)!='try-error' & class(dat)!='character'){
      ecov_sim   <- dat$truth$Ecov_obs
      ecov_slope <- summary(lm(ecov_sim ~ I(1:40)))$coefficients[2,1]
      ssb_cv    <- sd(dat$truth$SSB)/mean(dat$truth$SSB)
    }

    # get convergence ====
    conv  <- sapply(1:nrow(df.ems), function(em) {
      cc <- tryCatch(readRDS(file.path(res.path, paste0("om", om, '/','sim',sim,'_','em',em,'.RDS')) ),
                     error = function(e) conditionMessage(e))
      
      conv.info <- rep(NA, 9)  #length of vector returned from convergence_fn
      if(class(cc)!='try-error' & class(cc)!='character'){
        if(length(cc$fit)>0) conv.info = convergence_fn(cc) 
      }
      conv.info
    })
    
        em_match <- which(df.ems$ecov_how==df.oms$Ecov_how[om] &
                            df.ems$r_mod==df.oms$recruit_mod[om])                     #1 if ecov_how AND r_mod match for EM and OM
        sr_good <- ifelse(df.ems$r_mod==df.oms$recruit_mod[om], 1, 0)   #1 if r_mod matches for EM and OM
        ecov_good <- ifelse(df.ems$ecov_how==df.oms$Ecov_how[om], 1, 0) #1 if ecov_how matches for EM and OM
        
        aic_pick <- which(DAT==min(DAT,na.rm=TRUE))
        aic_order <- match(seq(1,nrow(df.ems)), order(DAT) )  #gives rank of each EM
        aic_next <- DAT[aic_order==1] - DAT[aic_order==2]  #what is dAIC between best and second best
        aic_diffs <-  DAT - DAT[aic_order==1]  # dAIC for all EMs
        aic_last <- max(aic_diffs)  #what is dAIC for worst model
        
        aic_weight <- exp(-0.5*aic_diffs)/sum(exp(-0.5*aic_diffs), na.rm=TRUE )
    
    AIC[k,]  <- data.frame(sim=sim, df.oms[om,],aic_pick=aic_pick,
                           dAIC_next = round(aic_next,3), dAIC_last=round(aic_last,3),
                           correct_form=ifelse(aic_pick==em_match,1,0),
                           correct_ecov = ecov_good[aic_pick],
                           correct_SR= sr_good[aic_pick],
                           EM_ecov=df.ems[aic_pick, 1], EM_SR = df.ems[aic_pick,2],
                           conv=t(conv[, aic_pick]),
                           ecov_slope=as.numeric(ecov_slope),
                           ssb_cv=as.numeric(ssb_cv))

    AIC_weight[(nrow(df.ems)*(k-1)+1):(nrow(df.ems)*k) , ] <- data.frame(sim=rep(sim,nrow(df.ems)), 
                                                                         OM=rep(om, nrow(df.ems)), 
                                                                         EM=seq(1,nrow(df.ems)), 
                                                                         EM_ecov_how=df.ems[,1], 
                                                                         EM_r_mod= df.ems[,2],
                                                                         AIC=DAT, 
                                                                         dAIC=aic_diffs, 
                                                                         AIC_rank=aic_order,                                             
                                                                         Model_prob=aic_weight, conv=t(conv) , 
                                                                         do.call("rbind", replicate( nrow(df.ems), 
                                                                                                     df.oms[om,], 
                                                                                                     simplify = FALSE)),
                                                                         ecov_slope=as.numeric(ecov_slope),
                                                                         ssb_cv=as.numeric(ssb_cv))
     k <- k + 1
     
    } # end sim loop
    if(save.rds==T){
      saveRDS(AIC,file.path(save.path, paste0('AIC', save.suffix, '.rds')) )
      saveRDS(AIC_weight,file.path(save.path, paste0('AIC_weight',  save.suffix,'.rds'))  )
    }
  } # end om loop
} # end function

