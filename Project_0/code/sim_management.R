verify_version = function(){
  required_wham_version <-  "1.0.6 / Github (timjmiller/wham@97577f1)" 
  required_commit = substr(strsplit(required_wham_version, "@", fixed = TRUE)[[1]][2], 1,7)
  ver <- sessioninfo::package_info() %>% as.data.frame %>% dplyr::filter(package=="wham") %>% dplyr::select(loadedversion, source) %>% unname
  wham_version <- paste0(ver, collapse=" / ")
  wham_commit = substr(strsplit(wham_version, "@", fixed = TRUE)[[1]][2], 1,7)

  if(wham_commit != required_commit) {
    stop(paste0("your wham version:", wham_version, "is not the required version:", 
    required_wham_version, ".\n", "Install the right version using \n",
    "devtools::install_github('timjmiller/wham', dependencies=TRUE, ref='97577f1') \n"))
  } else{
    cat("The right version (commit 97577f1) of wham is loaded! \n")
  }
}
#verify_version()

checkout_jobs = function(sims, member = "TJM", write.job.sheet = FALSE){
  job.sheet = readRDS(file.path(here::here(),"Project_0", "inputs", "naa.sim.jobs.RDS"))
  these.tasks = filter(job.sheet, sim %in% sims)
  if(any(!is.na(these.tasks$member))) stop("some of the specified sims have already been checked out")
  job.sheet$member[job.sheet$sim %in% sims] <- member
  if(write.job.sheet) saveRDS(job.sheet, file.path(here::here(),"Project_0", "inputs", "naa.sim.jobs.RDS"))
  return(job.sheet)
}

#job.sheet = checkout_jobs(sims = 6)

run_jobs = function(sims, n.cores = NULL){
  script.full.path = file.path(here::here(),"Project_0", "code", "naa_om_sim_fit_script.R")
  results.path = file.path(here::here(),"Project_0", "results", "naa_om")
  df.ems = readRDS(file.path(here::here(),"Project_0","inputs", "df.ems.RDS"))
  df.oms = readRDS(file.path(here::here(),"Project_0","inputs", "df.oms.RDS"))
  om_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "NAA_om_inputs.RDS"))
  em_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "em_inputs.RDS"))
  #######################################################
  #need to have matching assumptions about CVs for catch and indices, too
  obs_names = c("agg_catch","agg_catch_sigma", "agg_indices", "agg_index_sigma", "catch_paa", "index_paa", 
    "Ecov_obs", "obs", "obsvec")
  #######################################################
  seeds = readRDS(file.path(here::here(), "Project_0", "inputs","seeds.RDS"))
  n_oms = length(om_inputs)
  #library(snowfall) # used for parallel computing
  parallel::detectCores()
  if(is.null(n.cores)) n.cores = parallel::detectCores()/2
  snowfall::sfInit(parallel=TRUE, cpus=n.cores)
  oms = 1:n_oms
  temp = expand.grid(om = oms, sim = sims)
  snowfall::sfExportAll()
  snowfall::sfLapply(1:NROW(temp), function(row_i){
    this_om = temp$om[row_i]
    this_sim = temp$sim[row_i]
    write.dir <- file.path(here::here(),"Project_0", "results", "naa_om", paste0("om_", this_om))
    dir.create(write.dir, recursive = T, showWarnings = FALSE)
    rds.fn = file.path(write.dir, paste0("sim_", this_sim, ".RDS"))
    saveRDS(list(), rds.fn) #make list file that can be populated with the em fits.
    for(this_em in 1:length(em_inputs)){
      system(paste0("Rscript --vanilla ", script.full.path, " " , this_om, " ",  this_em, " ", this_sim, " \n"))
    }
  })
}

update_job_sheet_commits = function(sims, job.sheet){
  for(i in 1:24) for(j in sims){
    tfile = file.path(here::here(),"Project_0", "results", "naa_om", paste0("om_",i), paste0("sim_", j, ".RDS"))
    print(tfile)
    res = try(readRDS(tfile)) 
    if(is.character(re)) stop(paste0(tfile, "is not in the results directory"))
    temp = temp[!sapply(temp, is.null)]
    job.sheet$wham_commit[job.sheet$sim == j & job.sheet$om == paste0("om_",i)] <- 
      substr(strsplit(temp[[1]]$fit$wham_version, "@", fixed = TRUE)[[1]][2], 1,7)
  }
  return(job.sheet)
}

run_hpcc_jobs = function(this_sim, this_om, this_em){
  script.full.path = file.path(here::here(),"Project_0", "code", "naa_om_sim_fit_script_hpcc.R")
  results.path = file.path(here::here(),"Project_0", "results", "naa_om")
  df.ems = readRDS(file.path(here::here(),"Project_0","inputs", "df.ems.RDS"))
  df.oms = readRDS(file.path(here::here(),"Project_0","inputs", "df.oms.RDS"))
  om_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "NAA_om_inputs.RDS"))
  em_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "em_inputs.RDS"))
  seeds = readRDS(file.path(here::here(), "Project_0", "inputs","seeds.RDS"))
  #######################################################
  #need to have matching assumptions about CVs for catch and indices, too
  obs_names = c("agg_catch","agg_catch_sigma", "agg_indices", "agg_index_sigma", "catch_paa", "index_paa", 
    "Ecov_obs", "obs", "obsvec")
  #######################################################
  write.dir <- file.path(here::here(),"Project_0", "results", "naa_om")
  dir.create(write.dir, recursive = T, showWarnings = FALSE)
  rds.fn = file.path(write.dir, paste0("om", this_om, "_sim", this_sim, "_em", this_em, ".RDS"))
  saveRDS(lapply(ems, function(x) NULL), rds.fn) #make list file that can be populated with the em fits.
  system(paste0("Rscript --vanilla ", script.full.path, " " , this_om, " ",  this_em, " ", this_sim, " \n"))
}