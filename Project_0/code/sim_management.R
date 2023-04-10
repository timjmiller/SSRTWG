verify_version = function(commit = "77bbd94"){
  if(!"wham" %in% .packages()) stop("using verify_version: wham library is not loaded")
  desc = utils::packageDescription("wham")
  if(is.null(desc$GithubSHA1)) stop("using verify_version: wham library was not installed from Github")
  wham_commit = substr(desc$GithubSHA1,1,7)
  wham_version = paste0(desc$Version, ", commit: ", wham_commit)

  if(wham_commit != commit) {
    stop(paste0("your wham commit:", wham_version, "is not the required commit:", 
    commit, ".\n", "Install the right commit using \n",
    "devtools::install_github('timjmiller/wham', dependencies=TRUE, ref=", commit, ") \n"))
  } else{
    cat(paste0("The right commit: ",commit, " of wham is loaded! \n"))
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
  ems = 1:length(em_inputs)
  rds.fn = file.path(write.dir, paste0("om", this_om, "_sim", this_sim, "_em", this_em, ".RDS"))
  saveRDS(lapply(ems, function(x) NULL), rds.fn) #make list file that can be populated with the em fits.
  system(paste0("Rscript --vanilla ", script.full.path, " " , this_om, " ",  this_em, " ", this_sim, " \n"))
}

aggregate_hpcc_results = function(sim, oms, ems = 1:20, res_dir = file.path(here::here(),"Project_0", "results", "naa_om"))
{
  for(i in oms){
    print(paste0("om: ", i, ", sim: ", sim))
    write_dir = file.path(res_dir,paste0("om_",i))
    print(write_dir)
    dir.create(write_dir, recursive = T, showWarnings = FALSE)
    sim_aggregated = lapply(ems, function(x){
      em = readRDS(file.path(res_dir, paste0("om",i,"_sim", sim, "_em", x, ".RDS")))[[x]]
      return(em)
    })
    saveRDS(sim_aggregated, file.path(write_dir, paste0("sim_", sim,".RDS")))
  }
}

#this should be usable for all the operating models now
run_hpcc_jobs_rev = function(this_sim, this_om, this_em, 
  script.full.path = file.path(here::here(),"Project_0", "code", "naa_om_sim_fit_script_hpcc.R"), 
  df.ems = readRDS(file.path(here::here(),"Project_0","inputs", "df.ems.RDS")), 
  df.oms = readRDS(file.path(here::here(),"Project_0","inputs", "df.oms.RDS")), 
  om_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "NAA_om_inputs.RDS")),
  em_inputs = readRDS(file.path(here::here(),"Project_0","inputs", "em_inputs.RDS")),
  seeds = readRDS(file.path(here::here(), "Project_0", "inputs","naa_seeds.RDS")),
  write.dir) #= file.path(here::here(),"Project_0", "results", "naa_om")
{
  seeds = readRDS(file.path(here::here(), "Project_0", "inputs","seeds.RDS"))
  #######################################################
  #need to have matching assumptions about CVs for catch and indices, too
  obs_names = c("agg_catch","agg_catch_sigma", "agg_indices", "agg_index_sigma", "catch_paa", "index_paa", 
    "Ecov_obs", "obs", "obsvec")
  #######################################################
  if(is.null(write.dir)) stop("in run_hpcc_jobs_rev(), write.dir is not defined")
  dir.create(write.dir, recursive = T, showWarnings = FALSE)
  ems = 1:length(em_inputs)
  rds.fn = file.path(write.dir, paste0("om", this_om, "_sim", this_sim, "_em", this_em, ".RDS"))
  saveRDS(lapply(ems, function(x) NULL), rds.fn) #make list file that can be populated with the em fits.
  system(paste0("Rscript --vanilla ", script.full.path, " " , this_om, " ",  this_em, " ", this_sim, " \n"))
}

get_failed_jobs = function(){
  bad_logs = system('grep -rn  --include=logfile.* -L "Success" ~/logs', intern = TRUE)
  linebefore = "# LSBATCH: User input"
  jobs = sapply(bad_logs, function(y){
    x = readLines(y)
    bashline = which(x == linebefore)+1 
    if(length(bashline)) strsplit(x[bashline], "/")[[1]][5]
  })
  return(jobs)
}

get_failed_jobs = function(om_type = "naa", ems = 1:20){
  bad_logs <- system('grep -rn  --include=*.log -L "Success" ~/logs', intern = TRUE)
  fn <- sapply(strsplit(bad_logs, "/"), function(x) x[5])
  om <- as.integer(sapply(strsplit(fn, "_"), function(x) x[which(x == "om")+1]))
  sim <- as.integer(sapply(strsplit(fn, "_"), function(x) x[which(x == "sim")+1]))
  #em <- sapply(strsplit(fn, "_"), function(x) x[length(x)])
  #em <- as.integer(unlist(strsplit(em, ".log")))
  #first_em <- as.integer(sapply(strsplit(fn, "_"), function(x) x[which(x == "ems")+1]))
  #last_em <- as.integer(sapply(strsplit(fn, "_"), function(x) x[which(x == "to")+1]))
  em_files = lapply(1:length(om), function(i) system(paste0('find ~/SSRTWG/Project_0/results/', om_type, '_om -name "om', om[i], '_sim', sim[i], '*.RDS"'), intern = TRUE))
  #ems_complete <- lapply(1:length(om), function(i) sort(as.integer(sapply(strsplit(ems_complete[[i]], "em"), function(x) strsplit(x[2], ".RDS")[[1]][1]))))
  ems_complete <- lapply(em_files, function(i) sapply(i, function(j) {
    print(j)
    x <- readRDS(j)
    em <- as.integer(strsplit(strsplit(j, "em")[[1]][2], ".RDS")[[1]][1])
    print(em)
    out <- NULL
    print(length(x))
    print(length(x[[em]]))
    if(length(x[[em]])) out <- em
    return(out)
  }))
  ems_redo <- lapply(1:length(ems_complete), function(i) setdiff(ems, ems_complete[[i]]))

  return(list(bad_jobs= cbind(om = om, sim = sim), redo = ems_redo))
}

make_redo_failed_jobs_file <- function(failed_jobs, om_type = "naa", fn = paste0("~/SSRTWG/Project_0/code/redo_failed_jobs.txt")){

  redo = failed_jobs$redo
  bad_jobs <- failed_jobs$bad_jobs
  if(file.exists(fn)) stop("the file to write commands to already exists \n \n")
  cat("Redoing these failed jobs \n \n", file = fn, append = F)

  for(i in 1:length(redo)){
    ems <- redo[[i]]
    print(ems)
    om <- bad_jobs[i,"om"]
    sim <- bad_jobs[i,"sim"]
    for(j in 1:length(ems)){
      cat(paste0('bsub -n 1 -q long -W 24:00 -o ~/logs/sim_', sim, '_om_', om, '_em_', ems[j], 
        '.log -R "rusage[mem=5000]" -R "span[hosts=1]" -J ', sim, '_', om, "_", ems[j], 
        " bash ~/SSRTWG/Project_0/code/", om_type, "_om_hpcc_args.sh ", sim, " " , sim, " ", om, " ", om, " ", ems[j], " ", ems[j], "\n"), file = fn, append = TRUE)
    }
  }
}
#x <- get_failed_jobs()
#make_redo_failed_jobs_file(x)
#x <- get_failed_jobs("M",5:24)
#make_redo_failed_jobs_file(x, "M")
