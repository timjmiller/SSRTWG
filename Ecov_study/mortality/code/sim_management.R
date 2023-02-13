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
