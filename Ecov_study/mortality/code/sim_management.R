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
  bad_logs = system('grep -rn  --include=*.log -L "Success" ~/logs', intern = TRUE)
  fn <- sapply(strsplit(bad_logs, "/"), function(x) x[5])
  om <- as.integer(sapply(strsplit(fn, "_"), function(x) x[which(x == "om")+1]))
  em <- sapply(strsplit(fn, "_"), function(x) x[length(x)])
  em <- as.integer(unlist(strsplit(em, ".log")))
  first_sim <- as.integer(sapply(strsplit(fn, "_"), function(x) x[which(x == "sims")+1]))
  last_sim <- as.integer(sapply(strsplit(fn, "_"), function(x) x[which(x == "to")+1]))
  sims_complete = lapply(1:length(om), function(i) system(paste0('find ~/SSRTWG/Ecov_study/mortality/results/om', om[i], '/', ' -name "*em', em[i], '.RDS"'), intern = TRUE))
  sims_complete <- lapply(1:length(om), function(i) sort(as.integer(sapply(strsplit(sims_complete[[i]], "sim"), function(x) strsplit(x[2], "_")[[1]][1]))))
  sims_redo <- lapply(1:length(sims_complete), function(i) setdiff(first_sim[i]:last_sim[i], sims_complete[[i]]))

  return(list(bad_jobs= cbind(om = om, em = em, first = first_sim, last = last_sim), redo = sims_redo))
}

make_redo_failed_jobs_file <- function(failed_jobs, fn = paste0("~/SSRTWG/Ecov_study/mortality/code/redo_failed_jobs.txt")){

  redo = failed_jobs$redo
  bad_jobs <- failed_jobs$bad_jobs
  if(file.exists(fn)) stop("the file to write commands to already exists \n \n")
  cat("Redoing these failed jobs \n \n", file = fn, append = F)

  for(i in 1:length(redo)){
    sims <- redo[[i]]
    om <- bad_jobs[i,"om"]
    em <- bad_jobs[i,"em"]
    for(j in 1:length(sims)){
      cat(paste0('bsub -n 1 -q long -W 24:00 -o ~/logs/sim_', sims[j], '_om_', om, '_em_', em, 
        '.log -R "rusage[mem=5000]" -R "span[hosts=1]" -J ', sims[j], '_', om, "_", em, 
        " bash ~/SSRTWG/Ecov_study/mortality/code/M_Ecov_om_hpcc_args.sh ", sims[j], " " , sims[j], " ", om, " ", om, " ", em, " ", em, "\n"), file = fn, append = TRUE)
    }
  }
}
