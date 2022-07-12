source(file.path(here::here(), "Project_0","code", "sim_management.R")
library(wham)
library(snowfall)
verify_version()

job.sheet = checkout_jobs(sims = 6, member = "TJM")

run_jobs(sims = 6, n.cores = 4)

job.sheet = update_job_sheet_commits(sims = 6, job.sheet= job.sheet)
