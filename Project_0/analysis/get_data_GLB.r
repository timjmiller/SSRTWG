#load the libraries
library(tidyverse)
library(googledrive)
library(doParallel)
registerDoParallel(10)


#folder link to id
results <- "https://drive.google.com/drive/u/0/folders/1adYe52gSy9DS8NCGDwhviwhtwFckmzR7"
folder <- drive_get(as_id(results))$id

proc_ids <- drive_ls(folder)$id
proc_nms <- drive_ls(folder)$name

#setwd("~/dropbox/working/state_space_assessments/project0_google2/")
dir <- "~/dropbox/working/state_space_assessments/project0_google2"

for(i in 1:length(proc_ids)){
  dir_i <- file.path(dir,proc_nms[i])
  dir.create(dir_i)
  setwd(dir_i)
  
  om_ids <- drive_ls(proc_ids[i])$id
  om_nms <- drive_ls(proc_ids[i])$name
  for(j in 1:length(om_ids)){
    dir_j <- file.path(dir_i,om_nms[j])
    dir.create(dir_j)
    setwd(dir_j)
    
    sim_ids <- drive_ls(om_ids[j])$id
    sim_nms <- drive_ls(om_ids[j])$name
    
    foreach(k=1:length(sim_ids)) %dopar% {
      drive_download(sim_ids[k],overwrite=TRUE)
    }      
#    for(k in 1:length(sim_ids)){
#      drive_download(sim_ids[k],overwrite=TRUE)
#    }
  }
}

