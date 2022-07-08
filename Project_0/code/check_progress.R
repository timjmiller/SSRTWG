library(here)
library(wham)
oms = list.files(file.path(here(),"Project_0", "results", "naa_om"))
oms = grep("^om", oms, value = TRUE) #begins with "om"
n_ems_complete = sapply(oms, function(x){
  sims = list.files(file.path(here(),"Project_0", "results", "naa_om", x))
  print(sims)
  sapply(sims, function(y){
    ems = readRDS(file.path(here(),"Project_0", "results", "naa_om", x, y))
    return(length(ems))
  })})

n_ems_null = sapply(oms, function(x){
  sims = list.files(file.path(here(),"Project_0", "results", "naa_om", x))
  print(sims)
  em_out = sapply(sims, function(y){
    ems = readRDS(file.path(here(),"Project_0", "results", "naa_om", x, y))
    out = rep(NA,20)
    
    if(length(ems)) for(i in 1:length(ems)) out[i] = ifelse(is.null(ems[[i]]),1,0)
    return(out)
  })
  print(em_out)
  return(apply(em_out,1,sum, na.rm = TRUE))
})
n_ems_null = n_ems_null[,match(df.oms$Model,colnames(n_ems_null))]

apply(n_ems_null,1,sum)
cbind(df.oms, n_null = n_ems_null[8,match(df.oms$Model,colnames(n_ems_null))])

cbind(df.ems, n_null = apply(n_ems_null,1,sum))

cbind(df.oms, n_null = n_ems_null[8,match(df.oms$Model,colnames(n_ems_null))])


n_ems_sdrep = sapply(oms, function(x){
  sims = list.files(file.path(here(),"Project_0", "results", "naa_om", x))
  print(sims)
  em_out = sapply(sims, function(y){
    ems = readRDS(file.path(here(),"Project_0", "results", "naa_om", x, y))
    out = rep(0,20)
    if(length(ems)) for(i in 1:length(ems)) if(!is.null(ems[[i]])) {
      if(sum(length(ems[[i]]$fit$sdrep)==4)) out[i] = 1
    }
    return(out)
  })
  return(apply(em_out,1,sum, na.rm = TRUE))
})
n_ems_sdrep = n_ems_sdrep[,match(df.oms$Model,colnames(n_ems_sdrep))]
save(n_ems_sdrep, n_ems_null, 
  file =  file.path(here(),"Project_0", "results", "prelim_res.RData"))

load(file.path(here(),"Project_0", "results", "prelim_res.RData"))

all_ems_sdrep = n_ems_sdrep/4



sims = list.files((here(),"Project_0", "results", "naa_om", paste0("om_", this_om))

  write.dir <- file.path(here(),"Project_0", "results", "naa_om", paste0("om_", this_om))
  dir.create(write.dir, recursive = T, showWarnings = FALSE)
  rds.fn = file.path(write.dir, paste0("sim_", this_sim, ".RDS"))
  saveRDS(list(), rds.fn) #make list file that can be populated with the em fits.
