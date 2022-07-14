library(here)
library(wham)
oms = list.files(file.path(here(),"Project_0", "results", "naa_om"))
oms = grep("^om", oms, value = TRUE) #begins with "om"
df.ems = readRDS(file.path(here(),"Project_0","inputs", "df.ems.RDS"))
df.oms = readRDS(file.path(here(),"Project_0","inputs", "df.oms.RDS"))

for(i in 1:length(oms)){
  sim1 = lapply(1:20, function(x){
    em = try(readRDS(file.path(here(),"Project_0", "results", "naa_om", paste0("naa_om_results_om_",i,"_em_", x, "_sim_1_test.RDS"))))
    if(is.character(em)) em = NULL
    return(em)
  })
  saveRDS(sim1, file.path(here(),"Project_0", "results", "naa_om", paste0("om_",i), "sim_1.RDS"))
}
#om1_sim1 = readRDS(file.path(here(),"Project_0", "results", "naa_om", "om_1", "sim_1.RDS"))
substr(strsplit(om1_sim1[[1]]$fit$wham_version, "@", fixed = TRUE)[[1]][2], 1,7)

om1_sim2 =readRDS(file.path(here::here(),"Project_0", "results", "naa_om", "om_1", "sim_2.RDS"))
substr(strsplit(om1_sim2[[1]]$fit$wham_version, "@", fixed = TRUE)[[1]][2], 1,7)

job.sheet = readRDS(file.path(here::here(),"Project_0", "inputs", "naa.sim.jobs.RDS"))
job.sheet$member[job.sheet$sim %in% 1:5] = "TJM"
for(i in 1:24) for(j in 1:5){
  tfile = file.path(here::here(),"Project_0", "results", "naa_om", paste0("om_",i), paste0("sim_", j, ".RDS"))
  print(tfile)
  temp = readRDS(tfile)
  temp = temp[!sapply(temp, is.null)]
  job.sheet$wham_commit[job.sheet$sim == j & job.sheet$om == paste0("om_",i)] <- 
  substr(strsplit(temp[[1]]$fit$wham_version, "@", fixed = TRUE)[[1]][2], 1,7)
}
unique(job.sheet$wham_commit) #just 97577f1
saveRDS(job.sheet, file.path(here(),"Project_0", "inputs", "naa.sim.jobs.RDS"))

n_ems_complete = sapply(oms, function(x){
  sims = list.files(file.path(here::here(),"Project_0", "results", "naa_om", x))
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

#checking that the keep indicators are the same for all of the oms and ems
em_inputs = readRDS(file.path(here(),"Project_0","inputs", "em_inputs.RDS"))
om_inputs = readRDS(file.path(here(),"Project_0","inputs", "NAA_om_inputs.RDS"))
keep_names = c("keep_E", "keep_C", "keep_I", "keep_Ipaa", "keep_Cpaa")
keep_names = c("keep_C", "keep_I", "keep_Ipaa", "keep_Cpaa")
for(i in 2:length(om_inputs)) for(j in 1:(i-1)){
  for(k in keep_names){
    temp <- om_inputs[[i]]$data[[k]] == om_inputs[[j]]$data[[k]]
    if(sum(!temp)) print(paste0("for ", k, ", i: ", i, " j: ", j, " sum !=", sum(!temp)))
  }
}
for(i in 1:length(om_inputs)) for(j in 1:length(em_inputs)){
  for(k in keep_names){
    temp <- om_inputs[[i]]$data[[k]] == em_inputs[[j]]$data[[k]]
    if(sum(!temp)) print(paste0("for ", k, ", i: ", i, " j: ", j, " sum !=", sum(!temp)))
  }
}


# sims = list.files((here(),"Project_0", "results", "naa_om", paste0("om_", this_om))

#   write.dir <- file.path(here(),"Project_0", "results", "naa_om", paste0("om_", this_om))
#   dir.create(write.dir, recursive = T, showWarnings = FALSE)
#   rds.fn = file.path(write.dir, paste0("sim_", this_sim, ".RDS"))
#   saveRDS(list(), rds.fn) #make list file that can be populated with the em fits.
