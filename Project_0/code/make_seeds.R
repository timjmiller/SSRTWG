
set.seed(8675309)
#naa_oms
df.oms <- readRDS(file.path(here::here(),"Project_0", "inputs", "df.oms.RDS"))
noms <- NROW(df.oms)
seeds = sample(x = (-1e9):(1e9), size = noms*1000, replace = FALSE)
seeds <- lapply(1:noms, function(x) seeds[(1:1000) + 1000*(x-1)])
saveRDS(seeds, file.path(here::here(), "Project_0", "inputs","naa_seeds.RDS"))

#M_oms

df.oms <- readRDS(file.path(here::here(),"Project_0", "inputs", "df.M.oms.RDS"))
noms <- NROW(df.oms)
seeds = sample(x = (-1e9):(1e9), size = noms*1000, replace = FALSE)
seeds <- lapply(1:noms, function(x) seeds[(1:1000) + 1000*(x-1)])
saveRDS(seeds, file.path(here::here(), "Project_0", "inputs","M_seeds.RDS"))

#sel_oms

df.oms <- readRDS(file.path(here::here(),"Project_0", "inputs", "df.sel.oms.RDS"))
noms <- NROW(df.oms)
seeds = sample(x = (-1e9):(1e9), size = noms*1000, replace = FALSE)
seeds <- lapply(1:noms, function(x) seeds[(1:1000) + 1000*(x-1)])
saveRDS(seeds, file.path(here::here(), "Project_0", "inputs","sel_seeds.RDS"))

#q_oms

df.oms <- readRDS(file.path(here::here(),"Project_0", "inputs", "df.q.oms.RDS"))
noms <- NROW(df.oms)
seeds = sample(x = (-1e9):(1e9), size = noms*1000, replace = FALSE)
seeds <- lapply(1:noms, function(x) seeds[(1:1000) + 1000*(x-1)])
saveRDS(seeds, file.path(here::here(), "Project_0", "inputs","q_seeds.RDS"))


df.oms <- readRDS(file.path(here::here(),"Project_0", "inputs", "df.oms.RDS"))
