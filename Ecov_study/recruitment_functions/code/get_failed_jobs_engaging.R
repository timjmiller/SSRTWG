library(here)

df.ems    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.ems.RDS"))
df.oms    <- readRDS(file.path(here::here(),"Ecov_study", "recruitment_functions", "inputs", "df.oms.RDS"))

n.em <- nrow(df.ems)
n.om <- nrow(df.oms)
om.names <- paste0("om", seq(1,n.om))
iters <- seq(1,100)
n.iter <- length(iters)
n.iter.em <- n.iter*n.em

sim.names <- paste0(rep("sim",n.iter.em), rep(iters, each=n.em), "_em", rep(seq(1,n.em), n.iter), ".RDS")

om.fails <- c()
em.fails <- c()
iter.fails <- c()

fail.list <- list()

for (i in 1:n.om) {
print(i)
  sims.run <- list.files(file.path(here::here(),"Ecov_study", "recruitment_functions", "results", om.names[i]), pattern = ".RDS")
  fails <- which((sim.names %in% sims.run)==FALSE)
print(fails)  
  if(length(fails)==0) fail.list[[i]] <- NA
  
  if(length(fails>0)) {
    om.fails <- as.numeric(substr(om.names[i], 3, nchar(om.names[i]) ) )
    tmp <- strsplit(sim.names[fails], split="_")
    iter.fails <- as.numeric(c(unlist(lapply(tmp, function(x) {substr(x[1], 4, nchar(x[1]) ) }))) )
    em.fails <- as.numeric(c(unlist( lapply(tmp, function(x) {substr(x[2], 3, 3 ) }))  ) )
    
    fail.list[[i]] <- list(om=om.fails, nfails = length(tmp) , iter_em=cbind(iter.fails, em.fails))
    
  }
  om.fails <- c()
  em.fails <- c()
  iter.fails <- c()
  
}

###-make DF
om.fails.df <- c()  #
for(i in 1:nrow(df.oms)) {
print(i)
  if(is.na(fail.list[[i]][1])==FALSE) om.fails.df <- rbind(om.fails.df,
                                                           cbind(OM=    rep(i, fail.list[[i]]$nfails), 
                                                                 N.fail=rep(fail.list[[i]]$nfails, fail.list[[i]]$nfails), 
                                                                 fail.list[[i]]$iter_em) )
}

saveRDS(om.fails.df, file.path(here::here("Ecov_study", 'recruitment_functions', 'results'),'om.fails.df.RDS'))


