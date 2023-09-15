msy_res_fn <- function(df.oms, df.ems, om_type = "naa"){
  all_out <- lapply(1:NROW(df.oms), function(y){
    res <- lapply(1:100, function(x){
      print(paste0("om_", y, ", sim_",x))
      sim = readRDS(file.path(here::here(),"Project_0", "results", paste0(om_type,"_om"), paste0("om_", y), paste0("sim_",x,".RDS")))
      sim <- sim[which(df.ems$SR_model == 3)]
      msy_res <- lapply(sim,function(z) { #across ems
        out <- matrix(NA, 4,3) #a, b, Fmsy_static, SSBmsy_static
        if(length(z)) if(length(z$fit)) {
          out[,1] <- z$fit$rep$SSB/z$truth$SSB
        }
        if(!is.null(z$fit$sdrep)){
          #ind <- which(!is.na(z$fit$sdrep$SE_rep$log_SSB))
          if(res == "ssb") {
            out[,2] <- z$fit$sdrep$SE_rep$log_SSB #se
            ci <- matrix(z$fit$sdrep$Estimate_rep$log_SSB - log(z$truth$SSB),40,2) + qnorm(0.975) * cbind(-out[,2],out[,2])
          }
          if(res == "R") {
            out[,2] <- z$fit$sdrep$SE_rep$log_NAA_rep[,1] #se
            ci <- matrix(z$fit$sdrep$Estimate_rep$log_NAA_rep[,1] - log(z$truth$NAA[,1]), 40,2) + qnorm(0.975) * cbind(-out[,2],out[,2])
          }
          if(res == "F") {
            # print(z$fit$sdrep$SE_rep$log_F[,1])
            # print(length(z$fit$sdrep$SE_rep$log_F[,1]))
            # print(dim(out))
            out[,2] <- z$fit$sdrep$SE_rep$log_F[,1] #se
            ci <- matrix(z$fit$sdrep$Estimate_rep$log_F[,1] - log(z$truth$F[,1]), 40,2) + qnorm(0.975) * cbind(-out[,2],out[,2])
          }
          out[,3] <- 0 >= ci[,1] & 0 <= ci[,2]
        }
        return(out)
      })
      return(relres)
    })
    return(res)
  })
  return(all_out)
}

sim = readRDS(here("Project_0", "results", paste0("naa","_om"), paste0("om_", 1), paste0("sim_",50,".RDS")))
