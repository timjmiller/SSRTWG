get_SPR = function(F, M, sel, mat, waassb, fracyrssb, at.age = FALSE)
{
  n_ages = length(sel)
  SPR = numeric()
  n = 1
  F = F * sel
  Z = F + M
  for(a in 1:(n_ages-1))
  {
    SPR[a] = n[a] * mat[a] * waassb[a] * exp(-fracyrssb * Z[a])
    n[a+1] = n[a] * exp(-Z[a])
  }
  n[n_ages] = n[n_ages]/(1-exp(-Z[n_ages]))
  SPR[n_ages] = n[n_ages] * mat[n_ages] * waassb[n_ages] * exp(-fracyrssb * Z[n_ages])
  if(at.age) return(SPR)
  else return(sum(SPR))
}


#-------Y/R -----------------------------
get_YPR_fleets = function(FAA,M, waacatch, at.age = FALSE)
{
  n_ages = dim(FAA)[2]
  n_fleets <- dim(FAA)[1]
  YPR = matrix(0,n_fleets, n_ages)
  #F = F * sel
  Z = apply(FAA,2,sum) + M
  for(f in 1:n_fleets) {
    n = 1
    for(a in 1:(n_ages-1)){
      YPR[f,a] <- n[a] * FAA[f,a] * waacatch[f,a] * (1.0 - exp(-Z[a]))/Z[a]
      n[a+1] = n[a] * exp(-Z[a])
    }
    n[n_ages] = n[n_ages]/(1 - exp(-Z[n_ages]))
    YPR[f,n_ages] = n[n_ages] * FAA[f,n_ages] * waacatch[f,n_ages] * (1.0 - exp(-Z[n_ages]))/Z[n_ages]
  }
  if(at.age) return(YPR) #by fleet and age
  else return(apply(YPR,1,sum)) #by fleet
}

get_SPR_BRPS_fn <- function(mod, spr_yrs, percent){
  dat = mod$env$data
  if(missing(percent)) percent <- dat$percentSPR
  if(missing(spr_yrs)) spr_yrs <- dat$avg_years_ind+1 #c++
  R_yrs <- dat$XSPR_R_avg_yrs+1
  R_type <- dat$XSPR_R_opt
  if(R_type %in% c(1,3)) Rxspr <- mean(mod$rep$NAA[Ryrs,1])
  else Rxspr <- mean(mod$rep$pred_NAA[R_yrs,1])
  n_ages<- dat$n_ages
  years <- mod$years
  n_years <- length(years)
  n_fleets <- dat$n_fleets

  mat.age <- apply(dat$mature[spr_yrs,],2,mean)
  ssb.waa <- apply(dat$waa[dat$waa_pointer_ssb,,][spr_yrs,],2,mean)
  catch.waa <- apply(dat$waa[dat$waa_pointer_fleets,spr_yrs,,drop = FALSE],c(1,3),mean)
  M.age <- apply(mod$rep$MAA[spr_yrs,],2,mean)
  FAA <- apply(mod$rep$FAA[spr_yrs,,,drop=FALSE], 2:3, mean) #mean FAA by fleet
  FAAtot = apply(FAA,2,sum) #sum across fleet averages 
  seltot <- FAAtot/max(FAAtot)
  selAA <- FAA/max(FAAtot) 
  spawn.time <- mean(dat$fracyr_SSB[spr_yrs])
  spr0 = get_SPR(F=0, M=M.age, sel=seltot, mat=mat.age, waassb=ssb.waa, fracyrssb = spawn.time)
  F.start <- 0.11  # starting guess for optimization routine to find F_SPR%
  spr.f <- function(F.start) {
    spr = get_SPR(F=F.start, M=M.age, sel=seltot, mat=mat.age, waassb=ssb.waa, fracyrssb = spawn.time)
    abs(100*spr/spr0 - percent)
  }
  opt <- nlminb(start=F.start, objective=spr.f, lower=0, upper=10)
  Fxspr <- opt$par
  spr_Fxspr <- get_SPR(Fxspr, M=M.age, sel=seltot, mat=mat.age, waassb=ssb.waa, fracyrssb = spawn.time)
  FAA_xspr <- Fxspr * selAA
  ypr_Fxspr <- get_YPR_fleets(FAA = FAA_xspr, M=M.age, waacatch= catch.waa)
  Y_Fxspr <- Rxspr * ypr_Fxspr
  SSB_Fxspr <- Rxspr*spr_Fxspr
  return(list(Fxspr = Fxspr, FAA_xspr = FAA_xspr, SSB_Fxspr = SSB_Fxspr, Y_Fxspr = Y_Fxspr, spr_Fxspr = spr_Fxspr, Rxspr = Rxspr))
}
