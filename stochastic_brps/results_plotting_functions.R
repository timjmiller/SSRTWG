get_SSB <- function(res, median_NAA = NULL, total = TRUE, pred = FALSE){
    zaa <- res$FAA_tot + res$MAA
    S <- exp(-zaa * res$fracyr_SSB)
    if(is.null(res$n_ages_pop)) res$n_ages_pop <- res$n_ages
    ind <- c(1:res$n_ages, rep(res$n_ages, res$n_ages_pop-res$n_ages))
    ssbaa <- res$NAA * S *res$waa[res$waa_pointer_ssb,,ind] * res$mature[,ind]
    if(pred) ssbaa <- res$pred_NAA * S *res$waa[res$waa_pointer_ssb,,ind] * res$mature[,ind]
    if(!is.null(median_NAA)) ssbaa <- median_NAA * S *res$waa[res$waa_pointer_ssb,,ind] * res$mature[,ind]
    out <- ssbaa
    if(total) out <- apply(out,1,sum)
    return(out)
}

get_SPR = function(F, M, sel, mat, waassb, fracyrssb, at.age = FALSE, marg_sig = NULL)
{
  n_ages = length(sel)
  SPR = numeric()
  n = 1
  F = F * sel
  Z = F + M
  for(a in 1:(n_ages-1)) {
    SPR[a] = n[a] * mat[a] * waassb[a] * exp(-fracyrssb * Z[a])
    n[a+1] = n[a] * exp(-Z[a])
    if(!is.null(marg_sig)) n[a+1] <- n[a+1] * exp(-0.5*marg_sig[a+1]^2)
  }
  S_A <- exp(-Z[n_ages])
  if(!is.null(marg_sig)) S_A <- S_A * exp(-0.5*marg_sig[n_ages]^2)
  n[n_ages] = n[n_ages]/(1-S_A)
  SPR[n_ages] = n[n_ages] * mat[n_ages] * waassb[n_ages] * exp(-fracyrssb * Z[n_ages])
  if(!at.age) SPR <- sum(SPR)
  return(SPR)
}

get_SSBAA_eq <- function(res){
    fullF <-sapply(1:NROW(res$FAA_tot), function(x) res$FAA_tot[x,res$which_F_age[x]])
    sel <- res$FAA_tot/fullF
    if(is.null(res$n_ages_pop)) res$n_ages_pop <- res$n_ages
    ind <- c(1:res$n_ages, rep(res$n_ages, res$n_ages_pop-res$n_ages))
    mat <- res$mature[,ind]
    waa <- res$waa[res$waa_pointer_ssb,,ind]
    Rhat <- exp(res$mean_rec_pars)
    if(res$bias_correct_pe==1){
      print("bias-correcting")
      logit_rhos <- res$trans_NAA_rho
      #[1] 0.3474890 0.8860973
      rhos <- -1 + 2/(1 + exp(-2*logit_rhos))
      marg_sig <- exp(res$log_NAA_sigma)/sqrt(prod(1-rhos^2)) #marginal standard deviation or the 2dar1 process
      if(res$n_NAA_sigma > 1){
        marg_sig <- marg_sig[res$NAA_sigma_pointers]
      } else{
        marg_sig <- c(marg_sig, rep(0,res$n_ages_pop))
      }
      spr <- t(sapply(1:NROW(waa), function(x) get_SPR(fullF[x], res$MAA[x,], sel[x,], mat[x,],waa[x,], res$fracyr_SSB[x], at.age = TRUE, marg_sig = marg_sig)))
      Rhat <- Rhat * exp(-0.5*marg_sig[1]^2)
    } else{
      spr <- t(sapply(1:NROW(waa), function(x) get_SPR(fullF[x], res$MAA[x,], sel[x,], mat[x,],waa[x,], res$fracyr_SSB[x], at.age = TRUE)))
    }
    #spr <- t(sapply(1:NROW(waa), function(x) wham:::get_SPR(fullF[x], res$MAA[x,], sel[x,], mat[x,],waa[x,], res$fracyr_SSB[x], at.age = TRUE)))
    return(Rhat * spr)
    #return(res$NAA[,1] * spr)
}
