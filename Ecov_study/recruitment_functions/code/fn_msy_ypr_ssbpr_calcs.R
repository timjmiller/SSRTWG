# liz brooks
# dec 2025
# aggregate functions used for SSRTWG : SRR-Ecov ms
# selectivity; ypr; MSY; etc.

# functions =====
logistic.fn <-function(a50, k, ages) {
  tmp <- 1.0/(1.0 + exp(-(ages - a50)/k)) 
  return(tmp)
}

# ------------------------------
ypr<-function(nages, wgt.age, mm, F.mult, sel.age ) {
  
  
  if (length(mm)==nages) { M=mm  
  }    else  {  
    M=rep(mm , nages) 
    }
  
  yield=0.0
  cum.survive=1.0
  z=0.0
  
  
  for (i in 1:(nages-1)  ) {
    z=M[i] + F.mult*sel.age[i]
    yield=yield + wgt.age[i]*F.mult*sel.age[i]*(1-exp(-z) )*cum.survive/z
    cum.survive=cum.survive*exp(-z)  
  }
  
  cum.survive=cum.survive*exp(-M[nages]-F.mult*sel.age[nages])
  yield=yield + wgt.age[nages]*F.mult*sel.age[nages]*cum.survive/z
  
  return(yield)
  
}

# ------------------------------------
ssb.eq<-function(alpha.BH, R0.BH, spr, spr0 ) {
  
  sprF=spr/spr0
  
  ssb=spr0*R0.BH*(sprF*alpha.BH - 1.0)/(alpha.BH-1.0)
  
  
  return(ssb)
  
}

# ------------------------------------
s.per.recr<-function(nages,fec.age,mat.age,mm, F.mult, sel.age, spawn.time ) {
  
  spr=0.0
  cum.survive=1.0
  z=0.0
  if (length(mm)==nages) {     M.age=mm     
  }     else  {   
    M.age=rep(mm , nages) 
    }
  
  for (i in 1:(nages-1)  ) {
    z=M.age[i] + F.mult*sel.age[i]
    z.ts=(M.age[i]+F.mult*sel.age[i])*spawn.time
    spr=spr+cum.survive*fec.age[i]*mat.age[i]*exp(-z.ts)
    cum.survive=cum.survive*exp(-z )
    
  }
  
  z= M.age[nages] + F.mult*sel.age[nages]
  z.ts=(M.age[nages]+F.mult*sel.age[nages])*spawn.time
  spr=spr + fec.age[nages]*mat.age[nages]*cum.survive*exp(-z.ts)/( 1- exp(-z ) )
  
  return(spr)
  
}

# ------------------------------------
bev.holt.alpha<-function(S,R0,recr.par,spr0, is.steepness=T){
  
  if (is.steepness==T)  alpha.BH <- 4*recr.par/(1-recr.par)
  if (is.steepness==F)  alpha.BH <- recr.par
  
  
  y=rep(0,length(S))
  y=R0*S*alpha.BH/(R0*spr0 +(alpha.BH-1.0)*S)
  return(y)
  
} 

# ------------------------------------
# note: should typically run this more than once with different starting values to make sure it finds the min
get.MSY <- function(nages=NULL, steep=NULL, R0=NULL, fec.age=NULL, mat.age=NULL, M.age=NULL, 
                    fleet.waa=NULL, sel.age=NULL, spawn.time=NULL,
                    F.start=NULL) 
  {
  

   

    spr0 <- s.per.recr(nages=nages,fec.age=fec.age,mat.age=mat.age, mm=M.age, F.mult=0, sel.age=sel.age, spawn.time=spawn.time ) 
    
    ahat <- 4*steep/(1-steep)
    MSY.soln <- rep(NA, 9)
    names(MSY.soln) <- c("Fmsy", "MSY", "SPRmsy", "SSBmsy", "Rmsy", "YPRmsy",
                            "Rel.SSBmsy", "Rel.Rmsy", "Conv.Code")

 
      
      get.yield.f.min <- function(F.start) {
        temp.spr = s.per.recr(nages=nages, fec.age=fec.age, mat.age=mat.age, mm= M.age, F.mult=F.start, sel.age=sel.age, spawn.time=spawn.time)   
        temp.ypr = ypr(nages=nages, wgt.age=fleet.waa, mm=M.age,  F.mult=F.start, sel.age=sel.age )
        temp.SSB = ssb.eq( alpha.BH=ahat, R0.BH=R0, spr=temp.spr, spr0=spr0  )
        yield = temp.ypr*temp.SSB/temp.spr           #harvest in weight
        yy=-1*yield
        return(yy)  
      }   # end get.yield.f.min function
      
      F.nlmin <- nlminb( start=F.start , objective=get.yield.f.min, lower=0.01, upper=3.0,
                         control=list(eval.max=500, iter.max=200  )   )
      
      
      
      MSY.soln[1 ] <- F.nlmin$par  #Fmsy
      MSY.soln[2 ] <- -1*F.nlmin$objective #MSY
      
      # calculate other MSY quantities
      spr.nlmin <- s.per.recr(nages=nages, fec.age=fec.age, mat.age=mat.age, mm= M.age, F.mult=F.nlmin$par, 
                              sel.age=sel.age, spawn.time=spawn.time)   
      MSY.soln[3 ] <- spr.nlmin/spr0 #SPRmsy
      
      ssb.nlmin <- ssb.eq( alpha.BH=ahat, R0.BH=R0, spr=spr.nlmin, spr0=spr0 )
      MSY.soln[4 ] <-  ssb.nlmin #SSBmsy
      
      
      Rmsy.nlmin <- bev.holt.alpha(S=ssb.nlmin,R0=R0,recr.par=ahat,spr0=spr0, is.steepness=FALSE)
      MSY.soln[5  ] <- Rmsy.nlmin #Rmsy
      
      ypr.nlmin <- ypr(nages=nages, wgt.age=fleet.waa, mm=M.age,  F.mult=F.nlmin$par, sel.age=sel.age )
      MSY.soln[6  ] <- ypr.nlmin #YPRmsy
      
      rel.ssb.nlmin <- ssb.nlmin/(R0*spr0)
      MSY.soln[7 ] <-  rel.ssb.nlmin #SSBmsy/SSB0
      
      MSY.soln[8 ] <- Rmsy.nlmin/R0 #Rmsy/R0
      
      MSY.soln[9 ] <-  F.nlmin$convergence  #code=0 indicates convergence
 
      return(MSY.soln) 
}  # end function
# end functions =======
