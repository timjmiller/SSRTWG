
l  <- function(SSB,X,a,b,c) (a*SSB/(1 + exp(c*X)*b*SSB)) #limiting
m  <- function(SSB,X,a,b,c) (a*SSB/(exp(c*X) + b*SSB)) #masking
c  <- function(SSB,X,a,b,c) (exp(c*X)*a*SSB/(1 + b*SSB)) #controlling

#default SR parmaeters used throughout studies
a      <- 5.955694e-01 
b      <- 2.404283e-05
SSBMSY <- 132096.3 

dRdX <- function(beta,form,sdX){
  if(form==2) return((l(SSBMSY,sdX,a,b,beta) - l(SSBMSY,-sdX,a,b,beta))/sdX)
  if(form==4) return((m(SSBMSY,sdX,a,b,beta) - m(SSBMSY,-sdX,a,b,beta))/sdX)
  if(form==1) return((c(SSBMSY,sdX,a,b,beta) - c(SSBMSY,-sdX,a,b,beta))/sdX)
}

##--TEST--#####################
#set parameters
c0   <- -10000 #desired sensitivity
sdX  <- 1      #standard deviation of the ecov process
form <- 1

beta_cost <- function(beta){
  return((dRdX(beta,form=form,sdX=sdX) - c0)^2)
}

form=1; sdX=1
optimize(beta_cost,lower=-10,upper=10)

##--APPLY STANDARDIZATION--##############
df.oms$beta_std <- 
  sapply(1:nrow(df.oms), function(i){
  c0   <- -10000
  sdX  <<- sqrt(df.oms$Ecov_re_sig[i]^2/(1-df.oms$Ecov_re_cor[i]^2))
  form <<- df.oms$Ecov_how[i] 
  ifelse(form==0,df.oms$Ecov_effect[i],optimize(beta_cost,lower=-10,upper=10)$minimum)
})

