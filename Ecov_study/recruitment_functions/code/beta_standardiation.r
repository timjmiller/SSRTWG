
l   <- function(SSB,X,a,b,c) (a*SSB/(1 + exp(c*X)*b*SSB)) #limiting
m   <- function(SSB,X,a,b,c) (a*SSB/(exp(c*X) + b*SSB)) #masking
c   <- function(SSB,X,a,b,c) (exp(c*X)*a*SSB/(1 + b*SSB)) #controlling

#default SR parmaeters used throughout studies
a      <- 5.955694e-01 
b      <- 2.404283e-05
SSBMSY <- 132096.3 

dRdX <- function(beta,form,sdX){
  if(form==2) return((l(SSBMSY,sdX,a,b,beta)  - l(SSBMSY,-sdX,a,b,beta))/sdX)
  if(form==4) return((m(SSBMSY,sdX,a,b,beta)  - m(SSBMSY,-sdX,a,b,beta))/sdX)
  if(form==1) return((c(SSBMSY,sdX,a,b,beta) - c(SSBMSY,-sdX,a,b,beta))/sdX)
}

#set default sensitivities
c0_low  <- median(abs(base::c(dRdX(beta=0.1,form=1,sdX=1),
             dRdX(beta=0.1,form=2,sdX=1),
             dRdX(beta=0.1,form=4,sdX=1))))
c0_high <- median(abs(base::c(dRdX(beta=0.6,form=1,sdX=1),
                             dRdX(beta=0.6,form=2,sdX=1),
                             dRdX(beta=0.6,form=4,sdX=1))))

#function to minimize for necessary beta
beta_cost <- function(beta){
  return((dRdX(beta,form=form,sdX=sdX) - c0)^2)
}

##--TEST--#####################
#set parameters
#c0   <- -10000 #desired sensitivity
#form=1; sdX=1
#optimize(beta_cost,lower=-10,upper=10)

##--APPLY STANDARDIZATION--##############
df.oms$beta_std <- 
  sapply(1:nrow(df.oms), function(i){
  if(df$Ecov_effect[i]==0.1) c0 <<- c0_low
  if(df$Ecov_effect[i]==1.0) c0 <<- c0_high
  sdX  <<- sqrt(df.oms$Ecov_re_sig[i]^2/(1-df.oms$Ecov_re_cor[i]^2))
  form <<- df.oms$Ecov_how[i] 
  ifelse(form==0,df.oms$Ecov_effect[i],optimize(beta_cost,lower=-10,upper=10)$minimum)
})


dRdX_plot <- function(beta,form,sdX,SSBMSY){
  if(form==2) return((l(SSBMSY,sdX,a,b,beta)  - l(SSBMSY,-sdX,a,b,beta))/sdX)
  if(form==4) return((m(SSBMSY,sdX,a,b,beta)  - m(SSBMSY,-sdX,a,b,beta))/sdX)
  if(form==1) return((c(SSBMSY,sdX,a,b,beta) - c(SSBMSY,-sdX,a,b,beta))/sdX)
}


ssbs <- seq(0,SSBMSY*2,length.out=1000)
par(mfrow=base::c(3,2),mar=base::c(3,5,2,2),oma=base::c(4,4,2,2))
plot(ssbs,dRdX_plot(beta=0.07603968,form=1,sdX=1,SSBMSY=ssbs),xlab='SSB',ylab=expression(Delta*'R/'*Delta*'sd(X)'),type='l')
mtext(line=0,0.07603968,cex=0.7)
abline(v=SSBMSY,lty=2)
abline(h=c0_low,lty=2)
plot(ssbs,dRdX_plot(beta=0.45313471,form=1,sdX=1,SSBMSY=ssbs),xlab='SSB',ylab=expression(Delta*'R/'*Delta*'sd(X)'),type='l')
mtext(line=0,0.45313471,cex=0.7)
abline(v=SSBMSY,lty=2)
abline(h=c0_high,lty=2)


plot(ssbs,dRdX_plot(beta=-0.09998358,form=2,sdX=1,SSBMSY=ssbs),xlab='SSB',ylab=expression(Delta*'R/'*Delta*'sd(X)'),type='l')
abline(v=SSBMSY,lty=2)
abline(h=c0_low,lty=2)
mtext(line=0,-0.09998358,cex=0.7)
plot(ssbs,dRdX_plot(beta=-0.59645767,form=2,sdX=1,SSBMSY=ssbs),xlab='SSB',ylab=expression(Delta*'R/'*Delta*'sd(X)'),type='l')
abline(v=SSBMSY,lty=2)
abline(h=c0_high,lty=2)
mtext(line=0,-0.59645767,cex=0.7)


plot(ssbs,dRdX_plot(beta=-0.3175733,form=4,sdX=1,SSBMSY=ssbs),xlab='SSB',ylab=expression(Delta*'R/'*Delta*'sd(X)'),type='l')
abline(v=SSBMSY,lty=2)
abline(h=c0_low,lty=2)
mtext(line=0,-0.3175733,cex=0.7)
plot(ssbs,dRdX_plot(beta=-2.0999895,form=4,sdX=1,SSBMSY=ssbs),xlab='SSB',ylab=expression(Delta*'R/'*Delta*'sd(X)'),type='l')
abline(v=SSBMSY,lty=2)
abline(h=c0_high,lty=2)
mtext(line=0,-2.0999895,cex=0.7)




