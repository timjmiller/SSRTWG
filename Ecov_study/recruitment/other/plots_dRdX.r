library(viridis)

matplot(replicate(n=10,arima.sim(list(ar=0.95), 50, sd=0.1)),ylim=c(-5,5),type='l',col='black',lty=1)

matplot(replicate(n=10,arima.sim(list(ar=0.95), 50, sd=0.5)),ylim=c(-5,5),type='l',col='red',lty=1,add=TRUE)


matplot(replicate(n=10,arima.sim(list(ar=0.95), 50, sd=0.5)),ylim=c(-5,5),type='l',col='black',lty=1)
matplot(replicate(n=10,arima.sim(list(ar=0.5), 50, sd=0.5)),ylim=c(-5,5),type='l',col='red',lty=1,add=TRUE)

limiting     <- function(SSB,X,a,b,c) SSB/(b+a*SSB*exp(c*X))
masking      <- function(SSB,X,a,b,c) SSB/(b*exp(c*X)+a*SSB)
controlling  <- function(SSB,X,a,b,c) (SSB/(b+a*SSB))*exp(c*X)

dlimiting_dx     <- function(SSB,X,a,b,c) -((SSB^2)*a*c*exp(c*X))/(b+a*SSB*exp(c*X))^2
dmasking_dx      <- function(SSB,X,a,b,c) -(SSB*b*c*exp(c*X))/(b*exp(c*X)+a*SSB)^2
dcontrolling_dx  <- function(SSB,X,a,b,c) (SSB*c*exp(c*X))/(a*SSB+b)

limiting_ratio   <- function(SSB,X,a,b,c) dlimiting_dx(SSB,X,a,b,c)/limiting(SSB,X,a,b,c)
masking_ratio    <- function(SSB,X,a,b,c) dmasking_dx(SSB,X,a,b,c)/masking(SSB,X,a,b,c)
controlling_ratio <- function(SSB,X,a,b,c) dcontrolling_dx(SSB,X,a,b,c)/controlling(SSB,X,a,b,c)

ssbs <- seq(0,1,0.01)

xs   <- c(0,1,2)
ncol  <- 3
ylimr <- c(0,1.2)
ylimdr <- c(-0.12,0.12)


pdf('dRdX_v2.pdf',height=8,width=8)
par(mfrow=c(3,3),mar=c(2,2,1,2),oma=c(2,2,2,2))
matplot(ssbs, sapply(xs, function(x) limiting(ssbs,x,1,0.1,0.1)),type='l',lty=1,col=turbo(ncol),ylab='Recruitment', xlab='SSB',ylim=ylimr)
  legend('bottomright',bty='n',lty=1,col=turbo(ncol),legend=xs)
  mtext(side=2,'Recruitment',line=2.5)
  mtext('Limiting')
matplot(ssbs, sapply(xs, function(x) masking(ssbs,x,1,0.1,0.1)),type='l',lty=1,col=turbo(ncol),ylab='Recruitment', xlab='SSB',ylim=ylimr)
  mtext('Masking')
matplot(ssbs, sapply(xs, function(x) controlling(ssbs,x,1,0.1,0.1)),type='l',lty=1,col=turbo(ncol),ylab='Recruitment', xlab='SSB',ylim=ylimr)
  mtext('Controlling')

matplot(ssbs, sapply(xs, function(x) dlimiting_dx(ssbs,x,1,0.1,0.1)),type='l',lty=1,col=turbo(ncol),ylab='Recruitment', xlab='SSB',ylim=ylimdr)
  mtext(side=2,'dR/dX',line=2.5)
  abline(h=0,lty=2)
matplot(ssbs, sapply(xs, function(x) dmasking_dx(ssbs,x,1,0.1,0.1)),type='l',lty=1,col=turbo(ncol),ylab='Recruitment', xlab='SSB',ylim=ylimdr)
abline(h=0,lty=2)
matplot(ssbs, sapply(xs, function(x) dcontrolling_dx(ssbs,x,1,0.1,0.1)),type='l',lty=1,col=turbo(ncol),ylab='Recruitment', xlab='SSB',ylim=ylimdr)
abline(h=0,lty=2)
mtext(side=1,outer=TRUE,'SSB',line=0.5)  

matplot(ssbs, sapply(xs, function(x) limiting_ratio(ssbs,x,1,0.1,0.1)),type='l',lty=1,col=turbo(ncol),ylab='Recruitment', xlab='SSB',ylim=ylimdr)
mtext(side=2,'(1/R)dR/dX',line=2.5)
abline(h=0,lty=2)
matplot(ssbs, sapply(xs, function(x) masking_ratio(ssbs,x,1,0.1,0.1)),type='l',lty=1,col=turbo(ncol),ylab='Recruitment', xlab='SSB',ylim=ylimdr)
abline(h=0,lty=2)
matplot(ssbs, sapply(xs, function(x) controlling_ratio(ssbs,x,1,0.1,0.1)),type='l',lty=1,col=turbo(ncol),ylab='Recruitment', xlab='SSB',ylim=ylimdr)
abline(h=0,lty=2)
mtext(side=1,outer=TRUE,'SSB',line=0.5)  

dev.off()





cs   <- c(0.1,0.2,0.5)
ncol  <- 3
ylimr <- c(0,1.5)
ylimdr <- c(-0.8,0.8)
ylimdrr <- c(-0.6,0.6)


pdf('dRdX_v2.pdf',height=8,width=8)
par(mfrow=c(3,3),mar=c(2,2,1,2),oma=c(2,2,2,2))
matplot(ssbs, sapply(cs, function(x) controlling(ssbs,1,1,0.1,x)),type='l',lty=1,col=turbo(ncol),ylab='Recruitment', xlab='SSB',ylim=ylimr)
mtext(side=2,'Recruitment',line=2.5)
mtext('Controlling')
legend('bottomright',bty='n',lty=1,col=turbo(ncol),legend=cs)
matplot(ssbs, sapply(cs, function(x) limiting(ssbs,1,1,0.1,x)),type='l',lty=1,col=turbo(ncol),ylab='Recruitment', xlab='SSB',ylim=ylimr)
mtext('Limiting')
matplot(ssbs, sapply(cs, function(x) masking(ssbs,1,1,0.1,x)),type='l',lty=1,col=turbo(ncol),ylab='Recruitment', xlab='SSB',ylim=ylimr)
mtext('Masking')

matplot(ssbs, sapply(cs, function(x) dcontrolling_dx(ssbs,1,1,0.1,x)),type='l',lty=1,col=turbo(ncol),ylab='Recruitment', xlab='SSB',ylim=ylimdr)
abline(h=0,lty=2)
mtext(side=2,'dR/dX',line=2.5)
mtext(side=1,outer=TRUE,'SSB',line=0.5)  
matplot(ssbs, sapply(cs, function(x) dlimiting_dx(ssbs,1,1,0.1,x)),type='l',lty=1,col=turbo(ncol),ylab='Recruitment', xlab='SSB',ylim=ylimdr)
abline(h=0,lty=2)
matplot(ssbs, sapply(cs, function(x) dmasking_dx(ssbs,1,1,0.1,x)),type='l',lty=1,col=turbo(ncol),ylab='Recruitment', xlab='SSB',ylim=ylimdr)
abline(h=0,lty=2)

matplot(ssbs, sapply(cs, function(x) controlling_ratio(ssbs,1,1,0.1,x)),type='l',lty=1,col=turbo(ncol),ylab='Recruitment', xlab='SSB',ylim=ylimdrr)
mtext(side=2,'(1/R)dR/dX',line=2.5)
abline(h=0,lty=2)
matplot(ssbs, sapply(cs, function(x) limiting_ratio(ssbs,1,1,0.1,x)),type='l',lty=1,col=turbo(ncol),ylab='Recruitment', xlab='SSB',ylim=ylimdrr)
abline(h=0,lty=2)
matplot(ssbs, sapply(cs, function(x) masking_ratio(ssbs,1,1,0.1,x)),type='l',lty=1,col=turbo(ncol),ylab='Recruitment', xlab='SSB',ylim=ylimdrr)
abline(h=0,lty=2)
mtext(side=1,outer=TRUE,'SSB',line=0.5)  


