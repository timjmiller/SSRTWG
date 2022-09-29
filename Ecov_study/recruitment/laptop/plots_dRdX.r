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


ssbs <- seq(0,1,0.01)


xs <- c(-1,-0.1,0,0.1,1)

pdf('dRdX.pdf',height=5,width=8)
par(mfrow=c(2,3),mar=c(2,3,2,2),oma=c(2,2,2,2))
matplot(ssbs, sapply(xs, function(x) limiting(ssbs,x,1,0.1,0.1)),type='l',lty=1,col=turbo(5),ylab='Recruitment', xlab='SSB')
  legend('bottomright',bty='n',lty=1,col=turbo(5),legend=xs)
  mtext(side=2,'Recruitment',line=2.5)
  mtext('Limiting')
matplot(ssbs, sapply(xs, function(x) masking(ssbs,x,1,0.1,0.1)),type='l',lty=1,col=turbo(5),ylab='Recruitment', xlab='SSB')
  mtext('Masking')
matplot(ssbs, sapply(xs, function(x) controlling(ssbs,x,1,0.1,0.1)),type='l',lty=1,col=turbo(5),ylab='Recruitment', xlab='SSB')
  mtext('Controlling')

matplot(ssbs, sapply(xs, function(x) dlimiting_dx(ssbs,x,1,0.1,0.1)),type='l',lty=1,col=turbo(5),ylab='Recruitment', xlab='SSB')
  mtext(side=2,'dR/dX',line=2.5)
matplot(ssbs, sapply(xs, function(x) dmasking_dx(ssbs,x,1,0.1,0.1)),type='l',lty=1,col=turbo(5),ylab='Recruitment', xlab='SSB')
matplot(ssbs, sapply(xs, function(x) dcontrolling_dx(ssbs,x,1,0.1,0.1)),type='l',lty=1,col=turbo(5),ylab='Recruitment', xlab='SSB')
mtext(side=1,outer=TRUE,'SSB',line=0.5)  
dev.off()






