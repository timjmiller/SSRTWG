
library(fields)

basic <- make_basic_info()

pdf('basic_info.pdf',height=10,width=5)
par(mfrow=c(5,1),mar=c(2,4,2,2),oma=c(2,3,2,2))
  image.plot(y=basic$ages,x=basic$years,z=basic$waa[1,,],xlab='',ylab='Ages')
  mtext(adj=0.1,'Weight-at-age')
  
  image.plot(y=basic$ages,x=basic$years,z=basic$maturity,xlab='',ylab='Ages')
  mtext(adj=0.1,'Maturity-at-age')
  
  plot(basic$years,basic$F,ylab='F',xlab='')
  
  plot(basic$ages,NAA_re$N1_pars,ylab='N1_pars')
  
  plot(basic$ages,M$initial_means,ylab='M initial means')
  
dev.off()


