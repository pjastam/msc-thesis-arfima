figure4_1 <- function(title, y, from, to, by) {
  
  t <-  seq(as.Date(from), as.Date(to), by = by)
  dy <- diff(y)
  dt <- t[-1]
  
  library("TSA")
  
  par(mfrow = c(2, 2))
  plot(t,y,type='l',main=title,xlab="time",ylab="");
  legend("topright",legend=c(bquote(mean: .(round(mean(y),4))),bquote(variance: .(round(var(y),4)))), bty = "n")
  acf(y,lag.max=20,type="correlation",main="Autocorrelation Coefficients",xlab="lags 1-20",ylab="")
  periodogram(y,main="Periodogram",xlab=expression("Fractions of" ~ 2*pi),ylab="");  abline(h=0)
  acf(y,lag.max=20,type="partial",main="Partial Autocorrelation Coefficients",xlab="lags 1-20",ylab="")
  
  plot(dt,dy,type='l',main='First Differences',xlab="time",ylab="");
  legend("topright",legend=c(bquote(mean: .(round(mean(dy),4))),bquote(variance: .(round(var(dy),4)))), bty = "n")
  acf(dy,lag.max=20,type="correlation",main="Autocorrelation Coefficients",xlab="lags 1-20",ylab="")
  periodogram(dy,main="Periodogram",xlab=expression("Fractions of" ~ 2*pi),ylab="");  abline(h=0)
  acf(dy,lag.max=20,type="partial",main="Partial Autocorrelation Coefficients",xlab="lags 1-20",ylab="")
  par(mfrow = c(1, 1)) 
}