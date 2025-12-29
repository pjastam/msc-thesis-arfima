# The function figure2 generates 
# - three plots of times series with parameter values from the set (0.1, 0.5, 0.9)
# - one plot with the three theoretical autocorrelation functions for these time series
# - one plot with the three theoretical spectral densities for these time series
#
# nobs: number of observations generated for each of the three time series
# modelchoice: 'ar' or 'ma'
# pos_corr: indicator +1 or -1 to reflect positively or negative correlated observations, respectively

figure2_1_2_3_4 <- function(nobs,modelchoice,pos_corr) {
  
  # Set the seed for reproducibility
  set.seed(123)
  
  # Change the graphics parameters
  par(mfrow = c(4, 2))
  layout(matrix(c(1,2,3,4,5,5,6,6), 4, 2, byrow = TRUE))
  
  # Helper function
  if (modelchoice == 'ar') {
    greeksymbol <- '\u03B1'
  }
  else if (modelchoice == 'ma') {
    greeksymbol <- '\u03B8'
  }
  
  # Generate the realizations of three data generating processes
  # Source: https://stackoverflow.com/questions/21893165/assigning-names-in-a-list-using-variables
  model1 <- list(ar=NULL, ma=NULL); model1[[modelchoice]] <- pos_corr*0.1
  model2 <- list(ar=NULL, ma=NULL); model2[[modelchoice]] <- pos_corr*0.5
  model3 <- list(ar=NULL, ma=NULL); model3[[modelchoice]] <- pos_corr*0.9
  
  sim1 <- arima.sim(model1, n=nobs)
  sim2 <- arima.sim(model2, n=nobs)
  sim3 <- arima.sim(model3, n=nobs)
  
  # Plot the three graphs with realizations of the data generating processes 
  plot(sim1,main=paste0(toupper(modelchoice),'(1) (T=',nobs,', ',greeksymbol,'=',pos_corr*0.1,')'),xlab="time",ylab="",ylim=c(-4,4))
  plot(sim2,main=paste0(toupper(modelchoice),'(1) (T=',nobs,', ',greeksymbol,'=',pos_corr*0.5,')'),xlab="time",ylab="",ylim=c(-4,4))
  plot(sim3,main=paste0(toupper(modelchoice),'(1) (T=',nobs,', ',greeksymbol,'=',pos_corr*0.9,')'),xlab="time",ylab="",ylim=c(-4,4))
  plot.new()
  
  # Plot the graph with the three theoretical autocorrelation function
  acf1 <- do.call(stats::ARMAacf, within(Filter(Negate(is.null), model1), lag <- 100))
  acf2 <- do.call(stats::ARMAacf, within(Filter(Negate(is.null), model2), lag <- 100))
  acf3 <- do.call(stats::ARMAacf, within(Filter(Negate(is.null), model3), lag <- 100))
  
  plot(acf1,type="l",main='Theoretical autocorrelation functions',xlab="lags",ylab="",ylim=c(-1,1))
  lines(acf2,type="l",lty="dashed")
  lines(acf3,type="l",lty="dotted")
  legend("topright",legend=c(paste0(greeksymbol,'=',pos_corr*0.9),
                             paste0(greeksymbol,'=',pos_corr*0.5),
                             paste0(greeksymbol,'=',pos_corr*0.1)),
         lty=c("dotted","dashed","solid"), bty = "n")
  
  # Plot the graph with the three spectral densities
  sd1 <- do.call(astsa::arma.spec, within(Filter(Negate(is.null), model1), plot <- FALSE))
  sd2 <- do.call(astsa::arma.spec, within(Filter(Negate(is.null), model2), plot <- FALSE))
  sd3 <- do.call(astsa::arma.spec, within(Filter(Negate(is.null), model3), plot <- FALSE))
  
  plot(sd1$freq,sd1$spec,type="l",main='Spectral densities',xlab=expression("fractions of 2\u03C0"),ylab="",ylim=c(0,10))
  lines(sd2$freq,sd2$spec,type="l",lty="dashed");
  lines(sd3$freq,sd3$spec,type="l",lty="dotted");
  legend("topright",legend=c(paste0(greeksymbol,'=',pos_corr*0.9),
                             paste0(greeksymbol,'=',pos_corr*0.5),
                             paste0(greeksymbol,'=',pos_corr*0.1)),
         lty=c("dotted","dashed","solid"), bty = "n")
  
  # Reset the graphics parameters
  par(mfrow = c(1, 1))
}
