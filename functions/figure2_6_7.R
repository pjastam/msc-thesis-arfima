# The function figure2_6_7 generates 
# - three plots of times series with parameter values from the set (0.1, 0.5, 0.9)
# - one plot with the three theoretical autocorrelation functions for these time series
# - one plot with the three theoretical spectral densities for these time series
#
# nobs: number of observations generated for each of the three time series
# modelchoice: 'ar' or 'ma'
# pos_corr: indicator +1 or -1 to reflect positively or negative correlated observations, respectively

figure2_6_7 <- function(nobs, pos_corr) {
  # Set the seed for reproducibility
  set.seed(123)
  
  # Change the graphics parameters
  par(mfrow = c(4, 2))
  layout(matrix(c(1, 2, 3, 4, 5, 5, 6, 6), 4, 2, byrow = TRUE))
  
  # Helper function
#  if (modelchoice == 'ar') {
#    greeksymbol <- '\u03B1'
#  } else if (modelchoice == 'ma') {
#    greeksymbol <- '\u03B8'
#  }
  
  # Generate the realizations of three data generating processes
  model1 <- list(ar = NULL, ma = NULL, d = pos_corr * 0.25)
  model2 <- list(ar = NULL, ma = NULL, d = pos_corr * 1/3)
  model3 <- list(ar = NULL, ma = NULL, d = pos_corr * 0.45)
  
  sim1 <- arfima::arfima.sim(n = nobs, model = model1)
  sim2 <- arfima::arfima.sim(n = nobs, model = model2)
  sim3 <- arfima::arfima.sim(n = nobs, model = model3)
  
  # Plot the three graphs with realizations of the data generating processes
  plot(sim1, main = paste0('ARFIMA(0,d,0)', '(T=', nobs, ', d=', pos_corr * 0.25, ')'), xlab = "time", ylab = "", ylim = c(-4, 4))
  plot(sim2, main = paste0('ARFIMA(0,d,0)', '(T=', nobs, ', d=', MASS::fractions(pos_corr * 1/3), ')'), xlab = "time", ylab = "", ylim = c(-4, 4))
  plot(sim3, main = paste0('ARFIMA(0,d,0)', '(T=', nobs, ', d=', pos_corr * 0.45, ')'), xlab = "time", ylab = "", ylim = c(-4, 4))
  plot.new()
  
  # Plot the graph with the three theoretical autocorrelation function
  acf1 <- acf_arfima(n = 100, d = pos_corr * 0.25)
  acf2 <- acf_arfima(n = 100, d = pos_corr * 1/3)
  acf3 <- acf_arfima(n = 100, d = pos_corr * 0.45)
  
  plot(acf1, type = "l", main = 'Theoretical autocorrelation functions', xlab = "lags", ylab = "", ylim = c(0, 1))
  lines(acf2, type = "l", lty = "dashed")
  lines(acf3, type = "l", lty = "dotted")
  legend("topright", legend = c(paste0('d=', pos_corr * 0.45), paste0('d=', MASS::fractions(pos_corr * 1/3)), paste0('d=', pos_corr * 0.25)), lty = c("dotted", "dashed", "solid"), bty = "n")
  
  # Compute the theoretical pseudo-spectrum
  sd1 <- spectrum_arfima(d = pos_corr * 0.25)
  sd2 <- spectrum_arfima(d = pos_corr * 1/3)
  sd3 <- spectrum_arfima(d = pos_corr * 0.45)
  
  plot(sd1$freq, sd1$spec, type = "l", main = 'Spectral densities', xlab = expression("fractions of 2\u03C0"), ylab = "", ylim = c(0, 2.5))
  lines(sd2$freq, sd2$spec, type = "l", lty = "dashed")
  lines(sd3$freq, sd3$spec, type = "l", lty = "dotted")
  legend("topright", legend = c(paste0('d=', pos_corr * 0.45), paste0('d=', MASS::fractions(pos_corr * 1/3)), paste0('d=', pos_corr * 0.25)), lty = c("dotted", "dashed", "solid"), bty = "n")
  
  # Reset the graphics parameters
  par(mfrow = c(1, 1))
}

# Usage: 
# figure2_6_7(nobs = 100, modelchoice = 'ar', pos_corr = 0.8, modeltype = 'arfima', d = 0.4)
