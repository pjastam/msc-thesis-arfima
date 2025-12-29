library(urca)

loopadf <- function(endogenen, kk, trend = 1, constante = 1) {
  resultaat <- matrix(0, nrow = kk, ncol = 3)
  mend <- length(endogenen)
  
  if (trend == 0) {
    if (constante == 0) {
      # neither an intercept nor a trend is included in the test regression
      adf_type <- "none"
    } else {
      # an inteercept is included in the test regression
      adf_type <- "drift"
    }
  } else if (trend == 1) {
    # an intercept and a trend is included in the test regression
    adf_type <- "trend"
  }
  
  for (k in 1:kk) {
    res <- ur.df(endogenen, type = adf_type, lags = k, selectlags = c("Fixed"))
    resultaat[k, 1] <- res@teststat[1]
    resultaat[k, 2] <- -(log(sum(res@res^2) / mend) + log(mend) * k / mend)
    resultaat[k, 3] <- -(log(sum(res@res^2) / mend) + 2 * k / mend)
  }
  
  return(resultaat)
}

figure4_2 <- function(endogenen, kk, trend = 1, constante = 1) {
  
  dendogenen <- diff(endogenen)
  
  # Calculate the ADF test statistics
  # Type 3: with drift and trend
  library("aTSA")
  adf_endogenen <- loopadf(endogenen, kk, trend, constante)
  adf_dendogenen <- loopadf(dendogenen, kk, trend, constante)
  
  # Plot the test statistics for the endogenen
  par(mfrow = c(2, 2))
  
  plot(
    c(1:kk),
    adf_endogenen[, 1],
    main = "ADF t statistic",
    xlab = "lags",
    ylab = "",
    type = "l"
#    ylim = c(min(adf_endogenen[, 1]), 0)
  )
  
  plot(
    c(1:kk),
    adf_endogenen[, 2],
    main = "SIC criterion",
    xlab = "order AR polynomial",
    ylab = "",
    type = "l"
#    ylim = c(min(adf_endogenen[, 2]), 0)
  )

  plot(
    c(1:kk),
    adf_endogenen[, 3],
    main = "AIC criterion",
    xlab = "order AR polynomial",
    ylab = "",
    type = "l"
#    ylim = c(min(adf_endogenen[, 3]), 0)
  )
  
  par(mfrow = c(1, 1))
  par(mfrow = c(2, 2))
  
  # Plot the test statistics for the first differences
  plot(
    c(1:kk),
    adf_dendogenen[, 1],
    main = "ADF t statistic",
    xlab = "lags",
    ylab = "",
    type = "l"
#    ylim = c(min(adf_dendogenen[, 1]), 0)
  )
  
  plot(
    c(1:kk),
    adf_dendogenen[, 2],
    main = "SIC criterion",
    xlab = "order AR polynomial",
    ylab = "",
    type = "l"
#    ylim = c(min(adf_dendogenen[, 2]), 0)
  )

  plot(
    c(1:kk),
    adf_dendogenen[, 3],
    main = "AIC criterion",
    xlab = "order AR polynomial",
    ylab = "",
    type = "l"
#    ylim = c(min(adf_dendogenen[, 3]), 0)
  )
  
  par(mfrow = c(1, 1)) 
}

# Example usage:
#figure4_2(LNGNP82, 51, trend = 1, constante = 1)
