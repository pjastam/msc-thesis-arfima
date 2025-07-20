figure2_9_10_11 <- function(N, T, mu = 0, sigmasq = 1) {

  # Set the seed for reproducibility
  set.seed(123)
  
  # Parameters
  alpha <- 0.5
  d <- 1/3

  # Calculate the standard deviation
  sigma <- sqrt(sigmasq)
  
  # Function to calculate correlation based on the process type 1
  corrfrac <- function(T, d) {
    corr <- rep(1, T)
    for (i in 2:T) {
      corr[i] <- corr[i-1] * (i-2+d) / (i-1-d)
    }
    return(corr)
  }
  
  # Function to calculate correlation based on the process type 2
  corrar1 <- function(T, alpha) {
    corr <- rep(1, T)
    for (i in 2:T) {
      corr[i] <- corr[i-1] * alpha
    }
    return(corr)
  }
  
  m_matrix <- function(corr) {
    # Create Toeplitz matrix from the correlation vector
    corrmat <- toeplitz(corr)
    
    # Cholesky decomposition
    m <- chol(corrmat)
    
    return(m)
  }
  
  # Calculate the theoretical correlations based on the process type
  corr_frac <- corrfrac(T, d)
  corr_ar1 <- corrar1(T, alpha)
  
  # Calculate the Cholesky decomposition of the correlation matrix
  m_frac <- m_matrix(corr_frac)
  m_ar1 <- m_matrix(corr_ar1)

  # Storage for means, variances and autocovariances
  mean_est_arfima <- numeric(N)
  mean_est_arima <- numeric(N)
  
  var_pop_arfima <- numeric(N)
  var_est_arfima <- numeric(N)
  var_pop_arima <- numeric(N)
  var_est_arima <- numeric(N)
  
  acov_pop_arfima <- numeric(N)
  acov_est_arfima <- numeric(N)
  acov_pop_arima <- numeric(N)
  acov_est_arima <- numeric(N)
  
  # Population mean
  mean_pop_arfima <- 0  # Assumed known population mean for ARFIMA
  mean_pop_arima <- 0   # Assumed known population mean for ARIMA
  
  # Simulations
  for (i in 1:N) {
    # Simulate ARFIMA and ARIMA processes
    arfima_series <- mu + sigma * m_frac %*% rnorm(T)
    arima_series <- mu + sigma * m_ar1 %*% rnorm(T)
    
    # Compute sample mean
    mean_est_arfima[i] <- mean(arfima_series)
    mean_est_arima[i] <- mean(arima_series)
    
    # Compute variances
    var_pop_arfima[i] <- sum((arfima_series - mean_pop_arfima)^2) / (T - 1)
    var_est_arfima[i] <- sum((arfima_series - mean_est_arfima[i])^2) / (T - 1)
    var_pop_arima[i] <- sum((arima_series - mean_pop_arima)^2) / (T - 1)
    var_est_arima[i] <- sum((arima_series - mean_est_arima[i])^2) / (T - 1)
    
    # Compute autocovariances at lag 1
    acov_pop_arfima[i] <- sum((arfima_series[1:(T-1)] - mean_pop_arfima) * (arfima_series[2:T] - mean_pop_arfima)) / (T - 1)
    acov_est_arfima[i] <- sum((arfima_series[1:(T-1)] - mean_est_arfima[i]) * (arfima_series[2:T] - mean_est_arfima[i])) / (T - 1)
    acov_pop_arima[i] <- sum((arima_series[1:(T-1)] - mean_pop_arima) * (arima_series[2:T] - mean_pop_arima)) / (T - 1)
    acov_est_arima[i] <- sum((arima_series[1:(T-1)] - mean_est_arima[i]) * (arima_series[2:T] - mean_est_arima[i])) / (T - 1)    
  }
  
  return(list(mean_est_arfima = mean_est_arfima, mean_est_arima = mean_est_arima,
              var_pop_arfima = var_pop_arfima, var_est_arfima = var_est_arfima,
              var_pop_arima = var_pop_arima, var_est_arima = var_est_arima,
              acov_pop_arfima = acov_pop_arfima, acov_est_arfima = acov_est_arfima,
              acov_pop_arima = acov_pop_arima, acov_est_arima = acov_est_arima))  
}

# Usage:
#figure2_9_10_11(N = 1000, T=300, mu = 0, sigmasq = 1)
