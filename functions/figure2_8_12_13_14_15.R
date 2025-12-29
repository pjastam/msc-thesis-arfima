figure2_8_12_13_14_15 <- function(process, nrsimul, T, param, nk = T-1, mu = 0, sigmasq = 1, ymin = -1, ymax = 1, pop_mean = TRUE) {
  # This code has been translated from the original Matlab code by GitHub Co-Pilot.
  #
  # The original code can be found in the file named SIMUL.M in de directory scriptie_diskette1.
  # Note that this Matlab code is also available as Fortran code in the file GRANGERJ.F 
  # located in the directory scriptie_diskette2. The calculations differ somewhat, for example, 
  # 2 instead of 1.96 is used in the code for calculating the variables plusreeks and minreeks.

  #SIMUL(process,nrsimul,T,param,nk,mu,sigmasq)
  #	SIMUL geeft de theoretische en empirische autocorrelatiefunctie
  #	na nrsimul simulaties; verder worden +2*sigma en -2*sigma
  #	betrouwbaarheidsintervallen weggeschreven.
  #
  #	process=1: fractioneel geintegreerd proces
  #	process=2: ar(1) proces
  #	process=3: ma(1) proces
  
  #	P.J.A. Stam, 13-04-92.
  
  # Set the seed for reproducibility
  set.seed(123)
  
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
  
  # Function to calculate correlation based on the process type 3
  corrma1 <- function(T, gamma) {
    corr <- rep(0, T)
    corr[1] <- 1
    corr[2] <- gamma / (1 + gamma^2)
    return(corr)
  }
  
  # Calculate the theoretical correlations based on the process type
  if (process == 1) {
    corr <- corrfrac(T, param)
    label_param <- "d"
  } else if (process == 2) {
    corr <- corrar1(T, param)
    label_param <- "alpha"
  } else if (process == 3) {
    corr <- corrma1(T, param)
    label_param <- "gamma"
  } else {
    stop("Error! Process type not specified correctly!")
  }
  
  # Create Toeplitz matrix from the correlation vector
  corrmat <- toeplitz(corr)
  
  # Calculate the standard deviation
  sigma <- sqrt(sigmasq)
  
  # Cholesky decomposition
  m <- chol(corrmat)
  
  # Simulation
  gem <- numeric(nk)
  twmoment <- numeric(nk)
  
  # Loop over the number of simulations 
  for (j in 1:nrsimul) {
    # Generate random numbers
    y <- mu + sigma * m %*% rnorm(T)
    
    # Calculate empirical autocorrelation function of the generated time series
    empicorr <- numeric(T)
    
    # Calculate the population mean OR sample mean of the time series
    if (pop_mean == TRUE) {
      mn_y <- mu
      label_pop_mean <- "population"
    } else {
      mn_y <- mean(y)
      label_pop_mean <- "sample"
    }
    
    # Calculate the empirical autocorrelation function
    for (nkk in 1:nk) {
      empicorr[T-nkk] <- sum((y[1:(T-nkk)] - mn_y) * (y[(1+nkk):T] - mn_y)) / (T-nkk)
    }

    # Calculate means and second moments of the empirical autocorrelation functions    
    for (tt in 1:nk) {
      gem[tt] <- gem[tt] + (empicorr[T-tt] - gem[tt]) / j
      twmoment[tt] <- twmoment[tt] + empicorr[T-tt]^2
    }
  }
  
  # Calculate means and standard deviations of the empirical autocorrelation functions
  stdev <- numeric(nk)
  plusreeks <- numeric(nk)
  minreeks <- numeric(nk)
  
  # Calculate the standard deviation and confidence intervals of the empirical autocorrelation functions
  if (nrsimul > 1) {
    for (tt in 1:nk) {
      stdev[tt] <- sqrt((twmoment[tt] - nrsimul * gem[tt]^2) / (nrsimul - 1))
      plusreeks[tt] <- gem[tt] + 1.96 * stdev[tt]
      minreeks[tt] <- gem[tt] - 1.96 * stdev[tt]
    }
  }
  
  # Create the lag0 and lag1 matrices
  lag1 <- cbind(corr[2:(nk+1)], gem[1:nk], plusreeks[1:nk], minreeks[1:nk])
  lag0 <- rbind(c(1, 1, 1, 1), lag1)
  
  # Plot the theoretical and empirical autocorrelation functions, together with the confidence intervals
  plot(0:nk, lag0[,1], type = "l", ylim = c(ymin, ymax), main = paste0("simulation experiment GPH, T=", T,", ", label_pop_mean, " mean, ", label_param, "=", param), ylab = "", xlab = "")
  points(0:nk, lag0[,2], pch = 20)
  lines(0:nk, lag0[,3], lty = 2)
  lines(0:nk, lag0[,4], lty = 2)
}

# Usage:
#figure2_8_12_13_14_15(process = 1, nrsimul = 100, T=265, param = 0.25, nk = 50, mu = 0, sigmasq = 1, ymin = -1, ymax = 1, pop_mean = TRUE)
