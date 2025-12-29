# Function to calculate theoretical ACF for ARFIMA(0, 0.25, 0) with correct formula
acf_arfima <- function(n, d) {
  # Initialize the ACF vector
  acf_vals <- numeric(n + 1)
  acf_vals[1] <- 1  # ACF at lag 0 is always 1
  
  # Compute the ACF for each lag using the correct formula
  for (k in 1:n) {
    acf_vals[k + 1] <- gamma(1 - d) * gamma(k + d) / (gamma(d) * gamma(k + 1 - d))
  }
  
  return(acf_vals)
}

