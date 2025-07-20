# Function to compute the spectral density of an ARFIMA(0, d, 0) process
spectrum_arfima <- function(d, n.freq = 512) {
  freq <- seq(0, 0.5, length.out = n.freq)
  spec <- (2 * sin(pi * freq))^(-2 * d) / (2 * pi)
  
  return(list(freq = freq, spec = spec))
}
