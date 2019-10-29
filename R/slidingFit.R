slidingFit <- function(x) {
  
  raw_ts <- ts(x, start = 0, frequency = 44100)
  t_int <- seq.int(from = 1, to = 1.005, length.out = length(raw_ts))
  
  ssp <- spectrum(raw_ts, log = "no")
  ini_freq <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
  
  outres <-   optim(par = ini_freq,
                     t = t_int,
                     fn = fittedFit,
                     x = raw_ts,
                     method = "Brent",
                     lower = 1/((1/ini_freq)+100),
                     upper = 1/((1/ini_freq)-100),
                     control = list(fnscale = -1))

  
  return(outres$par) 
}
