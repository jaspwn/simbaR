#' @export

slidingFitR <- function(x, freq = 50000) {
  #message(ls())
  #x = dt[group_no == 456]$bpfiltered
  raw_ts <- ts(x, start = 0, frequency = freq)
  t_int <- seq.int(from = 1, to = 1 + length(raw_ts)/freq, length.out = length(raw_ts))

  ssp <- spectrum(raw_ts, log = "no", plot = FALSE)
  ini_freq <- 1/ssp$freq[ssp$spec==max(ssp$spec)]

  outres <-   optim(par = ini_freq,
                     t = t_int,
                     fn = fittedFit,
                     fitdat = raw_ts,
                     method = "Brent",
                     lower = 1/((1/ini_freq)+100),
                     upper = 1/((1/ini_freq)-100),
                     control = list(fnscale = -1))


  return(outres$value)
}
