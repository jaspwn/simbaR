#' @export

rollingFitfreq <- function(x, srate, stime) {

  raw_ts <- ts(x, start = stime, frequency = 1/srate)
  t_int <- seq.int(from = stime, to = stime + length(raw_ts)*srate, length.out = length(raw_ts))

  ssp <- spectrum(raw_ts, log = "no", plot = FALSE)
  ini_freq <- ssp$freq[ssp$spec==max(ssp$spec)]

  outres <-   optim(par = ini_freq,
                    t = t_int,
                    fn = fittedFit,
                    fitdat = raw_ts,
                    method = "Brent",
                    lower = ini_freq-100,
                    upper = ini_freq+100,
                    control = list(fnscale = -1))


  return(outres$par)
}
