#' @import data.table
#' @export

sineOptR <- function(x, t, samprate = dataSamprate, optMeth = "BFGS") {

  sampfreq <- 1/samprate

  fk <- fft(x)
  fk <- fk[2:length(fk)/2+1]
  fk <- 2*fk[seq(1, length(fk), by = 2)]/length(x)
  freq <- (1:(length(fk)))* sampfreq/(2*length(fk))
  fft_dt <- data.table(fur = fk, freq = freq, amp = Mod(fk))

  A <- max(fft_dt[, amp])
  f <- fft_dt[amp == A, freq]
  p <- 0
  o <- 0

  outres <-optim(par=c(A, f, p, o),
                 fn=sineFitR,
                 fitdat = x,
                 t = t,
                 method = optMeth,
                 control = list(fnscale = -1))


  return(list(amplitude = outres$par[1],
              freq = outres$par[2],
              rsqr = outres$value))

}
