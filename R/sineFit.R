#' @export

sineFitR <- function(par, fitdat, t) {

  A = par[1] #amplitude
  f = par[2] #frequency in Hz
  p = par[3] #phase in rads
  o = par[4] #offset

    reslm <- lm(fitdat ~ I(A*sin(2*pi*f*t + p) + o))
  summary(reslm)$adj.r.squared
}
