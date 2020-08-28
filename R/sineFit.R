#' @export

sineFitR <- function(par, fitdat, t) {

  A = rep(par[1], length(t)) #amplitude
  f = par[2] #frequency in Hz
  p = par[3] #phase in rads
  o = rep(par[4], length(t)) #offset

  reslm <- lm(fitdat ~ A*sin(2*pi*f*t + p) + o)
  summary(reslm)$adj.r.squared
}
