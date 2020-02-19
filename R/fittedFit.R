#' @export

fittedFit <- function(par, fitdat, t) {
    #sum(residuals(lm(x ~ sin(2*pi/par*test$t) + cos(2*pi/par*test$t))))
    #reslm <- lm(fitdat ~ sin(2*pi/(par*t)) + cos(2*pi/(par*t)))
    reslm <- lm(fitdat ~ sin(2*pi*par*t) + cos(2*pi*par*t))
    summary(reslm)$adj.r.squared
  }

