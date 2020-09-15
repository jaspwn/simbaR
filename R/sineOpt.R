#' @import data.table
#' @import nls2
#' @import TTR
#' @export

sineOptR <- function(y, t, samprate = dataSamprate, optMeth = "BFGS") {

  sampfreq <- 1/samprate

  #estimate offset
  o_est <- mean(y)

  #estimate amplitude
  #A_est <- max(y)-mean(y)
  # a1max <- max(y)-mean(y)
  # a1min <- mean(y)-min(y)
  # a1 <- (a1max+a1min)/2
  # A_est <- sqrt(o_est^2 + (a1^2/2))*1.41

  yoff <- y - o_est
  A_est <- sqrt((sum(abs((yoff))^2)/length(yoff)))*1.41

  #find indices where x goes from below offset to above offset
  ySMA <- TTR::SMA(y, n = 5)
  diffy <- which(diff(ySMA > o_est) == 1)

  #calculate time difference between first and last cycle
  tdiff <- t[tail(diffy, n=1)] - t[diffy[1]]

  #estimate frequencies from numbers of cycles/time difference
  f_est <- (length(diffy)-1)/tdiff

  #estimate phase by time difference between first midline crossing of data and non-phase shifted data
  y1 <- A_est*sin(f_est*2*t*pi) + o_est
  diffy1 <- which(diff(y1 > o_est) == 1)
  p_est <- (abs(t[diffy[1]]-t[diffy1[1]])) * f_est * 2*pi

  #check if all coefficients can be estimated

  if (any(is.na(c(A_est, o_est, f_est, p_est))) | length(c(A_est, o_est, f_est, p_est)) != 4 | f_est > 750) {

    print("coefficient estimation failed")

    return(list(fitStart = min(t),
                fitEnd = max(t),
                amplitude = A_est,
                freq = f_est,
                phase = p_est,
                offset = o_est,
                rsqr = NA,
                fit = "Coefficient estimation failed"))
  }

  if (grepl("singular convergence", nls2(formula = y ~ I(A*sin(2*pi*f*t + p) + o),
                                         start = list(A = A_est,
                                                      f = f_est,
                                                      p = p_est,
                                                      o = o_est),
                                         algorithm = "port",
                                         lower = c(-1000000, 0, -100, -1000000),
                                         upper = c(1000000, 1000000, 100, 1000000),
                                         control = nls.control(maxiter = 100,
                                                               warnOnly = TRUE))$message)) {

    A_coef <- A_est
    f_coef <- f_est
    p_coef <- p_est
    o_coef <- o_est

    reslm <- lm(y ~ I(A_coef*sin(2*pi*f_coef * t + p_coef) + o_coef))
    print("singular convergence")

    return(list(fitStart = min(t),
                fitEnd = max(t),
                amplitude = A_coef,
                freq = f_coef,
                phase = p_coef,
                offset = o_coef,
                rsqr = summary(reslm)$adj.r.squared,
                fit = "Singular convergence"))

  }

  if (grepl("Error", tryCatchWE(nls2(formula = y ~ I(A*sin(2*pi*f*t + p) + o),
                                     start = list(A = A_est,
                                                  f = f_est,
                                                  p = p_est,
                                                  o = o_est),
                                     algorithm = "port",
                                     lower = c(-1000000, 0, -100, -1000000),
                                     upper = c(1000000, 1000000, 100, 1000000),
                                     control = nls.control(maxiter = 100,
                                                           warnOnly = TRUE)))$value)) {

    A_coef <- A_est
    f_coef <- f_est
    p_coef <- p_est
    o_coef <- o_est

    reslm <- lm(y ~ I(A_coef*sin(2*pi*f_coef * t + p_coef) + o_coef))
    print("nls fit failed")

    return(list(fitStart = min(t),
                fitEnd = max(t),
                amplitude = A_coef,
                freq = f_coef,
                phase = p_coef,
                offset = o_coef,
                rsqr = summary(reslm)$adj.r.squared,
                fit = "FAILED"))

  } else {

    resnls <- nls2(formula = y ~ I(A*sin(2*pi*f*t + p) + o),
                   start = list(A = A_est,
                                f = f_est,
                                p = p_est,
                                o = o_est),
                   algorithm = "port",
                   lower = c(-1000000, 0, -100, -1000000),
                   upper = c(1000000, 1000000, 100, 1000000),
                   control = nls.control(maxiter = 100,
                                         warnOnly = TRUE))

  }


  if (grepl("Error", tryCatchWE(summary(resnls))$value)) {

    A_coef <- A_est
    f_coef <- f_est
    p_coef <- p_est
    o_coef <- o_est

    reslm <- lm(y ~ I(A_coef*sin(2*pi*f_coef * t + p_coef) + o_coef))
    print("nls fit failed")

    return(list(fitStart = min(t),
                fitEnd = max(t),
                amplitude = A_coef,
                freq = f_coef,
                phase = p_coef,
                offset = o_coef,
                rsqr = summary(reslm)$adj.r.squared,
                fit = "FAILED"))

  } else {

    A_coef <- abs(summary(resnls)$coefficients[1])
    f_coef <- summary(resnls)$coefficients[2]
    p_coef <- summary(resnls)$coefficients[3]
    o_coef <- summary(resnls)$coefficients[4]

    reslm <- lm(y ~ I(A_coef*sin(2*pi*f_coef * t + p_coef) + o_coef))
    print(summary(reslm)$adj.r.squared)

    return(list(fitStart = min(t),
                fitEnd = max(t),
                amplitude = A_coef,
                freq = f_coef,
                phase = p_coef,
                offset = o_coef,
                rsqr = summary(reslm)$adj.r.squared,
                fit = "nls2-port"))

  }




  #
  # fk <- fft(y)
  # fk <- fk[2:length(fk)/2+1]
  # fk <- 2*fk[seq(1, length(fk), by = 2)]/length(y)
  # freq <- (1:(length(fk)))* sampfreq/(2*length(fk))
  # fft_dt <- data.table(fur = fk, freq = freq, amp = Mod(fk))
  #
  # # A_est <- sqrt(mean(y)^2 + (((max(y)-min(y))/2)^2)/2) #RMS for amplitude estimate
  # # f_est <- fft_dt[amp == max(amp), freq] #freq estimate
  # # p_est <- asin(y[1]/A_est) #phase estimate
  # # o_est <- mean(y) #offset estimate
  #
  # A_est <- (max(y)-mean(y))
  # f_est <- fft_dt[amp == max(amp), freq] #freq estimate
  # p_est <- 100
  # o_est <- mean(y) #offset estimate

  # resnls <- nls(formula = y ~ I(A*sin(2*pi*f*t + p) + o),
  #               start = list(A = A_est,
  #                            f = f_est,
  #                            p = p_est,
  #                            o = o_est),
  #               algorithm = "port",
  #               lower = c(-1000000, 0, -100, -1000000),
  #               upper = c(1000000, 1000000, 100, 1000000),
  #               control = nls.control(printEval = TRUE,
  #                                     warnOnly = TRUE))

  ##old linear model optimisation
  # outres <-optim(par=c(A, f, p, o),
  #                fn=sineFitR,
  #                fitdat = x,
  #                t = t,
  #                method = optMeth,
  #                control = list(fnscale = -1))
  #
  #
  # return(list(amplitude = outres$par[1],
  #             freq = outres$par[2],
  #             phase = outres$par[3],
  #             offset = outres$par[4],
  #             rsqr = outres$value))

}
