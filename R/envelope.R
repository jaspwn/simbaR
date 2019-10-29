#' @export

##

envelopeR <- function(my_vector, samprate = samprate, samplingrate = 1/samprate, wsize = 1/samprate/100) {

  #full wave rectification of data
  rectdata <- abs(my_vector)

  #perform moving average over selected window size
  n <- length(rectdata)
  fvalues <- numeric()
  for (i in 1:n) {
    w_start <- max(1, i - wsize)
    w_end <- min(n, i + wsize)
    fvalues[i] <- mean(rectdata[w_start:w_end])
  }

  return(fvalues)

}
