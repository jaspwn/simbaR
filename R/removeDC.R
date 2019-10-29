#' @export

removeDC <- function(x, p) {

  z <- numeric(length = length(x))

  for (i in seq_along(x)) {


    if(i < p) {
      z[i] <- x[i] - mean(x[0:(i+p)])
    }

    if(i >= p) {
      z[i] <- x[i] - mean(x[(i-p):(i+p)])
    }

    if(i >= length(x) - p) {
      z[i] <- x[i] - mean(x[(i-p):length(x)])
    }

  }

  #raw_DC <<- z
  return(z)
}
