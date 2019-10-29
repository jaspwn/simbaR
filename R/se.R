#' @export

se <- function(x, na.rm = TRUE) {
  sd(x)/sqrt(length(x))
}
