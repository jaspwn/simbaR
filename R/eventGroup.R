#' @export

eventGroupR <- function(my_vector, threshold) {
  my_rle <- rle(my_vector > threshold)
  my_rle$values <- seq_along(my_rle$values)
  inverse.rle(my_rle)
}

