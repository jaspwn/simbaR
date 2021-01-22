#' @import data.table
#' @export

mattReadR <- function(outlist) {

  data <- fread(outlist)
  data[, filename := outlist]
  data[, c("species", "sex", "stime") := mattSplitR(outlist)]

  return(data)
}
