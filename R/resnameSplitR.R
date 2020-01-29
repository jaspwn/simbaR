#' @import data.table
#' @export

resnameSplitR <- function(x) {

  nameparts <- tstrsplit(x, "_")

  if (length(nameparts) > 5) {

    nameparts <- nameparts[1:5]

  }

  return(nameparts)

}
