#' @import data.table
#' @export

channameSplitR <- function(x) {

  nameparts <- tstrsplit(x, "_")

  if (length(nameparts) < 6) {

    nameparts <- append(nameparts, list(strftime(as.POSIXct(nameparts[[2]], format = "%H%M") + 3600, format = "%H%M")), 2)

  }

  return(nameparts)

}
