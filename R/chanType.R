#' @export
#' @import rhdf5
#' @import data.table

chanTypeR <- function(filename = filename, name = name) {

  chantitle <- intToUtf8(h5read(filename, paste0(name, "/title/")))

  return(data.table(name = name, chantitle = chantitle))

}
