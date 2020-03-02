#' @import data.table
#' @export

mattReadR <- function(outlist) {

  data <- fread(outlist)
  data[, filename := outlist]
  data[, c("species", "sex", "temp", "tod", "date") := resnameSplitR(tools::file_path_sans_ext(basename(outlist)))]

  return(data)
}
