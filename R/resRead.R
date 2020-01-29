#' @import data.table
#' @export

resReadR <- function(outlist) {

  data <- fread(outlist)
  data[, filename := outlist]
  data[, c("date", "stime", "tod", "treatment", "out") := resnameSplitR(tools::file_path_sans_ext(basename(outlist)))]

  return(data)
}
