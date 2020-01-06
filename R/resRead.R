#' @import data.table
#' @export

resReadR <- function(outlist) {

  data <- fread(outlist)
  data[, filename := tools::file_path_sans_ext(basename(outlist))]
  data[, c("date", "stime-etime", "tod", "treatment", "out") := tstrsplit(filename, "_")]

  return(data)
}
