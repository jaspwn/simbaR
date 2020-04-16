#' @import rhdf5
#' @export

## reads and saves metadate from .mat files exported from spike (.smrx) files

hdf5metaSaveR <- function(filename) {

  metalist <- hdf5metaReadR(filename)

  ## save output of hdf5metaReadR to .Rdata file
  if (!dir.exists(paste0(dirname(filename), "/metadata/"))) {
    dir.create(paste0(dirname(filename), "/metadata/"), recursive = TRUE)
  }

  save(metalist, file = paste0(paste0(dirname(filename), "/metadata/", gsub(".mat", "_metadata.Rdata", basename(filename)))))
  return(metalist)

}
