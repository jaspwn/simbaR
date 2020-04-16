#' @export

## looks for and loads metadata.Rdata, if it doesnt exist call hdf5metaSaveR to create it.

hdf5metaLoadR <- function(filename) {

  if (file.exists(paste0(paste0(dirname(filename), "/metadata/", gsub(".mat", "_metadata.Rdata", basename(filename)))))) {

    load(paste0(paste0(dirname(filename), "/metadata/", gsub(".mat", "_metadata.Rdata", basename(filename)))))

    return(metalist)

  } else {

    hdf5metaSaveR(filename)

    return(metalist)

  }

}
