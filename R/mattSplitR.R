#' @import data.table
#' @export

mattSplitR <- function(x) {

  bname <- gsub("\\..*", "", basename(x))

  nameparts <- tstrsplit(bname, "_")

  genotype <- nameparts[[1]]

  sex <- nameparts[[2]]

  ftime <- as.POSIXct(paste0(unlist(nameparts[4:8]), collapse = ""), format = "%m%d%Y%H%M%S")

  t_add <- as.numeric(tstrsplit(basename(x), "_")[[9]])

  if(t_add == 91) {

    t_add <-  10

  }

  stime <- ftime + (t_add - 1)*200

  return(list(genotype,
              sex,
              stime))

}
