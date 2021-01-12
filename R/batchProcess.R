#' @export
#' @import future.apply
#' @import rhdf5
#' @import data.table

batchProcessR <- function(filename, cfreqs = c(600, 1200), cores = availableCores()) {

  print(Sys.time())
  print(basename(filename))

  if (file.exists(paste0(tools::file_path_sans_ext(filename), "_out.txt"))) {

    return(message(paste(basename(filename), "has already been processed. Results are available in", paste0(tools::file_path_sans_ext(filename), "_out.txt"), sep = " ")))

  }



  #metalist <- hdf5metaReadR(filename)
  metalist <- hdf5metaLoadR(filename)
  filename <- filename
  assign("samprate", metalist$samprate, envir = .GlobalEnv)
  chandetails <- metalist$chandetails
  codedt <- metalist$codedt
  idx_lists <- metalist$idx_lists


  options(mc.cores = cores)
  options(future.globals.maxSize = 2097152000)
  plan(multisession)
  data_list <- future_lapply(X = idx_lists[1:60],
                             FUN = segmentProcessR,
                             cfreqs = cfreqs,
                             filename = filename,
                             samprate = samprate,
                             chandetails = chandetails,
                             codedt = codedt)

  dt_out <- rbindlist(data_list)
  fwrite(x = dt_out,
         file = paste0(tools::file_path_sans_ext(filename), "_out.txt"))

  return(print(paste(tools::file_path_sans_ext(filename), "finished processing @", Sys.time(), sep = " ")))

}
