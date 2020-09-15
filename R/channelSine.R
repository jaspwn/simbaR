#' @export
#' @import rhdf5
#' @import data.table

channelSineR <- function(filename, dataChan = "Ch1", markChan = "Ch4", offset = -200, width = 100, idxSkip = c()) {

  print(Sys.time())
  print(basename(filename))

  if (file.exists(paste0(tools::file_path_sans_ext(filename), "_out.txt"))) {

    return(message(paste(basename(filename), "has already been processed. Results are available in", paste0(tools::file_path_sans_ext(filename), "_out.txt"), sep = " ")))

  }

  names <- data.table(h5ls(filename))

  ## read in file channel information
  channames <- data.table(name = names[group == "/", name])
  channames[, "chan" := gsub("^.*_", "", name)]

  ## create data.table for data channel

  dataTitle <- intToUtf8(h5read(filename, paste0(channames[chan == dataChan, name], "/title/"))[1,])
  dataSamprate <- h5read(filename, paste0(channames[chan == dataChan, name], "/interval/"))[1,]
  dataStime <- h5read(filename, paste0(channames[chan == dataChan, name], "/start/"))[,1]
  dataLen <- h5read(filename, paste0(channames[chan == dataChan, name], "/length/"))[,1]

  assign("dataSamprate", dataSamprate, envir = .GlobalEnv)

  ## create time data.table for entire file
  t <- seq.int(from = dataStime, to = dataStime + dataLen*dataSamprate, length.out = dataLen)
  dt_time <- data.table(t)
  setattr(dt_time, "sorted", "t")
  rm(t)

  ## calculate import start/stop times based on offset and width
  startTime <- (offset-width*2)/1000
  endTime <- (offset+width*3)/1000

  ## calculate analysis start/stop times based on offset and width
  startFit <- (2*width/1000)/dataSamprate + 1
  endFit <- startFit + (width/1000)/dataSamprate - 1

  ## create data.table for marker channel

  markTitle <- intToUtf8(h5read(filename, paste0(channames[chan == markChan, name], "/title/"))[1,])
  markCodes <- h5read(filename, paste0(channames[chan == markChan, name], "/codes/"))
  markTimes <- h5read(filename, paste0(channames[chan == markChan, name], "/times/"))

  codedt <- data.table(time = markTimes[,1], ASCIIcode = markCodes[,1])
  rm(markTimes, markCodes)
  codedt[, lettercode := intToUtf8(ASCIIcode), by = time]

  ## get time index for the start of each marker code
  codedt[, t_idx := dt_time[.(time), roll = "nearest", which = TRUE]]
  codedt[, start_time := dt_time[.(time+startTime), roll = "nearest"]]
  codedt[, end_time := dt_time[.(time+endTime), roll = "nearest"]]
  codedt[, start_idx := dt_time[.(time+startTime), roll = "nearest", which = TRUE]]
  codedt[, end_idx := dt_time[.(time+endTime), roll = "nearest", which = TRUE]]
  rm(dt_time)

  idx_lists <- mapply(c, codedt[, start_idx], codedt[, end_idx], SIMPLIFY = FALSE, USE.NAMES = FALSE)

  if (length(idxSkip) > 0) {

    idx_lists[idxSkip] <- NULL
    print(paste0("Skipped idx pairs ", idxSkip))

  }

  data_list <- lapply(X = idx_lists,
                      FUN = sectionSineR,
                      filename = filename,
                      channames = channames,
                      codedt = codedt,
                      dataChan = dataChan,
                      dataSamprate = dataSamprate,
                      startFit = startFit,
                      endFit = endFit)

  dt_out <- rbindlist(data_list)
  dt_out[, file := basename(filename)]

  fwrite(x = dt_out,
         file = paste0(tools::file_path_sans_ext(filename), "_out.txt"))

  return(print(paste(tools::file_path_sans_ext(filename), "finished processing @", Sys.time(), sep = " ")))

  #return(dt_out)

}
