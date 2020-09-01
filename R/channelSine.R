#' @export
#' @import rhdf5
#' @import data.table

channelSineR <- function(filename, dataChan = "Ch1", markChan = "Ch4", offset = -200, width = 100) {

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
  endTime <- (offset+width*2)/1000

  ## calculate analysis start/stop times based on offset and width
  startFit <- abs((startTime-(offset/1000))/dataSamprate) + 1
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
  codedt[, start_idx := dt_time[.(time+startTime), roll = "nearest", which = TRUE]]
  codedt[, end_idx := dt_time[.(time+endTime), roll = "nearest", which = TRUE]]
  rm(dt_time)

  idx_lists <- mapply(c, codedt[, start_idx], codedt[, end_idx], SIMPLIFY = FALSE, USE.NAMES = FALSE)

  data_list <- future_lapply(X = idx_lists,
                             FUN = sectionSineR,
                             filename = filename,
                             channames = channames,
                             codedt = codedt,
                             dataChan = dataChan,
                             dataSamprate = dataSamprate,
                             startFit = startFit,
                             endFit = endFit)

  dt_out <- rbindlist(data_list)

  return(dt_out)

}
