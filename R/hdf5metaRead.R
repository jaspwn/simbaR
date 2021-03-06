#' @import rhdf5
#' @import data.table
#' @export

## reads appropriate meta data from hdf5 files and returns sample rate (samprate), channel details (chan details),
## data table with appropriate marker codes (codedt) and intervals to extract and analyses data over (idx_lists)

hdf5metaReadR <- function(filename) {

  names <- data.table(h5ls(filename))

  ## read in file channel information

  # channames <- data.table(name = names[group == "/", name])
  # samprate <<- h5read(filename, paste0(channames[1, name], "/interval/"))[1,]
  # stime <- h5read(filename, paste0(channames[1, name], "/start/"))[,1]
  # len <- h5read(filename, paste0(channames[1, name], "/length/"))[,1]
  # channames[, c("date", "sTod", "eTod", "tod", "treatment", "chan") := channameSplitR(name)]
  # chandetails <- channames[complete.cases(channames)]

  channames <- data.table(name = names[group == "/", name])
  samprate <<- h5read(filename, paste0(channames[1, name], "/interval/"))[1,]
  stime <- h5read(filename, paste0(channames[1, name], "/start/"))[,1]
  len <- h5read(filename, paste0(channames[1, name], "/length/"))[,1]
  channames[, chan := tail(tstrsplit(name, "_"), 1)]
  channames <- channames[name != "file"]

  ## find data channels
  chan_list <- lapply(X = channames[, name],
                      FUN = chanTypeR,
                      filename = filename)

  chan_list <- rbindlist(chan_list)
  chandetails <- merge(channames, chan_list, by = "name")
  chandetails[, chantype := ifelse(chantitle %in% c("DigMark", "Keyboard"), "marker", "data")]


  ## create data.table of marker channel
  marktimes <- h5read(filename, paste0(chandetails[chan == "Ch32", name], "/times/"))
  markcodes <- h5read(filename, paste0(chandetails[chan == "Ch32", name], "/codes/"))

  codedt <- data.table(time = marktimes[,1], ASCIIcode = markcodes[,1])
  rm(marktimes, markcodes)
  codedt[, lettercode := intToUtf8(ASCIIcode), by = time]

  ## create time data.table of file

  t <- seq.int(from = stime, to = stime + len*samprate, length.out = len)
  dt_time <- data.table(t)
  setattr(dt_time, "sorted", "t")
  rm(t)

  ## get time index for the start of each marker code
  codedt[, t_idx := dt_time[.(time), roll = "nearest", which = TRUE]]

  ## get time index of final data point

  end_idx <- dt_time[.(tail(codedt[, time], 1)+mean(diff(codedt[, time]))), roll = "nearest", which = TRUE]

  idx_pairs <- embed(codedt[, t_idx], 2)[, 2:1]
  idx_lists <- split(idx_pairs,
                     rep(1:nrow(idx_pairs), times = ncol(idx_pairs)))

  ## append final pair manually
  idx_lists <- append(idx_lists, list(c(tail(codedt[, t_idx], 1),
                                        end_idx)))

  ## analyse two minute segments
  # idx_trips <- embed(codedt$t_idx, 3)[seq.int(1,58,2), c(3,1)]
  # idx_trpl <- split(idx_trips,
  #                    rep(1:nrow(idx_trips), times = ncol(idx_trips)))
  # idx_lists <- append(idx_trpl, list(c(codedt$t_idx[length(codedt$t_idx)],
  #                                      end_idx)))


  return(list(samprate = samprate,
              chandetails = chandetails,
              codedt = codedt,
              idx_lists = idx_lists))

}

