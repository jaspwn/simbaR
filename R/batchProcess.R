
batchProcessR <- function(filename) {
  
  print(basename(filename))
  filename <- filename
  names <- data.table(h5ls(filename))
  
  ## read in file channel information
  channames <- data.table(name = names[group == "/", name])
  samprate <<- h5read(filename, paste0(channames$name[1], "/interval/"))[1,]
  stime <- h5read(filename, paste0(channames$name[1], "/start/"))[,1]
  len <- h5read(filename, paste0(channames$name[1], "/length/"))[,1]
  channames[, c("date", "sTod", "eTod", "tod", "treatment", "chan") := tstrsplit(name, "_")]
  chandetails <<- channames[complete.cases(channames)]
  
  ## create data.table of marker channel
  marktimes <- h5read(filename, paste0(chandetails[chan == "Ch32"]$name, "/times/"))
  markcodes <- h5read(filename, paste0(chandetails[chan == "Ch32"]$name, "/codes/"))
  
  codedt <<- data.table(time = marktimes[,1], ASCIIcode = markcodes[,1])
  rm(marktimes, markcodes)
  codedt[, lettercode := intToUtf8(as.vector(ASCIIcode)), by = "ASCIIcode"]
  
  ## create time data.table of file
  
  t <- seq.int(from = stime, to = stime + len*samprate, length.out = len)
  dt_time <- data.table(t)
  setattr(dt_time, "sorted", "t")
  rm(t)
  
  ## get time index for the start of each marker code
  codedt[, t_idx := dt_time[.(time), roll = "nearest", which = TRUE]]
  
  ## get time index of final data point
  end_idx <- dt_time[.(codedt$time[length(codedt$time)]+60), roll = "nearest", which = TRUE]
  
  idx_pairs <- embed(codedt$t_idx, 2)[, 2:1]
  idx_lists <- split(idx_pairs, 
                     rep(1:nrow(idx_pairs), times = ncol(idx_pairs)))
  ## append final pair manually
  idx_lists <<- append(idx_lists, list(c(codedt$t_idx[length(codedt$t_idx)],
                                         end_idx)))
  
  
  availableCores()
  options(mc.cores = 6)
  options(future.globals.maxSize = 2097152000)
  plan(multisession)
  data_list <- future_lapply(idx_lists, segmentProcessR, filename = filename)
  
  dt_out <- rbindlist(data_list)
  fwrite(x = dt_out,
         file = paste0(tools::file_path_sans_ext(filename), "_out.txt"))
  
  return(basename(filename))
  
}