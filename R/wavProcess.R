#' @import rhdf5
#' @import data.table
#' @import signal
#' @import TTR
#' @import zoo
#' @export

wavProcessR <- function(filename, cfreqs = c(600, 1200)) {

  if (file.exists(paste0(tools::file_path_sans_ext(filename), "_out.txt"))) {

    return(message(paste(basename(filename), "has already been processed. Results are available in", paste0(tools::file_path_sans_ext(filename), "_out.txt"), sep = " ")))

  }

  names <- data.table(h5ls(filename))

  ## read in file channel information
  channames <- data.table(name = names[group == "/", name])
  samprate <<- h5read(filename, paste0(channames[1, name], "/interval/"))[1,]
  stime <- h5read(filename, paste0(channames[1, name], "/start/"))[,1]
  len <- h5read(filename, paste0(channames[1, name], "/length/"))[,1]
  channames[, c("date", "recID",  "treatment", "chan") := tstrsplit(name, "_")]
  chandetails <- channames[complete.cases(channames)]

  bpfilt <- butter(n = 4, W = c(cfreqs[1]*samprate*2, cfreqs[2]*samprate*2), type = "pass", plane = "z")

  #message(idx_lists[1])

  chan01 <- as.vector(h5read(filename,
                             paste0(chandetails[chan == "Ch1", name], "/values/")))
  h5closeAll()
  chan02 <- as.vector(h5read(filename,
                             paste0(chandetails[chan == "Ch2", name], "/values/")))
  h5closeAll()
  chan03 <- as.vector(h5read(filename,
                             paste0(chandetails[chan == "Ch3", name], "/values/")))
  h5closeAll()
  chan04 <- as.vector(h5read(filename,
                             paste0(chandetails[chan == "Ch4", name], "/values/")))
  h5closeAll()
  chan05 <- as.vector(h5read(filename,
                             paste0(chandetails[chan == "Ch5", name], "/values/")))
  h5closeAll()
  chan06 <- as.vector(h5read(filename,
                             paste0(chandetails[chan == "Ch6", name], "/values/")))
  h5closeAll()

  t_seg <- seq.int(from = 0,
                   to = length(chan01)*samprate,
                   length.out = length(chan01))


  dt <- data.table(t = t_seg,
                   chan01 = chan01,
                   chan02 = chan02,
                   chan03 = chan03,
                   chan04 = chan04,
                   chan05 = chan05,
                   chan06 = chan06)

  rm(t_seg, chan01, chan02, chan03, chan04, chan05, chan06)

  moltendt <- melt(dt,

                   measure.vars = c("chan01",
                                    "chan02",
                                    "chan03",
                                    "chan04",
                                    "chan05",
                                    "chan06"),
                   variable.name = "channel")

  ## perform DC remove of each channel
  moltendt[, DC := removeDC(value, 50),
           by = "channel"]
  ## run band pass filter on DC remove signal
  moltendt[, bpfiltered := as.vector(filtfilt(bpfilt, x = DC)),
           by = "channel"]
  moltendt[, enved := envelopeR(bpfiltered, samprate = samprate),
           by = "channel"]
  moltendt[, envsd := mean(abs(bpfiltered), na.rm = TRUE),
           by = c("channel")]
  moltendt[, group_no := eventGroupR(enved, 2*envsd),
           by = "channel"]


  ## select events that meet smooth threshold
  dt_thresh <- moltendt[enved > 2*envsd]

  ## calculate length of each event but taking the length of time above threshold
  event_length <- dt_thresh[, .(n_dp = length(t)),
                            by = c("channel", "group_no")]


  ## merge thresholded events with event length
  data <- merge(dt_thresh, event_length, by = c("channel", "group_no"))
  ## select events longer than 10ms (0.01s) 0.01*50000 sample rate
  data <- data[data[,n_dp > 0.01/samprate]]

  ## save output of fly-by events that pass threshold
  if (!dir.exists(paste0(dirname(filename), "/raw_events/"))) {
    dir.create(paste0(dirname(filename), "/raw_events/"), recursive = TRUE)
  }

  fwrite(data, file = paste0(paste0(dirname(filename), "/raw_events/", gsub(".mat", ".txt", basename(filename)))))

  if (dim(data)[1] > 0) {

    data_sum <- data[, .(mint = min(t),
                         maxt = max(t),
                         minDC = min(DC),
                         maxDC = max(DC),
                         minfilt = min(bpfiltered),
                         maxfilt = max(bpfiltered)),
                     by = c("channel", "group_no")]

    data_sum <- unique(data_sum)


    ## use rollapply instead of winScan
    rol_win <- data[, .(fitfreq = rollapply(data = bpfiltered,
                                            width = 0.01/samprate,
                                            by = ceiling(0.01/samprate/2),
                                            FUN = rollingFitfreq,
                                            srate = samprate,
                                            stime = min(t),
                                            partial = FALSE,
                                            align = "left"),
                        rsqr = rollapply(data = bpfiltered,
                                         width = 0.01/samprate,
                                         by = ceiling(0.01/samprate/2),
                                         FUN = rollingFitR,
                                         srate = samprate,
                                         stime = min(t),
                                         partial = FALSE,
                                         align = "left"),
                        fitime = rollapply(data = t,
                                           width = 0.01/samprate,
                                           by = ceiling(0.01/samprate/2),
                                           FUN = min,
                                           partial = FALSE,
                                           align = "left")),
                    by = c("channel", "group_no")]

    rol_merge <- merge(rol_win, data_sum, by = c("channel", "group_no"))

    ## select fits over 0.9
    rol_sig <- rol_merge[rol_merge[, rsqr > 0.90]]

    if (!dir.exists(paste0(dirname(filename), "/fit_events/"))) {
      dir.create(paste0(dirname(filename), "/fit_events/"), recursive = TRUE)
    }

    fwrite(rol_sig, file = paste0(paste0(dirname(filename), "/fit_events/", gsub(".mat", "_fit.txt", basename(filename)))))

    rol_molten <- melt(data = rol_sig,
                       measure.vars = c("fitfreq",
                                        "rsqr"))


    sum_dt <- rol_molten[, .(mean = mean(value, na.rm = TRUE),
                             median = median(value, na.rm = TRUE),
                             stdev = sd(value, na.rm = TRUE),
                             se = se(value, na.rm = TRUE),
                             n_fits = length(value),
                             min_t = mint,
                             max_t = maxt,
                             evlength = maxt - mint,
                             min_DC = minDC,
                             max_DC = maxDC,
                             min_filt = minfilt,
                             max_filt = maxfilt),
                         by = c("channel", "group_no" ,"variable")]

    sum_dt <- unique(sum_dt)

    fwrite(x = sum_dt,
           file = paste0(tools::file_path_sans_ext(filename), "_out.txt"))

    return()

  } else {

    return(message(paste0(basename(filename), " contains no identifiable flyby events")))

  }


}


