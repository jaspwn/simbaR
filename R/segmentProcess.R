#' @import rhdf5
#' @import data.table
#' @import signal
#' @import TTR
#' @import windowscanr
#' @export

## takes a list of lists (idx_lists) containing indexes of segments to analyse in file (filename) exporting in HDF5
## (matlab export from spike2)

segmentProcessR <- function(idx_lists, filename, samprate, chandetails, codedt) {

  #message(ls())

  bpfilt <- butter(n = 4, W = c(600/(50000/2), 1200/(50000/2)), type = "pass", plane = "z")

  message(idx_lists[1])

  chan01 <- as.vector(h5read(filename,
                             paste0(chandetails[chan == "Ch1", name], "/values/"),
                             index = list(idx_lists[1]:idx_lists[2],1)))
  h5closeAll()
  chan02 <- as.vector(h5read(filename,
                             paste0(chandetails[chan == "Ch2", name], "/values/"),
                             index = list(idx_lists[1]:idx_lists[2],1)))
  h5closeAll()
  chan03 <- as.vector(h5read(filename,
                             paste0(chandetails[chan == "Ch3", name], "/values/"),
                             index = list(idx_lists[1]:idx_lists[2],1)))
  h5closeAll()
  chan04 <- as.vector(h5read(filename,
                             paste0(chandetails[chan == "Ch4", name], "/values/"),
                             index = list(idx_lists[1]:idx_lists[2],1)))
  h5closeAll()

  t_seg <- seq.int(from = codedt[t_idx == idx_lists[1], time],
                   to = codedt[t_idx == idx_lists[1], time] + length(chan01)*samprate,
                   length.out = length(chan01))

  dt <- data.table(t = t_seg,
                   chan01 = chan01,
                   chan02 = chan02,
                   chan03 = chan03,
                   chan04 = chan04)
  rm(t_seg, chan01, chan02, chan03, chan04)

  moltendt <- melt(dt,

                   measure.vars = c("chan01",
                                    "chan02",
                                    "chan03",
                                    "chan04"),
                   variable.name = "channel")

  ## rescale each channel
  moltendt[, code := codedt[t_idx == idx_lists[1], lettercode]]
  moltendt[, time := codedt[t_idx == idx_lists[1], time]]
  # moltendt[, scaled := scales::rescale(value, to = c(-1,1)),
  #          by = "channel"]
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
  data <- data[data[,n_dp > 500]]

  if (dim(data)[1] > 0) {

    data_sum <- data[, .(mint = min(t),
                         maxt = max(t),
                         time = time,
                         code = code),
                     by = c("channel", "group_no")]

    data_sum <- unique(data_sum)


    ## fit data over 10ms sliding window (50% slide)
    rol_win <- winScan(x = data,
                       groups = c("channel", "group_no"),
                       position = NULL,
                       values = c("bpfiltered"),
                       win_size = 500,
                       win_step = 250,
                       funs = c("slidingFitfreq", "slidingFitR"))

    rol_win <- data.table(rol_win)
    rol_merge <- merge(rol_win, data_sum, by = c("channel", "group_no"))

    ## select fits over 0.9
    rol_sig <- rol_merge[rol_merge[, bpfiltered_slidingFitR > 0.90]]

    rol_molten <- melt(data = rol_sig,
                       measure.vars = c("bpfiltered_slidingFitfreq",
                                        "bpfiltered_slidingFitR"))

    sum_dt <- rol_molten[, .(mean = mean(value, na.rm = TRUE),
                             median = median(value, na.rm = TRUE),
                             stdev = sd(value, na.rm = TRUE),
                             se = se(value, na.rm = TRUE),
                             n_fits = length(value),
                             min_t = mint,
                             max_t = maxt,
                             evlength = maxt - mint,
                             time = time,
                             code = code),
                         by = c("channel", "group_no" ,"variable")]
    #gc(TRUE)
    message(paste0("Segment commencing index ", idx_lists[1], " has finished processing"))
    return(unique(sum_dt))

  } else {

    return(message(paste0("Segment commencing index ", idx_lists[1], " contains no identifiable flyby events")))

  }


}
