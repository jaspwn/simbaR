#' @import tuneR
#' @import data.table
#' @import signal
#' @import windowscanr
#' @import TTR
#' @import ggplot2
#' @export

wavProcessR <- function() {
  filename <- file.choose()

  ## band pass filter for no playbacks
  bpfilt <- butter(n = 4, W = c(200/(44100/2), 1200/(44100/2)), type = "pass", plane = "z")

  ## read in wav file
  wav <- readWave(filename, toWaveMC = TRUE)

  ## read each channel into individual vector
  chan01 <- as.numeric(wav[,1]@.Data)
  chan02 <- as.numeric(wav[,2]@.Data)
  chan03 <- as.numeric(wav[,3]@.Data)
  chan04 <- as.numeric(wav[,4]@.Data)
  t <- seq.int(from = 0, to = length(wav[,1]@.Data)/wav@samp.rate, length.out = length(wav[,1]@.Data))

  ## create datatable of target channels and t vector
  dt <- data.table(t, chan01, chan02, chan03, chan04)

  ## melt datatable for datatable operations
  dt <- melt(dt,
             measure.vars = c("chan01",
                              "chan02",
                              "chan03",
                              "chan04"),
             variable.name = "channel")

  ## rescale each channel
  dt[, scaled := scales::rescale(value, to = c(-1,1)),
     by = "channel"]
  ## perform DC remove of each channel
  dt[, DC := removeDC(scaled, 44.1),
     by = "channel"]
  ## run band pass filter on DC remove signal
  dt[, bpfiltered := as.vector(filtfilt(bpfilt, x = DC)),
     by = "channel"]
  ## calculate sd of each channel
  dt[, stdev := sd(bpfiltered, na.rm = TRUE),
     by = "channel"]
  ## implement threshold to identify events
  dt[, thresh := ifelse(abs(bpfiltered) > (2*stdev), 1, 0),
     by = "channel"]
  ## smooth threshold
  dt[, smooth := TTR::SMA(thresh, n = 110),
     by = "channel"]
  ## create event groups
  dt[, group_no := eventGroupR(smooth, 0.2),
     by = "channel"]

  ## select events that meet smooth threshold
  dt_thresh <- dt[smooth > 0.2]

  ## calculate length of each event but taking the length of time above threshold
  event_length <- dt_thresh[, .(n_dp = length(t)),
                            by = c("channel", "group_no")]


  ## merge thresholded events with event length
  data <- merge(dt_thresh, event_length, by = c("channel", "group_no"))
  ## select events longer than 50ms (0.05s) 0.05*44100 sample rate
  data <- data[data[,n_dp > 2205]]

  ## fit data over 10ms sliding window (50% slide)
  rol_win <- winScan(x = data,
                     groups = c("channel", "group_no"),
                     position = NULL,
                     values = c("bpfiltered"),
                     win_size = 441,
                     win_step = 220,
                     funs = c("slidingFitfreq", "slidingFitR"))

  rol_win <- data.table(rol_win)

  ## select fits over 0.9
  rol_sig <- rol_win[bpfiltered_slidingFitR > 0.90]

  rol_molten <- melt(data = rol_sig,
                     measure.vars = c("bpfiltered_slidingFitfreq",
                                      "bpfiltered_slidingFitR"))
  ## calculate summary stats for each fit
  sum_dt <- rol_molten[, .(mean = mean(value, na.rm = TRUE),
                           median = median(value, na.rm = TRUE),
                           stdev = sd(value, na.rm = TRUE),
                           se = se(value, na.rm = TRUE),
                           n_fits = length(value)),
                       by = c("channel", "group_no" ,"variable")]

  ## add event_lengths to summary data
  sum_dt <- merge(sum_dt, event_length, by = c("channel", "group_no"))
  sum_dt[, ev_length := n_dp/44100]
  sum_dt[, file := tools::file_path_sans_ext(basename(filename))]
  sum_dt[, chamber := "semifield"]
  sum_dt[, sex := "males"]
  sum_dt[, temp := 24]
  sum_dt[, age := "3"]

  summary_dt <<- sum_dt
  out_dt <<- dt

  # write.csv(out_dt, paste0("./data/semi_field/output/", tools::file_path_sans_ext(basename(filename)), "_summary"))

}








