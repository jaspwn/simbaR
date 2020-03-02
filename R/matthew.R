#' @import rhdf5
#' @import data.table
#' @import signal
#' @import TTR
#' @import zoo
#' @export

matthew <- function(filename, cfreqs = c(600, 1200), srate = 50000) {

  if (file.exists(paste0(tools::file_path_sans_ext(filename), "_out.txt"))) {

    return(message(paste(basename(filename), "has already been processed. Results are available in", paste0(tools::file_path_sans_ext(filename), "_out.txt"), sep = " ")))

  }

  bpfilt <- butter(n = 4, W = c(cfreqs[1]/(50000/2), cfreqs[2]/(50000/2)), type = "pass", plane = "z")

  samprate <- 1/srate

  metadata <- data.table(h5ls(filename))

  chandt <- data.table( "chan01" = h5read(filename, metadata[, name])[,1])

  h5closeAll()

  t <- seq.int(from = 0,
               to = length(chandt[, chan01])*samprate,
               length.out = length(chandt[, chan01]))

  dt <- cbind(t, chandt)

  moltendt <- melt(dt,

                   measure.vars = c("chan01"),
                   variable.name = "channel")


  ## perform DC remove of each channel
  moltendt[, remsat := ifelse(value == -Inf, -1,
                              ifelse(value == Inf, 1, value))]

  moltendt[, DC := removeDC(remsat, 50),
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
                                          width = 500, by = 250,
                                          FUN = rollingFitfreq,
                                          srate = samprate,
                                          stime = min(t),
                                          partial = FALSE,
                                          align = "left"),
                      rsqr = rollapply(data = bpfiltered,
                                       width = 500, by = 250,
                                       FUN = rollingFitR,
                                       srate = samprate,
                                       stime = min(t),
                                       partial = FALSE,
                                       align = "left"),
                      fitime = rollapply(data = t,
                                         width = 500, by = 250,
                                         FUN = min,
                                         partial = FALSE,
                                         align = "left")),
                  by = c("channel", "group_no")]



  rol_merge <- merge(rol_win, data_sum, by = c("channel", "group_no"))

  ## select fits over 0.9
  rol_sig <- rol_merge[rol_merge[, rsqr > 0.90]]

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

  fwrite(x = unique(sum_dt),
         file = paste0(tools::file_path_sans_ext(filename), "_out.txt"))

  return(print(paste(tools::file_path_sans_ext(filename), "finished processing @", Sys.time(), sep = " ")))

}
