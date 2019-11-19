## merge data datatable with marker data table

dt[, code := cut(t, breaks = c(codedt$time, max(codedt$time)+60), labels = codedt$lettercode)]

#test <- cut(t, breaks = c(codedt$time, max(codedt$time)+60), labels = codedt$lettercode)

## remove unused data to save space
dt <- dt[complete.cases(dt)]

codepairs <- list(c("O", "o"),
                  c("B", "b"),
                  c("C", "c"),
                  c("D", "d"),
                  c("E", "e"),
                  c("G", "g"),
                  c("H", "h"),
                  c("I", "i"),
                  c("J", "j"),
                  c("K", "k"))

codepairs <- list(c("O", "o"))

markerProcessR <- function(codepair, dt) {

  dtsel <- dt[code %in% codepair]
  print(codepair)
  #dtsel <- dt[code %in% c("A", "a")]

  ## melt datatable for datatable operations
  moltendt <- melt(dt,

                   measure.vars = c("chan01",
                                    "chan02",
                                    "chan03",
                                    "chan04"),
                   variable.name = "channel")

  head(moltendt)
  ## rescale each channel
  moltendt[, scaled := scales::rescale(value, to = c(-1,1)),
           by = "channel"]
  ## perform DC remove of each channel
  moltendt[, DC := removeDC(scaled, 50),
           by = "channel"]
  ## run band pass filter on DC remove signal
  moltendt[, bpfiltered := as.vector(filtfilt(bpfilt, x = DC)),
           by = "channel"]
  ## calculate sd of each channel
  moltendt[, stdev := sd(bpfiltered, na.rm = TRUE),
           by = c("channel", "code")]
  ## implement threshold to identify events
  moltendt[, thresh := ifelse(abs(bpfiltered) > (2*stdev), 1, 0),
           by = "channel"]
  ## smooth threshold
  moltendt[, smooth := TTR::SMA(thresh, n = 110),
           by = "channel"]
  ## create event groups
  moltendt[, group_no := eventGroupR(smooth, 0.2),
           by = "channel"]

  ## select events that meet smooth threshold
  dt_thresh <- moltendt[smooth > 0.2]

  ## calculate length of each event but taking the length of time above threshold
  event_length <- dt_thresh[, .(n_dp = length(t)),
                            by = c("channel", "group_no")]


  ## merge thresholded events with event length
  data <- merge(dt_thresh, event_length, by = c("channel", "group_no"))
  ## select events longer than 10ms (0.01s) 0.01*50000 sample rate
  data <- data[data[,n_dp > 500]]



  data_sum <- data[, .(mint = min(t),
                       maxt = max(t),
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
                           evlength = maxt - mint,
                           code = code),
                       by = c("channel", "group_no" ,"variable")]

  return(unique(sum_dt))

}

data_list <- lapply(codepairs, markerProcessR, dt = dt)
