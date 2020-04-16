library(simbaR)
library(rhdf5)
library(data.table)
library(signal)
library(windowscanr)
library(ggplot2)
library(biosignalEMG)
library(cowplot)
library(zoo)

#bandpass filter design
bpfilt <- butter(n = 4, W = c(600/(50000/2), 1200/(50000/2)), type = "pass", plane = "z")

filename <- file.choose()
resfile <- file.choose()

## read in file channel information

meta <- hdf5metaReadR(filename)

samprate <- meta$samprate
chandetails <- meta$chandetails
codedt <- meta$codedt
idx_lists <- meta$idx_lists

## segment to analyse

idx <- idx_lists[[37]]

out_dt <- segmentProcessR(idx_lists = idx,
                          filename = filename,
                          samprate = samprate,
                          chandetails = chandetails,
                          codedt = codedt)


# create data.table out of data channels

chan01 <- as.vector(h5read(filename,
                           paste0(chandetails[chan == "Ch1", name], "/values/"),
                           index = list(idx[1]:idx[2],1)))
chan02 <- as.vector(h5read(filename,
                           paste0(chandetails[chan == "Ch2", name], "/values/"),
                           index = list(idx[1]:idx[2],1)))
chan03 <- as.vector(h5read(filename,
                           paste0(chandetails[chan == "Ch3", name], "/values/"),
                           index = list(idx[1]:idx[2],1)))
chan04 <- as.vector(h5read(filename,
                           paste0(chandetails[chan == "Ch4", name], "/values/"),
                           index = list(idx[1]:idx[2],1)))

t_seg <- seq.int(from = codedt[t_idx == idx[1], time],
                 to = codedt[t_idx == idx[1], time] + length(chan01)*samprate,
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

## add letter codes and start time to each segment
moltendt[, code := codedt[t_idx == idx[1], lettercode]]
moltendt[, time := codedt[t_idx == idx[1], time]]

## perform DC remove of each channel
moltendt[, DC := removeDC(value, 50),
         by = "channel"]
## run band pass filter on DC remove signal
moltendt[, bpfiltered := as.vector(filtfilt(bpfilt, x = DC)),
         by = "channel"]
## group using custom envelope function
moltendt[, enved := envelopeR(bpfiltered, samprate = samprate),
         by = "channel"]
moltendt[, envsd := mean(abs(bpfiltered), na.rm = TRUE),
         by = c("channel")]
moltendt[, group_no := eventGroupR(enved, 2*envsd),
         by = "channel"]

noise01 <- moltendt[t %between% c(2, 3), c("t", "channel", "value", "bpfiltered")]

save(noise01, file = "./data/noise01.rda", compress = TRUE)

uscaled <- round(scales::rescale(u, to = c(-7500000,7500000)))


flyby <- moltendt[channel == "chan01"][group_no == 44, value]
white_n <- noise(kind = "white", duration = 1, samp.rate = 50000, stereo = FALSE, xunit = "time")@left
pink_n <- noise(kind = "pink", duration = 1, samp.rate = 50000, stereo = FALSE, xunit = "time")@left
red_n <- noise(kind = "red", duration = 1, samp.rate = 50000, stereo = FALSE, xunit = "time")@left
chan01_n <- noise01[channel == "chan01", value]
chan02_n <- noise01[channel == "chan02", value]
chan03_n <- noise01[channel == "chan03", value]
chan04_n <- noise01[channel == "chan04", value]



uscaled <- c(round(scales::rescale(chan01_n, to = c(-7500000,7500000))),
             round(scales::rescale(flyby, to = c(-7500000,7500000))),
             round(scales::rescale(chan01_n, to = c(-7500000,7500000))))

uscaled <- round(scales::rescale(c(chan01_n, flyby, chan01_n), to = c(-7500000,7500000)))

w  <-  Wave(uscaled, samp.rate = 50000, bit = 24) #make the wave variable
x <- stereo(w,w)
play(x, player = "/usr/bin/vlc", "--play-and-stop --audio-visual visual")



dt_thresh <- moltendt[enved > 2*envsd]

## calculate length of each event but taking the length of time above threshold
event_length <- dt_thresh[, .(n_dp = length(t)),
                          by = c("channel", "group_no")]


## merge thresholded events with event length
data <- merge(dt_thresh, event_length, by = c("channel", "group_no"))
## select events longer than 10ms (0.01s) 0.01*50000 sample rate
data <- data[data[,n_dp > 500]]

if (!dir.exists(paste0(dirname(filename), "/raw_events/"))) {
  dir.create(paste0(dirname(filename), "/raw_events/"), recursive = TRUE)
}

fwrite(data, file = paste0(paste0(dirname(filename), "/raw_events/", gsub(".mat", paste0("_", idx[1], ".txt"), basename(filename)))))

  data_sum <- data[, .(mint = min(t),
                       maxt = max(t),
                       minDC = min(DC),
                       maxDC = max(DC),
                       minfilt = min(bpfiltered),
                       maxfilt = max(bpfiltered),
                       time = time,
                       code = code),
                   by = c("channel", "group_no")]

  data_sum <- unique(data_sum)


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

  if (!dir.exists(paste0(dirname(filename), "/fit_events/"))) {
    dir.create(paste0(dirname(filename), "/fit_events/"), recursive = TRUE)
  }

  fwrite(rol_sig, file = paste0(paste0(dirname(filename), "/fit_events/", gsub(".mat", paste0("_fit_", idx[1], ".txt"), basename(filename)))))

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
                         time = time,
                         code = code),
                     by = c("channel", "group_no" ,"variable")]

sum_dt <- unique(sum_dt)



  orig_rec <- moltendt
  chan <- "chan01"
  grp_no <- 21

  event_start <- min(orig_rec[channel == chan][group_no == grp_no, t]) - 0.1
  event_end <- max(orig_rec[channel == chan][group_no == grp_no, t]) + 0.1



  ggplot(data = orig_rec[t %between% c(event_start, event_end)], aes(x = t)) +
    geom_rect(inherit.aes = FALSE, aes(xmin = min(orig_rec[channel == chan][group_no == grp_no, t]),
                                       xmax = max(orig_rec[channel == chan][group_no == grp_no, t]),
                                       ymin = -Inf, ymax = Inf, fill = "group"), fill = "#ffeda0") +
    #geom_line(aes(y = DC, colour = "DC removed"), size = 0.2) +
    geom_line(aes(y = bpfiltered, colour = "bp_filtered"), size = 0.5) +
    geom_line(aes(y = enved, colour = "envelope")) +
    geom_line(aes(y = 2*envsd, colour = "threshold")) +
    facet_wrap(. ~ channel) +
    ggtitle(paste(chan, grp_no, sep = "-"))

  ggplot(data = orig_rec[t %between% c(1420,1420.5)], aes(x = t)) +
    geom_line(aes(y = bpfiltered, colour = "bp_filtered"), size = 0.5) +
    geom_line(aes(y = enved, colour = "envelope")) +
    geom_line(aes(y = 2*envsd, colour = "threshold")) +
    facet_wrap(. ~ channel) +
    ggtitle(paste(chan, grp_no, sep = "-"))
