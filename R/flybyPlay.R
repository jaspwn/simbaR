#' @import ggplot2
#' @import tuneR
#' @import data.table
#' @export

flybyPlayR <- function(grp_no, chan, orig_rec, filtered = TRUE) {

  event_start <- min(orig_rec[channel == chan][group_no == grp_no, t]) - 0.1
  event_end <- max(orig_rec[channel == chan][group_no == grp_no, t]) + 0.1

  plot <- ggplot(data = orig_rec[t %between% c(event_start, event_end)], aes(x = t)) +
    geom_rect(inherit.aes = FALSE, aes(xmin = min(orig_rec[channel == chan][group_no == grp_no, t]),
                                       xmax = max(orig_rec[channel == chan][group_no == grp_no, t]),
                                       ymin = -Inf, ymax = Inf, fill = "group"), fill = "#ffeda0") +
    geom_line(aes(y = DC, colour = "DC removed"), size = 0.2) +
    #geom_line(aes(y = smooth2, colour = "threshold_2std")) +
    geom_line(aes(y = bpfiltered, colour = "bp_filtered"), size = 0.5) +
    #geom_line(aes(y = 0.2, colour = "0.2 level")) +
    facet_wrap(. ~ channel) +
    ggtitle(paste(chan, grp_no, sep = "-"))

  print(plot)

  play_start <- min(orig_rec[channel == chan][group_no == grp_no, t]) - 1
  play_end <- max(orig_rec[channel == chan][group_no == grp_no, t]) + 1

  ## play filtered data by default, if not play raw
  if (filtered == TRUE) {
    u <- orig_rec[channel == chan][t %between% c(play_start, play_end), bpfiltered]
  } else {
    u <- orig_rec[channel == chan][t %between% c(play_start, play_end), value * 50]
  }

  uscaled <- scales::rescale(u, to = c(-7500000,7500000))

  w  <-  Wave(uscaled, samp.rate = 44100, bit = 24) #make the wave variable
  x <- stereo(w,w)
  play(x, player = "/usr/bin/vlc", "--play-and-exit --audio-visual visual")

  #writeWave(w, "femaleflyby_filtered.wav")


}



