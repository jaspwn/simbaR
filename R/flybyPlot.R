#' @import ggplot2
#' @import data.table
#' @export

flybyPlotR <- function(grp_no, chan, orig_rec, plot_thresh = FALSE, t_interval) {

  if (missing(t_interval)) {

    event_start <- min(orig_rec[channel == chan][group_no == grp_no, t]) - 0.1
    event_end <- max(orig_rec[channel == chan][group_no == grp_no, t]) + 0.1

    if (plot_thresh == TRUE) {

      plot <- ggplot(data = orig_rec[t %between% c(event_start, event_end) &
                                       channel == chan], aes(x = t)) +
        geom_rect(inherit.aes = FALSE, aes(xmin = min(orig_rec[channel == chan][group_no == grp_no, t]),
                                           xmax = max(orig_rec[channel == chan][group_no == grp_no, t]),
                                           ymin = -Inf, ymax = Inf, fill = "group"), fill = "#ffeda0") +
        geom_line(aes(y = DC, colour = "DC removed"), size = 0.2) +
        geom_line(aes(y = bpfiltered, colour = "bp_filtered"), size = 0.5) +
        geom_line(aes(y = enved, colour = "envelope")) +
        geom_line(aes(y = 2*envsd, colour = "threshold")) +
        ggtitle(paste(chan, grp_no, sep = "-"))

      print(plot)

    } else {

      plot <- ggplot(data = orig_rec[t %between% c(event_start, event_end) &
                                       channel == chan], aes(x = t)) +
        geom_rect(inherit.aes = FALSE, aes(xmin = min(orig_rec[channel == chan][group_no == grp_no, t]),
                                           xmax = max(orig_rec[channel == chan][group_no == grp_no, t]),
                                           ymin = -Inf, ymax = Inf, fill = "group"), fill = "#ffeda0") +
        geom_line(aes(y = scales::rescale(DC, to = c(-1,1)), colour = "DC removed"), size = 0.2) +
        geom_line(aes(y = scales::rescale(bpfiltered, to = c(-1,1)), colour = "bp_filtered"), size = 0.5) +
        ggtitle(paste(chan, grp_no, sep = "-"))

      print(plot)

    }

  }

  if (missing(grp_no)) {

    event_start <- t_interval[1]
    event_end <- t_interval[2]

    grp_no = paste0("Time interval ", event_start, "-", event_end)

    if (plot_thresh == TRUE) {

      plot <- ggplot(data = orig_rec[t %between% c(event_start, event_end) &
                                       channel == chan], aes(x = t)) +
        geom_line(aes(y = bpfiltered, colour = "bp_filtered"), size = 0.5) +
        geom_line(aes(y = enved, colour = "envelope")) +
        geom_line(aes(y = 2*envsd, colour = "threshold")) +
        ggtitle(paste(chan, grp_no, sep = "-"))

      print(plot)

    } else {

      plot <- ggplot(data = orig_rec[t %between% c(event_start, event_end) &
                                       channel == chan], aes(x = t)) +
        geom_line(aes(y = scales::rescale(DC, to = c(-1,1)), colour = "DC removed"), size = 0.2) +
        geom_line(aes(y = scales::rescale(bpfiltered, to = c(-1,1)), colour = "bp_filtered"), size = 0.5) +
        ggtitle(paste(chan, grp_no, sep = "-"))

      print(plot)

    }

  }

  # event_start <- min(orig_rec[channel == chan][group_no == grp_no, t]) - 0.1
  # event_end <- max(orig_rec[channel == chan][group_no == grp_no, t]) + 0.1



}



