#' @import ggplot2
#' @import tuneR
#' @import data.table
#' @import rhdf5
#' @export

flybyPlayR <- function(resdt, row, filtered = FALSE, printPlot = TRUE, filext = ".mat") {

  #bandpass filter design
  bpfilt <- butter(n = 4, W = c(600/(50000/2), 1200/(50000/2)), type = "pass", plane = "z")

  #data <- resReadR(resfile)

  flyby <- resdt[row]

  filename <- paste0(gsub("_out.txt", "", flyby[, filename]), filext)

  ## read in file channel information

  meta <- hdf5metaLoadR(filename)
  #h5closeAll()

  samprate <- meta$samprate
  chandetails <- meta$chandetails
  codedt <- meta$codedt
  idx_lists <- meta$idx_lists

  ## find where flyby occurs

  #find indexes of raw data file corresponding to lettercode of flyby
  #idx <- unlist(idx_lists[which(lapply(idx_lists, function(x) grep(codedt[lettercode == flyby[,code], t_idx], x)) == 1)])
  #idx <- unlist(idx_lists[which.min(abs(codedt[, time] - flyby[, min_t]))])
  idx <- unlist(idx_lists[which.max(codedt[, time][codedt[, time]<=flyby[, min_t]])])

  if(file.exists(paste0(paste0(dirname(filename), "/raw_events/", gsub(".mat", paste0("_", idx[1], ".txt"), basename(filename)))))) {

    moltendt <- fread(paste0(paste0(dirname(filename), "/raw_events/", gsub(".mat", paste0("_", idx[1], ".txt"), basename(filename)))))

  } else {

    # create data.table out of data channels
    if(!exists("moltendt", envir = .GlobalEnv) || flyby[, time] != moltendt[1, time]) {

      chan01 <- as.vector(h5read(filename,
                                 paste0(chandetails[chan == "Ch1", name], "/values/"),
                                 index = list(idx[1]:idx[2],1)))
      h5closeAll()
      chan02 <- as.vector(h5read(filename,
                                 paste0(chandetails[chan == "Ch2", name], "/values/"),
                                 index = list(idx[1]:idx[2],1)))
      h5closeAll()
      chan03 <- as.vector(h5read(filename,
                                 paste0(chandetails[chan == "Ch3", name], "/values/"),
                                 index = list(idx[1]:idx[2],1)))
      h5closeAll()
      chan04 <- as.vector(h5read(filename,
                                 paste0(chandetails[chan == "Ch4", name], "/values/"),
                                 index = list(idx[1]:idx[2],1)))
      h5closeAll()

      t_seg <- seq.int(from = codedt[t_idx == idx[1], time],
                       to = codedt[t_idx == idx[1], time] + length(chan01)*samprate,
                       length.out = length(chan01))

      dt <- data.table(t = t_seg,
                       chan01 = chan01,
                       chan02 = chan02,
                       chan03 = chan03,
                       chan04 = chan04)
      rm(t_seg, chan01, chan02, chan03, chan04)

      moltendt <<- melt(dt,

                        measure.vars = c("chan01",
                                         "chan02",
                                         "chan03",
                                         "chan04"),
                        variable.name = "channel")

      ## rescale each channel
      moltendt[, code := codedt[t_idx == idx[1], lettercode]]
      moltendt[, time := codedt[t_idx == idx[1], time]]
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

    }

  }






  ## print plot

  if (printPlot == TRUE) {

    ## select data for plotting and playing
    event_start <- flyby[, min_t]
    event_end <- flyby[, max_t]


    plot <- ggplot(data = moltendt[t %between% c(event_start, event_end)], aes(x = t)) +
      geom_rect(inherit.aes = FALSE, aes(xmin = flyby[, min_t],
                                         xmax = flyby[, max_t],
                                         ymin = -Inf, ymax = Inf, fill = "group"), fill = "#ffeda0") +
      geom_line(aes(y = DC, colour = "DC removed"), size = 0.2) +
      geom_line(aes(y = bpfiltered, colour = "bp_filtered"), size = 0.5) +
      geom_line(aes(y = enved, colour = "envelope"), size = 1) +
      scale_color_brewer(type = "qual", palette = "Set1") +
      facet_wrap(. ~ channel) +
      cowplot::theme_minimal_grid() +
      theme(legend.position = "bottom")

    flytable <- gridExtra::tableGrob(flyby[,c(1,2,5,6,8,9,10,11,12,13)], theme = gridExtra::ttheme_minimal())

    comboplot <- cowplot::plot_grid(flytable, plot, ncol = 1, rel_heights = c(1,4), axis = "l", align = "h")

    print(comboplot)

  }

  play_start <- flyby[, min_t]
  play_end <- flyby[, max_t]

  ## play filtered data by default, if not play raw
  if (filtered == TRUE) {
    u <- moltendt[channel == flyby[, channel]][t %between% c(play_start, play_end), bpfiltered]
  } else {
    u <- moltendt[channel == flyby[, channel]][t %between% c(play_start, play_end), value]
  }

  uscaled <- round(scales::rescale(c(noise01[channel == flyby[, channel], value],
                                     u,
                                     noise01[channel == flyby[, channel], value]),
                                   to = c(-7500000,7500000)))

  #uscaled <- round(scales::rescale(u, to = c(-7500000,7500000)))

  w  <-  Wave(uscaled, samp.rate = 1/samprate, bit = 24) #make the wave variable
  x <- stereo(w,w)
  play(x, player = "/usr/bin/vlc", "--play-and-stop --audio-visual visual")

  #return(moltendt)

}
