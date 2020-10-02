#' @import ggplot2
#' @import tuneR
#' @import data.table
#' @import rhdf5
#' @export

wavPlayR <- function(resdt, row, cfreqs = c(600, 1200), filtered = FALSE, printPlot = TRUE, filext = ".mat") {


  filename <- paste0(gsub("_out.txt", "", resdt[row, filename]), filext)
  names <- data.table(h5ls(filename))

  ## read in file channel information
  channames <- data.table(name = names[group == "/", name])
  samprate <<- h5read(filename, paste0(channames[1, name], "/interval/"))[1,]
  stime <- h5read(filename, paste0(channames[1, name], "/start/"))[,1]
  len <- h5read(filename, paste0(channames[1, name], "/length/"))[,1]
  channames[, c("date", "recID",  "treatment", "chan") := tstrsplit(name, "_")]
  chandetails <- channames[complete.cases(channames)]

  bpfilt <- butter(n = 4, W = c(cfreqs[1]*samprate*2, cfreqs[2]*samprate*2), type = "pass", plane = "z")

  if(file.exists(paste0(paste0(dirname(filename), "/raw_events/", gsub(".mat", ".txt", basename(filename)))))) {

    moltendt <- fread(paste0(paste0(dirname(filename), "/raw_events/", gsub(".mat", ".txt", basename(filename)))))

  } else {

    # create data.table out of data channels
    if(!exists("moltendt", envir = .GlobalEnv) || resdt[row, time] != moltendt[1, time]) {

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

    }

  }






  ## print plot

  if (printPlot == TRUE) {

    ## select data for plotting and playing
    event_start <- resdt[row, min_t]
    event_end <- resdt[row, max_t]


    plot <- ggplot(data = moltendt[t %between% c(event_start, event_end) &
                                     channel == resdt[, channel]], aes(x = t)) +
      geom_rect(inherit.aes = FALSE, aes(xmin = resdt[row, min_t],
                                         xmax = resdt[row, max_t],
                                         ymin = -Inf, ymax = Inf, fill = "group"), fill = "#ffeda0") +
      geom_line(aes(y = DC, colour = "DC removed"), size = 0.2) +
      geom_line(aes(y = bpfiltered, colour = "bp_filtered"), size = 0.5) +
      geom_line(aes(y = enved, colour = "envelope"), size = 1) +
      scale_color_brewer(type = "qual", palette = "Set1") +
      cowplot::theme_minimal_grid() +
      theme(legend.position = "bottom")

    flytable <- gridExtra::tableGrob(resdt[row,c(1,2,5,6,8,9,10,11,12,13)], theme = gridExtra::ttheme_minimal())

    comboplot <- cowplot::plot_grid(flytable, plot, ncol = 1, rel_heights = c(1,4), axis = "l", align = "h")

    print(comboplot)

  }

  play_start <- resdt[row, min_t]
  play_end <- resdt[row, max_t]

  ## play filtered data by default, if not play raw
  if (filtered == TRUE) {
    u <- moltendt[channel == resdt[row, channel]][t %between% c(play_start, play_end), bpfiltered]
  } else {
    u <- moltendt[channel == resdt[row, channel]][t %between% c(play_start, play_end), value]
  }

  uscaled <- round(scales::rescale(c(noise01[channel == resdt[row, channel], value],
                                     u,
                                     noise01[channel == resdt[row, channel], value]),
                                   to = c(-7500000,7500000)))

  #uscaled <- round(scales::rescale(u, to = c(-7500000,7500000)))

  w  <-  Wave(uscaled, samp.rate = 1/samprate, bit = 24) #make the wave variable
  x <- stereo(w,w)
  play(x, player = "/usr/bin/vlc", "--play-and-exit --audio-visual visual")

  #return(moltendt)

}
