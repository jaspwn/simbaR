#' @import tuneR
#' @import data.table
#' @export

## takes single channel exported .txt files from spike2 and converts them to wav file. volfac range 100-500

txt2wavR <- function(samprate, volfac) {
  filename <- file.choose()

  chan <- fread(filename, skip = 18, col.names = "v", blank.lines.skip = TRUE)

  t <- seq.int(from = 0, to = (length(chan$v))/samprate, length.out = length(chan$v))

  ## create datatable of target channel and t vector
  dt <- data.table(t, chan*100000)
  dt[, scaled := v*volfac]

  w  <-  Wave(as.integer(dt$scaled), samp.rate = samprate, bit = 24) #make the wave variable
  x <- stereo(w,w)
  play(x, player = "/usr/bin/vlc", "--play-and-exit --audio-visual visual")

  writeWave(x, paste0(tools::file_path_sans_ext(filename), ".wav"))

  message(paste0("wav ouput to ", tools::file_path_sans_ext(filename), ".wav"))

}
