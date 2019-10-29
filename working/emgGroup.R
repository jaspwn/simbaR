#' @import biosignalEMG


##

emgGroupR <- function(my_vector, samprate = samprate, samplingrate = 1/samprate, window = 1/samprate/100) {

  #perform full wave rectification
  emg_vector <- emg(my_vector, samplingrate = samprate)
  rect_vector <- rectification(emg_vector, rtype = "fullwave")
  #calculate envelope
  MAlope <- envelope(rect_vector, method = "MA", wsize = window)
  #return envelope values
  return(MAlope$values)

}
