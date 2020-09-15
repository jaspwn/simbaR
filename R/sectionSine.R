#' @export
#' @import rhdf5
#' @import data.table


sectionSineR <- function(idx_lists, filename, channames, codedt, dataChan, dataSamprate, startFit, endFit) {

  assign("dataSamprate", dataSamprate, envir = .GlobalEnv)

  h5closeAll()

  data <- as.vector(h5read(filename,
                           paste0(channames[chan == dataChan, name], "/values/"),
                           index = list(idx_lists[1]:idx_lists[2],1)))
  h5closeAll()

  t_seg <- seq.int(from = codedt[start_idx == idx_lists[1], start_time],
                   to = codedt[start_idx == idx_lists[1], end_time],
                   length.out = length(data))

  dt <- data.table(t = t_seg,
                   channel = data)
  rm(t_seg, data)

  ## add time and code markers
  dt[, code := codedt[start_idx == idx_lists[1], lettercode]]
  dt[, time := codedt[start_idx == idx_lists[1], time]]

  print(paste0("idx number ", codedt[start_idx == idx_lists[1], which = TRUE]))

  ## perform DC remove of each channel

  #dt[, DC := removeDC(channel, 1000)]

  dt_sum <- dt[startFit:endFit, sineOptR(y = channel, t = t, samprate = dataSamprate), by = c("code", "time")]



  return(dt_sum)

}
