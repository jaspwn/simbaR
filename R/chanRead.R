#' @export
#' @import rhdf5
#' @import data.table

chanReadeR <- function(filename = filename,
                       chandetails = chandetails,
                       idx_lists = idx_list,
                       codedt = codedt,
                       samprate = samprate,
                       channum) {

  chandata <- as.vector(h5read(filename,
                               paste0(chandetails[chan == channum, name], "/values/"),
                               index = list(idx_lists[1]:idx_lists[2],1)))
  h5closeAll()

  t_seg <- seq.int(from = codedt[t_idx == idx_lists[1], time],
                   to = codedt[t_idx == idx_lists[1], time] + length(chandata)*samprate,
                   length.out = length(chandata))

  dt <- data.table(t = t_seg, value = chandata)
  dt[, channel := channum]
  dt[, code := codedt[t_idx == idx_lists[1], lettercode]]
  dt[, time := codedt[t_idx == idx_lists[1], time]]

  return(dt)

}
