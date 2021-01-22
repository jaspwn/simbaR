#' @import data.table
#' @import behavr
#' @export

## Creates swarm window by setting dayPhase to peaktime +/- sWidth for each expDay (experimental day)
## x is time of event, y is the experimentdal day, and z (if provided) is the individual channel.

dayPhaseR <- function(x, y, z, sWidth = 30, channelWise = FALSE) {

  if(channelWise == FALSE) {

    if(x %between% c(peak_dt[expDay == y, peaktime] - behavr::mins(sWidth),
                     peak_dt[expDay == y, peaktime] + behavr::mins(sWidth))) {

      res <- paste0(y, " swarm")

      return(res)

    } else {

      res <- paste0(y, " other")

      return(res)

    }
  } else {

    if(x %between% c(peak_dt[expDay == y & channel == z, peaktime] - behavr::mins(sWidth),
                     peak_dt[expDay == y & channel == z, peaktime] + behavr::mins(sWidth))) {

      res <- paste0(y, " swarm")

      return(res)

    } else {

      res <- paste0(y, " other")

      return(res)

    }

  }

}
