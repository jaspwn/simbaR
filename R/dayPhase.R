#' @import data.table
#' @import behavr
#' @export

## Creates swarm window by setting dayPhase to peaktime +/- sWidth for each expDay (experimental day)

dayPhaseR <- function(x, y, sWidth = 30) {

  if(x %between% c(peak_dt[expDay == y, peaktime] - behavr::mins(sWidth),
                   peak_dt[expDay == y, peaktime] + behavr::mins(sWidth))) {

    res <- paste0(y, " swarm")

    return(res)

  } else {

    res <- paste0(y, " other")

    return(res)

  }


}
