
#' Calculate the joint power
#'
#' @param power01 power01
#'
#' @return ???
#'
#'
#' @examples
Cal_JointPower <- function(power01) {
  power02 <- 1 - power01[1]

  if (length(power01) < 2) {
  } else {
    for (i in 2:length(power01)) {
      power02 <- power02 * (1 - power01[i])
    }
  }

  return(1 - power02)
}
