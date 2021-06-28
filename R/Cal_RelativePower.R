###

#' Calculate the relative power
#'
#' @param power01 power
#' @importFrom psych geometric.mean
#'
#'
#' @examples
Cal_RelativePower <- function(power01) {
  power02 <- c()
  if (length(power01) < 2) {
    power02 <- power01
  } else {
    for (i in 1:(length(power01) - 1)) {
      power01g <- psych::geometric.mean(power01[1:i])
      if ((power01g - power01[i + 1]) < 0 & (power01g - power01[i + 1]) > -10e-10) {
        power02[i] <- 0
      } else {
        power02[i] <- (power01g * (power01g - power01[i + 1]))^0.5
      }
    }
  }

  return(power02)
}
