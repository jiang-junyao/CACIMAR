
#' Calculate the joint power
#'
#' @param MarkerRoc01 MarkerRoc01
#'
#' @examples
Cal_JointPower2 <- function(MarkerRoc01) {


  ### revise power to negative if avg_diff < 0
  MarkerRoc02 <- apply(MarkerRoc01, 1, function(x01) {
    for (i0 in seq(2, length(x01), 3)) {
      if (x01[i0] < 0) {
        x01[i0 + 1] <- -x01[i0 + 1]
      }
    }
    return(x01)
  })

  MarkerRoc021 <- as.matrix(MarkerRoc02[seq(3, dim(MarkerRoc02)[1], 3), ])
  if (dim(MarkerRoc02)[1] == 1) {
    MarkerRoc021 <- t(MarkerRoc021)
  }

  ### calculate joint power for each cluster
  MarkerRoc03 <- apply(MarkerRoc021, 1, function(x01) {
    x02 <- x01[x01 > 0]
    x03 <- x01[x01 < 0]
    if (length(x02) > 0) {
      power01 <- Cal_JointPower(x02)
    } else {
      power01 <- 0
    }
    if (length(x03) > 0) {
      power02 <- Cal_JointPower(-x03)
    } else {
      power02 <- 0
    }
    power03 <- power01 - power02
    return(power03)
  })

  return(MarkerRoc03)
}
