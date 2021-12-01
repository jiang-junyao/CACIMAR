Sort_MarkersPower <- function(MarkerRoc01, Thr01 = 0) {
  MarkerRoc02 <- sort.int(as.numeric(MarkerRoc01), index.return = T, decreasing = T)
  Name01 <- colnames(MarkerRoc01)
  if (length(MarkerRoc01) == 0 | max(MarkerRoc02$x) <= Thr01) {
    MarkerRoc03 <- c("None", "None", "None")
  } else if (length(MarkerRoc02$x) == 1 & length(MarkerRoc02$x[MarkerRoc02$x > Thr01]) == 1 |
             length(MarkerRoc02$x) != 1 & length(MarkerRoc02$x[MarkerRoc02$x > Thr01]) == 1 & MarkerRoc02$x[2] < 0) {
    MarkerRoc03 <- c(paste(Name01[MarkerRoc02$ix[1]], sep = ","),
                     paste(MarkerRoc02$x[1], sep = ","), Name01[MarkerRoc02$ix[1]])
  } else {
    if (length(MarkerRoc02$x) != 1 & length(MarkerRoc02$x[MarkerRoc02$x > Thr01]) == 1) {
      MarkerRoc02x <- MarkerRoc02$x[1:2]
      MarkerRoc02ix <- MarkerRoc02$ix[1:2]
    } else {
      MarkerRoc02x <- c(MarkerRoc02$x[MarkerRoc02$x > Thr01], Thr01)
      MarkerRoc02ix <- c(MarkerRoc02$ix[MarkerRoc02$x > Thr01], length(MarkerRoc02$ix) + 1)
    }
    Power01 <- Cal_RelativePower(MarkerRoc02x)
    maxPower01 <- which.max(Power01)
    if (maxPower01 == length(MarkerRoc02x)) {
      maxPower01 <- length(Power01) - 1
    }
    MarkerRoc03 <- c(paste(Name01[MarkerRoc02ix[1:(maxPower01 + 1)]], collapse = ","),
                     paste(MarkerRoc02x[1:(maxPower01 + 1)], collapse = ","), Name01[MarkerRoc02ix[1]])
  }

  return(MarkerRoc03)
}
