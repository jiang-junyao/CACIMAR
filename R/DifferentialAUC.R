#' Title internal from seurat
#'
#' @param x x
#' @param y y
#'
#'
#'
#'
#'
DifferentialAUC <- function(x, y) {
  prediction.use <- ROCR::prediction(
    predictions = c(x, y),
    labels = c(rep(x = 1, length(x = x)), rep(x = 0, length(x = y))),
    label.ordering = 0:1
  )
  perf.use <- ROCR::performance(prediction.obj = prediction.use, measure = "auc")
  auc.use <- round(x = perf.use@y.values[[1]], digits = 3)
  return(auc.use)
}
