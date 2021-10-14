

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

Cal_MarkersRoc <- function(pbmc01, marker01) {
  clus01 <- levels(pbmc01@active.ident)
  for (i in 1:length(clus01)) {
    print(clus01[i])
    markerRoc01 <- Markertest(pbmc01, cells.1 = Seurat::WhichCells(pbmc01, idents = clus01[i]), cells.2 = Seurat::WhichCells(pbmc01, idents = clus01[i], invert = TRUE), genes.use = marker01)
    if (i == 1) {
      markerRoc02 <- markerRoc01[order(rownames(markerRoc01)), ]
    } else {
      markerRoc02 <- cbind(markerRoc02, markerRoc01[order(rownames(markerRoc01)), ])
    }
    colnames(markerRoc02)[(ncol(markerRoc02) - 2):ncol(markerRoc02)] <- paste0(paste(paste0("C", clus01[i]), collapse = ""), "_", colnames(markerRoc02)[(dim(markerRoc02)[2] - 2):dim(markerRoc02)[2]])
  }

  return(markerRoc02)
}

AUCMarkerTest <- function(data1, data2, mygenes, print.bar = TRUE) {
  myAUC <- unlist(x = lapply(
    X = mygenes,
    FUN = function(x) {
      return(DifferentialAUC(
        x = as.numeric(x = data1[x, ]),
        y = as.numeric(x = data2[x, ])
      ))
    }
  ))
  myAUC[is.na(x = myAUC)] <- 0
  iterate.fxn <- ifelse(test = print.bar, yes = pbapply::pblapply, no = lapply)
  avg_diff <- unlist(x = iterate.fxn(
    X = mygenes,
    FUN = function(x) {
      return(
        Seurat::ExpMean(
          x = as.numeric(x = data1[x, ])
        ) - Seurat::ExpMean(
          x = as.numeric(x = data2[x, ])
        )
      )
    }
  ))
  toRet <- data.frame(cbind(myAUC, avg_diff), row.names = mygenes)
  toRet <- toRet[rev(x = order(toRet$myAUC)), ]
  return(toRet)
}

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
