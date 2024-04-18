Cal_MarkersRoc <- function(pbmc01, marker01) {
  cluster01 <- levels(pbmc01@active.ident)
  for (i in 1:length(cluster01)) {
    print(cluster01[i])
    markerRoc01 <- Markertest(pbmc01, cells.1 = Seurat::WhichCells(pbmc01, idents = cluster01[i]),
                              cells.2 = Seurat::WhichCells(pbmc01, idents = cluster01[i], invert = TRUE), genes.use = marker01)
    if (i == 1) {
      markerRoc02 <- markerRoc01[order(rownames(markerRoc01)), ]
    } else {
      markerRoc02 <- cbind(markerRoc02, markerRoc01[order(rownames(markerRoc01)), ])
    }
    colnames(markerRoc02)[(ncol(markerRoc02)-2):ncol(markerRoc02)] <- paste0(paste(paste0("", cluster01[i]), collapse = ""),
                                                                             "_", colnames(markerRoc02)[(ncol(markerRoc02)-2):ncol(markerRoc02)])
  }

  return(markerRoc02)
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

  ### calculate the joint power for each cluster
  MarkerRoc03 <- apply(MarkerRoc021, 1, function(x01) {
    x02 <- x01[x01 > 0]; x03 <- x01[x01 < 0]
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


AUCMarkerTest <- function(data1, data2, mygenes, print.bar = TRUE) {
  myAUC <- unlist(x = lapply(X = mygenes,
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


Markertest <- function(object, cells.1, cells.2, genes.use = NULL, print.bar = TRUE,
                       assay.type = "RNA") {
  data.test <- Seurat::GetAssayData(
    object = object,
    slot = "data"
  )
  if (is.null(genes.use)) {
    genes.use <- rownames(x = data.test)
  }
  to.return <- AUCMarkerTest(
    data1 = data.test[, cells.1],
    data2 = data.test[, cells.2], mygenes = genes.use, print.bar = print.bar
  )
  to.return$power <- abs(x = to.return$myAUC - 0.5) * 2
  return(to.return)
}

