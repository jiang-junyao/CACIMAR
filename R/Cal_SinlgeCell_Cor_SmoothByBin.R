Cal_SinlgeCell_Cor_SmoothByBin<-function(){
  source('/users/jwang3/RetReg/Programs/Seurat_SubsetData.R')
  source('/users/jwang3/RetReg/Programs/SmoothByBin.R')

  TotalBin1 <- 50; NumRand1 <- 5; FixPseudotimeRange1 <- T; CtrBin1 <- 0
  ByBin1 <- c('Equal.Pseudotime'); SubG1 <- 'Treatment';
  if(Spec1=='Mm'){ Ctr1 <- 'mmCtrP'; BranchState1 <- NULL
  Group1 <- list(c('mmNMDA'), c('mmLD'))
  }else if(Spec1=='Zf'){ Ctr1 <- 'zfAd'; BranchState1 <- NULL
  Group1 <- list(c('zfNMDA'), c('zfLD'), c('zfTR'))
  }else if(Spec1=='Ch'){ Ctr1 <- 'chP';
  Group1 <- list(c('chNMDA'), c('chFI','chNMFI')); BranchState1 <- NULL
  }

  if(!is.null(BranchState1)){ Part1 <- 2 }else{ Part1 <- 1 }
  Meta1 <- pbmc@meta.data; NumGroup1 <- combn(length(Group1), 2)
  for(i in 1:ncol(NumGroup1)){
    Group2 <- Group1[NumGroup1[,i]]; print(Group2)

    ###separate two parts according to BranchState1
    for(i1 in 1:Part1){
      if(i1==1){
        State1 <- unique(Meta1$State); OtherS1 <- setdiff(State1, BranchState1)
        if(!is.null(BranchState1)){
          pbmc1 <- Seurat_SubsetData(pbmc, SubG1='State', SubS1=OtherS1)
        }else{ pbmc1 <- pbmc }
      }else{
        pbmc1 <- Seurat_SubsetData(pbmc, SubG1='State', SubS1=BranchState1)
      }
      Meta1 <- pbmc1@meta.data

      ###Identify minimal and maximal pseudotime across the models
      Meta12 <- list(); Pseudotime1 <- c();
      for(j in 1:length(Group2)){
        SubS1 <- c(Ctr1, Group2[[j]]); Meta11 <- c()
        for(k in 1:length(SubS1)){
          Meta11 <- rbind(Meta11, Meta1[Meta1[, SubG1]==SubS1[k], ])
        }
        Pseudotime1 <- rbind(Pseudotime1, c(min(Meta11$Pseudotime), max(Meta11$Pseudotime)))
        Meta12[[j]] <- Meta11
      }
      Pseudotime12 <- c(max(Pseudotime1[,1]), min(Pseudotime1[,2]))

      Exp1 <- list(); Group3 <- c();
      for(j in 1:length(Group2)){ Group3 <- c(Group3, Group2[[j]][1])
      SubS1 <- c(Ctr1, Group2[[j]]);
      pbmc2 <- Seurat_SubsetData(pbmc1, SubG1=SubG1, SubS1=SubS1)
      print(table(pbmc2@meta.data$Condition))
      if(FixPseudotimeRange1==T){
        Exp1[[j]] <- SmoothByBin(as.matrix(pbmc2@assays$RNA@data), pbmc2@meta.data, PseudotimeRange1=Pseudotime12, SmoothLength1=TotalBin1, ByBin1=ByBin1[1])
      }else{ Exp1[[j]] <- SmoothByBin(as.matrix(pbmc2@assays$RNA@data), pbmc2@meta.data, SmoothLength1=TotalBin1, ByBin1=ByBin1[1]) }
      if(i1==1 & j==1){ Exp2 <- Exp1[[j]]
      }else{ Exp2 <- cbind(Exp2, Exp1[[j]]) }
      }
    }

    RandCor1 <- c()
    for(j in 1:(NumRand1+1)){ print(j)
      if(j==1){ Exp3 <- Exp2
      }else{ Ind1 <- c(sample(TotalBin1, TotalBin1), sample((TotalBin1+1):(TotalBin1*2), TotalBin1))
      Exp3 <- Exp2[, Ind1]
      }
      Cor1 <- t(apply(Exp3, 1, function(x1){
        x11 <- rbind(x1[1:TotalBin1], x1[(TotalBin1+1):(TotalBin1*2)])
        sx11 <- colSums(x11)
        if(length(sx11[sx11!=0]) < ceiling(TotalBin1/10)){ x12 <- x11
        }else{ x12 <- x11[, colSums(x11)!=0] }
        cx12 <- cor.test(x12[1, ], x12[2, ])
        if(is.na(cx12$est) & is.na(cx12$p.val) ){ cx12$est <- 0; cx12$p.val <- 1 }
        return(c(cx12$est, cx12$p.val))
      }))
      if(j==1){ Cor12 <- Cor1
      }else{
        RandCor1 <- rbind(RandCor1, t(Cor1[, 1]))
      }
    }

    colnames(Cor12) <- paste0(c('Cor_','Pcor_'), paste(Group3, collapse='_'))
    if(NumRand1>1){ RandCor12 <- as.matrix(colMeans(RandCor1))
    }else{ RandCor12 <- t(RandCor1) }

    if(i==1){ Cor2 <- cbind(Cor12, RandCor12)
    }else{ Cor2 <- cbind(Cor2, Cor12, RandCor12) }
    colnames(Cor2)[ncol(Cor2)] <- paste0('RandomCor_',  paste(Group3, collapse='_'))
  }
  Gene1 <- Converse_GeneIDSymbol(rownames(Cor2), Symb1[[1]])
  Cor3 <- cbind(Gene1, Cor2)
  write.table(Cor3, paste0(File1,'_Bin',TotalBin1,'_R',NumRand1,'_GeneCor.txt'), quote=F, sep='\t', row.names=F)
}

#' Internal function
smoothByBin <- function(Exp1, Pseudotime1, PseudotimeRange1 = NULL, SmoothLength1 = 100, ByBin1 = c("Equal.Pseudotime", "Equal.Cells")) {
  if (ncol(Exp1) != nrow(Pseudotime1)) {
    stop("The length of pseudotime is not equal to the number of cells")
  } else if (!is.element("Pseudotime", colnames(Pseudotime1))) {
    stop("No pseudotime inforamtion in variable Pseudotime1")
  }

  if (ByBin1[1] == "Equal.Pseudotime") {
    if (is.null(PseudotimeRange1)) {
      PseudotimeBin1 <- seq(min(Pseudotime1$Pseudotime), max(Pseudotime1$Pseudotime), length.out = SmoothLength1 + 1)
    } else {
      PseudotimeBin1 <- seq(PseudotimeRange1[1], PseudotimeRange1[2], length.out = SmoothLength1 + 1)
    }
  } else {
    Pseudotime1 <- Pseudotime1[order(Pseudotime1$Pseudotime), ]
    Bin1 <- ceiling(nrow(Pseudotime1) / SmoothLength1)
    PseudotimeBin1 <- seq(1, nrow(Pseudotime1), by = Bin1)
    PseudotimeBin1[length(PseudotimeBin1) + 1] <- nrow(Pseudotime1)
  }

  Exp2 <- array(0, dim = c(nrow(Exp1), SmoothLength1))
  for (i in 1:(length(PseudotimeBin1) - 1)) {
    if (ByBin1[1] == "Equal.Pseudotime") {
      if (i == length(PseudotimeBin1) - 1) {
        Cells1 <- Pseudotime1[Pseudotime1$Pseudotime >= PseudotimeBin1[i] & Pseudotime1$Pseudotime <= PseudotimeBin1[i + 1], ]
      } else {
        Cells1 <- Pseudotime1[Pseudotime1$Pseudotime >= PseudotimeBin1[i] & Pseudotime1$Pseudotime < PseudotimeBin1[i + 1], ]
      }
    } else {
      if (i == length(PseudotimeBin1) - 1) {
        Cells1 <- Pseudotime1[PseudotimeBin1[i]:PseudotimeBin1[i + 1], ]
      } else {
        Cells1 <- Pseudotime1[PseudotimeBin1[i]:(PseudotimeBin1[i + 1] - 1), ]
      }
    }

    if (nrow(Cells1) > 1) {
      Exp2[, i] <- rowMeans(Exp1[, match(rownames(Cells1), colnames(Exp1))])
    } else if (nrow(Cells1) == 1) {
      Exp2[, i] <- Exp1[, match(rownames(Cells1), colnames(Exp1))]
    }
  }
  rownames(Exp2) <- rownames(Exp1)

  return(Exp2)
}
