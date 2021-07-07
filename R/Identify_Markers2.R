
#' Inner function for indentifying markers
#'

Identify_Markers2 <- function(pbmc, Marker, GeneSymb1=NULL, PowerThr1=1/3){

  TotoalCluster <- length(unique(pbmc@active.ident))
  MarkerRoc1 <- Marker
  NumCluster <- apply(MarkerRoc1,1,function(X1){length(strsplit(as.character(X1[1]),',')[[1]])})
  MarkerRoc2 <- cbind(MarkerRoc1,NumCluster)
  MarkerRoc3 <- MarkerRoc2[MarkerRoc2$Diff>PowerThr1 & MarkerRoc2$NumCluster>1 & MarkerRoc2$NumCluster<=(ceiling(TotoalCluster/4)+1),]

  Cluster1 <- strsplit(levels(as.factor(MarkerRoc3$Cluster)),',')
  Cluster3 <- c()
  for(i in 1:length(Cluster1)){
    Cluster2 <- Cluster1[[i]]
    Cluster3 <- c(Cluster3,Cluster2[1:(length(Cluster2)-1)])
  }
  uCluster3 <- unique(Cluster3)

  Num1 <- c(); MarkerRoc5 <- c()
  for(i in 1:length(uCluster3)){ print(uCluster3[i])
    Cluster4 <- gsub('C','',uCluster3[i])
    MarkerRoc4 <- MarkerRoc3[grep(paste0(uCluster3[i],','),MarkerRoc3$Cluster),]
    Num1 <- rbind(Num1,c(uCluster3[i],dim(MarkerRoc4)[1]))
    if(dim(MarkerRoc4)[1]>0){
      MarkerRoc5 <- rbind(MarkerRoc5,cbind(MarkerRoc4,Cluster4,rownames(MarkerRoc4)))
    }
  }
  colnames(MarkerRoc5)[c(1,dim(MarkerRoc5)[2]-1,dim(MarkerRoc5)[2])] <- c('AllCluster','cluster','gene')
  colnames(Num1) <- c('Cluster','NumMarkers')
  Num1 <- Num1[order(Num1[,'Cluster']),]
  print(Num1)

  MarkerRoc51 <- apply(MarkerRoc5,1,function(x1){ x11 <- strsplit(as.character(x1[1]),',')[[1]]; x12 <- paste(x11[1:(length(x11)-1)],collapse=',')
  x21 <- strsplit(as.character(x1[2]),',')[[1]]; x22 <- paste(x21[1:(length(x21)-1)],collapse=',')
  return(c(x12,x22))
  })
  MarkerRoc52 <- cbind(t(MarkerRoc51),MarkerRoc5[,5:6])
  colnames(MarkerRoc52)[1:2] <- colnames(MarkerRoc5)[1:2]

  MarkerRoc53 <- cbind(MarkerRoc52, GeneSymb1[match(MarkerRoc52$gene,GeneSymb1[,1]),4])
  colnames(MarkerRoc53)[ncol(MarkerRoc53)] <- 'Symbol'

  return(MarkerRoc53)
}
