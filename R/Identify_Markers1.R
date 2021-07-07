#' Inner function for markers identifying

Identify_Markers1<-function(Seurat_object, PowerThr1=1/3){
  Marker0 <- FindAllMarkers(object = Seurat_object, test.use ='roc', return.thresh = PowerThr1)
  Marker1 <- as.data.frame(Marker0)
  Marker2 <- Marker1[Marker1[,'avg_diff']>0 & Marker1[,'power']>PowerThr1,]
  uMarker2 <- unique(Marker2[,'gene'])
  MarkerRoc1 <- Cal_MarkersRoc(Seurat_object, uMarker2)
  MarkerRoc2 <- Select_MarkersPower(MarkerRoc1)
  return(MarkerRoc2)
}



#' Inner function for markers identifying
#'

Select_MarkersPower <- function(MarkerRoc01){

  MarkerRoc02 <- t(Revise_MarkersPower(MarkerRoc01))
  colnames(MarkerRoc02) <- gsub('_power','',colnames(MarkerRoc02))
  Power02 <- c()
  for(i in 1:dim(MarkerRoc02)[1]){
    MarkerRoc03 <- t(as.matrix(MarkerRoc02[i,seq(3,dim(MarkerRoc02)[2],3)]))
    Power01 <- Sort_MarkersPower(MarkerRoc03,0)
    Power02 <- rbind(Power02,Power01)
  }
  rownames(Power02) <- rownames(MarkerRoc02)
  colnames(Power02) <- c('Cluster','Power','Diff')

  return(Power02)
}


#' Inner function for markers identifying
Revise_MarkersPower <- function(MarkerRoc01){
  MarkerRoc02 <- apply(MarkerRoc01,1,function(x01){
    for(i0 in seq(2,length(x01),3)){
      if(x01[i0]<0){ x01[i0+1] <- -x01[i0+1] }
    }
    return(x01)
  } )

  return(MarkerRoc02)
}



#' Inner function for markers identifying

Sort_MarkersPower <- function(MarkerRoc01,Thr01=0){

  MarkerRoc02 <- sort.int(as.numeric(MarkerRoc01),index.return=T,decreasing=T)
  Name01 <- colnames(MarkerRoc01)
  if(length(MarkerRoc01)==0|max(MarkerRoc02$x)<=Thr01){
    MarkerRoc03 <- c('', '', '');
  }else if(length(MarkerRoc02$x)==1 & length(MarkerRoc02$x[MarkerRoc02$x>Thr01])==1 |
           length(MarkerRoc02$x)!=1 & length(MarkerRoc02$x[MarkerRoc02$x>Thr01])==1 & MarkerRoc02$x[2]<0){
    MarkerRoc03 <- c(paste(Name01[MarkerRoc02$ix[1]],paste0('Ctr',Thr01),sep=','), paste(MarkerRoc02$x[1],Thr01,sep=','), MarkerRoc02$x[1]);
  }else{
    if(length(MarkerRoc02$x)!=1 & length(MarkerRoc02$x[MarkerRoc02$x>Thr01])==1){
      MarkerRoc02x <- MarkerRoc02$x[1:2]
      MarkerRoc02ix <- MarkerRoc02$ix[1:2]
    }else{
      MarkerRoc02x <- c(MarkerRoc02$x[MarkerRoc02$x>Thr01],Thr01)
      MarkerRoc02ix <- c(MarkerRoc02$ix[MarkerRoc02$x>Thr01],length(MarkerRoc02$ix)+1)
      Name01 <- c(Name01,paste0('Ctr',Thr01))
    }
    Power01 <- Cal_RelativePower(MarkerRoc02x)
    maxPower01 <- which.max(Power01)
    if(maxPower01==length(MarkerRoc02x)){ maxPower01=length(Power01)-1; }

    MarkerRoc03 <- c(paste(Name01[MarkerRoc02ix[1:(maxPower01+1)]],collapse=','), paste(MarkerRoc02x[1:(maxPower01+1)],collapse=','), max(Power01));
  }

  return(MarkerRoc03)
}
