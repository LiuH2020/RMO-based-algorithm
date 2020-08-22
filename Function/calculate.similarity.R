calculate.similarity<-function(reversal_pair,Beta,reversal_pair_cancer=NULL){
  similarity<-matrix(NA,ncol = length(reversal_pair),nrow = ncol(Beta))
  existence_pair<-list()
  for(i in 1:length(reversal_pair)){
    pair_existence <- reversal_pair[[i]][,1] %in% rownames(Beta) & reversal_pair[[i]][,2] %in% rownames(Beta)
    count <- colMeans((Beta[reversal_pair[[i]][pair_existence,2],]-Beta[reversal_pair[[i]][pair_existence,1],])>0)
    similarity[,i]<-count
    existence_pair[[i]]<-reversal_pair[[i]][pair_existence,]
  }
  names(existence_pair)<-names(reversal_pair)
  colnames(similarity)<-sapply(strsplit(names(reversal_pair),'-'),function(x) x[2])
  rownames(similarity)<-colnames(Beta)
  
  cancer_highest_similarity_index<-apply(similarity,1,order)[length(reversal_pair),]
  predction<-colnames(similarity)[cancer_highest_similarity_index]
  names(predction)<-rownames(similarity)
  
  max_similarity<-apply(similarity,1,max)
  if(any(max_similarity<0.5)){
    predction[max_similarity<0.5] <- reversal_pair_cancer
  }
  result<-list(similarity=similarity,prediction_label=predction,existence_reversal_pair=existence_pair)
  return(result)
}