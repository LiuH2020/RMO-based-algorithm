calculate.similarity<-function(reversal_pair,Beta,reversal_pair_cancer=NULL){

  if (!is.list(reversal_pair)|is.data.frame(reversal_pair)){
    stop("Warning: reversal_pair is not a list\n")
  }
  
  for(j in 1:length(reversal_pair)){
    name <- unlist(strsplit(names(reversal_pair)[j],'[^a-zA-Z0-9]'))
    index <- match(reversal_pair_cancer,name)
    if(index!=1){
      reversal_pair[[j]] <- as.data.frame(cbind(reversal_pair[[j]][,2],reversal_pair[[j]][,1]))
    }
  }
  similarity <- matrix(NA,ncol = length(reversal_pair),nrow = ncol(Beta))
  existence_pair <- list()
  for(i in 1:length(reversal_pair)){
    pair_existence <- reversal_pair[[i]][,1] %in% rownames(Beta) & reversal_pair[[i]][,2] %in% rownames(Beta)
    count <- colMeans((Beta[reversal_pair[[i]][pair_existence,2],]-Beta[reversal_pair[[i]][pair_existence,1],])>0)
    similarity[,i] <- count
    existence_pair[[i]] <- reversal_pair[[i]][pair_existence,]
  }
  names(existence_pair) <- names(reversal_pair)
  colnames(similarity) <- gsub('[^a-zA-Z]','',gsub('pair','',gsub(reversal_pair_cancer,'',names(reversal_pair))))
  rownames(similarity) <- colnames(Beta)
  
  cancer_highest_similarity_index <- apply(similarity,1,order)[length(reversal_pair),]
  predction <- colnames(similarity)[cancer_highest_similarity_index]
  names(predction) <- rownames(similarity)
  
  max_similarity <- apply(similarity,1,max)
  if(any(max_similarity<0.5)){
    predction[max_similarity<0.5] <- reversal_pair_cancer
  }
  result <- list(similarity=similarity,prediction_label=predction,existence_reversal_pair=existence_pair)
  return(result)
}
