classification<-function(reversal_pair,Beta,reversal_pair_cancer){

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
  prediction_matrix <- matrix(NA,ncol = length(reversal_pair),nrow = ncol(Beta))
  colnames(prediction_matrix) <- names(reversal_pair)
  rownames(prediction_matrix) <- colnames(Beta)
  existence_pair <- list()
  if(length(reversal_pair)==1){
    pair_existence <- reversal_pair[[1]][,1] %in% rownames(Beta) & reversal_pair[[1]][,2] %in% rownames(Beta)
    count <- colMeans((Beta[reversal_pair[[1]][pair_existence,1],]-Beta[reversal_pair[[1]][pair_existence,2],])>0)
    prediction_matrix[count>0.5,1] <- reversal_pair_cancer
    prediction_matrix[count<=0.5,1] <- gsub('[^a-zA-Z]','',gsub('pair','',gsub(reversal_pair_cancer,'',names(reversal_pair))))
    existence_pair[[1]] <- reversal_pair[[1]][pair_existence,]
    prediction_label <- as.vector(prediction_matrix)
    names(prediction_label) <- rownames(prediction_matrix)
  } else {
    for(i in 1:length(reversal_pair)){
      pair_existence <- reversal_pair[[i]][,1] %in% rownames(Beta) & reversal_pair[[i]][,2] %in% rownames(Beta)
      count <- colMeans((Beta[reversal_pair[[i]][pair_existence,1],]-Beta[reversal_pair[[i]][pair_existence,2],])>0)
      prediction_matrix[count>0.5,i] <- reversal_pair_cancer
      prediction_matrix[count<=0.5,i] <- paste('non-',reversal_pair_cancer,sep='')
      existence_pair[[i]] <- reversal_pair[[i]][pair_existence,]
    } 
    prediction_label <- rowSums(prediction_matrix==reversal_pair_cancer)
    GM_index <- prediction_label==length(reversal_pair)
    prediction_label[GM_index] <- reversal_pair_cancer
    prediction_label[!GM_index] <- paste('non',reversal_pair_cancer,sep='-')
  }
  
  result <- list(prediction_matrix=prediction_matrix,prediction_label=prediction_label,existence_reversal_pair=existence_pair)
  return(result)
}