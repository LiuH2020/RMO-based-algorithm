######step 1:Select the stable CpG site pairs for eight cancers######
#BRC, CRC, GC, GM, LC, LIC, MA, PC
source('F:/The origin of brain metastases/code/selec.stable.pair.R')
input <- 'F:/The origin of brain metastases/profile/training'
output <- 'F:/The origin of brain metastases/output/stable_pair'
cancer_type <- c('BRC','CRC','GC','GM','LC','LIC','MA','PC')
for(i in 1:length(cancer_type)){
  file_name <- paste(input,'/',cancer_type[i],'.RData',sep='')
  load(file_name)
  pair<-selec.stable.pair(Beta = Beta,cut_off = 0.95)
  file_name <- paste(output,'/',cancer_type[i],'_pair.RData',sep='')
  save(pair,file=file_name)
  rm(pair,Beta,annotation)
}
rm(list=ls())

#################step 2:Select the reverse CpG site pairs between GM and seven other cancer types####################
library(dplyr)
input <- 'F:/The origin of brain metastases/output/stable_pair'
cancer <- 'GM' #select the reverse CpG site pairs of GM
output <- 'F:/The origin of brain metastases/output/reversal_pair'

file <- list.files(input,pattern = '*_pair.RData')
filecase <- paste(input,file,sep='/')
cancer_index <- grep(cancer,filecase)
load(filecase[cancer_index])
cancer_pair <- as.data.frame(pair)
rm(pair)
filecase <- filecase[-cancer_index]
file <- file[-cancer_index]
reversal_pair <- list()
for(i in 1:length(filecase)){
  load(filecase[i])
  pair <- as.data.frame(cbind(pair[,2],pair[,1]))
  reversal_pair[[i]] <- dplyr::intersect(cancer_pair,pair)
  names(reversal_pair)[i] <- paste(cancer,gsub('_pair.RData','',file[i]),sep='-')
  print(dim(reversal_pair[[i]]))
  rm(pair)
}
output_name <- paste(output,'/',cancer,'_reversal_pair.RData',sep='')
save(reversal_pair,file=output_name)

############step 3:The performance of reversal pair in Train group and reclassify################# 
###Distinguishing GMs from other cancers
source('F:/The origin of brain metastases/code/calculate.similarity.R')
source('F:/The origin of brain metastases/code/classification.R')
input1 <- 'F:/The origin of brain metastases/output/reversal_pair'
input2 <- 'F:/The origin of brain metastases/profile'
GM_reversal_file_name <- paste(input1,'GM_reversal_pair.RData',sep='/')
train_beta_file_name <- paste(input2,'training.RData',sep='/')
output <- 'F:/The origin of brain metastases/output'

load(GM_reversal_file_name)#The reversal pair of GM
load(train_beta_file_name)
train_beta <- Beta
train_annotation <- annotation
GM_classification<-classification(reversal_pair = reversal_pair,Beta = train_beta,reversal_pair_cancer = 'GM')

table(train_annotation$Annotation,GM_classification$prediction_label)
GM_prediction_index <- GM_classification$prediction_label=='GM'
prediction_label<-GM_classification$prediction_label
###Calculate the similaity to identify the other cancers
similarity<-calculate.similarity(reversal_pair = reversal_pair,Beta = train_beta[,!GM_prediction_index],reversal_pair_cancer = 'GM')
prediction_label[!GM_prediction_index]<-similarity$prediction_label
table(train_annotation$Annotation,prediction_label)
train_no_optimize_prediction_name <- paste(output,'/','train_no_optimize_prediction.RData',sep='')
save(train_annotation,train_beta,prediction_label,file=train_no_optimize_prediction_name)
rm(list=ls())
###reclassify the above identified LC samples
#select the reversal CpG site pairs for LC
source('F:/The origin of brain metastases/code/selec.stable.pair.R')
input <- 'F:/The origin of brain metastases/profile/training'
output <- 'F:/The origin of brain metastases/output/stable_pair_90'
cancer_type <- c('BRC','CRC','GC','LC','LIC','MA','PC')
for(i in 1:length(cancer_type)){
  file_name <- paste(input,'/',cancer_type[i],'.RData',sep='')
  load(file_name)
  pair<-selec.stable.pair(Beta = Beta,cut_off = 0.9)
  file_name <- paste(output,'/',cancer_type[i],'_pair.RData',sep='')
  save(pair,file=file_name)
  rm(pair,Beta,annotation)
}
rm(list=ls())


library(dplyr)
input <- 'F:/The origin of brain metastases/output/stable_pair_90'
cancer <- 'LC' #select the reverse CpG site pairs of GM
output <- 'F:/The origin of brain metastases/output/reversal_pair'

file <- list.files(input,pattern = '*_pair.RData')
filecase <- paste(input,file,sep='/')
cancer_index <- grep(cancer,filecase)
load(filecase[cancer_index])
cancer_pair <- as.data.frame(pair)
rm(pair)
filecase <- filecase[-cancer_index]
file <- file[-cancer_index]
reversal_pair <- list()
for(i in 1:length(filecase)){
  load(filecase[i])
  pair <- as.data.frame(cbind(pair[,2],pair[,1]))
  reversal_pair[[i]] <- dplyr::intersect(cancer_pair,pair)
  names(reversal_pair)[i] <- paste(cancer,gsub('_pair.RData','',file[i]),sep='-')
  print(dim(reversal_pair[[i]]))
  rm(pair)
}
output_name <- paste(output,'/',cancer,'_reversal_pair.RData',sep='')
save(reversal_pair,file=output_name)
rm(list=ls())
#reclassify the above identified LC samples
source('F:/The origin of brain metastases/code/calculate.similarity.R')
load('F:/The origin of brain metastases/output/train_no_optimize_prediction.RData')
load('F:/The origin of brain metastases/output/reversal_pair/LC_reversal_pair.RData')
output <- 'F:/The origin of brain metastases/output'

LC_predict_index <- prediction_label=='LC'
similarity <- calculate.similarity(reversal_pair = reversal_pair,Beta = train_beta[,LC_predict_index],reversal_pair_cancer = 'LC')
prediction_label[LC_predict_index] <-similarity$prediction_label
train_LC_reclassify_prediction_name <- paste(output,'/','train_LC_reclassify_prediction.RData',sep='')
save(train_annotation,train_beta,prediction_label,file=train_LC_reclassify_prediction_name)

###reclassify the above identified GC and CRC samples
#select CRC-GC reversal CpG site pair
library(dplyr)
input <- 'F:/The origin of brain metastases/output/stable_pair_90'
output <- 'F:/The origin of brain metastases/output/reversal_pair'
CRC_pair_data_name <- paste(input,'CRC_pair.RData',sep='/')
GC_pair_data_name <- paste(input,'GC_pair.RData',sep='/')
load(CRC_pair_data_name)
CRC_pair <- as.data.frame(pair)
rm(pair)
load(GC_pair_data_name)
GC_pair <- as.data.frame(pair)
rm(pair)
CRC_GC_pair <- dplyr::intersect(CRC_pair,as.data.frame(cbind(GC_pair[,2],GC_pair[,1])))
CRC_GC_pair_file_name <-paste(output,'CRC_GC_pair.RData',sep='/')
save(CRC_GC_pair,file=CRC_GC_pair_file_name)
rm(list=ls())

#reclassify the above identified GC and CRC samples
input1 <-'F:/The origin of brain metastases/output/reversal_pair'
input2 <- 'F:/The origin of brain metastases/output'
CRC_GC_pair_file_name <-paste(input1,'CRC_GC_pair.RData',sep='/')
train_LC_reclassify_prediction_name <-paste(input2,'train_LC_reclassify_prediction.RData',sep='/')
source('F:/The origin of brain metastases/code/classification.R')
load(CRC_GC_pair_file_name)
load(train_LC_reclassify_prediction_name)
CRC_and_GC_predict_index <- prediction_label=='GC'|prediction_label=='CRC'

classfict<-classification(reversal_pair = list(CRC_GC_pair=CRC_GC_pair),Beta = train_beta[,CRC_and_GC_predict_index],reversal_pair_cancer = 'CRC')
prediction_label[CRC_and_GC_predict_index] <- classfict$prediction_label
table(train_annotation$Annotation,prediction_label)

################step 4:The performance of RMO-based method in validation group##############
##identify the origin of BMs with GM reversal pair 
source('F:/The origin of brain metastases/code/calculate.similarity.R')
source('F:/The origin of brain metastases/code/classification.R')
input1 <- 'F:/The origin of brain metastases/output/reversal_pair'
input2 <-'F:/The origin of brain metastases/profile'
GM_reversal_pair_file_name <- paste(input1,'GM_reversal_pair.RData',sep='/')
LC_reversal_pair_file_name <- paste(input1,'LC_reversal_pair.RData',sep='/')
CRC_GC_pair_file_name <- paste(input1,'CRC_GC_pair.RData',sep='/')
validation_beta_file_name <- paste(input2,'validation.RData',sep='/')
output <- 'F:/The origin of brain metastases/output'

load(validation_beta_file_name)
test_beta <- Beta
test_annotation <- annotation
load(GM_reversal_pair_file_name)#The reversal pair of GM
classfict <- classification(reversal_pair = reversal_pair,Beta = test_beta,reversal_pair_cancer = 'GM')
prediction_label <- classfict$prediction_label
table(test_annotation$Annotation,prediction_label)

GM_predict_index <- prediction_label=='GM'
similarity <- calculate.similarity(reversal_pair = reversal_pair,Beta = test_beta[,!GM_predict_index],reversal_pair_cancer = 'GM')
prediction_label[!GM_predict_index] <- similarity$prediction_label
rm(reversal_pair)
##reclassify the above identified LC samples
load(LC_reversal_pair_file_name)
LC_predict_index <- prediction_label=='LC'
similarity <- calculate.similarity(reversal_pair = reversal_pair,Beta = test_beta[,LC_predict_index], reversal_pair_cancer = 'LC')
prediction_label[LC_predict_index] <- similarity$prediction_label

##reclassify the above identified GC and CRC samples
load(CRC_GC_pair_file_name)
CRC_or_GC_predict_index <- prediction_label=='CRC'|prediction_label=='GC'
classify <- classification(reversal_pair = list(CRC_GC_pair=CRC_GC_pair),Beta = test_beta[,CRC_or_GC_predict_index],reversal_pair_cancer = 'CRC')
prediction_label[CRC_or_GC_predict_index] <- classify$prediction_label
table(test_annotation$Annotation,prediction_label)
