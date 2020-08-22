selec.stable.pair<-function(Beta,cut_off=0.99){
  library(data.table)
  path<-getwd()
  time<-as.numeric(Sys.time())
  tmp<-paste(path,'/',time,'_tmp.txt',sep='')
  n_col<-dim(Beta)[2]
  n_row<-dim(Beta)[1]
  row_name<-rownames(Beta)
  exp2<-Beta
  exp_name<-rownames(exp2)
  Stat<-Sys.time()
  for(i in 1:(n_row-1)){
    #print(i)
    
    exp1<-matrix(1,ncol=1,nrow=n_row-i)%*%Beta[i,]
    exp2<-exp2[-1,]
    count<-exp1-exp2
    exp_name<-exp_name[-1]
    re1<-exp_name[rowMeans(count>0)>cut_off]
    re2<-exp_name[rowMeans(count<0)>cut_off]
    len1<-length(re1)
    len2<-length(re2)
    CpG_i<-c(rep(row_name[i],len1),re2)
    CpG_j<-c(re1,rep(row_name[i],len2))
    CpG<-data.frame(cbind(CpG_i,CpG_j))
    fwrite(CpG,file=tmp,append=T,sep='\t',row.names = F,col.names = F)
    rm(CpG,CpG_i,CpG_j,exp1,count,re1,re2,len1,len2)
  }
    end1<-Sys.time()
    print(end1-Stat)
  pair<-fread(tmp,header = F)
  pair<-as.data.frame(pair)
  unlink(tmp)
  end2<-Sys.time()
  print(end2-Stat)
  return(pair)
}