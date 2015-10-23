
dist.KLD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=KLD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  resultsMatrix
  return(resultsMatrix) 
}

.libPaths(c(.libPaths(),'/groups/churchman/ue4/Rlibrary/'))
library(data.table)

source('getParameters.r')
out<-getParameters()

inDir<-file.path(out$baseDir,out$projectName,"corpus")
saveDir<-file.path(out$baseDir,out$projectName,"KLdivergence")

KLdistMat<-list()
for(scl in 1:as.numeric(out$scales)){

  load(file.path(inDir,paste('corpusDataTableAllMutantsScale',scl,'.rdata',sep='')))
  allWords<-NULL
  
  for(mut in out$samples){
    allWords<-union(allWords,get(paste('corDT_',mut,sep=''))[,Word[log(Enrichment)>-2&log10(Observed)>-5]])
  }
  inMat<-NULL
  allWords<-allWords[!is.na(allWords)]
  for(mut in out$samples){
    inMat<-as.data.frame(cbind(inMat,get(paste('corDT_',mut,sep=''))[allWords,Observed]))
  }
  inMat[inMat<0]<-0
  inMat[is.na(inMat)]<-0
  
  rownames(inMat)<-allWords
  colnames(inMat)<-out$samples
  cs<-colSums(inMat)
  inMat<-t(apply(inMat,MAR=1,FUN=function(x) x/cs))
  
  KLdistMat<-c(KLdistMat,list(dist.KLD(inMat)))

}
save(KLdistMat,file=file.path(saveDir,'KLdistanceForAllScales.rdata'))