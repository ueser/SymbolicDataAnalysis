.libPaths(c(.libPaths(),'/groups/churchman/ue4/Rlibrary/'))
library(data.table)

source('getParameters.r')
out<-getParameters()

#load annotation
# load(out$annotationFile)
# print("Annotation loaded ")



projectName<-out$projectName
r_n<-out$resolution
inDirCov<-out$initialDataDir
saveDir<-file.path(out$baseDir,out$projectName,"symbolicData")


args = commandArgs(trailingOnly = T)

scl<-as.numeric(args[[1]])
mut<-args[[2]]
wn<-exp(scl)
k2 <- kernel("daniell", round(wn)) 


print(paste("Proccessing ",mut,'...'))


SYMDATA<-data.table()





for(ch in 1:16)
{
  covFile = file.path(inDirCov,paste(mut,'chr',ch,'.wig',sep=''))

  dat<-read.table(covFile,colClasses = "numeric",skip=1)
  readSense<-dat[,2]
  if(scl>0){
      ### -- For scale-space
        
        x2<-kernapply(readSense,k2)
        q<-cut(1:length(x2),seq(0,length(readSense)+1,wn))
        mean.mag <- tapply(x2, q, mean,na.rm=T)
      ### --- 
      readSense<-max(readSense)*(mean.mag-min(mean.mag,na.rm=T))/max(mean.mag-min(mean.mag,na.rm=T),na.rm=T)
     
      
  }
  # take log of the reads to classify the pauses based on their order of magnitudes. Underlying assumption is
  # that different order of magnitudes come from different mechanisms
  # define words by normalizing the reads between 0 and (r_n-1) within each gene
  qr<-log(round(readSense)+.001)
  sym<-round(r_n*qr/max(qr,na.rm=T))
  sym[round(readSense)==0]<-0
  sym<-sym[!is.na(sym)]
  


  nDT<-data.table(chrName=paste('chr',ch,sep=""),symSignal=paste(sym,collapse=''))
  
  SYMDATA<-rbindlist(list(SYMDATA,nDT))
  


}


setkey(SYMDATA,'chrName')

ngram<-list(unigram="0",bigram="00",trigram="000",fourgram="0000",fivegram="00000",sixgram="000000")
for(i in 1:dim(SYMDATA)[1])
{
  vec<-SYMDATA[i,symSignal]
  for(m_n in 1:6){
#        print(paste("Proccessing ",mut,'...',m_n,'grams'))
    idx<-1:(nchar(vec)-m_n)
    newtext<-substring(vec,idx,idx+m_n-1)
    jdx<-which(newtext%in% paste(rep(0,m_n),collapse=""))
    ngram[[m_n]]<-unique(c(ngram[[m_n]],newtext))
  }
}



lngt<-SYMDATA[,nchar(symSignal),by=chrName]

save(list=c('SYMDATA','ngram','lngt'),file=file.path(saveDir,paste(projectName,'_',mut,'_ScaleExp',scl,'.rdata',sep='')))



