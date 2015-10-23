
envToDT<-function(qword,envir=lexicon){
  mns<-get(qword,envir=lexicon)[[2]][,mean(V1,na.rm=T),by=chrName]
  
  DT<-data.table(Word=qword,
      			Expected=get(qword,envir=lexicon)[[3]],
  				Observed=get(qword,envir=lexicon)[[1]],
  				Enrichment=get(qword,envir=lexicon)[[4]],
  				PositionMN=mean(mns$V1/lngt[mns[,chrName],]$V1,na.rm=T),
  				PositionSD=sd(mns$V1/lngt[mns[,chrName],]$V1),
  				PositionSEM=sd(mns$V1/lngt[mns[,chrName],]$V1,na.rm=T)/sqrt(sum(!is.na(mns$V1))),
  				Ngram=nchar(qword))
  return(DT)
}

source('getParameters.r')
out<-getParameters()

.libPaths(c(.libPaths(),'/groups/churchman/ue4/Rlibrary/'))
library(data.table)

inDir<-file.path(out$baseDir,out$projectName,"lexicon")
saveDir<-file.path(out$baseDir,out$projectName,"corpus")


args = commandArgs(trailingOnly = T)

scl<-as.numeric(args[[1]])

saveList<-NULL


for(mut in out$samples){
	load(file.path(inDir,paste('lexicon_',mut,'_','Scale',scl,'.rda',sep='')))
	nameSpace<-ls(lexicon)
	nameSpace<-nameSpace[!is.na(as.numeric(nameSpace))]



	corDT<-list()
	for(i in 1:length(nameSpace))
	  {
	  	tryCatch({
	    nDT<-envToDT(nameSpace[i],envir=lexicon)
	    corDT<-rbindlist(list(corDT,nDT))
	    }, error=function(e) {})
	  }
       
	setkey(corDT,'Word','Ngram')
	assign(paste('corDT_',mut,sep=''),corDT)
	saveList<-c(saveList,paste('corDT_',mut,sep=''))
}

save(list=saveList,file=file.path(saveDir,paste('corpusDataTableAllMutantsScale',scl,'.rdata',sep='')))