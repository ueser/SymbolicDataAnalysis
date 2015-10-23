.libPaths(c(.libPaths(),'/groups/churchman/ue4/Rlibrary/'))
library(data.table)

source('getParameters.r')
out<-getParameters()
source(file.path(out$baseDir,"Codes","SymbolicAnalysis","sourceFunctionsNuc.R"))



r_n<-out$resolution
projectName<-out$projectName

inDir<-file.path(out$baseDir,out$projectName,"symbolicData")
saveDir<-file.path(out$baseDir,out$projectName,"lexicon")



args = commandArgs(trailingOnly = T)

scl<-as.numeric(args[[1]])
mut<-args[[2]]
wn<-exp(scl)

t0<-proc.time()

load(file.path(inDir,paste(projectName,'_',mut,'_ScaleExp',scl,'.rdata',sep='')))

print(paste('Loading the previous data took ', (proc.time()-t0)[[3]], ' seconds',sep=""))

SYMDATA<<-SYMDATA
lngt<<-lngt
lexicon<<-new.env()
initializeLexicon(r_n,envir=lexicon)

print(paste('Lexicon is initialized for:', mut, ' -at scale:',scl,' -with time:', (proc.time()-t0)[[3]], ' seconds',sep=""))

# Record the observed frequency of each word
# Calculate the frequencies for 6 ngram categories
ngramFreq<-data.table(Ngram=1:r_n,Frequency=rep(0,r_n))
for(i in 1:length(ngram))
{
  print(paste('Evaluating ',i,'gram Observed freq.s', '-with time:', (proc.time()-t0)[[3]], ' seconds',sep=""))

  Png<-0
  for(j in 1:length(ngram[[i]]))
  {
    qword<-ngram[[i]][j]
    if(nchar(qword)==0){
      print(paste('Skipping empty word in ',i,'gram'))
    }
    # print(paste('Evaluating ',qword,'and',i,'gram'))    
    

    Pobs<-getObservedFreq(qword,envir=lexicon)
    Png<-Png+Pobs$Observed
    
  }


  ngramFreq[i,Frequency:=Png]
  
  # print(paste('Evaluating ',i,'gram'))
  
}

ngramFreq<<-ngramFreq

# Calculate the expected frequency for unigrams
for(i in 1:r_n){
  obs<-get(ngram[[1]][i],envir=lexicon)
  assign(ngram[[1]][i],c(obs,Expected=1/r_n,Enrichment=obs$Observed*r_n),envir=lexicon)
}

# Calculate the expected frequency for the rest of the ngrams
for(i in 2:length(ngram))
{
    print(paste('Evaluating ',i,'gram Expected freq.s', '-with time:', (proc.time()-t0)[[3]], ' seconds',sep=""))

  for(j in 1:length(ngram[[i]]))
  {
    qword<-ngram[[i]][j]
    if(nchar(qword)==0){
      print(paste('Skipping empty word in ',i,'gram'))
    }
    # print(paste('Evaluating ',qword,'and',i,'gram'))    
    
    Pexp<-getExpectedFreq(qword,envir=lexicon)
    
  }
  # print(paste('Evaluating ',i,'gram'))
  
}


save(list=c('lexicon','SYMDATA','ngram','lngt'),file=file.path(saveDir,paste('lexicon_',mut,'_','Scale',log(wn),'.rda',sep='')))

  print(paste('succesfully saved -with time:', (proc.time()-t0)[[3]], ' seconds',sep=""))

