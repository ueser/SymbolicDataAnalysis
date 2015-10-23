s1<-readWave('~/Downloads/HMjQygwPI1c.wav')
LR=(s1@left+s1@right)/2


scl<-8

r_n<-6
wn<-exp(scl)
k2 <- kernel("daniell", round(wn)) 
SYMDATA<-data.table()



x2<-kernapply(LR,k2)
q<-cut(1:length(x2),seq(0,length(LR)+1,wn))
mean.mag <- tapply(x2, q, mean,na.rm=T)

readSense<-max(LR)*(mean.mag-min(mean.mag,na.rm=T))/max(mean.mag-min(mean.mag,na.rm=T),na.rm=T)

qr<-round(readSense)
sym<-round(r_n*qr/max(qr,na.rm=T))
sym[round(readSense)==0]<-0
sym<-sym[!is.na(sym)]

nDT<-data.table(geneName='Mozart',symPausesSense=paste(sym,collapse=''))
SYMDATA<-rbindlist(list(SYMDATA,nDT))
setkey(SYMDATA,'geneName')

corpus<-new.env()
SYMDATA<<-SYMDATA
lngt<<-SYMDATA[,nchar(symPausesSense),by=geneName]

ngram<-list(unigram="0",bigram="00",trigram="000",fourgram="0000",fivegram="00000",sixgram="000000")

vec<-SYMDATA[1,symPausesSense]
  for(m_n in 1:6){
    #        print(paste("Proccessing ",mut,'...',m_n,'grams'))
    idx<-1:(nchar(vec)-m_n)
    newtext<-substring(vec,idx,idx+m_n-1)
    jdx<-which(newtext%in% paste(rep(0,m_n),collapse=""))
    ngram[[m_n]]<-unique(c(ngram[[m_n]],newtext))
  }





source('~/Google Drive/Projects/MutantScreening/buildCorpus.R')

initializeCorpus(r_n,envir=corpus)
for(i in 1:length(ngram))
{
  for(j in 1:length(ngram[[i]]))
  {
    qword<-ngram[[i]][j]
    if(nchar(qword)==0){
      print(paste('Skipping empty word in ',i,'gram'))
    }
    print(paste('Evaluating ',qword,'and',i,'gram'))    
    
    Pexp<-getExpectedFreq(qword,envir=corpus)
    
  }
  print(paste('Evaluating ',i,'gram'))
  
}



nameSpace<-ls(corpus)
nameSpace<-nameSpace[!is.na(as.numeric(nameSpace))]



corDT<-data.table(Word=rep(0,0),Ngram=rep(0,0))
for(i in 1:length(nameSpace))
{
  nDT<-envToDT(nameSpace[i],envir=corpus)
  corDT<-rbindlist(list(corDT,nDT))
}

setkey(corDT,'Word','Ngram')
assign(paste('corDT_Mozart',sep=''),corDT)



##########

load("/Users/umuteser/Google Drive/Projects/MutantScreening/CorpusWT_Nov27.rda")
nameSpace<-corDT_WT[,Word[log(Enrichment)>.2]]


idx<-jdx<-vdx<-NULL
dimnames<-list(RowName=NULL,ColName=NULL)
for(i in 1:length(nameSpace))
{
  if(i%%100==1){print(i/length(nameSpace))}
  obj<-get(nameSpace[i],envir=corpus)
  tbl<-table(match(obj$Position$geneName,geneNames))
  jdx<-c(jdx,as.numeric(names(tbl)))
  vdx<-c(vdx,tbl)
  idx<-c(idx,rep(i,length(tbl)))
  dimnames$RowName<-c(dimnames$RowName,rep(nameSpace[i],length(tbl)))
  dimnames$ColName<-c(dimnames$ColName,names(table(obj$Position$geneName)))
}

stm<-simple_triplet_matrix(idx, jdx, vdx, dimnames = dimnames)
save(stm,'simpleTripletMatrix.rdata')


###

load('CorpusWT_Nov27.rda')
vocab<-corDT_WT[,Word[log(Enrichment)>.2]]
# vocab<-corDT_RCO1D[,Word[log(Enrichment)>.2]]
documents<-list()
for(i in 1:length(vocab))
{
  if(i%%100==1){print(i/length(vocab))}
  obj<-get(vocab[i],envir=corpus)
  tbl<-table(obj$Position$geneName)
    
  for(gene in names(tbl)){
    
    documents[[gene]]<-cbind(documents[[gene]],as.integer(c(i-1,tbl[[gene]])))
  }
 
}


ldam<-lda.collapsed.gibbs.sampler(documents, K=5, vocab, num.iterations = 100, alpha = 1, eta = 10)
topDocs<-top.topic.documents(ldam$document_sums)
topWords<-top.topic.words(ldam$topics)

#ldamWT<-ldam
#documentsWT<-documents
reads<-readNETseq(geneID = names(documents)[topDocs[,2]][1],annot = ScerAnnot,covFwdFile = )
####


load('CorpusDST1D_Nov27.rda')
corpusDST1D<-corpus
load('CorpusWT_Nov27.rda')

load('corpusDataTableAllMutantsScale0.rdata')
vocabDST1D<-corDT_DST1D[,Word[log(Enrichment)>.2]]
vocabWT<-corDT_WT[,Word[log(Enrichment)>.2]]


vocab<-corDT_DST1D[,Word[log(Enrichment)>.2]]
wordTransformDST1D<-data.table()
for(i in 1:length(vocab))
{
  if(i%%100==1){print(i/length(vocab))}
  obj<-get(vocab[i],envir=corpusDST1D)
 
 nDT<-cbind(rep(vocab[i],length(obj$Position$geneName)),obj$Position) 
  wordTransformDST1D<-rbindlist(list(wordTransformDST1D,nDT))
 
}
setnames(wordTransformWT,c('V1','geneName'),c('Word','geneName'))
setkey(wordTransformWT,'geneName')



####

newDat<-rawReads[1:100000,]

library(ggplot2)

ggplot(newDat,aes(x=WTSense,y=DST1DSense))+
  geom_point(alpha=0.3)+
  geom_density2d()+
  theme_bw()+
  scale_x_log10()+
  scale_y_log10()+
  labs(x='WT Reads ',y='dst1∆ Reads')
