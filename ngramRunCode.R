#set working directory
setwd('~/mutantAnalysis/')

#load All's gene annotation workspace
load('~/mutantAnalysis/genesAnnotationBurak.RData')

options("scipen"=100, "digits"=22)


geneAll = as.data.frame(genes)
trxAll  = as.data.frame(transcripts, stringsAsFactors =FALSE)
rownames(geneAll) = geneAll$gene
rownames(trxAll)  = trxAll$gene



Allannot = geneAll[,c(2,4,5,3,6)]
# to make the annotation 0 index based (we don't do it for the end position because it is not comprised)
Allannot$start  = Allannot$start - 1
Allannot$strand = replace(Allannot$strand, Allannot$strand==1, '+') # to make the strand written the same way as in Stirling's
Allannot$strand = replace(Allannot$strand, Allannot$strand==-1,'-') # annotation


# bash script to extract gene locus coverage:
#/usr/local/bin/bedextract $1 $2
# where $1 is the bed coverage file and $2 is the gene Coordinates (format : chr\tStart\tEnd)
# Note: the coverage file is 0 base Indexed!!
#bashCommand='bash /hms/scratch1/jd187/GenomeChromosome/sacCer3/Annotation/scripts/geneExtraction2.sh '
# bashCommand='bash ~/mutantAnalysis/geneExtraction2.sh '
bashCommand='bash geneExtraction2.sh '

inDirCov='inputData/MagdalenaData/all.files/'

lst<-dir(path=inDirCov,pattern='*.bed$')

mutList1<-c('WT','SET1D','SET2D','RCO1D','EAF3D','DST1D')

for(i in 3:length(mutList1)){
mt<-mutList1[i]
scl<-2

covFwdFile = paste("churchman_",mt,"_minus.bed",sep="")
covRevFile = paste("churchman_",mt,"_plus.bed",sep="")
s<-mt


source('~/Documents/Codes/ngramCore.R')

}
#####

recoderFunc <- function(data, oldvalue, newvalue) {
  
  # convert any factors to characters
  
  if (is.factor(data))     data     <- as.character(data)
  if (is.factor(oldvalue)) oldvalue <- as.character(oldvalue)
  if (is.factor(newvalue)) newvalue <- as.character(newvalue)
  
  # create the return vector
  
  newvec <- data
  
  # put recoded values into the correct position in the return vector
  
  for (i in unique(oldvalue)) newvec[data == i] <- newvalue[oldvalue == i]
  
  newvec
  
}

library(seqinr) 
