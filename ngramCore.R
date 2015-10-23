

covFwdFile = paste("churchman_",mt,"_minus.bed",sep="")
covRevFile = paste("churchman_",mt,"_plus.bed",sep="")
s<-mt




tmpSense<-tmpASense<-vecS<-NULL
for (j in 1:length(rownames(Allannot))){
  # j is printed in order to have an idea of how slow/fast the algorythm is.
  if (j%%100 == 1){print(j)}
  
  geneID = rownames(Allannot)[j]
  
  #Write the gene coordinates in a temporary file to be used in piping the bash command [UE]
    write.table(Allannot[geneID,c(1:3)],paste('tmpCoord',s,scl,'.txt',sep=''), sep="\t", quote=F, row.names=F, col.names=F)
  
  cr<-Allannot[geneID,c(1:3)]
  
  
  
  if(Allannot[geneID,'strand']=='+'){
    covSense<-covFwdFile
    covASense<-covRevFile
 
  }else
  {

    covSense<-covRevFile
    covASense<-covFwdFile
  }
  
  
  #Custom coordinates
#   write.table(cr,paste('tmpDivCoord',s,scl,'.txt',sep=''), sep="\t", quote=F, row.names=F, col.names=F)
  
  
  #bashCmdSense<-paste('cat tmpCoord6.txt | /hms/scratch1/jd187/Program/bedops/bin/bedextract inputData/',covSense,' -',sep='')
  # bashCmdAntisense<-paste('cat tmpCoord6.txt | /hms/scratch1/jd187/Program/bedops/bin/bedextract inputData/',covASense,' -',sep='')
  
  # Pipe runs the bedextract command on coverage file for the designated gene coordinates [UE]
  
  temp_sense = read.table(pipe(paste(bashCommand,inDirCov, covSense, paste(' tmpCoord',s,scl,'.txt',sep=''), sep='')), sep="\t", header=FALSE,col.names=c('chr', 'start', 'end', 'Cov'), colClasses=c('character', rep('numeric',3)))
  
#   temp_antisense = read.table(pipe(paste(bashCommand, inDirCov,covASense,paste(' tmpCoord',s,scl,'.txt',sep=''), sep='')), sep="\t", header=FALSE,col.names=c('chr', 'start', 'end','Cov'), colClasses=c('character', rep('numeric',3)))
  
  
  if (nrow(temp_sense) > 0){
    temp<-temp_sense
    temp[1,'start'] = cr[2]
    temp[nrow(temp),'end'] = cr[3]
    vv<-temp[,'end'] - temp[,'start']
    nbStarts = rep(temp[vv>0,'Cov'], vv[vv>0])
    posStarts= c(temp[1,'start']:(temp[nrow(temp),'end']-1 ))
    
    pos<-posStarts-posStarts[1]+1
    
    if(Allannot[geneID,'strand']=='-'){
      
      tmpSense<-rev(nbStarts)
    }else{
      
      tmpSense<-nbStarts
    }
    
    
  }
  if(sum(tmpSense)/(sum(tmpSense>0)+.0001)<2){next}

wn<-50
QNT<-5
vecT<-rep(0,length(tmpSense)-wn)
for(qq in seq(1,length(posStarts)-wn,wn))
{
  vecT[qq:(qq+wn-1)]<-round(tmpSense[qq:(qq+wn-1)]*QNT/max(tmpSense[qq:(qq+wn-1)]))
}
vecT<-vecT[!is.na(vecT)]

vecS<-paste(vecS,paste(recoderFunc(vecT,0:QNT,letters[1:(QNT+1)]),collapse=''),collapse='')

#   
#   if (nrow(temp_antisense) > 0){
#     temp<-temp_antisense
#     temp[1,'start'] = cr[2]
#     temp[nrow(temp),'end'] = cr[3]
#     vv<-temp[,'end'] - temp[,'start']
#     nbStarts = rep(temp[vv>0,'Cov'], vv[vv>0])
#     posStarts= c(temp[1,'start']:(temp[nrow(temp),'end']-1 ))
#     
#     pos<-posStarts-posStarts[1]+1
#     if(Allannot[geneID,'strand']=='-'){
#       tmpASense<-rev(nbStarts)
#     }else{
#       tmpASense<-nbStarts
#     }
#     
#     
#   }
  
}

write.table(vecS,file=paste('Symbols_',mt,'.txt',sep=''),sep="\t", quote=F, row.names=F, col.names=F)


r<-textcnt(vecS, method="ngram",n=5L,split=NULL, decreasing=TRUE)
a<-data.frame(counts = unclass(r), size = nchar(names(r)))
frq<-split(a,a$size)

write.table(frq[[5]],file=paste('~/Documents/Data/mutantAnalysis/ngram/',s,'.txt',sep=''), sep="\t", quote=F, row.names=T, col.names=T)

