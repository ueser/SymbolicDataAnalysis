Y<-matrix(0,nrow=length(scl),ncol=length(nbStarts))
Y<-sapply(scl,function(x) g_conv(nbStarts,x))

Y<-matrix(0,nrow=length(scl),ncol=length(Tdst1))
Y<-sapply(scl,function(x) g_conv(Tdst1,x))

Ynorm<-apply(Y,MARGIN=2,FUN=function(x) (x-min(x))/max(x-min(x)))
  
image.plot(x=pos,y=log(scl),Ynorm)


##

s<-'DST1D'

frqN<-read.table(file=paste('~/Documents/Data/mutantAnalysis/ngram/',s,'.txt',sep=''),sep="\t",row.names=1, header=T)
QV<-0:QNT
x<-strsplit(rownames(frqN),'*')
newX<-matrix(0,nrow=length(rownames(frqN)),ncol=QNT)
newF<-vector()
k<-0
for(i in 1:length(rownames(frqN))){

  if(any(x[[i]]=="f")){
    k<-k+1
  for(qq in 1:length(QV)){
   
      newX[k,x[[i]]==letters[qq]]<-QV[qq]}
  newF[k]<-frqN[i,2]
  }
}
newX<-newX[1:k,]


P11<-c(0,0,5,0,0)
P12<-c(5,0,0,0,5)

P21<-c(5,4,3,2,1)
P22<-c(1,2,3,4,5)

n1<-(P11-P12)/sqrt(sum((P11-P12)^2))
n2<-(P21-P22)/sqrt(sum((P21-P22)^2))
d1<-vector()
d2<-vector()
for(k in 1:length(newX[,1])){
  P<-newX[k,]
  d2[k]<-sqrt(sum(P11-P - (sum((P11-P)*n1)*n1))^2) *(sum(P*c(1,1,0,-1,-1))+.01)/(abs(sum(P*c(1,1,0,-1,-1))+.01))
  d1[k]<-sqrt(sum(P21-P - (sum((P21-P)*n2)*n2))^2) *(sum(P*c(-1,0,1,0,-1))+.01)/(abs(sum(P*c(-1,0,1,0,-1))+.01))
}


dfx = data.frame(ev1=d1, ev2=d2, ev3=sqrt(newF/sum(newF)))

symbols(x=dfx$ev1, y=dfx$ev2, circles=dfx$ev3, inches=1/30, ann=F, bg="dark green", fg=NULL,add=F)
abline(v=0,h=0)
title(main=st, xlab="Symmetry",ylab="Skewness")


