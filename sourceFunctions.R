
library(R.utils) #used in getExpectedFreq function

library(data.table) #commonly used

# library(Biostrings) # used in getObservedFreq function for gregexpr2 


gregexpr2<-function (pattern, text) 
{
  if (!is.character(pattern) || length(pattern) != 1 || is.na(pattern) || 
        nchar(pattern) == 0) 
    stop("invalid pattern")
  matches <- gregexpr(pattern, text, fixed = TRUE)
  nP <- nchar(pattern)
  for (i in 1:length(text)) {
    mi <- matches[[i]]
    if (length(mi) == 1 && mi == -1) {
      attr(matches[[i]], "match.length") <- NULL
    }
    else {
      subtexts <- substring(text[i], mi + 1, mi + 2 * nP -2)
      missing_matches <- gregexpr2(pattern, subtexts)
      for (j in 1:length(mi)) {
        mj <- missing_matches[[j]]
        if (length(mj) != 1 || mj != -1) 
          matches[[i]] <- c(matches[[i]], mi[j] + mj)
      }
      matches[[i]] <- sort(matches[[i]])
    }
  }
  matches
}

##################################################
# initialization of the lexicon                   #
##################################################
# r_n is the resolution, i.e. the number of symbols
# r_n is used for naive probability 1/r_n for the expected of 1-grams

initializeLexicon<-function(r_n,envir=lexicon)
{
  
  alph <- list(0:(r_n))
  for(i in 1:length(alph[[1]]))
  {
    getObservedFreq(as.character(alph[[1]][i]),envir=lexicon)
    obs<-get(as.character(alph[[1]][i]),envir=lexicon)
    assign(as.character(alph[[1]][i]),c(obs,Expected=1/r_n,Enrichment=obs$Observed*r_n),envir=lexicon)
  }
 
}


##################################################
# Find the observed frequency                    #
##################################################
# SYMDATA is a global variable



getObservedFreq<-function(qword,envir=lexicon)
{
  inst<-SYMDATA[,gregexpr2(qword,symPausesSense),by=geneName]
  inst<-inst[V1!=-1]
  
  Pobs<-dim(inst)[1]/lngt[,sum(V1)-length(V1)*nchar(qword)]
  assign(qword,list(Observed=Pobs,Position=inst),envir=lexicon)
}



##################################################
# Find the partitions of a word                  #
##################################################
getPart<-function(qword)
{
  Parts<-list()
  for(i in 1:(2^(nchar(qword)-1)-1))
  {
    a<-gregexpr('1',intToBin(i+2^(nchar(qword)-1)))
    sp<-unlist(a)[-1]-1
    part<-substring(qword,c(1,sp+1),c(sp,nchar(qword)))
    Parts[[i]]<-part
  }
  return(Parts)
}



##################################################
# Find the expected frequency                    #
##################################################

getExpectedFreq<-function(qword,envir=lexicon){


  if(exists(qword,envir=lexicon)){
        Pexp<-get(qword,envir=lexicon)$Expected
        if(length(Pexp)>0){
          return(Pexp)
        }
  }
    
  Parts<-getPart(qword)
  
  Pexp<-0
  for(j in 1:length(Parts)){
    
    tbl<-table(Parts[[j]])
    
    Ppart<-Pp<-1
    for(k in 1:length(tbl))
    {
      qq<-names(tbl)[k]
      Pobs<-vocab[qq,observed]$observed
      Ppart<-Ppart*Pobs^tbl[[k]]
      Pp<-Pp*ngramFreq[,frequency[ngram==nchar(qq)]]$frequency
    }
    Pexp<-Pexp+Ppart*Pp*lngt^length(Parts[[j]])
  }
  
  if(!exists(qword,envir=lexicon)){
        getObservedFreq(qword,envir=lexicon)
  }
  obs<-get(qword,envir=lexicon)
  assign(qword,c(obs,Expected=Pexp,Enrichment=obs$Observed/Pexp),envir=lexicon)
  return(Pexp)
}

# getExpectedFreq<-function(qword,envir=lexicon){
#   
#   if(exists(qword,envir=lexicon)){
#     Pfin<-get(qword,envir=lexicon)$Expected
#     if(length(Pfin)>0){
#       return(Pfin)
#     }
#   }
#   Pfin<-0
#   
#   for(i in 1:(2^(nchar(qword)-1)-1))
#   {
#     a<-gregexpr('1',intToBin(i+2^(nchar(qword)-1)))
#     sp<-unlist(a)[-1]-1
#     part<-substring(qword,c(1,sp+1),c(sp,nchar(qword)))
#     
#     PpartTot<-1
#     for(j in 1:length(part)){
#       
#       if(!exists(part[j],envir=lexicon)){
#        
#         print(paste('new calculation for ',part[j]))
#         getObservedFreq(part[j],envir=lexicon)
# 
#       }
#       Ppart<-get(part[j],envir=lexicon)$Observed
#       PpartTot<-PpartTot*Ppart
#       
#     }
#     
#     Pfin<-Pfin+PpartTot
#     
#   }
#  
#   if(!exists(qword,envir=lexicon)){
#     getObservedFreq(qword,envir=lexicon)
#   }
#   obs<-get(qword,envir=lexicon)
#   assign(qword,c(obs,Expected=Pfin,Enrichment=obs$Observed/Pfin),envir=lexicon)
#   return(Pfin)
# }
