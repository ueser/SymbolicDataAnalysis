setDirectories<-function()
{
	source('getParameters.r')
	out<-getParameters()
	
	dir.create(file.path(out$baseDir,out$projectName), showWarnings = FALSE)
	setwd(file.path(out$baseDir,out$projectName))
	dir.create("symbolicData")
	dir.create("lexicon")
	dir.create("corpus")
	dir.create("KLdivergence")
	f<-file("Parameters.txt","w")
	cat("Date:",Sys.Date(),"\n",file=f,sep="\t")
	for( i in 1:length(out))
	{
		cat(names(out[i]),"-",out[[i]],"\n",file=f,sep="\t")
	}
	close(f)
}

setDirectories()