getParameters<-function()
{
	out<-list(
		samples= c("WT_delR1","WT_delR2","WT_delR3","chd1","cyc8","msi1","rlf2","rpb2","set2","snf2","spt16","swr1","tbf1","tup1","WT_rpb2","WT_spt16"),
		notification="umuteser@gmail.com",
		user="ue4",
		scales=6,
		resolution=6,
		initialDataDir="/groups/churchman/ue4/Projects/MutantScreen/Nucleosome/initialData",
		posStrandExt="_pos.bg",
		negStrandExt="_neg.bg",
		annotationFile="/groups/churchman/ue4/Projects/MutantScreen/Core/stringentAnnotation.rdata",
		projectName="Perturbation_Nucleosome",
		baseDir="/groups/churchman/ue4")
	return(out)
}