# multiscale symbolic analysis pipeline 


# put all the scripts that pipeline uses into a folder and cd into it
# cd Codes/SymbolicAnalysis/
# Rename the initial data files to the sample names

Samples= "WT\nDST1D\nRCO1D\nSET1D\nSET2D\nEAF3D"
Notification="umuteser@gmail.com"
User="ue4"
scales=6
initialDataDir="/groups/churchman/ue4/Projects/MutantScreen/NETseq"
posStrandExt="_pos.bg"
negStrandExt="_neg.bg"
annotationFile="/groups/churchman/ue4/Projects/MutantScreen/Core/stringentAnnotation.rdata"
projectName="Perturbation_NETseq"


baseDir="/groups/churchman/ue4/${projectName}"



Rscript setDirectories.r

mutList=("WT_delR1"
"WT_delR2"
"WT_delR3"
"chd1"
"cyc8"
"msi1"
"rlf2"
"rpb2"
"set2"
"snf2"
"spt16"
"swr1"
"tbf1"
"tup1"
"WT_rpb2"
"WT_spt16")



mutList=("WT_delR1")



# to split the wig files into chromosomes
count=0
while [ "x${mutList[count]}" != "x" ]
do
	   
	   f=${mutList[count]}
	awk -v f="$f" 'match($0,"chrom.*"){x=substr($0,RSTART+6);}{print > (f x ".wig");}' ${f}.wig
	count=$(( $count + 1 ))
done

### Symbolic Data Build

for i in `seq 0 6`;
do
	count=0
	while [ "x${mutList[count]}" != "x" ]
	do
	   
	   mut=${mutList[count]}


	   bsub  -R "rusage[mem=8000]"  -J convert_${mut}_${i} -W 11:57 -q short -u umuteser@gmail.com -N \
	        "bash runRcode.sh convertToSymbolicNuc.r $i $mut"
	   count=$(( $count + 1 ))
    done

done



### Lexicon Build

for i in `seq 0 6`;
do
	count=0
	while [ "x${mutList[count]}" != "x" ]
	do
	   
	   mut=${mutList[count]}


	   bsub  -R "rusage[mem=8000]"  -J buildLexicon_${mut}_${i} -W 11:57 -q short -u umuteser@gmail.com -N \
	        "bash runRcode.sh buildLexicon.r $i $mut"
	   count=$(( $count + 1 ))
    done

done


### Make the data table

for i in `seq 1 6`;
do

	   bsub  -R "rusage[mem=8000]"  -J makeDataTable_${i} -W 11:57 -q short -u umuteser@gmail.com -N \
	        "bash runRcode.sh makeDataTable.r $i"

done


### Get KL divergence

bsub  -R "rusage[mem=8000]"  -J KLdivergence -W 11:57 -q short -u umuteser@gmail.com -N \
	        "Rscript getKLdivergence.r"

### Get JS divergence
bsub  -R "rusage[mem=8000]"  -J JSdivergence -W 11:57 -q short -u umuteser@gmail.com -N \
	        "Rscript getJSdivergence.r"
