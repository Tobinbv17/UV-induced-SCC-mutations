#!/bin/bash
#BSUB -J starrd1905
#BSUB -o output/starrd1905-%J.out
#BSUB -e output/starrd1905-%J.err


####----- filename array

UNMAPPED='/PHShome/yz843/czlabwork/STAR0923/fastq'
declare -A filenames
pushd ${UNMAPPED}
for filename in *
do
	seq=${filename:12:17}
	filenames[$seq]+="${filename},"
	# echo ${filenames[$seq]}
done
popd

for seq in "${!filenames[@]}"
do
	filenamesBySeq=${filenames[$seq]}
	readFilesIn=
	for part in 1 2
	do
		tmp=
		for filename in $(echo ${filenamesBySeq} | sed "s/,/ /g")
		do
			if [[ ${filename:39:1} = "${part}" ]] ; then 
		    	tmp=${tmp},${UNMAPPED}/${filename}
			fi
		done
	readFilesIn="${readFilesIn} ${tmp:1}"
  	done
	./starnormal.sh "${seq}" "${readFilesIn:1}" #"${outSAMattrRGline}" #STAR alignment
done

####----- star alignment
#!/bin/bash


echo '---Running STAR NOW---'
module load star/2.5.3
echo '----Module loaded running now----'

OUT='./output_star/'
OUTBAM='./output_bam/'
TMP='tmp'
STAR_OUT='./star_output/'

mkdir ${OUT}
mkdir ${OUTBAM}
mkdir ${STAR_OUT}

COUNT='./counts' 
SUM='./summary'
mkdir ${COUNT}
mkdir ${SUM}


SMP=$1


# Run the actual star aligner
STAR --outTmpDir ${TMP} --genomeDir /PHShome/yz843/czlabwork/STAR0923/genomedir \
--outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${STAR_OUT} \
--readFilesCommand zcat --readFilesIn $2 --runThreadN 32 \
--sjdbGTFfile /PHShome/yz843/czlabwork/annotation/gencode.vM20.annotation.gtf


cp ${STAR_OUT}Aligned.sortedByCoord.out.bam ${OUTBAM}${SMP}.sortedByCoord.bam
cp ${STAR_OUT}SJ.out.tab ${OUT}${SMP}.SJ.out.tab
cp ${STAR_OUT}Log.final.out ${OUT}${SMP}.Log.final.out
cp ${STAR_OUT}Log.out ${OUT}${SMP}.Log.out
cp ${STAR_OUT}Log.progress.out ${OUT}${SMP}.Log.progress.out
rm -r ${STAR_OUT}
rm -r ${TMP}







