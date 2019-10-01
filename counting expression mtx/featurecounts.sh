echo '---Running subread NOW---'
module load subread/1.5.0
echo '----Module loaded running now----'

featureCounts -T 5 -t exon -g gene_id -a /PHShome/yz843/czlabwork/annotation/gencode.vM20.annotation.gtf \
-o countmatrix.txt /PHShome/yz843/czlabwork/STAR0923/output_bam/*
