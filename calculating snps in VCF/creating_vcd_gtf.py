file = '/PHShome/yz843/czlabwork/vcfczold/20191007_py/vcf.dat'
#vcf_reader = vcf.Reader(open(path + 'SK-30.snps.vcf', 'r'))

#vcf_reader = vcf.Reader(filename = file)

#NOTE: this version differs from uv_snp_strand.py in that it inclludes both transcript_id and gene_id in the output....

dict_gtf={}

with open("/PHShome/yz843/czlabwork/vcfczold/20191007_py/gtf1.dat", 'r') as gtf_in:

    read_gtf = gtf_in.readline()
    read_gtf = gtf_in.readlines()
    for line in read_gtf:
        line = line.strip().split('\t')
        chrom, start, end, strand, transcript_id, gene_id, gene_name = line[0], line[1], line[2], line[3], line[4], line[5], line[6]


        if chrom not in dict_gtf:
             dict_gtf[chrom] = {}
             if gene_id not in dict_gtf[chrom]:
                 dict_gtf[chrom][gene_id] = {}
                 dict_gtf[chrom][gene_id][transcript_id] = [start, end, strand, gene_name]
             else:
                 dict_gtf[chrom][gene_id][transcript_id] = [start, end, strand, gene_name]
        
        else:
            if gene_id not in dict_gtf[chrom]:
                 dict_gtf[chrom][gene_id] = {}
                 dict_gtf[chrom][gene_id][transcript_id] = [start, end, strand, gene_name]
            else:
                 dict_gtf[chrom][gene_id][transcript_id] = [start, end, strand, gene_name]


###################

with open("/PHShome/yz843/czlabwork/vcfczold/20191007_py/vcf_gtf_gname.dat", 'w') as fout:
            
    with open(file, 'r') as vcf_in:

    	read_vcf = vcf_in.readline()
    	read_vcf = vcf_in.readlines()
    	for line in read_vcf:
            line = line.strip().split('\t')
            chrom, pos, ref, alt, sample, GT, AD = line[0], line[1], line[2], line[3], line[4], line[5], line[6]

            for gene_id in dict_gtf[chrom]:
                for transcript_id in dict_gtf[chrom][gene_id]:
                    s = dict_gtf[chrom][gene_id][transcript_id][0]
                    e = dict_gtf[chrom][gene_id][transcript_id][1]
                    strand = dict_gtf[chrom][gene_id][transcript_id][2]

                    if int(pos) in range(int(s),int(e)+1):
                        if ref == 'C' and alt == 'T':
                            if strand == '+':
                                #print('yes')
                                #print(chrom, int(pos), ref, alt, strand, sample, GT, AD, gene_id, transcript_id, 'non-transcribed')
                                fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, int(pos), ref, alt, strand, sample, GT, AD, gene_id, transcript_id, 'non-transcribed', gene_name))

                            elif strand == '-':
                                #print(chrom, pos, ref, alt, strand, 'transcribed')
                                fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, int(pos), ref, alt, strand, sample, GT, AD, gene_id, transcript_id,'transcribed',gene_name))

                        elif ref == 'G' and alt == 'A': 
                            if strand == '+':
                                fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, int(pos), ref, alt, strand, sample, GT, AD, gene_id, transcript_id,'transcribed',gene_name))
                            elif strand == '-':
                                fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, int(pos), ref, alt, strand, sample, GT, AD, gene_id, transcript_id,'non-transcribed',gene_name))
                            


