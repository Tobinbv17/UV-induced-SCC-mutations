########## Creating all snp map to gtf



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
        chrom, start, end, strand, transcript_id, gene_id = line[0], line[1], line[2], line[3], line[4], line[5]


        if chrom not in dict_gtf:
             dict_gtf[chrom] = {}
             if gene_id not in dict_gtf[chrom]:
                 dict_gtf[chrom][gene_id] = {}
                 dict_gtf[chrom][gene_id][transcript_id] = [start, end, strand]
             else:
                 dict_gtf[chrom][gene_id][transcript_id] = [start, end, strand]
        
        else:
            if gene_id not in dict_gtf[chrom]:
                 dict_gtf[chrom][gene_id] = {}
                 dict_gtf[chrom][gene_id][transcript_id] = [start, end, strand]
            else:
                 dict_gtf[chrom][gene_id][transcript_id] = [start, end, strand]


###################

with open("/PHShome/yz843/czlabwork/vcfczold/20191007_py/vcf_gtf_all.dat", 'w') as fout:
            
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
                        fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, int(pos), ref, alt, strand, sample, GT, AD, gene_id, transcript_id))
                                                                                     
                                                                                     
                                                                                     
