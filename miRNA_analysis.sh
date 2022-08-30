# Prepare software:
# blastall: v2.2.26
# bowtie: v1.1.2
# miRDeep2: v2.0.0.7
# edgeR: v3.12.1
# targetscan: v7.0
# miranda: v3.3a

# Prepare Database:
# Rfam: 11.0
# genome.fa: Ensembl release 94(Bos taurus)
# rRNA.fa: Ensembl release 94(Bos taurus)
# 3UTR.fa: Ensembl release 94(Bos taurus)
# mirbase_animal.fa: miRBase Release 22.1
# mirbase_animal.txt: miRBase Release 22.1
# hairpin_bta.fa: miRBase Release 22.1
# mirbase_bta.txt: miRBase Release 22.1

# Step1: Filter raw data
# Perform QC on the raw data to get clean reads
perl clean_fastq.pl -i NC-1.fq.gz  -z no -b CAGATCGGAAGAG -p 33 -m 0 -c 0 -d clean_fq_dir/NC-1
perl clean_fastq.pl -i NC-2.fq.gz  -z no -b CAGATCGGAAGAG -p 33 -m 0 -c 0 -d clean_fq_dir/NC-2
perl clean_fastq.pl -i NC-3.fq.gz  -z no -b CAGATCGGAAGAG -p 33 -m 0 -c 0 -d clean_fq_dir/NC-3
perl clean_fastq.pl -i si-1.fq.gz  -z no -b CAGATCGGAAGAG -p 33 -m 0 -c 0 -d clean_fq_dir/si-1
perl clean_fastq.pl -i si-2.fq.gz  -z no -b CAGATCGGAAGAG -p 33 -m 0 -c 0 -d clean_fq_dir/si-2
perl clean_fastq.pl -i si-3.fq.gz  -z no -b CAGATCGGAAGAG -p 33 -m 0 -c 0 -d clean_fq_dir/si-3

# Merge tags 
perl merge_all_sample_clean_txt.pl clean_fq_dir sample.list clean.txt clean.fa yes animal

# Step2: Align tags to Rfam database and rRNA with blastall(v2.2.26)
blastall -p blastn -a 3 -i clean.fa -d Rfam.fasta -F F -e 0.01 -m 8 -o rfam.out
blastall -p blastn -a 3 -i clean.fa -d rRNA.fa -F F -e 0.01 -m 8 -o rRNA.out

# Step3: Align tags to genome using bowtie(v1.1.2)
bowtie genome.fa -f clean.fa -v 0 --best --strata -a -S | samtools view -bS - > genome.bam

# Step4: Identify existing mirna
# Align tags to precursor of Bos taurus in mirbase using bowtie(v1.1.2)
bowtie hairpin_bta.fa -f clean.fa -v 0 --best --strata -a --norc -S | samtools view -bS - > exist.bam
samtools view exist.bam  | perl bam2mirna.pl mirbase_bta.txt - clean.fa exist_mirna exist.out
perl sam2exist_mirna.pl exist.out > exist_mirna.list
perl split_mirna.pl clean.txt exist_mirna.list sample.list exist_mirna.txt exist_mirna.fa

# Step5: Identify known mirna
# Remove tags that can be aligned to rfam, rRNA and exist
perl remove_clean_fa_some_tag.pl clean.fa aling_dir.list known_mirna 0 > known_clean.fa

# Align the rest tags to precursor in mirbase using bowtie(v1.1.2)
bowtie mirbase_animal.fa -f known_clean.fa -v 2 --best --strata -a --norc -S | samtools view -bS - > known.bam
samtools view known.bam  | perl bam2mirna.pl mirbase_animal.txt - clean.fa known_mirna known.out
perl sam2known_mirna_miFam.pl mirbase_animal.txt known.out known_mirna.list
perl split_mirna.pl clean.txt known_mirna.list sample.list known_mirna.txt known_mirna.fa

# Step6: Identify novel mirna
# Remove tags that can be aligned to rfam, rRNA, exist and known
perl remove_clean_fa_some_tag.pl clean.fa align_dir.list known_mirna 21 > novel_clean.fa

# Use miRDeep2(v2.0.0.7) to identify novel mirna with the rest tags
samtools view genome.bam | perl bam2txt.pl clean.fa - > genome.out
perl miRDeep2.pl novel_clean.fa genome.fa genome.out none none none novel_outdir -v -g 50000
perl mirdeep_result.pl novel_outdir/output.mrd novel_outdir/result.csv pdf_dir novel_outdir
perl split_mirna.pl clean.txt novel_outdir/novel_mirna.list sample.list novel_mirna.txt novel_mirna.fa

# Step7: Calculate expression(TPM)
perl tag2ann.pl clean.txt align_dir.list annotation.xls
perl miRNA.pl align_stat -clean_txt clean.txt -annotation annotation.xls -sample_list sample.list -type annotation -outdir ann_dir
perl get_sample_count_for_exp.pl ann_dir/match_annotation.stat.xls exist_mirna,known_mirna,novel_mirna > sample_count.txt
cat exist_mirna.txt <(sed 1d known_mirna.txt) <(sed 1d novel_mirna.txt) > all_mirna.txt
perl sample_count2tpm.pl all_mirna.txt sample_count.txt > all_mirna.exp_profile.xls

# Step8: Diff analysis using edgeR(v3.12.1)
cut -f1,5-10 all_mirna.txt > NC-vs-si.count.txt
cut -f1,11-16 all_mirna.txt > NC-vs-si.fpkm.txt
Rscript edgeR.r NC-vs-si.count.txt 3,3 NC,si NC-vs-si no 2 0.05 2 yes NC-vs-si.fpkm.txt

# Step9: Predict target gene
# Intersect the result of targetscan(v7.0) and miranda(v3.3a)
# targetscan
perl -e '$/ = ">"; while (<>){chomp; next unless $_; my ($id, $seq) = split /\n/,$_,2; print "$id\t9913\t$seq"}' 3UTR.fa > utr.input
cat exist_mirna.fa known_mirna.fa novel_mirna.fa  > all_mirna.fa
perl -e '$/ = ">"; while (<>){chomp; next unless $_; my ($id, $seq) = split /\n/,$_,2; $id =~ s/ .+//; $seq =~ s/\n//; print "$id\t$seq\t9913\n"}' all_mirna.fa > mirna_input
perl targetscan_70.v2.pl mirna_input utr.input targetscan.out
# miranda
miranda all_mirna.fa 3UTR.fa -en -10 -strict > miranda.aln
# merge
perl merge_3_soft_result_to_intersection.pl none miranda.aln targetscan.out all_target.aln.xls taget_dir

#Step10: Enrichment analysis on omicshare
kegg: https://www.omicshare.com/tools/Home/Soft/pathwaygseasenior
go: https://www.omicshare.com/tools/Home/Soft/gogseasenior
