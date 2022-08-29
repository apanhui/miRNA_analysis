#!/bin/sh
#------------------------------------
# Prepare software
# blastall: v2.2.26
# bowtie: v1.1.2
# miRDeep2: v2.0.0.7
# edgeR: v3.12.1
# targetscan: v7.0
# miranda: v3.3a
#------------------------------------

#Step1: Filter raw data
#Perform QC on the raw data to get clean reads
perl clean_fastq.pl -i NC-1.fq.gz  -z no -b CAGATCGGAAGAG -p 33 -m 0 -c 0 -d NC-1
perl clean_fastq.pl -i NC-2.fq.gz  -z no -b CAGATCGGAAGAG -p 33 -m 0 -c 0 -d NC-2
perl clean_fastq.pl -i NC-3.fq.gz  -z no -b CAGATCGGAAGAG -p 33 -m 0 -c 0 -d NC-3
perl clean_fastq.pl -i si-1.fq.gz  -z no -b CAGATCGGAAGAG -p 33 -m 0 -c 0 -d si-1
perl clean_fastq.pl -i si-2.fq.gz  -z no -b CAGATCGGAAGAG -p 33 -m 0 -c 0 -d si-2
perl clean_fastq.pl -i si-3.fq.gz  -z no -b CAGATCGGAAGAG -p 33 -m 0 -c 0 -d si-3

#Step2: Remove noncoding RNA
#Align GeneBank and Rfam database with blastall(v2.2.26) to remove rRNA, scRNA, snoRNA, snRNA and tRNA
blastall -p blastn -a 3 -i clean.fa -d Rfam.fasta -F F -e 0.01 -m 8 -o rfam.out
blastall -p blastn -a 3 -i clean.fa -d ncgb.fa -F F -e 0.01 -m 8 -o ncgb.out

#Step3: Align tags to genome using bowtie(v1.1.2)
bowtie genome.fa -f clean.fa -v 0 --best --strata -a -S | samtools view -bS - > genome.bam

#Step4: Identify existing mirna
#Align tags to precursor of Bos taurus in mirbase using bowtie(v1.1.2)
bowtie hairpin_bta.fa -f clean.fa -v 0 --best --strata -a --norc -S | samtools view -bS - > exist.bam

#Step5: Identify known mirna
#Align the rest tags to precursor in mirbase using bowtie(v1.1.2)
bowtie mirbase_animal.fa -f clean.fa -v 2 --best --strata -a --norc -S | samtools view -bS - > known.bam

#Step6: Identify novel mirna
#Use miRDeep2(v2.0.0.7) to identify novel mirna with the rest tags
perl miRDeep2.pl clean.fa genome.fa aln.txt none none none novel_mirna -v -g 50000

#Step7: Calculate expression(TPM)
perl sample_count2tpm.pl count.matrix sample_count.txt > fpkm.matrix

#Step8: Diff analysis using edgeR(v3.12.1)
library(edgeR)
count = read.table("count.txt", header = T, sep = "\t")
group = factor(c(1,1,1,2,2,2 ))
y = DGEList(counts, group)
y = calcNormFactors(y)
y = estimateCommonDisp(y)
y = estimateTagwiseDisp(y)
et = exactTest(y)
top <- topTags(et)
write.table(top, file = "edgeR.xls", quote = F, sep ="\t")

#Step9: Predict target gene
#Intersect the result of targetscan(v7.0) and miranda(v3.3a)
perl targetscan_70.v2.pl mirna_input utr.input targetscan.out
miranda mirna.fa 3UTR.fa -en -10 -strict > miranda.aln

#Step10: Enrichment analysis on omicshare
# kegg: https://www.omicshare.com/tools/Home/Soft/pathwaygseasenior
# go: https://www.omicshare.com/tools/Home/Soft/gogseasenior
