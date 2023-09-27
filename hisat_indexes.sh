#!/bin/bash

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/snp151Common.txt.gz
gzip -d snp151Common.txt.gz
awk 'BEGIN{OFS="\t"} {if($2 ~ /^chr/) {$2 = substr($2, 4)}; if($2 == "M") {$2 = "MT"} print}' snp151Common.txt > snp151Common.txt.ensembl

wget http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.chr.gtf.gz
gzip -d Homo_sapiens.GRCh38.107.chr.gtf.gz

hisat2_extract_splice_sites.py Homo_sapiens.GRCh38.107.chr.gtf > genome.ss
hisat2_extract_exons.py Homo_sapiens.GRCh38.107.chr.gtf > genome.exon

wget http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

hisat2_extract_snps_haplotypes_UCSC.py Homo_sapiens.GRCh38.dna.primary_assembly.fa snp151Common.txt.ensembl genome
hisat2-build  -p 16 --snp genome.snp --haplotype genome.haplotype --exon genome.exon --ss genome.ss Homo_sapiens.GRCh38.dna.primary_assembly.fa genome_snp_tran
