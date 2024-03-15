# ASTRA bulk pipeline

Snakemake workflow for the analysis of allele specific expression and/or chromatin accessibility from sequencing data.

<img src="./pipeline.svg">

## Input files needed

As input are needed the gzipped fastq file of *one* sample.
It is possible to use either RNA-seq and/or ATAC-seq and eventual WES fastq files and/or additional known SNPs (vcf file format).

## Reference files preparation

### Genome gtf and fasta files
Gtf and fasta reference files will be used several times thoughout the pipeline, for the alignment, variant calling, feature counting etc.

Files available on [Ensembl web site](https://www.ensembl.org/Homo_sapiens/Info/Index).
GRCh38 download example:

```bash 
wget http://ftp.ensembl.org/pub/release-107/gtf/homo_sapiens/Homo_sapiens.GRCh38.107.chr.gtf.gz
gzip -d Homo_sapiens.GRCh38.107.chr.gtf.gz
wget http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```
Add file paths to [config file](config/config.yaml) in `genome_gtf` and `genome_fa`

### Hisat2
The alignment is performed using Hisat2. Hisat2 specific index file should be created.

Prepare [hisat2](https://www.nature.com/articles/s41587-019-0201-4) index files available [here](http://daehwankimlab.github.io/hisat2/download/) for download or preparation instructions with custom reference available [here](http://daehwankimlab.github.io/hisat2/howto/#build-hgfm-index-with-snps-and-transcripts). 

For an example look at this [file](hisat_indexes.sh).

Add file paths to [config file](config/config.yaml) in `hisat_index`

### Haptree-X
Haplotype phasing is performed using 3 different tools, WhatsHap, Shapeit and Haptree-X. While the other methods are automaticaly included in the pipeline through conda environments or docker images, Haptree-X should be installed locally. In particular, you would just need to download the last released executable file, an example follows.

The Haptree-X excecutable file can be downloaded from [here](https://github.com/0xTCG/haptreex/releases).

Add path to [config file](config/config.yaml) in `haptreex_exe`

### Other reference files 

Reference file of phased variants for variant calling with Shapeit4 available [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working).

**Be careful to have the same chromosome notation in all the reference and target files (chr1 or 1).**

<details><summary>GRCh38 download example </summary>
<p> 

```bash 
for i in {1..22} X;do wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr$i.filtered.SNV_INDEL_SV_phased_panel.vcf.gz; done

for i in {1..22} X; do wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr$i.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi; done
```

</p>
</details>

## How to run

```bash
snakemake --cores [cores_number] --use-conda
```

Once the config file is fully completed you can run the entire pipeline by just typing the command above. 

## Results architecture
<details><summary> Output architecture </summary>
<p> 
  
```bash
results/
├── exome
│   ├── filtration
│   │   ├── filtered.vcf.gz
│   │   ├── filtered.vcf.gz.tbi
│   │   ├── indels.vcf.gz
│   │   ├── indels.vcf.gz.tbi
│   │   ├── indels_filtered.vcf.gz
│   │   ├── indels_filtered.vcf.gz.tbi
│   │   ├── metrics.variant_calling_detail_metrics
│   │   ├── metrics.variant_calling_summary_metrics
│   │   ├── sample.txt
│   │   ├── snps.vcf.gz
│   │   ├── snps.vcf.gz.tbi
│   │   ├── snps_filtered.vcf.gz
│   │   ├── snps_filtered.vcf.gz.tbi
│   │   ├── snps_het.vcf.gz
│   │   └── snps_het.vcf.gz.tbi
│   ├── haplotypeCaller
│   │   ├── exome.g.vcf.gz.tbi
│   │   ├── exome.vcf.gz
│   │   └── exome.vcf.gz.tbi
│   ├── prephasing
│   │   ├── pre_phased.vcf.gz
│   │   └── pre_phased.vcf.gz.tbi
│   └── recalibration
│       ├── exome.recal.bai
│       ├── exome.recal.bam
│       └── recal_data.table
├── rna
│   ├── ASEX
│   │   └── rna.table
│   ├── alignment
│   │   └── rna.splitted.bai
│   ├── filtration
│   │   ├── filtered.vcf.gz
│   │   ├── filtered.vcf.gz.tbi
│   │   ├── indels.vcf.gz
│   │   ├── indels.vcf.gz.tbi
│   │   ├── indels_filtered.vcf.gz
│   │   ├── indels_filtered.vcf.gz.tbi
│   │   ├── metrics.variant_calling_detail_metrics
│   │   ├── metrics.variant_calling_summary_metrics
│   │   ├── sample.txt
│   │   ├── snps.vcf.gz
│   │   ├── snps.vcf.gz.tbi
│   │   ├── snps_filtered.vcf.gz
│   │   ├── snps_filtered.vcf.gz.tbi
│   │   ├── snps_het.vcf.gz
│   │   └── snps_het.vcf.gz.tbi
│   ├── haplotypeCaller
│   │   ├── rna.g.vcf.gz.tbi
│   │   ├── rna.vcf.gz
│   │   └── rna.vcf.gz.tbi
│   ├── prephasing
│   │   ├── pre_phased.vcf.gz
│   │   └── pre_phased.vcf.gz.tbi
│   ├── recalibration
│   │   ├── recal_data.table
│   │   ├── rna.recal.bai
│   │   └── rna.recal.bam
│   └── transcripts_quant
│       ├── logs
│       │   └── salmon_quant.log
│       └── quant.sf
├── merged_vcf
│   ├── snps_het.vcf.gz
│   └── snps_het.vcf.gz.tbi
├── phased
│   ├── haptreex.tsv
│   ├── manual_refinment.vcf.gz
│   ├── manual_refinment.vcf.gz.tbi
│   ├── pre_phased.vcf.gz
│   ├── pre_phased.vcf.gz.tbi
│   ├── shapeit_whatshap.vcf.gz
│   └── shapeit_whatshap.vcf.gz.tbi
└── seesaw
    ├── salmon
    │   ├── logs
    │   │   └── salmon_quant.log
    │   └── quant.sf
    └── transcripts.fa
```

</p>
</details>
