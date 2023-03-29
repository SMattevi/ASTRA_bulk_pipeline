#Allelic specific expression: gatk ASEReadCounter over the phased vcf for each cluster bam
rule ASEReadCount_a:
    input:
        vcf= "results/merged_vcf/snps_het.vcf.gz",
        bam = "results/atac/mapping_result/atac.final.bam"
    output:
        directory("results/atac/ASE{chrom}")
    conda: "../envs/gatk.yml"
    params:
        fa=config["genome_fa"]
    shell:
        """ mkdir -p "results/atac/ASE{wildcards.chrom}"
        gatk ASEReadCounter -R {params.fa}  -I {input.bam} -V {input.vcf}  -O {output}/atac.table -L {wildcards.chrom}"""

#Allelic specific expression: gatk ASEReadCounter over the phased vcf for each cluster bam
rule ASEReadCount_e:
    input:
        vcf= "results/merged_vcf/snps_het.vcf.gz",
        bam = "results/exome/recalibration/exome.recal.bam"
    output:
        directory("results/exome/ASE{chrom}")
    conda: "../envs/gatk.yml"
    params:
        fa=config["genome_fa"]
    shell:
        """ mkdir -p "results/exome/ASE{wildcards.chrom}"
        gatk ASEReadCounter -R {params.fa}  -I {input.bam} -V {input.vcf}  -O {output}/exome.table -L {wildcards.chrom}"""

rule ASEReadCount_r:
    input:
        vcf= "results/merged_vcf/snps_het.vcf.gz",
        bam = "results/rna/recalibration/rna.recal.bam"
    output:
        directory("results/rna/ASE{chrom}")
    conda: "../envs/gatk.yml"
    params:
        fa=config["genome_fa"]
    shell:
        """ mkdir -p "results/rna/ASE{wildcards.chrom}"
        gatk ASEReadCounter -R {params.fa}  -I {input.bam} -V {input.vcf}  -O {output}/rna.table -L {wildcards.chrom}"""
