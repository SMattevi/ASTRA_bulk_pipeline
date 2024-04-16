#Allelic specific expression: gatk ASEReadCounter over the phased vcf for each cluster bam
rule ASEReadCount_a:
    input:
        vcf= "results_{sample_id}/merged_vcf/snps_het.vcf.gz",
        bam = "results_{sample_id}/atac/mapping_result/atac.final.bam"
    output:
        directory("results_{sample_id}/atac/ASE{chrom}")
    conda: "../envs/gatk.yml"
    params:
        fa=config["genome_fa"],
        sampleid="{sample_id}"
    shell:
        """ mkdir -p "results_{params.sampleid}/atac/ASE{wildcards.chrom}"
        gatk ASEReadCounter -R {params.fa}  -I {input.bam} -V {input.vcf}  -O {output}/atac.table -L {wildcards.chrom}"""

#Allelic specific expression: gatk ASEReadCounter over the phased vcf for each cluster bam
rule ASEReadCount_e:
    input:
        vcf= "results_{sample_id}/merged_vcf/snps_het.vcf.gz",
        bam = "results_{sample_id}/exome/recalibration/exome.recal.bam"
    output:
        directory("results_{sample_id}/exome/ASE{chrom}")
    conda: "../envs/gatk.yml"
    params:
        fa=config["genome_fa"],
        sampleid="{sample_id}"
    shell:
        """ mkdir -p "results_{params.sampleid}/exome/ASE{wildcards.chrom}"
        gatk ASEReadCounter -R {params.fa}  -I {input.bam} -V {input.vcf}  -O {output}/exome.table -L {wildcards.chrom}"""

rule ASEReadCount_r:
    input:
        vcf= "results_{sample_id}/merged_vcf/snps_het.vcf.gz",
        bam = "results_{sample_id}/rna/recalibration/rna.recal.bam"
    output:
        directory("results_{sample_id}/rna/ASE{chrom}")
    conda: "../envs/gatk.yml"
    params:
        fa=config["genome_fa"],
        sampleid="{sample_id}"
    shell:
        """ mkdir -p "results_{params.sampleid}/rna/ASE{wildcards.chrom}"
        gatk ASEReadCounter -R {params.fa}  -I {input.bam} -V {input.vcf}  -O {output}/rna.table -L {wildcards.chrom}"""
