rule prephasing_WHATSHAP_er:
    input:
        vcf="results/{tec}/filtration/snps_het.vcf.gz",
        bam="results/{tec}/recalibration/{tec}.recal.bam"
    output:
        v="results/{tec}/prephasing/pre_phased.vcf.gz",
        t="results/{tec}/prephasing/pre_phased.vcf.gz.tbi"
    params: 
        config["genome_fa"]
    conda:
        "../envs/whatshap.yml"
    shell:
        """ mkdir -p results/{wildcards.tec}/prephasing
        whatshap phase -o results/{wildcards.tec}/prephasing/pre_phased.vcf --reference={params} {input.vcf} {input.bam} --ignore-read-groups
        bgzip results/{wildcards.tec}/prephasing/pre_phased.vcf
        tabix {output.v} """

rule prephasing_WHATSHAP_atac:
    input:
        vcf="results/atac/filtration/snps_het.vcf.gz",
        bam="results/atac/mapping_result/atac.final.bam"
    output:
        v="results/atac/prephasing/pre_phased.vcf.gz",
        t="results/atac/prephasing/pre_phased.vcf.gz.tbi"
    params: 
        config["genome_fa"]
    conda:
        "../envs/whatshap.yml"
    shell:
        """ mkdir -p results/atac/prephasing
        whatshap phase -o results/atac/prephasing/pre_phased.vcf --reference={params} {input.vcf} {input.bam} --ignore-read-groups
        bgzip results/atac/prephasing/pre_phased.vcf
        tabix {output.v} """

rule merge_tec_vcf:
    input:
        pp=expand("results/{tec}/filtration/snps_het.vcf.gz",tec=config["tech"])
    output:
        v="results/merged_vcf/snps_het.vcf.gz",
        t="results/merged_vcf/snps_het.vcf.gz.tbi"
    conda:
        "../envs/samtools.yml"
    shell:
        """ bcftools concat {input.pp} -a | bcftools norm -d all -O z -o {output.v} 
        tabix {output.v} """

rule merge_tec_vcf_prephased:
    input:
        np=expand("results/{tec}/prephasing/pre_phased.vcf.gz",tec=config["tech"])
    output:
        vp="results/phased/pre_phased.vcf.gz",
        tp="results/phased/pre_phased.vcf.gz.tbi"
    conda:
        "../envs/samtools.yml"
    shell:
        """ bcftools concat {input.np} -a | bcftools norm -d all -O z -o {output.vp} 
        tabix {output.vp} """

rule divide_chr_VCF:
    input:
        "results/merged_vcf/variantsQC.vcf.gz"
    output:
        cv="results/merged_vcf/chr{chrom}QC.vcf.gz",
        ct="results/merged_vcf/chr{chrom}QC.vcf.gz.tbi"
    conda:
        "../envs/samtools.yml"
    shell:
        """ bcftools view {input} --regions {wildcards.chrom} -O z -o {output.cv}
        tabix {output.cv} """

##########################
######## phasing ########
##########################
        
#phase the quality controlled called SNPs with shapeit4
rule phasing_SHAPEIT4:
    input:
        vcftbi="results/phased/pre_phased.vcf.gz.tbi",
        vcf="results/phased/pre_phased.vcf.gz",
        ref_vcf=expand("{path}/{prefix}{chrom}.{extension}", path=config["path_ALLvcf"],prefix=config["prefix_ALLvcf"],extension=config["extension_ALLvcf"],  allow_missing=True)
    output: 
        temp("results/phased/chr{chrom}_phased.vcf")
    conda:
        "../envs/shapeit.yml"
    threads: config["threads_num"]
    params: 
        sex=config["sex"]
    shell:
        """  if [ {wildcards.chrom} == "X" ]
        then
            if [ {params.sex} == "male" ]
            then
                shapeit4 --input {input.vcf} --region X:10000-2781479,X:155701382-156030895 --output {output} --thread {threads} --reference {input.ref_vcf}
            else
                shapeit4 --input {input.vcf} --region X --output {output} --thread {threads} --reference {input.ref_vcf}
            fi
        else 
            shapeit4 --input {input.vcf} --region {wildcards.chrom} --output {output} --thread {threads} --reference {input.ref_vcf}
        fi """

rule bgzip_and_indexing:
    input:
        expand("results/phased/chr{chrom}_phased.vcf",chrom=config["chromosomes_to_phase"])
    output: 
        vcf="results/phased/shapeit_whatshap.vcf.gz",
        vcftbi="results/phased/shapeit_whatshap.vcf.gz.tbi"
    conda:
        "../envs/samtools.yml"
    shell:
        """ bcftools concat {input} -Oz -o {output.vcf}
            tabix {output.vcf} """