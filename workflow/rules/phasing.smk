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

rule phasing_haptreex:
    input:
        vcf="results/merged_vcf/snps_het.vcf.gz",
        bam="results/rna/recalibration/rna.recal.bam"
    output:
        "results/phased/haptreex.tsv"
    params:
        gtffile=config["genome_gtf"],
        haptreex=config["haptreex_exe"],
        tec=config["tech"]
    shell:
        """ bcftools view {input.vcf} -Ov -o temp.vcf
        if [[ "{params.tec}" == *exome* ]]
        then
            {params.haptreex} -v temp.vcf -r results/rna/recalibration/rna.recal.bam -g {params.gtffile} -o {output} -d results/exome/recalibration/exome.recal.bam
        elif [[ "{params.tec}" == *atac* ]]
        then
            {params.haptreex} -v temp.vcf -r results/rna/recalibration/rna.recal.bam -g {params.gtffile} -o {output} -d results/atac/recalibration/atac.recal.bam
        else
            {params.haptreex} -v temp.vcf -r results/rna/recalibration/rna.recal.bam -g {params.gtffile} -o {output}
        fi
        rm temp.vcf """

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

rule manual_phasing:
    input:
        haptreex="results/phased/haptreex.tsv",
        shapeit="results/phased/shapeit_whatshap.vcf.gz",
        not_phased="results/merged_vcf/snps_het.vcf.gz",
        ase="results/rna/ASE{chrom}"
    output:
        temp("results/phased/manual_phasing{chrom}.tsv")
    params:
        sample=config["sample_name"]
    shell:
        """ Rscript --vanilla workflow/scripts/iterative.R -c {wildcards.chrom} -v {input.not_phased} \
        -s {input.shapeit} \
        -x {input.haptreex} \
        -a {input.ase}/rna.table \
        -n {params.sample} \
        -o {output}
        """

rule tsv_to_vcf:
    input:
        "results/phased/manual_phasing{chrom}.tsv"
    output: 
        temp("results/phased/manual_refinment{chrom}.vcf")
    conda:
        "../envs/samtools.yml"
    params:
        sample=config["sample_name"]
    shell:
        """ day=$(date "+%d/%m/%4Y %T")
        echo "##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##fileDate=$day
##source=manual_phasing
##contig=<ID={wildcards.chrom}>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count">
##INFO=<ID=CM,Number=A,Type=Float,Description="Interpolated cM position">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotypes">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{params.sample}" >  {output}

        less {input} | awk '{{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t.\t.\tAC=1;AF=0.5\tGT\t"$6"|"$7}}'>> {output}
        """

rule bgzip_and_indexing_man:
    input:
        expand("results/phased/manual_refinment{chrom}.vcf",chrom=config["chromosomes_to_phase"])
    output: 
        vcf="results/phased/manual_refinment.vcf.gz",
        vcftbi="results/phased/manual_refinment.vcf.gz.tbi"
    conda:
        "../envs/samtools.yml"
    shell:
        """ bcftools concat {input} -Oz -o {output.vcf}
            tabix {output.vcf} """