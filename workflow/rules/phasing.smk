rule prephasing_WHATSHAP_e:
    input:
        vcf="results_{sample_id}/rna/filtration/snps_het.vcf.gz",
        bam="results_{sample_id}/rna/recalibration/rna.recal.bam"
    output:
        v="results_{sample_id}/rna/prephasing/pre_phased.vcf.gz",
        t="results_{sample_id}/rna/prephasing/pre_phased.vcf.gz.tbi"
    params: 
        genome=config["genome_fa"],
        outdir="results_{sample_id}"
    conda:
        "../envs/whatshap.yml"
    shell:
        """ mkdir -p {params.outdir}/rna/prephasing
        whatshap phase -o {params.outdir}/rna/prephasing/pre_phased.vcf --reference={params.genome} {input.vcf} {input.bam} --ignore-read-groups
        bgzip {params.outdir}/rna/prephasing/pre_phased.vcf
        tabix {output.v} """

rule prephasing_WHATSHAP_r:
    input:
        vcf="results_{sample_id}/exome/filtration/snps_het.vcf.gz",
        bam="results_{sample_id}/exome/recalibration/exome.recal.bam"
    output:
        v="results_{sample_id}/exome/prephasing/pre_phased.vcf.gz",
        t="results_{sample_id}/exome/prephasing/pre_phased.vcf.gz.tbi"
    params: 
        genome=config["genome_fa"],
        outdir="results_{sample_id}"
    conda:
        "../envs/whatshap.yml"
    shell:
        """ mkdir -p {params.outdir}/exome/prephasing
        whatshap phase -o {params.outdir}/exome/prephasing/pre_phased.vcf --reference={params.genome} {input.vcf} {input.bam} --ignore-read-groups
        bgzip {params.outdir}/exome/prephasing/pre_phased.vcf
        tabix {output.v} """

rule prephasing_WHATSHAP_atac:
    input:
        vcf="results_{sample_id}/atac/filtration/snps_het.vcf.gz",
        bam="results_{sample_id}/atac/mapping_result/atac.final.bam"
    output:
        v="results_{sample_id}/atac/prephasing/pre_phased.vcf.gz",
        t="results_{sample_id}/atac/prephasing/pre_phased.vcf.gz.tbi"
    params: 
        fa=config["genome_fa"],
        sampleid="{sample_id}"
    conda:
        "../envs/whatshap.yml"
    shell:
        """ mkdir -p results_{params.sampleid}/atac/prephasing
        whatshap phase -o results_{params.sampleid}/atac/prephasing/pre_phased.vcf --reference={params.fa} {input.vcf} {input.bam} --ignore-read-groups
        bgzip results_{params.sampleid}/atac/prephasing/pre_phased.vcf
        tabix {output.v} """

rule merge_tec_vcf:
    input:
        pp=expand("results_{sample_id}/{tec}/filtration/snps_het.vcf.gz",tec=config["tech"],sample_id=config["sample_name"])
    output:
        v="results_{sample_id}/merged_vcf/snps_het.vcf.gz",
        t="results_{sample_id}/merged_vcf/snps_het.vcf.gz.tbi"
    conda:
        "../envs/samtools.yml"
    shell:
        """ bcftools concat {input.pp} -a | bcftools norm -d all -O z -o {output.v} 
        tabix {output.v} """

rule merge_tec_vcf_prephased:
    input:
        np=expand("results_{sample_id}/{tec}/prephasing/pre_phased.vcf.gz",tec=config["tech"],sample_id=config["sample_name"])
    output:
        vp="results_{sample_id}/phased/pre_phased.vcf.gz",
        tp="results_{sample_id}/phased/pre_phased.vcf.gz.tbi"
    conda:
        "../envs/samtools.yml"
    shell:
        """ bcftools concat {input.np} -a | bcftools norm -d all -O z -o {output.vp} 
        tabix {output.vp} """

rule divide_chr_VCF:
    input:
        "results_{sample_id}/merged_vcf/variantsQC.vcf.gz"
    output:
        cv="results_{sample_id}/merged_vcf/{chrom}QC.vcf.gz",
        ct="results_{sample_id}/merged_vcf/{chrom}QC.vcf.gz.tbi"
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
        vcftbi="results_{sample_id}/phased/pre_phased.vcf.gz.tbi",
        vcf="results_{sample_id}/phased/pre_phased.vcf.gz",
        ref_vcf=expand("{path}/{prefix}{chrom}.{extension}", path=config["path_ALLvcf"],prefix=config["prefix_ALLvcf"],extension=config["extension_ALLvcf"],  allow_missing=True)
    output: 
        temp("results_{sample_id}/phased/{chrom}_phased.vcf")
    conda:
        "../envs/shapeit.yml"
    threads: config["threads_num"]
    params: 
        sex=config["sex"] #config["SAMPLES"][wildcards.sample_id]["sex"]
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
        vcf="results_{sample_id}/merged_vcf/snps_het.vcf.gz",
        bam="results_{sample_id}/rna/recalibration/rna.recal.bam"
    output:
        "results_{sample_id}/phased/haptreex.tsv"
    params:
        gtffile=config["genome_gtf"],
        haptreex=config["haptreex_exe"],
        tec=config["tech"], #config["SAMPLES"][wildcards.sample_id]["tech"]
        htslib=config["htslib_path"],
        sampleid="{sample_id}"
    conda: "../envs/samtools.yml"
    shell:
        """ bcftools view {input.vcf} -Ov -o temp.vcf
        export LD_LIBRARY_PATH={params.htslib}
        if [[ "{params.tec}" == *exome* ]]
        then
            {params.haptreex} -v temp.vcf -r results_{params.sampleid}/rna/recalibration/rna.recal.bam -g {params.gtffile} -o {output} -d results_{params.sampleid}/exome/recalibration/exome.recal.bam
        elif [[ "{params.tec}" == *atac* ]]
        then
            {params.haptreex} -v temp.vcf -r results_{params.sampleid}/rna/recalibration/rna.recal.bam -g {params.gtffile} -o {output} -d results_{params.sampleid}/atac/recalibration/atac.recal.bam
        else
            {params.haptreex} -v temp.vcf -r results_{params.sampleid}/rna/recalibration/rna.recal.bam -g {params.gtffile} -o {output}
        fi
        rm temp.vcf """

rule bgzip_and_indexing:
    input:
        expand("results_{sample_id}/phased/{chrom}_phased.vcf",chrom=config["chromosomes_to_phase"],sample_id=config["sample_name"])
    output: 
        vcf="results_{sample_id}/phased/shapeit_whatshap.vcf.gz",
        vcftbi="results_{sample_id}/phased/shapeit_whatshap.vcf.gz.tbi"
    conda:
        "../envs/samtools.yml"
    shell:
        """ bcftools concat {input} -Oz -o {output.vcf}
            tabix {output.vcf} """

rule manual_phasing:
    input:
        haptreex="results_{sample_id}/phased/haptreex.tsv",
        shapeit="results_{sample_id}/phased/shapeit_whatshap.vcf.gz",
        not_phased="results_{sample_id}/merged_vcf/snps_het.vcf.gz",
        ase="results_{sample_id}/rna/ASE{chrom}"
    output:
        temp("results_{sample_id}/phased/manual_phasing{chrom}.tsv")
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
        "results_{sample_id}/phased/manual_phasing{chrom}.tsv"
    output: 
        temp("results_{sample_id}/phased/manual_refinment{chrom}.vcf")
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
        expand("results_{sample_id}/phased/manual_refinment{chrom}.vcf",chrom=config["chromosomes_to_phase"],sample_id=config["sample_name"])
        #lambda wildcards: 
            # expand("results_{sample_id}/phased/manual_refinment{chrom}.vcf",
            #        chrom=config["chromosomes_to_phase"], # GLOBAL CHROMS LIST
            #        sample_id=wildcards.sample_id)
    output: 
        vcf="results_{sample_id}/phased/manual_refinment.vcf.gz",
        vcftbi="results_{sample_id}/phased/manual_refinment.vcf.gz.tbi"
    conda:
        "../envs/samtools.yml"
    shell:
        """ bcftools concat {input} -Oz -o {output.vcf}
            tabix {output.vcf} """