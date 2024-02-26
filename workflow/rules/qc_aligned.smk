#Sort alignment result and mark duplicate 
rule mark_dup:
    input:
        "results/{tec}/alignment/{tec}.sam"
    output:
        final_markdup=temp("results/{tec}/alignment/{tec}.positionsort.bam"),
        final_markdup_index=temp("results/{tec}/alignment/{tec}.positionsort.bam.bai"),
        tmpsort=temp("results/{tec}/{tec}.sorted.bam"),
        out0=temp("results/{tec}/{tec}.positionsort0.bam"),
        fixmate=temp("results/{tec}/{tec}.fixmate.bam")
    conda: "../envs/samtools.yml"
    params: config["memory"]
    threads: config["threads_num"]
    shell: 
        """ mkdir -p tmp
        samtools sort -m {params} -T tmp/ -@ {threads} -n -o {output.tmpsort} {input}
        
        samtools fixmate -@ {threads} -m {output.tmpsort} {output.fixmate}
        
        samtools sort -m {params} -@ {threads} -T tmp/ -o {output.out0} {output.fixmate}
       
        samtools markdup -@ {threads} {output.out0} {output.final_markdup}
        samtools index -@ {threads} {output.final_markdup} 
        
        rmdir tmp """

rule SplitNCigarReads:
    input:
        "results/rna/alignment/rna.positionsort.bam"
    output:
        temp("results/rna/alignment/rna.splitted.bam")
    conda: "../envs/gatk.yml"
    params: config["genome_fa"]
    shell:
        """ gatk SplitNCigarReads \
        -R {params} \
        -I {input} \
        -O {output} """
    
rule BaseRecalibrator:
    input:
        myinput 
    output: 
        bam=temp("results/{tec}/recalibration/{tec}.positionsort.gatkgroup.bam"),
        table="results/{tec}/recalibration/recal_data.table"
    conda: "../envs/gatk.yml"
    params: 
        sites=config["known_sites"], 
        fa=config["genome_fa"],
        samp=config["sample_name"]
    threads: config["threads_num"]
    shell:
        """ if [[ {wildcards.tec} == "exome" ]]
        then
            input_real="results/exome/alignment/exome.positionsort.bam"
        else
            input_real="results/rna/alignment/rna.splitted.bam"
        fi
        echo $input_real
        gatk AddOrReplaceReadGroups -I $input_real -O {output.bam} -RGLB DNA -RGPL ILLUMINA -RGPU {wildcards.tec} -RGSM {wildcards.tec}_{params.samp} -VALIDATION_STRINGENCY SILENT
        gatk BaseRecalibrator -I {output.bam} --known-sites {params.sites} -R {params.fa} -O {output.table} """

rule index:
    input:
        "results/{tec}/recalibration/{tec}.positionsort.gatkgroup.bam"
    output:
        temp("results/{tec}/recalibration/{tec}.positionsort.gatkgroup.bam.bai")
    conda: "../envs/samtools.yml"
    threads: config["threads_num"]
    shell: "samtools index -@ {threads} {input}"

rule BQSR:
    input:
        bam="results/{tec}/recalibration/{tec}.positionsort.gatkgroup.bam",
        bai="results/{tec}/recalibration/{tec}.positionsort.gatkgroup.bam.bai",
        table="results/{tec}/recalibration/recal_data.table"
    output:
        "results/{tec}/recalibration/{tec}.recal.bam"
    conda: "../envs/gatk.yml"
    shell: 
        "gatk ApplyBQSR --bqsr {input.table} -I {input.bam} -O {output}"

rule featureCount:
    input:
        R1=expand("{path}/{sample}_R1.fastq.gz", path=config["path_rna"],sample=config["fastqs_rna"],sep=","),
        R2=expand("{path}/{sample}_R2.fastq.gz", path=config["path_rna"],sample=config["fastqs_rna"],sep=",")
    output:
        "results/rna/transcripts_quant/quant.sf"
    conda:
        "../envs/salmon.yml"
    threads: config["threads_num"]
    params:
        index=config["salmon_index"],
        outdir="results/rna/transcripts_quant"
    threads: config["threads_num"]
    shell:
        """ R1=$(echo {input.R1})
        #R1new=$(echo $R1 | sed 's/ /,/g ')
        R2=$(echo {input.R2})
        #R2new=$(echo $R2 | sed 's/ /,/g ')
        salmon quant -i {params.index} \
            -l A -1 $R1 \
            -2  $R2 \
            --validateMappings -o {params.outdir} \
            -p {threads}"""

#QC of the aligned bam file-> params -q 30 -> used to extract cells
rule QC_bam_atac:
    input:
        "results/atac/alignment/atac.positionsort.bam"
    threads: config["threads_num"]
    conda: "../envs/samtools.yml"
    output:
        temp("results/atac/mapping_result/atac.positionsort.MAPQ20.bam")
    shell:
        """ samtools view -f 0x2 -b -h -q 20 -@ {threads} {input} -o {output}
        samtools index -@ {threads} {output} 
        bash workflow/scripts/createsummary.sh {input} {output} {threads} atac """

#modify bam header-> add read group needed for ASEReadCounter
rule GATK_AddorRep:
    input:
        "results/atac/mapping_result/atac.positionsort.MAPQ20.bam"
    output:
        "results/atac/mapping_result/atac.final.bam"
    conda: "../envs/gatk.yml"
    params:
        samp=config["sample_name"]
    shell:
        """ gatk AddOrReplaceReadGroups -I {input} -O {output} -RGLB DNA -RGPL ILLUMINA -RGPU atac -RGSM atac_{params.samp} -VALIDATION_STRINGENCY SILENT """

rule index_bam:
    input:
        bam="results/atac/mapping_result/atac.final.bam"
    output:
        "results/atac/mapping_result/atac.final.bam.bai"
    conda: "../envs/samtools.yml"
    shell:
        """ samtools index {input}"""
