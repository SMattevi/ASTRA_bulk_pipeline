#Sort alignment result and mark duplicate 
rule mark_dup:
    input:
        "results_{sample_id}/{tec}/alignment/{tec}.sam"
    output:
        final_markdup=temp("results_{sample_id}/{tec}/alignment/{tec}.positionsort.bam"),
        final_markdup_index=temp("results_{sample_id}/{tec}/alignment/{tec}.positionsort.bam.bai"),
        tmpsort=temp("results_{sample_id}/{tec}/{tec}.sorted.bam"),
        out0=temp("results_{sample_id}/{tec}/{tec}.positionsort0.bam"),
        fixmate=temp("results_{sample_id}/{tec}/{tec}.fixmate.bam")
    conda: "../envs/samtools.yml"
    params: config["memory"]
    threads: config["threads_num"]
    log: "logs/{sample_id}/{tec}_alignmentQC.log"
    shell: 
        """ mkdir -p tmp
        samtools sort -m {params} -T tmp/ -@ {threads} -n -o {output.tmpsort} {input} 2> {log}
        
        samtools fixmate -@ {threads} -m {output.tmpsort} {output.fixmate} 2>> {log}
        
        samtools sort -m {params} -@ {threads} -T tmp/ -o {output.out0} {output.fixmate} 2>> {log}
       
        samtools markdup -@ {threads} {output.out0} {output.final_markdup} 2>> {log}
        samtools index -@ {threads} {output.final_markdup}  2>> {log}
        
        rmdir tmp """

rule SplitNCigarReads:
    input:
        "results_{sample_id}/rna/alignment/rna.positionsort.bam"
    output:
        temp("results_{sample_id}/rna/alignment/rna.splitted.bam")
    conda: "../envs/gatk.yml"
    params: config["genome_fa"]
    log: "logs/{sample_id}/rna_SplitNCigarReads.log"
    shell:
        """ gatk SplitNCigarReads \
        -R {params} \
        -I {input} \
        -O {output} \
        2> {log}"""
    
rule BaseRecalibrator:
    input:
        myinput 
    output: 
        bam=temp("results_{sample_id}/{tec}/recalibration/{tec}.positionsort.gatkgroup.bam"),
        table="results_{sample_id}/{tec}/recalibration/recal_data.table"
    conda: "../envs/gatk.yml"
    params: 
        sites=config["known_sites"], 
        fa=config["genome_fa"],
        samp="{sample_id}"
    threads: config["threads_num"]
    log: "logs/{sample_id}/{tec}_BaseRecalibrator.log"
    shell:
        """ if [[ {wildcards.tec} == "exome" ]]
        then
            input_real="results_{params.samp}/exome/alignment/exome.positionsort.bam"
        else
            input_real="results_{params.samp}/rna/alignment/rna.splitted.bam"
        fi
        echo $input_real
        gatk AddOrReplaceReadGroups -I $input_real -O {output.bam} -RGLB DNA -RGPL ILLUMINA -RGPU {wildcards.tec} -RGSM {wildcards.tec}_{params.samp} -VALIDATION_STRINGENCY SILENT 2> {log}
        gatk BaseRecalibrator -I {output.bam} --known-sites {params.sites} -R {params.fa} -O {output.table} 2>> {log} """

rule index:
    input:
        "results_{sample_id}/{tec}/recalibration/{tec}.positionsort.gatkgroup.bam"
    output:
        temp("results_{sample_id}/{tec}/recalibration/{tec}.positionsort.gatkgroup.bam.bai")
    conda: "../envs/samtools.yml"
    threads: config["threads_num"]
    log: "logs/{sample_id}/{tec}_BaseRecalibrator.log"
    shell: "samtools index -@ {threads} {input} 2>> {log}"

rule BQSR:
    input:
        bam="results_{sample_id}/{tec}/recalibration/{tec}.positionsort.gatkgroup.bam",
        bai="results_{sample_id}/{tec}/recalibration/{tec}.positionsort.gatkgroup.bam.bai",
        table="results_{sample_id}/{tec}/recalibration/recal_data.table"
    output:
        "results_{sample_id}/{tec}/recalibration/{tec}.recal.bam"
    conda: "../envs/gatk.yml"
    log: "logs/{sample_id}/{tec}_BQSR.log"
    shell: 
        "gatk ApplyBQSR --bqsr {input.table} -I {input.bam} -O {output} 2> {log}"

rule featureCount:
    input:
        R1=expand("{path}/{sample}_R1.fastq.gz", path=config["path_rna"],sample=config["fastqs_rna"],sep=","),
        R2=expand("{path}/{sample}_R2.fastq.gz", path=config["path_rna"],sample=config["fastqs_rna"],sep=",")
        # R1=lambda wildcards: find_fastqs(
        #     path=config["SAMPLES"][wildcards.sample_id]["path_rna"],
        #     prefixes=config["SAMPLES"][wildcards.sample_id]["fastqs_rna"],
        #     read_num=1
        # ),
        # R2=lambda wildcards: find_fastqs(
        #     path=config["SAMPLES"][wildcards.sample_id]["path_rna"],
        #     prefixes=config["SAMPLES"][wildcards.sample_id]["fastqs_rna"],
        #     read_num=2
        # )
    output:
        outfile="results_{sample_id}/rna/transcripts_quant/quant.sf"
    conda:
        "../envs/salmon.yml"
    threads: config["threads_num"]
    params:
        index=config["salmon_index"],
        outdir=subpath(output.outfile, parent=True)
    threads: config["threads_num"]
    log: "logs/{sample_id}/rna_salmonQuant.log"
    shell:
        """ R1=$(echo {input.R1})
        #R1new=$(echo $R1 | sed 's/ /,/g ')
        R2=$(echo {input.R2})
        #R2new=$(echo $R2 | sed 's/ /,/g ')
        salmon quant -i {params.index} \
            -l A -1 $R1 \
            -2  $R2 \
            --validateMappings -o {params.outdir} \
            -p {threads} 2> {log} """

#QC of the aligned bam file-> params -q 30 -> used to extract cells
rule QC_bam_atac:
    input:
        "results_{sample_id}/atac/alignment/atac.positionsort.bam"
    threads: config["threads_num"]
    conda: "../envs/samtools.yml"
    output:
        temp("results_{sample_id}/atac/mapping_result/atac.positionsort.MAPQ20.bam")
    log: "logs/{sample_id}/atac_alignmentQC_summ.log"
    shell:
        """ samtools view -f 0x2 -b -h -q 20 -@ {threads} {input} -o {output} 2> {log}
        samtools index -@ {threads} {output} 2>> {log}
        bash workflow/scripts/createsummary.sh {input} {output} {threads} atac {wildcards.sample_id} 2>> {log}"""

#modify bam header-> add read group needed for ASEReadCounter
rule GATK_AddorRep:
    input:
        "results_{sample_id}/atac/mapping_result/atac.positionsort.MAPQ20.bam"
    output:
        "results_{sample_id}/atac/mapping_result/atac.final.bam"
    conda: "../envs/gatk.yml"
    params:
        samp="{sample_id}"
    log: "logs/{sample_id}/atac_addorrep.log"
    shell:
        """ gatk AddOrReplaceReadGroups -I {input} -O {output} -RGLB DNA -RGPL ILLUMINA -RGPU atac -RGSM atac_{params.samp} -VALIDATION_STRINGENCY SILENT 2> {log}"""

rule index_bam:
    input:
        bam="results_{sample_id}/atac/mapping_result/atac.final.bam"
    output:
        "results_{sample_id}/atac/mapping_result/atac.final.bam.bai"
    conda: "../envs/samtools.yml"
    log: "logs/{sample_id}/atac_addorrep.log"
    shell:
        """ samtools index {input} 2>> {log}"""