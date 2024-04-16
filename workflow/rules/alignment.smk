rule alignment_exome:
    input:
        R1=expand("{path}/{sample}_1.fastq.gz", path=config["path_exome"],sample=config["fastqs_exome"],sep=","),
        R2=expand("{path}/{sample}_2.fastq.gz", path=config["path_exome"],sample=config["fastqs_exome"],sep=",")
    threads: 
        config["threads_num"]
    params: 
        index_gen = config["hisat_index"]
    output:
        temp("results_{sample_id}/exome/alignment/exome.sam")
    log: "logs/exome_{sample_id}.log"
    conda:
        "../envs/hisat.yml"
    shell:
        """ R1=$(echo {input.R1})
            R1new=$(echo $R1 | sed 's/ /,/g ')
            R2=$(echo {input.R2})
            R2new=$(echo $R2 | sed 's/ /,/g ')
            hisat2 -x {params.index_gen} \
            -1 $R1new \
            -2 $R2new \
            -S {output} \
            --no-spliced-alignment \
            -p {threads} \
            --summary-file {log} """

rule alignment_rna:
    input:
        R1=expand("{path}/{sample}_R1.fastq.gz", path=config["path_rna"],sample=config["fastqs_rna"],sep=","),
        R2=expand("{path}/{sample}_R2.fastq.gz", path=config["path_rna"],sample=config["fastqs_rna"],sep=",")
    threads: 
        config["threads_num"]
    params: 
        index_gen = config["hisat_index"],
    output:
        temp("results_{sample_id}/rna/alignment/rna.sam")
    log: "logs/rna_{sample_id}.log"
    conda:
        "../envs/hisat.yml"
    shell:
        """ R1=$(echo {input.R1})
            R1new=$(echo $R1 | sed 's/ /,/g ')
            R2=$(echo {input.R2})
            R2new=$(echo $R2 | sed 's/ /,/g ')
            hisat2 -x {params.index_gen} \
            -1 $R1new \
            -2 $R2new \
            -S {output} \
            -p {threads} \
            --summary-file {log} """

rule alignment_atac:
    input:
        R1=expand("{path}/{sample}_1.fastq.gz", path=config["path_atac"], sample=config["fastqs_atac_lanes"],sep=","),
        R2=expand("{path}/{sample}_2.fastq.gz", path=config["path_atac"], sample=config["fastqs_atac_lanes"],sep=",")
    threads: 
        config["threads_num"]
    params: 
        index_gen = config["hisat_index"]
    output:
        temp("results_{sample_id}/atac/alignment/atac.sam")
    log: "logs/atac_hisat_{sample_id}.log"
    conda:
        "../envs/hisat.yml"
    shell:
        """ R1=$(echo {input.R1})
            R1new=$(echo $R1 | sed 's/ /,/g ')
            R2=$(echo {input.R2})
            R2new=$(echo $R2 | sed 's/ /,/g ')
            hisat2 -x {params.index_gen} \
            -1 $R1new \
            -2 $R2new \
            -S {output} \
            --no-spliced-alignment \
            -p {threads} \
            --summary-file {log}"""
