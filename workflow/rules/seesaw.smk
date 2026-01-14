rule transcripts_extraction:
    input: 
        "results_{sample_id}/phased/manual_refinment.vcf.gz"
    output:
        gtf=temp("results_{sample_id}/seesaw/g2gtools.gtf"),
        db=temp("results_{sample_id}/seesaw/g2gtools.db"),
        fa=temp("results_{sample_id}/seesaw/g2gtools.fa"),
        vcitbi=temp("results_{sample_id}/seesaw/g2gtools.vci.gz.tbi"),
        vcigz=temp("results_{sample_id}/seesaw/g2gtools.vci.gz"),
        gtfunm=temp("results_{sample_id}/seesaw/g2gtools.gtf.unmapped"),
        fafai=temp("results_{sample_id}/seesaw/g2gtools.fa.fai"),
        transcript="results_{sample_id}/seesaw/transcripts.fa"
    params:
        genomefa=config["genome_fa"],
        gtf=config["genome_gtf"],
        sample="{sample_id}"
    conda: "../envs/g2gtools.yml"
    threads: config["threads_num"]
    shell:
        """ vci={output.vcigz}
        vci=${{vci%.gz}}
        mkdir -p results_{params.sample}/seesaw/

        g2gtools vcf2vci -f {params.genomefa} -o $vci -s {params.sample} --diploid -p {threads} -i {input}

        # Incorporate SNP
        g2gtools patch -i {params.genomefa} -c {output.vcigz} -o {output.fa} -p {threads}

        # Liftover gene annotation onto sample coordinates.
        g2gtools convert -f gtf -i {params.gtf} -c {output.vcigz} -o {output.gtf}
                            
        # Parse custom gene annotation into a G2G database file.
        g2gtools gtf2db -i {output.gtf} -o {output.db}
                            
        # Extract exons, transcripts, and gene regions from the custom sample genome.
        g2gtools extract --transcripts -i {output.fa} -db {output.db} >{output.transcript}
        """

rule isoform_quantification:
    input: 
        transcripts="results_{sample_id}/seesaw/transcripts.fa",
        #fq1=expand("{path}/{sample}_R1.fastq.gz", path=config["path_rna"],sample=config["fastqs_rna"],sep=","),
        #fq2=expand("{path}/{sample}_R2.fastq.gz", path=config["path_rna"],sample=config["fastqs_rna"],sep=",")
        fq1=lambda wildcards: find_fastqs(
            path=config["SAMPLES"][wildcards.sample_id]["path_rna"],
            prefixes=config["SAMPLES"][wildcards.sample_id]["fastqs_rna"],
            read_num=1
        ),
        fq2=lambda wildcards: find_fastqs(
            path=config["SAMPLES"][wildcards.sample_id]["path_rna"],
            prefixes=config["SAMPLES"][wildcards.sample_id]["fastqs_rna"],
            read_num=2
        )
    output:
        final="results_{sample_id}/seesaw/salmon/quant.sf",
        txome=directory("results_{sample_id}/seesaw/salmon/diploid_txome")
    conda:
        "../envs/salmon.yml"
    threads: config["threads_num"]
    params: 
        outdir="results_{sample_id}/seesaw/salmon"
    shell:
        """
        R1=$(echo {input.fq1})
        #R1new=$(echo $R1 | sed 's/ /,/g ')
        R2=$(echo {input.fq2})
        #R2new=$(echo $R2 | sed 's/ /,/g ')

        salmon index -p {threads} -i {output.txome} -t {input.transcripts} --keepDuplicates

        salmon quant -i {output.txome} \
        -l A -p {threads} --numBootstraps 30 \
        -o {params.outdir} -1 $R1 -2 $R2
        """