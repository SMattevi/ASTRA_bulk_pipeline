rule alignment_LR:
    input:
        config["fastqs_LR"]
    threads: 
        config["threads_num"]
    params: 
        genome_ref = config["genome_fa"]
    output:
        bam="results_{sample_id}/LR/alignment/{sample_id}_sorted.bam",
    conda:
        "../envs/longreads.yml"
    shell:
        """ 
        mkdir -p tmp/ref_LR_alignment
        # Add if missing "chr" to reference - needed for following steps
        less {params.genome_ref} | sed -e '/^>chr/!s/^>/>chr/' > tmp/ref_LR_alignment/ref.fa

        mkdir -p results_{wildcards.sample_id}/LR/alignment
        # Align with minimap
        minimap2 -ax splice -uf -k14 tmp/ref_LR_alignment/ref.fa {input} > results_{wildcards.sample_id}/LR/alignment/BJ.sam

        # Sort alignment results with samtools
        samtools sort results_{wildcards.sample_id}/LR/alignment/{wildcards.sample_id}.sam -o {output.bam}
        samtools index {output.bam}

        rm tmp/ref_LR_alignment/ref.fa
        """

rule add_chr:
    input: 
        "results_{sample_id}/phased/manual_refinment.vcf.gz"
    output:
        "results_{sample_id}/phased/manual_refinment_chr.vcf.gz"
    shell:
        """
        for i in {{1..22}}; do echo "$i chr$i"; done > chr_map.txt
        echo "X chrX" >> chr_map.txt
        echo "Y chrY" >> chr_map.txt
        echo "MT chrMT" >> chr_map.txt
        bcftools annotate --rename-chrs chr_map.txt {input} -Ov -o {output} 
        rm chr_map.txt
        """

rule make_new_vcf:
    input: 
        bam="results_{sample_id}/LR/alignment/{sample_id}_sorted.bam",
        vcf="results_{sample_id}/phased/manual_refinment_chr.vcf.gz"
    params:
        reffa=config["genome_fa"],
        newdir=subpath(output.hap1, parent=True)
    output:
        hap1="results_{sample_id}/LR/lorals/make_new_vcf/{sample_id}_sorted_hap1.fa",
        hap2="results_{sample_id}/LR/lorals/make_new_vcf/{sample_id}_sorted_hap2.fa",
        vcf="results_{sample_id}/LR/lorals/make_new_vcf/{sample_id}_sorted.vcf"
    conda:
        "../envs/lorals.yml"
    shell:
        """
        pwd=$PWD
        cp /mnt/projects-2025/iPSC_lr/lorals/lorals/HapCUT2VCF.py .
        export PATH="$PATH:$pwd" 
        
        mkdir -p {params.newdir}

        make_new_vcf.sh -b {input.bam} -G {params.reffa} -V {input.vcf} -o {params.newdir}
        
        cd {params.newdir}
        mv {wildcards.sample_id}_sorted.hap2.fa {wildcards.sample_id}_sorted_hap2.fa
        mv {wildcards.sample_id}_sorted.hap1.fa {wildcards.sample_id}_sorted_hap1.fa
        """

rule hap_aligner:
    input:
        fastq=config["fastqs_LR"],
        new_vcf_hap1="results_{sample_id}/LR/lorals/make_new_vcf/{sample_id}_sorted_hap1.fa",
        new_vcf_hap2="results_{sample_id}/LR/lorals/make_new_vcf/{sample_id}_sorted_hap2.fa"
    output:
        "results_{sample_id}/LR/lorals/hap_aligner/{sample_id}_reads_aln_sorted.merged.bam"
    params:
        newdir="results_{sample_id}/LR/lorals/hap_aligner/",
        input_suff="results_{sample_id}/LR/lorals/make_new_vcf/{sample_id}_sorted_"
    conda:
        "../envs/lorals.yml"
    shell:
        """
        mkdir -p {params.newdir} 
        hap_aligner.sh \
            -f {input.fastq} \
            -G {params.input_suff} \
            -o {params.newdir}
        """

rule calc_annotate_ase:
    input:
        bam="results_{sample_id}/LR/lorals/hap_aligner/{sample_id}_reads_aln_sorted.merged.bam",
        vcf="results_{sample_id}/LR/lorals/make_new_vcf/{sample_id}_sorted.vcf"
    output:
        ase="results_{sample_id}/LR/lorals/ase/ase.tsv",
        ase_ann="results_{sample_id}/LR/lorals/ase/ase_annotated.tsv"
    params:
        outdir=subpath(output.ase, parent=True),
    conda:
        "../envs/lorals.yml"
    shell:
        """
        mkdir -p {params.outdir}
        calc_ase -b {input.bam} -f {input.vcf} -o {params.outdir}
        annotate_ase -i {output.ase} -b /mnt/projects-2025/iPSC_lr/resources/ref.bed -o {params.outdir}
        """

rule alignment_transcript:
    input:
        config["fastqs_LR"]
    threads: 
        config["threads_num"]
    params: 
        transcript_ref= config["transcriptome_fa"]
    output:
        sam=temp("results_{sample_id}/LR/alignment/{sample_id}_transcript.sam"),
        bam="results_{sample_id}/LR/alignment/{sample_id}_transcript.bam"
    conda:
        "../envs/longreads.yml"
    shell:
        """ 
        minimap2 -ax map-ont {params.transcript_ref} {input} > {output.sam}
        # Sort alignment results with samtools
        samtools sort {output.sam} -o {output.bam}
        samtools index {output.bam}
        """

rule calc_annotate_ast:
    input:
        bam_trs="results_{sample_id}/LR/alignment/{sample_id}_transcript.bam",
        bam_gen="results_{sample_id}/LR/alignment/{sample_id}_sorted.bam",
        ase="results_{sample_id}/LR/lorals/ase/ase.tsv"
    output:
        ast="results_{sample_id}/LR/lorals/ast/asts.tsv_asts_quant.tsv"
    params:
        outdir=subpath(output.ast, parent=True),
        t2gene=config["t2gene"]
    conda:
        "../envs/lorals.yml"
    threads: 
        config["threads_num"]
    shell:
        """
        mkdir -p {params.outdir}
        calc_asts -m quant -b {input.bam_trs} -x {input.bam_gen} -i {input.ase} -o {params.outdir} 
        process_asts -i {output.ast} -g {params.t2gene} -t {threads} -o {params.outdir}
        """