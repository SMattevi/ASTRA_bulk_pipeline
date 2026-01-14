rule het_selection_genome:
    input: 
        # config["genome_vcf"]
        lambda wildcards: config["SAMPLES"][wildcards.sample_id]["genome_vcf"]
    output:
        fin="results_{sample_id}/genome/filtration/snps_het.vcf.gz",
        inter=temp("results_{sample_id}/genome/filtration/snps_het_1.vcf.gz"),
        prephased=temp("results_{sample_id}/genome/prephasing/pre_phased.vcf.gz")
    params:
        samp="{sample_id}",
        samplefile="results_{sample_id}/genome/filtration/sample.txt",
        outdir="results_{sample_id}/genome/prephasing"
    conda: "../envs/samtools.yml"
    shell:
        """ SAMPLE_OR=$(bcftools query -l {input})
        echo "$SAMPLE_OR {params.samp}"> {params.samplefile}
        bcftools view {input} -i 'GT=="het"' -m2 -M2 --types snps -O z -o {output.inter}
        bcftools reheader -s {params.samplefile} {output.inter} -o {output.fin}
        tabix {output.fin} 
        mkdir -p {params.outdir}
        cp {output.fin} {output.prephased}
        tabix {output.prephased}
        """