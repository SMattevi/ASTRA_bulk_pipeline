rule het_selection_genome:
    input: 
        config["genome_vcf"]
    output:
        fin="results/genome/filtration/snps_het.vcf.gz",
        inter=temp("results/genome/filtration/snps_het_1.vcf.gz"),
        prephased=temp("results/genome/prephasing/pre_phased.vcf.gz")
    params:
        samp=config["sample_name"],
        samplefile="results/genome/filtration/sample.txt"
    conda: "../envs/samtools.yml"
    shell:
        """ SAMPLE_OR=$(bcftools query -l {input})
        echo "$SAMPLE_OR {params.samp}"> {params.samplefile}
        bcftools view {input} -i 'GT=="het"' -m2 -M2 -O z -o {output.inter}
        bcftools reheader -s {params.samplefile} {output.inter} -o {output.fin}
        tabix {output.fin} 
        mkdir -p results/genome/prephasing
        cp {output.fin} {output.prephased}
        tabix {output.prephased}
        """