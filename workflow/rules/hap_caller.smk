rule HaplotypeCaller_e:
    input:
        "results/exome/recalibration/exome.recal.bam"
    output:
        initial=temp("results/exome/haplotypeCaller/exome.g.vcf.gz"),
        final="results/exome/haplotypeCaller/exome.vcf.gz"
    conda: "../envs/gatk.yml"
    params: 
        fa=config["genome_fa"]
    shell: 
        """ gatk --java-options "-Xmx4g" HaplotypeCaller -R {params.fa} -I {input} -O {output.initial} -ERC GVCF
        gatk --java-options "-Xmx4g" GenotypeGVCFs -R {params.fa} -V {output.initial} -O {output.final} """

rule HaplotypeCaller:
    input:
        "results/rna/recalibration/rna.recal.bam"
    output:
        initial=temp("results/rna/haplotypeCaller/rna.g.vcf.gz"),
        final="results/rna/haplotypeCaller/rna.vcf.gz"
    conda: "../envs/gatk.yml"
    params: 
        fa=config["genome_fa"]
    shell: 
        """ gatk --java-options "-Xmx4g" HaplotypeCaller -R {params.fa} -I {input} -O {output.initial} -ERC GVCF
        gatk --java-options "-Xmx4g" GenotypeGVCFs -R {params.fa} -V {output.initial} -O {output.final} """

rule HardFiltering:
    input:
        "results/{tec}/haplotypeCaller/{tec}.vcf.gz"
    output:
        final="results/{tec}/filtration/filtered.vcf.gz",
        ind_in="results/{tec}/filtration/indels.vcf.gz",
        snps_in="results/{tec}/filtration/snps.vcf.gz",
        ind_fin="results/{tec}/filtration/indels_filtered.vcf.gz",
        snps_fin="results/{tec}/filtration/snps_filtered.vcf.gz"
    conda: "../envs/gatk.yml"
    params: 
        dbsnp=config["known_sites"],
        dict_gen=config["genome_dict"]
    shell:
        """gatk SelectVariants -V {input} -select-type SNP -O {output.snps_in}
        gatk SelectVariants -V {input} -select-type INDEL -O {output.ind_in}

        gatk VariantFiltration -V {output.snps_in} -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O {output.snps_fin}

        gatk VariantFiltration  -V {output.ind_in} -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O {output.ind_fin}

        gatk MergeVcfs I={output.snps_fin} I={output.ind_fin} O={output.final}

        gatk CollectVariantCallingMetrics -I {output.final} --DBSNP {params.dbsnp} -SD {params.dict_gen} -O results/{wildcards.tec}/filtration/metrics """
    
rule het_selection_rna:
    input: 
        "results/rna/filtration/snps_filtered.vcf.gz"
    output:
        fin="results/rna/filtration/snps_het.vcf.gz",
        inter=temp("results/rna/filtration/snps_het_1.vcf.gz")
    params:
        ad=config["AD"],
        samp=config["sample_name"],
        samplefile="results/rna/filtration/sample.txt"
    conda: "../envs/samtools.yml"
    shell:
        """ echo "rna_{params.samp} {params.samp}"> {params.samplefile}
        bcftools view {input} -i 'GT=="het" & sMIN(AD)>{params.ad}' -m2 -M2 -O z -o {output.inter}
        bcftools reheader -s {params.samplefile} {output.inter} -o {output.fin}
        tabix {output.fin} """

rule het_selection_exome:
    input: 
        "results/exome/filtration/snps_filtered.vcf.gz"
    output:
        fin="results/exome/filtration/snps_het.vcf.gz",
        inter=temp("results/exome/filtration/snps_het_1.vcf.gz")
    params:
        samp=config["sample_name"],
        samplefile="results/exome/filtration/sample.txt"
    conda: "../envs/samtools.yml"
    shell:
        """ echo "exome_{params.samp} {params.samp}"> {params.samplefile}
        bcftools view {input} -i 'GT=="het"' -m2 -M2 -O z -o {output.inter}
        bcftools reheader -s {params.samplefile} {output.inter} -o {output.fin}
        tabix {output.fin} """  

rule het_selection_atac:
    input: 
        "results/atac/filtration/snps_filtered.vcf.gz"
    output:
        fin="results/atac/filtration/snps_het.vcf.gz",
        inter=temp("results/atac/filtration/snps_het_1.vcf.gz")
    params:
        samp=config["sample_name"],
        samplefile="results/atac/filtration/sample.txt"
    conda: "../envs/samtools.yml"
    shell:
        """ echo "atac_{params.samp} {params.samp}"> {params.samplefile}
        bcftools view {input} -i 'GT=="het"' -m2 -M2 -O z -o {output.inter}
        bcftools reheader -s {params.samplefile} {output.inter} -o {output.fin}
        tabix {output.fin} """    

# rule concatenate:
#     input: 
#         expand("results/{tec_both}/filtration/snps_het.vcf.gz",tec_both=['exome','rna'])
#     output:
#         "results/merged_vcf/rna_exome.vcf.gz"
#     conda: "../envs/samtools.yml"
#     shell:
#        """ bcftools concat {input} -O z -a -o {output}
#        tabix {output} """

#variant calling performed with strelka using the "GermlineWorkflow" and the "--exome" option and selection of only peaks (+-200 bp) regions

rule GATK_haplotypeCall:
    input:
        peaks= "results/atac/peaks/MACS2/atac_features_BlacklistRemoved.bed.gz",
        bam= "results/atac/mapping_result/atac.final.bam",
        bai= "results/atac/mapping_result/atac.final.bam.bai"
    output:
        initial= "results/atac/haplotypeCaller/atac.g.vcf.gz",
        vcf= "results/atac/haplotypeCaller/atac.vcf.gz",
        vcf_biall= "results/atac/haplotypeCaller/atac.bial.vcf.gz"
    params:
        genome= config["genome_fa"]
    conda: "../envs/gatk.yml"
    shell:
        """
        # run gatk haplotype caller
        gatk --java-options "-Xmx4g" HaplotypeCaller \
        -L {input.peaks} \
        -R {params.genome} \
        -I {input.bam} \
        -O {output.initial} \
        -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation

        # correct the genotypes that come out of haplotype caller
        gatk --java-options "-Xmx4g" GenotypeGVCFs \
        --include-non-variant-sites \
        -L {input.peaks} \
        -R {params.genome} \
        -V {output.initial} \
        -O {output.vcf}

        gatk SelectVariants \
        -R {params.genome} \
        --variant {output.vcf} \
        --restrict-alleles-to BIALLELIC \
        -select 'vc.getHetCount()==1' --select-type-to-include SNP \
        -O {output.vcf_biall}
        """

#QC: extract only het variants, with AF and DP as specified in config file
rule QC_VCF:
    input:
        "results/atac/haplotypeCaller/atac.bial.vcf.gz"
    output:
        final="results/atac/filtration/snps_filtered.vcf.gz",
	    initial="results/atac/haplotypeCaller/variantsPASS.vcf.gz",
	    vcftbi="results/atac/filtration/snps_filtered.vcf.gz.tbi"
    conda:
        "../envs/samtools.yml"
    params:
        DP=config["DP"]
    shell:
        """ 
        bcftools norm --rm-dup all {input} | bgzip > {output.initial}
        
	    bcftools view {input} -i 'FORMAT/DP>{params.DP}' -O z -o {output.final}

	    tabix {output.final} """