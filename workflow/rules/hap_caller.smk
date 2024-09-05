rule HaplotypeCaller_e:
    input:
        "results_{sample_id}/exome/recalibration/exome.recal.bam"
    output:
        initial=temp("results_{sample_id}/exome/haplotypeCaller/exome.g.vcf.gz"),
        final="results_{sample_id}/exome/haplotypeCaller/exome.vcf.gz"
    conda: "../envs/gatk.yml"
    params: 
        fa=config["genome_fa"]
    shell: 
        """ gatk --java-options "-Xmx4g" HaplotypeCaller -R {params.fa} -I {input} -O {output.initial} -ERC GVCF
        gatk --java-options "-Xmx4g" GenotypeGVCFs -R {params.fa} -V {output.initial} -O {output.final} """

rule HaplotypeCaller:
    input:
        "results_{sample_id}/rna/recalibration/rna.recal.bam"
    output:
        initial=temp("results_{sample_id}/rna/haplotypeCaller/rna.g.vcf.gz"),
        final="results_{sample_id}/rna/haplotypeCaller/rna.vcf.gz"
    conda: "../envs/gatk.yml"
    params: 
        fa=config["genome_fa"]
    shell: 
        """ gatk --java-options "-Xmx4g" HaplotypeCaller -R {params.fa} -I {input} -O {output.initial} -ERC GVCF
        gatk --java-options "-Xmx4g" GenotypeGVCFs -R {params.fa} -V {output.initial} -O {output.final} """

rule HardFiltering_r:
    input:
        "results_{sample_id}/rna/haplotypeCaller/rna.vcf.gz"
    output:
        final="results_{sample_id}/rna/filtration/filtered.vcf.gz",
        ind_in="results_{sample_id}/rna/filtration/indels.vcf.gz",
        snps_in="results_{sample_id}/rna/filtration/snps.vcf.gz",
        ind_fin="results_{sample_id}/rna/filtration/indels_filtered.vcf.gz",
        snps_fin="results_{sample_id}/rna/filtration/snps_filtered.vcf.gz"
    conda: "../envs/gatk.yml"
    params: 
        dbsnp=config["known_sites"],
        dict_gen=config["genome_dict"],
        outdir="results_{sample_id}",
        tec="rna"
    shell:
        """gatk SelectVariants -V {input} -select-type SNP -O {output.snps_in}
        gatk SelectVariants -V {input} -select-type INDEL -O {output.ind_in}

        gatk VariantFiltration -V {output.snps_in} -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O {output.snps_fin}

        gatk VariantFiltration  -V {output.ind_in} -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O {output.ind_fin}

        gatk MergeVcfs I={output.snps_fin} I={output.ind_fin} O={output.final}

        gatk CollectVariantCallingMetrics -I {output.final} --DBSNP {params.dbsnp} -SD {params.dict_gen} -O {params.outdir}/{params.tec}/filtration/metrics """
    
rule HardFiltering_e:
    input:
        "results_{sample_id}/exome/haplotypeCaller/exome.vcf.gz"
    output:
        final="results_{sample_id}/exome/filtration/filtered.vcf.gz",
        ind_in="results_{sample_id}/exome/filtration/indels.vcf.gz",
        snps_in="results_{sample_id}/exome/filtration/snps.vcf.gz",
        ind_fin="results_{sample_id}/exome/filtration/indels_filtered.vcf.gz",
        snps_fin="results_{sample_id}/exome/filtration/snps_filtered.vcf.gz"
    conda: "../envs/gatk.yml"
    params: 
        dbsnp=config["known_sites"],
        dict_gen=config["genome_dict"],
        outdir="results_{sample_id}",
        tec="exome"
    shell:
        """gatk SelectVariants -V {input} -select-type SNP -O {output.snps_in}
        gatk SelectVariants -V {input} -select-type INDEL -O {output.ind_in}

        gatk VariantFiltration -V {output.snps_in} -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O {output.snps_fin}

        gatk VariantFiltration  -V {output.ind_in} -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O {output.ind_fin}

        gatk MergeVcfs I={output.snps_fin} I={output.ind_fin} O={output.final}

        gatk CollectVariantCallingMetrics -I {output.final} --DBSNP {params.dbsnp} -SD {params.dict_gen} -O {params.outdir}/{params.tec}/filtration/metrics """
    
rule het_selection_rna:
    input: 
        "results_{sample_id}/rna/filtration/snps_filtered.vcf.gz"
    output:
        fin="results_{sample_id}/rna/filtration/snps_het.vcf.gz",
        inter=temp("results_{sample_id}/rna/filtration/snps_het_1.vcf.gz"),
        noed=temp("results_{sample_id}/rna/filtration/snps_noediting.vcf.gz")
    params:
        ad=config["AD"],
        samp=config["sample_name"],
        samplefile="results_{sample_id}/rna/filtration/sample.txt",
        editingsites=config["editing_sites"]
    conda: "../envs/samtools.yml"
    shell:
        """ echo "rna_{params.samp} {params.samp}"> {params.samplefile}
        bcftools view {input} -i 'GT=="het" & sMIN(AD)>{params.ad}' -m2 -M2 -O z -o {output.inter}
        bcftools view -T ^{params.editingsites} {output.inter} -Oz -o {output.noed}
        bcftools reheader -s {params.samplefile} {output.inter} -o {output.fin}
        tabix {output.fin} """

rule het_selection_exome:
    input: 
        "results_{sample_id}/exome/filtration/snps_filtered.vcf.gz"
    output:
        fin="results_{sample_id}/exome/filtration/snps_het.vcf.gz",
        inter=temp("results_{sample_id}/exome/filtration/snps_het_1.vcf.gz")
    params:
        samp=config["sample_name"],
        samplefile="results_{sample_id}/exome/filtration/sample.txt"
    conda: "../envs/samtools.yml"
    shell:
        """ echo "exome_{params.samp} {params.samp}"> {params.samplefile}
        bcftools view {input} -i 'GT=="het"' -m2 -M2 -O z -o {output.inter}
        bcftools reheader -s {params.samplefile} {output.inter} -o {output.fin}
        tabix {output.fin} """  

rule het_selection_atac:
    input: 
        "results_{sample_id}/atac/filtration/snps_filtered.vcf.gz"
    output:
        fin="results_{sample_id}/atac/filtration/snps_het.vcf.gz",
        inter=temp("results_{sample_id}/atac/filtration/snps_het_1.vcf.gz")
    params:
        samp=config["sample_name"],
        samplefile="results_{sample_id}/atac/filtration/sample.txt"
    conda: "../envs/samtools.yml"
    shell:
        """ echo "atac_{params.samp} {params.samp}"> {params.samplefile}
        bcftools view {input} -i 'GT=="het"' -m2 -M2 -O z -o {output.inter}
        bcftools reheader -s {params.samplefile} {output.inter} -o {output.fin}
        tabix {output.fin} """    

# rule concatenate:
#     input: 
#         expand("results_{sample_id}/{tec_both}/filtration/snps_het.vcf.gz",tec_both=['exome','rna'])
#     output:
#         "results_{sample_id}/merged_vcf/rna_exome.vcf.gz"
#     conda: "../envs/samtools.yml"
#     shell:
#        """ bcftools concat {input} -O z -a -o {output}
#        tabix {output} """

#variant calling performed with strelka using the "GermlineWorkflow" and the "--exome" option and selection of only peaks (+-200 bp) regions

rule GATK_haplotypeCall:
    input:
        peaks= "results_{sample_id}/atac/peaks/MACS2/atac_features_BlacklistRemoved.bed.gz",
        bam= "results_{sample_id}/atac/mapping_result/atac.final.bam",
        bai= "results_{sample_id}/atac/mapping_result/atac.final.bam.bai"
    output:
        initial= "results_{sample_id}/atac/haplotypeCaller/atac.g.vcf.gz",
        vcf= "results_{sample_id}/atac/haplotypeCaller/atac.vcf.gz",
        vcf_biall= "results_{sample_id}/atac/haplotypeCaller/atac.bial.vcf.gz"
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
        "results_{sample_id}/atac/haplotypeCaller/atac.bial.vcf.gz"
    output:
        final="results_{sample_id}/atac/filtration/snps_filtered.vcf.gz",
        initial="results_{sample_id}/atac/haplotypeCaller/variantsPASS.vcf.gz",
        vcftbi="results_{sample_id}/atac/filtration/snps_filtered.vcf.gz.tbi"
    conda:
        "../envs/samtools.yml"
    params:
        DP=config["DP"]
    shell:
        """ 
        bcftools norm --rm-dup all {input} | bgzip > {output.initial}
        
	    bcftools view {input} -i 'FORMAT/DP>{params.DP}' -O z -o {output.final}

	    tabix {output.final} """
