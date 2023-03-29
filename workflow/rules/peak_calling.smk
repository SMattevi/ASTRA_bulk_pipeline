#peak calling with MACS2 
rule peak_calling:
    input: 
        "results/atac/mapping_result/atac.positionsort.MAPQ20.bam"
    output:
        "results/atac/peaks/MACS2/atac_features_BlacklistRemoved.bed"
    params:
        config["bl_file"] #blacklist regions file
    conda:
        "../envs/peakcall.yml"
    shell:
        """ 
        mkdir -p results/atac/peaks/MACS2/
        macs2 callpeak -t {input} --outdir results/atac/peaks/MACS2/ -n atac -f BAMPE -q 0.05 -g hs --nomodel --extsize 200 --slocal 1000
        
        awk '{{ if ($1>=1 && $1<=22 || $1=="X" || $1=="Y" || $1=="M") {{print $0}}}}' results/atac/peaks/MACS2/atac_peaks.narrowPeak\
            >results/atac/peaks/MACS2/tmp.narrowPeak
        mv results/atac/peaks/MACS2/tmp.narrowPeak results/atac/peaks/MACS2/atac_peaks.narrowPeak

        bedtools intersect -a results/atac/peaks/MACS2/atac_peaks.narrowPeak -b {params} -v \
            > results/atac/peaks/MACS2/atac_features_BlacklistRemoved.bed
        """

rule bgzip_peakfile:
    input:
        "results/atac/peaks/MACS2/atac_features_BlacklistRemoved.bed"
    output:
        "results/atac/peaks/MACS2/atac_features_BlacklistRemoved.bed.gz"
    conda:
        "../envs/samtools.yml"
    shell:
        """ bgzip {input}
        tabix {output} """

rule annotate_peaks:
    input:
        "results/atac/peaks/MACS2/atac_features_BlacklistRemoved.bed.gz"
    output:
        "results/atac/peaks/annotated.tsv",
        "results/atac/peaks/annotated.pdf"
    params: config["TSS_region"]
    shell:
        "Rscript workflow/scripts/annotate_peaks.R {params}"