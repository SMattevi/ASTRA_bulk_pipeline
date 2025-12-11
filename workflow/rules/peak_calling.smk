#peak calling with MACS2 
rule peak_calling:
    input: 
        "results_{sample_id}/atac/mapping_result/atac.positionsort.MAPQ20.bam"
    output:
        "results_{sample_id}/atac/peaks/MACS2/atac_features_BlacklistRemoved.bed"
    params:
        bl= config["bl_file"], #blacklist regions file
        sampleid="{sample_id}"
    conda:
        "../envs/peakcall.yml"
    shell:
        """ 
        mkdir -p results_{params.sampleid}/atac/peaks/MACS2/
        macs2 callpeak -t {input} --outdir results_{params.sampleid}/atac/peaks/MACS2/ -n atac -f BAMPE -q 0.05 -g hs --nomodel --extsize 200 --slocal 1000
        
        awk '{{ if ($1>=1 && $1<=22 || $1=="X" || $1=="Y" || $1=="M") {{print $0}}}}' results_{params.sampleid}/atac/peaks/MACS2/atac_peaks.narrowPeak\
            >results_{params.sampleid}/atac/peaks/MACS2/tmp.narrowPeak
        mv results_{params.sampleid}/atac/peaks/MACS2/tmp.narrowPeak results_{params.sampleid}/atac/peaks/MACS2/atac_peaks.narrowPeak

        bedtools intersect -a results_{params.sampleid}/atac/peaks/MACS2/atac_peaks.narrowPeak -b {params.bl} -v \
            > results_{params.sampleid}/atac/peaks/MACS2/atac_features_BlacklistRemoved.bed
        """

rule bgzip_peakfile:
    input:
        "results_{sample_id}/atac/peaks/MACS2/atac_features_BlacklistRemoved.bed"
    output:
        "results_{sample_id}/atac/peaks/MACS2/atac_features_BlacklistRemoved.bed.gz"
    conda:
        "../envs/samtools.yml"
    shell:
        """ bgzip {input}
        tabix {output} """

rule annotate_peaks:
    input:
        "results_{sample_id}/atac/peaks/MACS2/atac_features_BlacklistRemoved.bed.gz"
    output:
        "results_{sample_id}/atac/peaks/annotated.tsv",
        "results_{sample_id}/atac/peaks/annotated.pdf"
    params: config["TSS_region"]
    shell:
        "Rscript workflow/scripts/annotate_peaks.R {params} {wildcards.sample_id}"