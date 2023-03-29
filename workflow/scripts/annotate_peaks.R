args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

library(tidyr)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
tss_range<-as.integer(args[1])

path="results/atac/peaks/MACS2/atac_features_BlacklistRemoved.bed.gz"

peak_file <- readPeakFile(path)
peak_file <- diffloop::addchr(peak_file)

peakAnno <- annotatePeak(peak = peak_file, tssRegion=c(-tss_range,tss_range),
                         TxDb=txdb, annoDb="org.Hs.eg.db", addFlankGeneInfo=TRUE, flankDistance=5000)

flank_gene_IDS<-as.data.frame(peakAnno@anno$flank_geneIds)
annotat<-as.data.frame(peakAnno@anno)

annotat$seqnames<-gsub("chr","",as.character(annotat$seqnames))
write.table(annotat,paste0("results/atac/peaks/annotated.tsv"),sep="\t",quote=F,row.names = F)

gene_list<-unique(separate_rows(flank_gene_IDS, 'peakAnno@anno$flank_geneIds', convert = TRUE, sep = ";"))

pdf(paste0("results/atac/peaks/annotated.pdf"))
plotAnnoPie(peakAnno)
plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")
dev.off()
