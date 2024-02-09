#!/usr/bin/env Rscript
library("optparse")
library("ASTRA")

option_list = list(
  make_option(c("-c", "--chr"), type="character", default=NULL,
              help="chromosome number", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-v","--vcf"), type="character", default=NULL,
              help="VCF not phased", metavar="character"),
  make_option(c("-s","--shapeit"), type="character", default=NULL,
              help="VCF phased shapeit format", metavar="character"),
  make_option(c("-x","--haptreex"), type="character", default=NULL,
              help="VCF phased haptreex format", metavar="character"),
  make_option(c("-a","--ase"), type="character", default=NULL,
              help="ASE table", metavar="character"),
  make_option(c("-n","--sample"), type="character", default=NULL,
              help="ASE table", metavar="character")

);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt)

BM<-biomart_df()

chromosome=opt$chr
BMsel<-BM[BM$chromosome_name==chromosome,]

varcall<-opt$vcf
phased<-opt$shapeit
phasedHP<-opt$haptreex

ASE<-opt$ase

df_sh<-read_shapeit(varcall,phased,sample_name = opt$sample)
df_sh<-df_sh[df_sh$CHROM==chromosome,]

df_hp<-read_haptreex(varcall,phasedHP,sample_name = opt$sample)
df_hp<-df_hp[df_hp$CHROM==chromosome,]


olaps<-conshap(df1=df_sh,df2 = df_hp,BMsel = BMsel,DF1="WShapeit",DF2="HapTreeX")

df_compl<-addASE(olaps,ASE,BMsel)
df_tmp<-df_compl

df_compl<-biallmonoall(df_compl)

df_ph<-manualP(df_compl)

final_df<-prep_outcome(df_initial = olaps,df_phased = df_ph[df_ph$Ref!="NP",])
final_df$Tool<-ifelse(is.na(final_df$Tool),"Manual",final_df$Tool)

final_df$Tool<-as.factor(final_df$Tool)

write.table(final_df[,c(1:7)],opt$out, sep="\t",quote=F, row.names = F,col.names = F)
