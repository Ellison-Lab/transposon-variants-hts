library(tidyverse)
library(rtracklayer)
library(VariantAnnotation)

pileups <- snakemake@input[['pileups']] # pileups <- Sys.glob('results/pileups/*pileups.tsv.gz')
te_fa <- snakemake@input[['fasta']] #te_fa <- '/media/mlawlor/work/data/genome/dmel_repbase_lib.fasta'
min_snp_depth <- snakemake@params[['min_snp_depth']] # min_snp_depth <- 15
sample_2_fingerprint <- snakemake@params[['sample_2_fingerprint']] # sample_2_fingerprint <- 'w1118_male'


te_dss <- import(te_fa)

tes <- names(te_dss) %>% str_extract('.+(?=#)')

names(te_dss) <- tes

pileups.df2 <- map_df(pileups, read_csv) %>% filter(seqnames %in% tes)

sample_names <- pileups.df2$sample %>% unique

samples_w_coverage <- pileups.df2 %>%
  group_by(sample,seqnames, pos) %>%
  summarise(total_depth=sum(count)) %>%
  group_by(seqnames,pos) %>%
  summarize(samps.w.cov = sum(total_depth > 1)) %>%
  ungroup()

pileups.df2 <- pileups.df2 %>%
  spread(sample, count, fill = 0)

# remove positions that are completely uncovered in a sample. This prevents wide stretches of positions from being used as sample=specific
pileups.df2 <- left_join(pileups.df2, samples_w_coverage, by=c('seqnames','pos')) %>%
  filter(samps.w.cov == length(sample_names))

# take only snps with 0 depth in all but 1 sample
pileups.df2 <- pileups.df2 %>%
  rowwise() %>%
  filter(sum(c_across(any_of(sample_names)) > 0) == 1) %>%
  ungroup()

# make sure variant supported by more than X reads
pileups.df2 <-pileups.df2 %>%
  rowwise() %>%
  filter(sum(c_across(any_of(sample_names))) > min_snp_depth) %>%
  ungroup()

# add col with snp specificity, ie the name of the sample that exclusively has that snp
pileups.df2 <- pileups.df2 %>%
  mutate(specificity= colnames(dplyr::select(.,any_of(sample_names)))[max.col(dplyr::select(.,any_of(sample_names)))])

pileups.df2 <- pileups.df2 %>%
  filter(specificity == sample_2_fingerprint)

# arbitrarily add genotype calls - this is bc star needs this to run allelic expression
pileups.df2$GT <- c("0/1")

gr <- pileups.df2 %>%
  mutate(start = pos, end=pos, name=paste(specificity,nucleotide, sep = ',')) %>%
  mutate(ref = "", alt=nucleotide) %>%
  GRanges()

refs <- getSeq(te_dss, gr) %>% as.character()

gr$ref <- refs

# store N where the male specific allele is also ref
gr[(gr$ref == gr$alt)]$ref <- "N"

message('making vr')
vr <- VariantAnnotation::makeVRangesFromGRanges(gr)[,c('GT','specificity')]

sampleNames(vr) <- vr$specificity

output_vcf <- snakemake@output[['vcf']] #output_vcfs <- c('~/Downloads/snps/w1118_male-snps.vcf')

message('writing vr')

writeVcf(obj = vr, filename = output_vcf, index=F)

export(gr, snakemake@output[['bed']])

write_csv(pileups.df2, snakemake@output[['csv']])
