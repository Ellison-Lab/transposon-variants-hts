library(VariantAnnotation)
library(Rsamtools)
library(tidyverse)

bam <- snakemake@input[['bam']] 
# bam <- 'results/merged/w1118_male.bam'
vcf <- snakemake@input[['vcf']] 
# vcf <- 'results/snps/snps.vcf'
sample_name <- snakemake@params[['sample_2_fingerprint']]
# sample_name <- 'aa'

vr <- VariantAnnotation::readVcfAsVRanges(vcf)

gr <- vr %>% GRanges()

allele.lookup <- vr %>% as.data.frame() %>% as_tibble() %>%
  dplyr::select(seqnames, pos=start, ref, alt, specificity) %>%
  filter(specificity == 'w1118_male')

pup <- PileupParam(max_depth = 10e6,
                   distinguish_nucleotides = T,
                   distinguish_strands = F,
                   min_minor_allele_depth = 0,
                   min_nucleotide_depth = 1)

scbp <- ScanBamParam(which = gr)

pileup <- pileup(bam, scanBamParam = scbp, pileupParam = pup)

pileups.df <- pileup %>%
  as_tibble(.) %>%
  mutate(sample = sample_name)

pileups.df <- pileups.df %>%
  left_join(allele.lookup, by=c('seqnames','pos')) %>%
  mutate(sex = ifelse(nucleotide == alt, 'male','unknown')) %>%
  group_by(seqnames, pos, sex, sample) %>%
  summarise(depth = sum(`count`),.groups = 'drop') %>%
  dplyr::select(sample, seqnames,pos,sex, depth) %>%
  arrange(seqnames, pos)

write_tsv(pileups.df, snakemake@output[['tsv']])
