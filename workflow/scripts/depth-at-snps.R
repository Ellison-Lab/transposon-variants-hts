library(VariantAnnotation)
library(Rsamtools)
library(tidyverse)

bam <- snakemake@input[['bam']] #bam <- 'results/merged/w1118_male.bam'
vcf <- snakemake@input[['vcf']] #'results/snps/snps.vcf'
sample_name <- snakemake@wildcards[['sample']]# sample_name <- 'aa'

gr <-VariantAnnotation::readVcfAsVRanges(vcf) %>% GRanges()

pup <- PileupParam(max_depth = 10e6, 
                   distinguish_nucleotides = F,
                   distinguish_strands = F, 
                   min_minor_allele_depth = 0, 
                   min_nucleotide_depth = 1)

scbp <- ScanBamParam(which = gr)

pileup <- pileup(bam, scanBamParam = scbp, pileupParam = pup)

pileups.df <- as_tibble(pileup) %>%
  mutate(sample = sample_name) %>%
  dplyr::select(seqnames,pos,total.depth=count) %>%
  arrange(seqnames, pos)

write_tsv(pileups.df, snakemake@output[['tsv']])
