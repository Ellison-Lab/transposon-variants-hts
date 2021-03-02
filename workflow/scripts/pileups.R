library(tidyverse)
library(Rsamtools)
library(rtracklayer)

pu.param <- PileupParam(max_depth=10e6,
            min_base_quality=snakemake@params[['bq']],
            min_nucleotide_depth=1,
            min_minor_allele_depth=0, distinguish_strands=F,
            distinguish_nucleotides=TRUE, ignore_query_Ns=TRUE,
            include_deletions=F, include_insertions=F)

te_fa <- snakemake@input[['fasta']] #'/media/mlawlor/work/data/genome/dmel_repbase_lib.fasta'

te_dss <- import(te_fa)

sb.param <- ScanBamParam(which = GRanges(seqnames =  str_extract(names(te_dss),'.+(?=#)'), IRanges(start=1, width=width(te_dss))))


pileups.df <- pileup(snakemake@input[['bam']], pileupParam = pu.param, scanBamParam = sb.param) %>%
  as_tibble() %>%
  mutate(sample = snakemake@wildcards[['sample']]) %>%
  arrange(seqnames, pos)

write_csv(pileups.df, snakemake@output[['csv']])


# # --------------------------
#
# library(GenomicRanges)
# library(arrow)
#
# te.lookup <- read_tsv("~/work/TestisTEs2021/resources/te_id_lookup.curated.tsv.txt")
#
# geps <- open_dataset('~/finalized/larval-w1118-testes/optimal_gep_membership/', format='arrow') %>%
#   collect()
#
# top_mods <- geps %>%
#   filter(qval < 0.005) %>%
#   group_by(module) %>% summarise(n_tes = sum(!str_detect(X1,'FBgn'))) %>%
#   arrange(-n_tes)
#
# tep_tes <- geps %>%
#   filter(qval < 0.005) %>%
#   filter(module == top_mods$module[1]) %>%
#   left_join(te.lookup, by=c(X1='merged_te')) %>%
#   pull(gene_id) %>%
#   unique()
#
