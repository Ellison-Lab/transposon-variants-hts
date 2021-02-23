library(tidyverse)
library(rtracklayer)

#coverage_fl <- "results/depth/w1118_female/w1118_female.mosdepth.summary.txt"
#te_fa <- '/media/mlawlor/work/data/genome/dmel_repbase_lib.fasta'

coverage_fl <- snakemake@input[['cov']]
te_fa <- snakemake@input[['fasta']]

coverage <- read_tsv(coverage_fl) %>%
  filter(!str_detect(chrom,'_region'))

tes <- import(te_fa) %>%
  names() %>% str_extract('.+(?=#)')

autosomes <- c('2L','2R','3L','3R','4')

sex.chr <- c('X','Y')

diploid_thresh <- coverage %>%
  filter(chrom %in% autosomes) %>%
  pull(mean) %>%
  mean()
  
copies <- coverage %>%
  filter(chrom %in% c(tes, autosomes, sex.chr)) %>%
  mutate(est.copies = 2*mean/diploid_thresh) %>%
  dplyr::select(sequence = chrom, length, bases, median.cov = mean, est.copies) %>%
  arrange(-est.copies)
  
write_tsv(copies,snakemake@output[['tsv']])


# ---------------------------------------------------------------------------------
# library(ggpubr)
# library(plotly)
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
#   filter(!str_detect(X1,'FBgn')) %>%
#   pull(gene_id) %>%
#   unique()
# 
# copies.male <- read_tsv('resuts/copies/w1118_male.tsv')
# copies.female <- read_tsv('resuts/copies/w1118_female.tsv')
# 
# (left_join(copies.male, copies.female, by="sequence", suffix=c('.male','.female')) %>%
#   mutate(tep = ifelse(sequence %in% tep_tes,'TEP','other')) %>%
#   #filter(tep!='TEP') %>%
#   filter(sequence %in% tes) %>%
#   filter(!str_detect(sequence,'LTR')) %>%
#   filter(!sequence %in% excl.seqs) %>%
#   ggplot(aes(est.copies.male,est.copies.female, label=sequence)) +
#   geom_point(aes(color=tep)) +
#   theme(aspect.ratio = 1) +
#   geom_abline(intercept = 0, slope = 1)) %>%
#   ggplotly
# 
# left_join(copies.male, copies.female, by="sequence", suffix=c('.male','.female')) %>%
#   filter(!sequence %in% excl.seqs) %>%
#     filter(sequence %in% tes) %>%
#   filter(!str_detect(sequence,'[-_]LTR')) %>%
#   mutate(tep = ifelse(sequence %in% tep_tes,'TEP','other')) %>%
#   ggplot(aes(tep,est.copies.male/est.copies.female, label=sequence))+
#   geom_boxplot() +
#   geom_jitter(width = 0.1) +
#   theme(aspect.ratio = 1) +
#   stat_compare_means() +
#   xlab('')



