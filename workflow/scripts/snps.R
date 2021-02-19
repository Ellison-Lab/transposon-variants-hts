library(tidyverse)
library(rtracklayer)

pileups <- snakemake@input[['pileups']] # Sys.glob('results/pileups/*pileups.tsv.gz')
te_fa <- snakemake@input[['fasta']] #'/media/mlawlor/work/data/genome/dmel_repbase_lib.fasta'
min_snp_depth <- snakemake@params[['min_snp_depth']]

tes <- import(te_fa) %>%
  names() %>% str_extract('.+(?=#)')

pileups.df2 <- map_df(pileups, read_tsv) %>% filter(seqnames %in% tes)

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
  mutate(specificity= colnames(select(.,any_of(sample_names)))[max.col(select(.,any_of(sample_names)))])

pileups.df2 %>%
  mutate(start = pos, end=pos, name=paste(specificity,nucleotide, sep = ',')) %>%
  GRanges() %>% 
  export(snakemake@output[['bed']])

write_tsv(pileups.df2, snakemake@output[['tsv']])


# ----------------------

# 
# 
# pileups.df2 %>%
#   ggplot(aes(specificity)) +
#   geom_bar()
# 
# 
# pileups.df2 %>%
#   gather(sex,count,w1118_female, w1118_male) %>%
#   filter(count > 0) %>%
#   ggplot(aes(count)) +
#   geom_histogram() +
#   facet_wrap(~sex) +
#   scale_x_log10() +
#   ggtitle('distribution of read counts supporting sex-specific variants') +
#   xlab('read depth supporting variant')
# 
# 
# pileups.df2 %>%
#   filter(w1118_male > 10 | w1118_female > 10) %>%
#   gather(sex,count,w1118_female, w1118_male) %>%
#   filter(count > 0) %>%
#   dplyr::select(seqnames, pos, sex) %>%
#   distinct() %>%
#   #mutate(seqnames = str_remove(seqnames,'#.+')) %>%
#   ggplot(aes(forcats::fct_infreq(seqnames))) +
#   geom_bar() +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
#   facet_wrap(~sex, ncol=1) +
#   scale_y_log10() +
#   ggtitle('Number of sex-specific variant sites detected on each TE') +
#   ylab('number of variant sites detected') +
#   xlab("")
