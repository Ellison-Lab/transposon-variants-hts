library(tidyverse)
library(Rsamtools)

pu.param <- PileupParam(max_depth=1e6, 
            min_base_quality=30, min_mapq=30, 
            min_nucleotide_depth=10, 
            min_minor_allele_depth=3, distinguish_strands=F,
            distinguish_nucleotides=TRUE, ignore_query_Ns=TRUE, 
            include_deletions=F, include_insertions=F)

pileups.df <- list(male = 'results/merged/w1118_male.bam',
                   female = 'results/merged/w1118_female.bam') %>% 
  map(~pileup(., pileupParam = pu.param)) %>%
  map_df(as_tibble, .id='sex') %>%
  spread(sex, count, fill = 0) %>%
  arrange(seqnames, pos) %>%
  group_by(seqnames, pos) %>%
  filter(n() > 1) %>%
  ungroup()

# get male and female specific minor alleles
pileups.df2 <- pileups.df %>% 
  filter(female == 0 | male == 0) %>%
  mutate(sex.specific = ifelse(male == 0, 'female','male'))


# --------------------------

library(GenomicRanges)
library(arrow)

te.lookup <- read_tsv("~/work/TestisTEs2021/resources/te_id_lookup.curated.tsv.txt")

geps <- open_dataset('~/finalized/larval-w1118-testes/optimal_gep_membership/', format='arrow') %>%
  collect()

top_mods <- geps %>% 
  filter(qval < 0.005) %>%
  group_by(module) %>% summarise(n_tes = sum(!str_detect(X1,'FBgn'))) %>%
  arrange(-n_tes)

tep_tes <- geps %>%
  filter(module == top_mods$module[1]) %>%
  left_join(te.lookup, by=c(X1='merged_te')) %>%
  pull(gene_id) %>%
  unique()


pileups.df2 %>%
  gather(sex,count,female, male) %>%
  filter(count > 0) %>%
  ggplot(aes(count)) +
  geom_histogram(bins = 100) +
  facet_wrap(~sex.specific) +
  scale_x_log10() +
  ggtitle('distribution of read counts supporting sex-specific variants') +
  xlab('read depth supporting variant')


pileups.df2 %>%
  filter(male + female > 100) %>%
  dplyr::select(seqnames, pos, sex.specific) %>%
  distinct() %>%
  mutate(seqnames = str_remove(seqnames,'#.+')) %>%
  mutate(GEP=ifelse(seqnames %in% tep_tes,'TEP','other')) %>%
  ggplot(aes(forcats::fct_infreq(seqnames), fill=GEP)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  facet_wrap(~sex.specific, ncol=1) +
  scale_y_log10() +
  ggtitle('Number of sex-specific variant sites detected on each TE') +
  ylab('number of variant sites detected') +
  xlab("")


