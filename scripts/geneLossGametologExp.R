library(tidyverse)

###### raw data
reads_pg_pvals<-read_csv("data/rna_ase_results_eqtl_sept12.csv.gz")


###### gene loss
source("scripts/geneloss.R")
###### chisquare
# are old NC genes NOT shared in texas enriched for decreased y/x??
# 3 bins
# Y/X oe?
# Y/X equal ish?
# Y/X ue

## wait does that even make sense?
# merged the data

 

# mutate - which of these genes are not lost in nc? whihc are?
# this excludes shared homomz. or should!

reads_pg_pvals2 <- reads_pg_pvals %>% mutate(sigLFC = case_when(
  res.allele.padj < 0.1 &
    abs(res.allele.log2FoldChange) > 1 ~ "True",
  TRUE ~"False" )) %>%
  mutate(readRatio = male_mean.pat/male_mean.mat) %>% 
  mutate(shared = case_when(
    id_tx_mat %nin% ncGenesToKeep2$id_tx_mat ~ "True",
    TRUE ~ "False"
  ))

