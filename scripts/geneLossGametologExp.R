library(tidyverse)

###### raw data
reads_pg_pvals<-read_csv("data/rna_ase_results_eqtl_sept12.csv.gz")
## replace all instances of -RA with nothing
pgMatOrths<-read_csv("data/pangenes_tx/tx_mat_pg_synorths.txt.gz") %>% 
  mutate(across(all_of(colnames(.)), ~str_replace(.,"-RA$","")))


reads_pg_pvals_join<-left_join(reads_pg_pvals, pgMatOrths, by = c("id_tx_mat" = "tx_mat"))
###### gene loss
source("scripts/geneloss.R")
###### chisquare
# TX genes uniquely deleted from NC X and ones that still remain in NC have diff x/y exp?

# 3 bins
# Y/X oe?
# Y/X equal ish?
# Y/X ue

## wait does that even make sense?
# merged the data


reads_pg_pvals2 <- reads_pg_pvals_join %>% mutate(sigLFC = case_when(
  res.allele.padj < 0.1 &
    abs(res.allele.log2FoldChange) > 1 ~ "True",
  TRUE ~"False" )) %>%
  mutate(readRatio = male_mean.pat/male_mean.mat) %>% 
  mutate(log2readRatio = log2(readRatio)) %>%
  filter(!is.na(nc_hap1)) %>%
  mutate(shared = case_when(
    is.na(nc_hap2) ~ "NC_lost",
    TRUE ~ "NC_kept"
  ))

## explorator plots before chisqr/summarize into expression bins
ggplot(reads_pg_pvals2, aes(x = shared, y = log2readRatio, fill = shared)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Expression of TX genes uniquely deleted from NC X and ones that still remain in NC",
       x = "Gene status",
       y = "log2(readRatio)")
## geom hist
ggplot(reads_pg_pvals2, aes(x = log2readRatio, fill = shared)) +
  geom_histogram(binwidth = 0.5, position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Expression of TX genes uniquely deleted from NC X and ones that still remain in NC",
       x = "log2(readRatio)",
       y = "Count")

# looks weird without the log, naturally...
ggplot(reads_pg_pvals2, aes(x = readRatio, fill = shared)) +
  geom_histogram(binwidth = 0.5, position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Expression of TX genes uniquely deleted from NC X and ones that still remain in NC",
       x = "log2(readRatio)",
       y = "Count")
#### ok bin bin bin bin time!

reads_pg_pvals3 <- reads_pg_pvals2 %>% mutate(readBin = case_when(
  log2readRatio < 0 ~ "Y>X",
  #log2readRatio >= -1 & log2readRatio <= 1 ~ "Y=X",
  log2readRatio > 0 ~ "Y<X"
)) %>% filter(sigTest == "sig")

## chisqr
summaryTab<- reads_pg_pvals3 %>% group_by(readBin, shared) %>% summarize(n = n()) %>% ungroup()
summaryTab
contingencyTable<-xtabs(n~readBin+shared, data = summaryTab)
print(contingencyTable)
chisq.test(contingencyTable) # not sigificant!
