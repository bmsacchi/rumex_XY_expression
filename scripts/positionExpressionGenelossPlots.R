library(ggplot2)
library(tidyverse)
library(cowplot)
library(RColorBrewer)
#library(slider)


## read count and significance data from yx gametolog expression analysis
## and positions
reads_pg_pvals<-read_csv("data/rna_ase_results_eqtl_sept12.csv.gz")
#pvals_genes_norm <- read_csv("data/norm_rna_ase_resultsOnlySept2024.csv")
#normCounts<-read_csv("data/normalizedCounts_xyExp.csv")
## data frame with gene positions
#pgMatOrths<-read_csv("data/pgMatOrths.csv")

## plot it baby!
# genomic position on the x chromosome - x axis
# y/x ratio from reads - y axis
## unnorm'd

reads_pg_pvals2 <- reads_pg_pvals %>% mutate(sigLFC = case_when(
    res.allele.padj < 0.1 &
      abs(res.allele.log2FoldChange) > 1 ~ "True",
    TRUE ~"False" )) %>%
  mutate(readRatio = male_mean.pat/male_mean.mat) %>% 
  select(!c("39dF.mat":"43bM.pat"))


## x axis is position
## y axis is read ratio
#figure1b<-
ggplot(reads_pg_pvals2) + 
  geom_point(aes(x = start_tx_mat, y = log2(readRatio))) + 
  ylim(-20,20) +
  theme_bw() +
  xlab("Genomic position on X chromosome") +
  ylab("Read ratio (Y/X)") 
#ggsave("figures/positionExpressionPlots.png", width = 10, height = 10, units = "in", dpi = 300)

## window it!
# plot proportion of y/x underexp sig sites along the genome
# window size = 100
# using ggplot and 

# first need to combine the data with a larger file containing more of the genes (e.g. a gfF)
# or the orthologs file

tx_mat<-data.table::fread("data/pangenes_tx/tx_mat_pg_all.txt.gz") %>% 
  filter(genome == "tx_mat") %>%
  filter(chr =="X") %>% # trim the -RA from all gene names
  mutate(tx_mat = gsub("-RA", "", tx_mat)) %>%
  select(chr, tx_mat, start, end) %>%
  rename(id_tx_mat = tx_mat) %>%
  full_join(reads_pg_pvals2, by = "id_tx_mat") %>% # turn NA in sigLFC to FALSE
  mutate(sigLFC = ifelse(is.na(sigLFC), "False", sigLFC))

reads_pg_pvals3 <- tx_mat %>% mutate(window = floor(start/1000000)*1000000) %>%
  group_by(window) %>%
  summarize(propSigYoe = 
              sum(sigLFC == "True" & (res.allele.log2FoldChange > 1))/n(),
            propSigXoe = 
              sum(sigLFC == "True" & (res.allele.log2FoldChange < 1))/n()) %>%
  pivot_longer(cols = c(propSigYoe, propSigXoe), names_to = "propSig", values_to = "propSigValue")
# plot this
ggplot(reads_pg_pvals3) +
  geom_point(aes(x = window, y = propSigValue, color = propSig)) +
  theme_bw() +
  xlab("Genomic position on X chromosome") +
  ylab("Proportion of genes with significant gametolog-specific expression") +   scale_color_discrete(labels = c('log2FC < 1 \nY/X < 1','log2FC >1 \nY/X > 1' )) +
  labs(color = '')  +# make the text bigger
  theme(legend.text = element_text(size = 12)) + theme_classic()

ggsave("figures/genomicPositionYXProportion.png", width = 10, height = 8, units = "in", dpi = 300)
