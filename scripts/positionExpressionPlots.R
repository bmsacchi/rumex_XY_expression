library(ggplot2)
library(tidyverse)
library(cowplot)
library(RColorBrewer)

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
  mutate(readRatio = male_mean.pat/male_mean.mat)

## x axis is position
## y axis is read ratio
#figure1b<-
ggplot(reads_pg_pvals2) + 
  geom_point(aes(x = start_tx_mat, y = log2(readRatio))) + 
  ylim(-20,20) +
  theme_bw() +
  xlab("Genomic position on X chromosome") +
  ylab("Read ratio (Y/X)") 
ggsave("figures/positionExpressionPlots.png", width = 10, height = 10, units = "in", dpi = 300)

##
