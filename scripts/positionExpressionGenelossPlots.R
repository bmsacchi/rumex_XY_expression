library(ggplot2)
library(tidyverse)
library(cowplot)
#library(RColorBrewer)
#library(slider)
library("khroma")

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
  select(!c("39dF.mat":"43bM.pat")) %>%
  mutate(log2male_mean.mat = log2(male_mean.mat)) %>%
  mutate(log2male_mean.pat = log2(male_mean.pat)) %>%
  #filter(male_mean.mat > 10 & male_mean.pat >10) #%>% 
  filter(log2male_mean.pat != "-Inf" & log2male_mean.mat != "-Inf") %>%
  mutate(start_mb = start_tx_mat/1000000)
reads_pg_pvals2$sigLFC<-factor(reads_pg_pvals2$sigLFC, levels = c("True","False"))

vibrant<- color("vibrant")

## x axis is position
## y axis is read ratio
#figure1b<-
colors <- vibrant(7)[c(4,7)]

#colors
ggplot(reads_pg_pvals2) + 
  geom_point(aes(x = start_mb, y = log2(readRatio), color = sigLFC)) + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.7) +
  ylim(-13,13) +
  theme_classic() +
  scale_color_manual(values = colors)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.position = "none") +
  xlab("Genomic position on X chromosome (Mb)") +
  ylab("Read ratio (Y/X)") 
#ggsave("figures/Figure1B_2025Jan.png", width = 8, height = 6, units = "in", dpi = 300)

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
  mutate(sigLFC_all = case_when(is.na(sigLFC) ~ "False", .default = sigLFC ))

reads_pg_pvals3 <- tx_mat %>% mutate(window = floor(start/1000000)) %>%
  group_by(window) %>%
  summarize(propSigYoe = 
              sum(sigLFC_all == "True" & (res.allele.log2FoldChange > 1))/n(),
            propSigXoe = 
              sum(sigLFC_all == "True" & (res.allele.log2FoldChange < 1))/n()) %>%
  pivot_longer(cols = c(propSigYoe, propSigXoe), names_to = "propSig", values_to = "propSigValue") %>%
  filter(propSigValue >0)
# plot this

colors2<- vibrant(7)[c(2,5)]
ggplot(reads_pg_pvals3) +
  geom_point(aes(x = window, y = propSigValue, color = propSig)) +
  #theme_bw() +
  xlab("Genomic position on X chromosome (Mb)") +
  ylab("Proportion of genes with significant \ngametolog-specific expression") +  
  #scale_color_discrete(labels = c('log2FC < 1 \nY/X < 1','log2FC >1 \nY/X > 1' )) +
  #labs(color = '')  +# make the text bigger
  scale_x_continuous(breaks = seq(0,350, 50)) +
  theme_classic() +
  scale_color_manual(values = colors2) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.position = "none") 


#ggsave("figures/Figure1B_windows_2025Jan.png", width = 8, height = 4, units = "in", dpi = 300)

####### Repeat with Y position #############
# just first plot
reads_pg_pvals2 <- reads_pg_pvals %>% mutate(sigLFC = case_when(
  res.allele.padj < 0.1 &
    abs(res.allele.log2FoldChange) > 1 ~ "True",
  TRUE ~"False" )) %>%
  mutate(readRatio = male_mean.pat/male_mean.mat) %>% 
  select(!c("39dF.mat":"43bM.pat")) %>%
  mutate(log2male_mean.mat = log2(male_mean.mat)) %>%
  mutate(log2male_mean.pat = log2(male_mean.pat)) %>%
  #filter(male_mean.mat > 10 & male_mean.pat >10) #%>% 
  filter(log2male_mean.pat != "-Inf" & log2male_mean.mat != "-Inf") %>%
  mutate(start_mb = start_tx_pat/1000000)
reads_pg_pvals2$sigLFC<-factor(reads_pg_pvals2$sigLFC, levels = c("True","False"))

vibrant<- color("vibrant")

## x axis is position
## y axis is read ratio
#figure1b<-
colors <- vibrant(7)[c(4,7)]

#colors
ggplot(reads_pg_pvals2) + 
  geom_point(aes(x = start_mb, y = log2(readRatio), color = sigLFC)) + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.7) +
  ylim(-13,13) +
  theme_classic() +
  #geom_vline(xintercept = 45, linetype = 2) +
  scale_color_manual(values = colors)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.position = "none") +
  xlab("Genomic position on Y chromosome (Mb)") +
  ylab("Read ratio (Y/X)") 
ggsave("figures/Figure1B_3C_Ypos_2025Jan.png", width = 8, height = 6, units = "in", dpi = 300)


##### median sliding windows
## UGH

windowmaker<- function(df,win){df %>%
    arrange(start_mb) %>%
    mutate(medianRatio = slider::slide_index_dbl(.x = readRatio, .i = start_mb, .f = median, .before = win/2, .after = win/2, .complete = TRUE))
}

reads_pg_pvals4 <- reads_pg_pvals2 %>% windowmaker(win=10) #%>%


ggplot(reads_pg_pvals4) + 
  geom_point(aes(x = start_mb, y = log2(medianRatio))) + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.7) +
  ylim(-3,3) +
  theme_classic() +
  #geom_vline(xintercept = 45, linetype = 2) +
  scale_color_manual(values = colors)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.position = "none") +
  xlab("Genomic position on Y chromosome \n(100 Mb sliding windows)") +
  ylab("Median Y/X read ratio (log2)") 

### calc median in _ gene windows? 1 mb windows?
## combine?

ggplot(reads_pg_pvals4) + 
  geom_point(aes(x = start_mb, y = log2(readRatio), color = sigLFC)) + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.7) +
  geom_line(aes(x=start_mb,y=log2(medianRatio))) +
  ylim(-13,13) +
  theme_classic() +
  #geom_vline(xintercept = 45, linetype = 2) +
  scale_color_manual(values = colors)+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.position = "none") +
  xlab("Genomic position on Y chromosome (Mb)") +
  ylab("Read ratio (Y/X)") 
#ggsave("figures/Figure3C_Ypos_windows50Mb_2025Jan.png", width = 8, height = 6, units = "in", dpi = 300)



