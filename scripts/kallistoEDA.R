library(tidyverse)

## read in kallisto output
kallisto_df <- read_delim("expression_03_2026/kallisto_out/merged_est_counts.tsv") %>%
  mutate(target_id = str_replace(target_id, "-RA$", ""))


### clean up colnames
kallisto_colnames <- colnames(kallisto_df)

## regex to fix colnames
kallisto_colnames_repaired <- str_extract(kallisto_colnames, "target_id|\\d+[a-z][FM](?=LR)")

### change colnames
colnames(kallisto_df) <- kallisto_colnames_repaired
## prev samples that did not pass QC
rnacolsremove<-c("67eF","14bF","45aF","40aF","24fM","35aM")

# remove col if name matches rnacolsremove
kallisto_df_clean <- kallisto_df %>% select(-(all_of(rnacolsremove))) 

kallisto_df_sum<-kallisto_df_clean %>%rowwise() %>% mutate( # Calculate mean counts
  femaleCountsTotal = sum(c_across(contains("F", ignore.case = FALSE))),  # Calculate total female counts
  female_mean = (femaleCountsTotal / num_female_cols), # Calculate mean female counts dynamically
  female_SD = stats::sd(c_across(contains("F", ignore.case = FALSE))),  # Calculate standard deviation
  female_SE = plotrix::std.error(c_across(contains("F", ignore.case = FALSE))),  # Calculate standard error
  maleCountsTotal = sum(c_across(contains("M", ignore.case = FALSE))),  # Calculate total male counts
  male_mean = (maleCountsTotal / num_male_cols),  # Calculate mean male counts dynamically
  male_SE = plotrix::std.error(c_across(contains("M", ignore.case = FALSE))),
  totalCountsGene = (femaleCountsTotal+maleCountsTotal)) %>%
  relocate(totalCountsGene, 
           femaleCountsTotal, female_mean, female_SE,
           maleCountsTotal, male_mean, male_SE)

#dplyr remove ### read in orthology information
pgMatOrthsRaw<- read_csv("data/pgMatOrths.csv.gz") 

pgMatOrths<-pgMatOrths %>% 
  select(pairID:chr_tx_pat) %>%
  unique()
  # filter(flag_tx_mat == "PASS" & flag_tx_pat == "PASS") 
### read in annotation information

## merge pgmatorths and kallisto

#Y genes

num_male_cols <- ncol(x_exp %>% select(contains("M", ignore.case = FALSE)))
num_female_cols<- ncol(x_exp %>% select(contains("F", ignore.case = FALSE)))
#ncol(x_exp)


y_exp <- kallisto_df_sum %>% filter(target_id %in% pgMatOrths$id_tx_pat)
x_exp <- kallisto_df_sum %>% filter(target_id %in% pgMatOrths$id_tx_mat)
df1<-left_join(pgMatOrths,x_exp, 
                             by = c("id_tx_mat"="target_id")) %>%
  left_join(.,y_exp, by = c("id_tx_pat" = "target_id"),
          suffix = c(".mat",".pat")) 

ase_data<-df1 %>% 
  select(c(pairID,contains("M", ignore.case = FALSE ))) %>% 
  column_to_rownames(var ="pairID")
IDs<-colnames(ase_data)
sorted_col_names <- gtools::mixedsort(IDs)
ase_data <- ase_data[, sorted_col_names]
#colnames(ase_data)
# yay!!
sample<-str_remove(sorted_col_names,"..at") # tm1 tm1 etc.
allele<-str_extract(sorted_col_names,".at$") # mat pat mat pat

colInfo<-data.frame(sorted_col_names,sample,allele) # wooo

ase_matrix<-as.matrix(ase_data)
ase_data <- round(ase_data)
design<- ~sample+allele

m <- 75
library(DESeq2)

dds <- DESeqDataSetFromMatrix(ase_data, colInfo, design)
sizeFactors(dds) <- rep(1, 2*m)
dds <- DESeq(dds, fitType="local")

resultsNames(dds)

res.allele <- results(dds, name="allele_pat_vs_mat")
head(res.allele$log2FoldChange)
adj_pvals<-res.allele$padj


pvals_genes<-data.frame(res.allele@rownames,res.allele$padj,
                        res.allele$log2FoldChange,
                        res.allele$lfcSE) 


#### compare crossmapping with kallisto vs. previous tools

## how many female reads map to Y?
badgenes<-df1 %>% 
  filter(female_mean.pat > 20) #%>% select(ta) # remove genes with improper female Y mapping
# some may be on Y PAR - which is fromm 0-45MB. Disregard and remove all for now!
badgenes_x <- pgMatOrths %>% select(id_tx_mat,id_tx_pat) %>% 
  filter(id_tx_pat%in%badgenes$id_tx_pat) %>% select(id_tx_mat)

pargenes <- pgMatOrthsRaw %>% filter(start_tx_mat < 450000)
#1500... thats a lot oh no  

`%nin%` <- purrr::negate(`%in%`)

df1_nomismap<-df1 %>% filter(id_tx_mat%nin%badgenes_x$id_tx_mat)  %>% 
  filter(id_tx_pat%nin%badgenes) %>%
  filter(id_tx_mat%nin%pargenes$id_tx_mat)
## 1267 from 2778 genesbadgenes_x## 1267 from 2778 genes
reads_pg_pvals<-inner_join(df1_nomismap, pvals_genes, 
                           by = c("pairID" ="res.allele.rownames")) %>%
  mutate(
    sigTest = # just do significance, leave the log2foldchane to be filtered later
      case_when(res.allele.padj < 0.05 ~ "sig", 
                TRUE ~"nonsig" )) 

  
#write_csv(reads_pg_pvals,"data/rna_ase_results_eqtl_sept12.csv.gz")

library(ggplot2)
library(tidyverse)
library(cowplot)
#library(RColorBrewer)
library(khroma)
df2 <- reads_pg_pvals %>% mutate(sigLFC = case_when(
  res.allele.padj < 0.05 &
    abs(res.allele.log2FoldChange) > 1 ~ "True",
  TRUE ~ "False")) %>%
  mutate(log2male_mean.mat = log2(male_mean.mat)) %>%
  mutate(log2male_mean.pat = log2(male_mean.pat)) %>%
  #filter(male_mean.mat > 10 & male_mean.pat >10) #%>% 
  filter(log2male_mean.pat != "-Inf" & log2male_mean.mat != "-Inf") 


df2$sigLFC<-factor(df2$sigLFC, levels = c("True","False"))

colnames(df2)


figure1a <- ggplot(df2_old,aes(x=male_mean.mat,y=male_mean.pat)) + 
  geom_point(alpha = 0.5, size =1,aes(color = sigLFC)) +
  geom_errorbar(aes(color = sigLFC,
                    ymin = ((male_mean.pat)-(male_SE.pat)), 
                    ymax = ((male_mean.pat)+(male_SE.pat)))) +
  geom_errorbarh(aes(color = sigLFC,
                     xmin = male_mean.mat-male_SE.mat,
                     xmax = male_mean.mat+male_SE.mat)) +
  geom_abline(slope = 1, intercept = 0, linetype=1) +
  theme_classic() +
  geom_smooth(method ="lm") #+ # all points regression line
  #geom_smooth(method = "lm", 
#              data = function(x) { filter(x, sigLFC == "True") },
 #             se = T, fullrange = TRUE) +
  xlab("mean expression of X\n gametolog in males") +
  ylab("mean expression of Y\n gametolog in males") 

figure1a
## we've just removed X genes that mismap to Y only, so the other situation is not covered...
## keeping the cutoff of 20 read count or greater

  #scale_x_continuous(transform = "log2",limits = c(NA,4000),breaks = c(0.0625,2,64,2048)) +
  #scale_y_continuous(transform = "log2",limits = c(NA,4000),breaks = c(0.0625,2,64,2048)) +
  #labs(color= "p-val < 0.1 & |log2FC| >1 ") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.position = "none") +
  scale_color_manual(values= colors) # from the vibran

  
#### compare biases?
reads_pg_pvals_old<-read_csv("data/rna_ase_results_eqtl_sept12.csv.gz")

### intersect:

df2_old <- reads_pg_pvals_old %>% mutate(sigLFC = case_when(
  res.allele.padj < 0.05 &
    abs(res.allele.log2FoldChange) > 1 ~ "True",
  TRUE ~ "False"))

## compare:

## sig lfc >1
pos_lfc<-filter(df2,res.allele.log2FoldChange > 1 & sigLFC =="True") %>% 
  select(pairID)
neg_lfc<-filter(df2, res.allele.log2FoldChange < -1 & sigLFC =="True") %>% 
  select(pairID)

pos_lfc_old<- filter(df2_old, res.allele.log2FoldChange > 1 & sigLFC =="True")
neg_lfc_old <- filter(df2_old, res.allele.log2FoldChange < -1 & sigLFC =="True")

match_pos <- pos_lfc %>% filter(pairID %in% pos_lfc_old$pairID)
match_pos_old <- pos_lfc_old %>% filter(pairID %in% pos_lfc$pairID)
match_neg <- neg_lfc %>% filter(pairID %in% neg_lfc_old$pairID)

cross_match<- pos_lfc %>% filter(pairID %in% neg_lfc_old$pairID)
cross_match2<- neg_lfc %>% filter(pairID %in% pos_lfc_old$pairID)
