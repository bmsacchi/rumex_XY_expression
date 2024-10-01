library(tidyverse)

## parsing relaxed selection results
relax<-read_csv("hyphy/RELAX_xy.csv.gz") %>% mutate(OG = str_remove(OG, "_relax.json"))
## rename so less wordy
colnames(relax)<-c("OG","branch","LRT","pval","param")
#head(relax)
## reduce and get rid of the repetitive stuff

## extract OG and branch info only
## mutate case when branch contains Rsag, new variable is "Rsag" for example
relax_OGs<- relax %>% 
  select(OG,branch) %>% 
  mutate(spp = case_when(
  grepl("Rsag", branch) ~ "Rsag",
  grepl("TX_maternal",branch) ~ "rhast_tx_mat",
  grepl("TX_paternal", branch) ~ "rhast_tx_pat",
  grepl("buc", branch) ~ "bucephalophorus",
  grepl("Node2", branch) ~ "node",
  grepl("Node3", branch) ~ "node",
  TRUE ~ "NA"
)) %>%  
  mutate(branch = gsub("_RA_1$", "", branch))


#relax_OGs$branch %>% gsub("_RA_1","",.)

## pivot wider so that branch is the value for the corresponding column from spp
## remove RA from all gene names
relax_OGs_wide<-relax_OGs %>% pivot_wider(id_cols = OG, names_from = spp, values_from = branch) %>% 
  unite("pairID", c("rhast_tx_mat", "rhast_tx_pat"), sep = "_", remove = FALSE)

## extract relax results and OG name only
relax_results <- relax %>% select(-branch) %>% unique() #%>% mutate(OG = str_remove(OG, "\\.json"))
# these are just the unique ones. original raw data had repeats for every branch in an OG
# oopsie from parsing the jsons. but easier to fix here.

## join the two dataframes
## remove duplciate OGs (randomly toss one)

relax_df<-left_join(relax_OGs_wide, relax_results, by = "OG") #%>% 

relax_df_dedup<- relax_df %>%
  group_by(.,across(c(node, pairID, Rsag, rhast_tx_mat, rhast_tx_pat, bucephalophorus))) %>%  
  # Group by all columns (you can also specify specific columns)
  slice_sample(n = 1) %>%             # Randomly keep one row from each group
  ungroup()                           # Ungroup the dataframe
## yay
## combine

## how many X-Y gametologs where tested? this is for multiple correction
#orths<-read_csv("hyphy/orthsUpdatedAnnos.csv.gz")
reads_pg<-read_csv("data/rna_ase_results_eqtl_sept12.csv.gz",show_col_types =FALSE) %>% 
  select(1:27) 
relax_reads_pg<-left_join(reads_pg, relax_df_dedup, by = "pairID") %>% 
  drop_na() %>%
  #mutate(padj = p.adjust(pval, method = 'fdr', n = length(pval))) %>%
  mutate(sigRelax = ifelse(pval < 0.05, "True", "False")) %>%
  mutate(sigLFC = case_when(
    res.allele.padj < 0.1 &
      abs(res.allele.log2FoldChange) > 1 ~ "True",
    TRUE ~"False" )) %>%
  mutate(log2male_mean.mat = log2(male_mean.mat)) %>%
  mutate(log2male_mean.pat = log2(male_mean.pat)) %>%
  mutate(type = case_when(
    param > 1 ~ "intensified",
    param < 1 ~ "relaxed",
  ))

### summary table
## sigTest, sigRelax, LogFC > or < 0
summarize_relax<- relax_reads_pg %>% 
  group_by(sigRelax, sigLFC) %>% dplyr::summarise(n = n())

tab1<-table(relax_reads_pg$sigRelax, relax_reads_pg$sigLFC)
print(tab1)
## fishers test
test1<-fisher.test(tab1)
test1 # yeah not significant no surprise there!!!
#### direction
summarize_relax2<- relax_reads_pg %>% filter(sigTest =="sig" & sigRelax == "True") #%>%
  group_by(sigRelax, res.allele.log2FoldChange>0,type) %>% dplyr::summarise(n = n()) %>% ungroup() 

# transform into contigency table format
tab2<-table(summarize_relax2$type, summarize_relax2$res.allele.log2FoldChange>0)
tab2

test2<-fisher.test(tab2)
test2 #NOP  
