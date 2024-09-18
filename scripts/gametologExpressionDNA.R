# cleaning and producing useable data for the gametolog expression analysis
# no plots just data wrangling
library(tidyverse)
library(DESeq2)
library(gtools)
source("scripts/generalFunctions.R") # load relevant data cleaning and processing functions

pgMatOrths<-read_csv("data/pgMatOrths.csv")
dnaCountsClean<-read_csv("data/dnaCountsClean.csv.gz")
#rnaCountsClean<-read_csv("data/rnaCountsClean.csv.gz")

#### data wrangling

setdiff(colnames(rnaCountsClean),colnames(dnaCountsClean))
# same samples!
dnaSamplenamesClean<-colnames(select(dnaCountsClean, !c("Geneid":"male_mean"))) # save sample names for later
sampleInfo<-as_tibble_col(dnaSamplenamesClean, column_name="Sample") %>% mutate(.,Sex = if_else(grepl("*M", Sample),"M", "F")) # compile sample information for deseq2

### DNA read ratio estimation
#do this first. For filtering reasons! We want to remove sites with high mapping bias
#This will also help eliminate genes where there's weird amounts of mapping on Y for females. 


eqtl_x_dna<- filter(dnaCountsClean, Chr =="X") %>% select(c(Geneid,totalCountsGene:161)) ## extract X gene counts
eqtl_y_dna<- filter(dnaCountsClean, Chr =="Y") %>% select(c(Geneid,totalCountsGene:161)) ## extract Y gene counts

reads_pg_dna<-left_join(pgMatOrths,eqtl_x_dna, by = c("id_tx_mat"="Geneid")) %>%
  left_join(.,eqtl_y_dna, by = c("id_tx_pat" = "Geneid"),suffix = c(".mat",".pat")) %>% 
  drop_na() %>% select(-flag_tx_pat,-flag_tx_mat) %>% 
  unique() %>%
  rowwise() %>% 
  mutate(totalCountAll=(totalCountsGene.mat+totalCountsGene.pat)) %>%
  ungroup() %>% filter(totalCountAll >20) %>%
  relocate(c(pairID:male_mean.mat, totalCountsGene.pat:male_mean.pat,totalCountAll)) 
reads_pg_dna_noPAR<-reads_pg_dna %>% filter(start_tx_pat > 45000000) # remove pseudoautosomal region genes

badgenes_dna<-reads_pg_dna_noPAR %>% # remove genes with female Y gene mapping (mismapping)
  filter(female_mean.pat >= 20 & chr_tx_pat =="Y") %>% select(pairID,id_tx_mat,id_tx_pat)

goodgenes_dna <-reads_pg_dna_noPAR %>% filter(female_mean.pat < 20 & chr_tx_pat =="Y") # genes that don't have this mismapping

## filter the dna data
ase_data_dna<-reads_pg_dna_noPAR %>% 
  filter(pairID%nin%badgenes_dna$pairID) %>% # remove bad genes 1450 down to 647
  select(c(pairID,contains("M", ignore.case = FALSE ))) %>% column_to_rownames(var ="pairID")


## ase matrix for DNA data
IDs<-colnames(ase_data_dna) # colnames
sorted_col_names <- mixedsort(IDs) # sort colnames
ase_data_dna <- ase_data_dna[, sorted_col_names] # turn to matrix
#colnames(ase_data)
# yay!!
sample<-str_remove(sorted_col_names,"..at") # tm1 tm1 etc. remove suffices
allele<-str_extract(sorted_col_names,".at$") # mat pat mat pat extract suffices

colInfo<-data.frame(sorted_col_names,sample,allele) # create column info df

file_name<-"data/dna_ase_results_eqtl_noPAR.csv"


### DNA ase counts
if (!file.exists(file_name)) {
 
  print("Running the command because the file does not exist.")

 
  ase_matrix_dna<-as.matrix(ase_data_dna) # convert to matrix
  design<- ~sample+allele # design formula
  
  m <- 75 # number of samples

  dds <- DESeqDataSetFromMatrix(ase_data_dna, colInfo, design) # create deseq object
  sizeFactors(dds) <- rep(1, 2*m) # set size factors to unity
  dds <- DESeq(dds, fitType="local") # fit model
  
  resultsNames(dds) # get results names

  res.allele <- results(dds, name="allele_pat_vs_mat") # get results
  head(res.allele$log2FoldChange) # check log2foldchange


  adj_pvals_dna<-res.allele$padj # get adjusted pvals


  pvals_genes_dna<-data.frame(res.allele@rownames,res.allele$padj,res.allele$log2FoldChange,res.allele$lfcSE) # create df with pvals
  reads_pg_pvals_dna<-inner_join(reads_pg_dna_noPAR, pvals_genes_dna, by =c("pairID" ="res.allele.rownames")) %>% # join with read data
    mutate(sigTest = 
             case_when((res.allele.padj < 0.1) ~ "sig",
                       TRUE ~"nonsig" )) %>% relocate(c(pairID:totalCountAll,res.allele.padj:sigTest))
  #write_csv(reads_pg_pvals_dna, "data/dna_ase_results_eqtl_noPAR.csv")
} else { 
  print(paste("Skipping command because", file_name, "exists."))
  #head(reads_pg_pvals)
  reads_pg_pvals_dna<- read_csv("data/dna_ase_results_eqtl_noPAR.csv") # read in data
}

fold_change_cutoffs <- c(0.1,0.5, 1, 1.5, 2) # fc cutoffs for summary table
obs_tables_dna <- lapply(fold_change_cutoffs, function(cutoff) { 
  get_obs_table(reads_pg_pvals_dna, cutoff)
})
combined_table_dna <- do.call(rbind, obs_tables_dna) # combine tables
print(combined_table_dna) # print table
#write_csv(combined_table_dna, "data/n_xy_oe_DNA_eQTL_lfcRanges.csv")
# genes with mapping bias. 
biasedgenes_dna<-filter(reads_pg_pvals_dna, sigTest == "sig" & abs(res.allele.log2FoldChange) > 0.5 ) #%>% select(pairID, id_tx_mat, id_tx_pat)





