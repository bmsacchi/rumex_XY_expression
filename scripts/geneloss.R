### gene loss
library(tidyverse, quietly = TRUE)
library(vroom, quietly = TRUE)
library(data.table, quietly = TRUE)

### gene loss and hemizyogus genes
# same approach to NC paper

mat_synorths_raw<-data.table::fread("data/pangenes_tx/tx_mat_pg_synorths.txt")
mat_nsorths_raw<-data.table::fread("data/pangenes_tx/tx_mat_pg_all.txt") 
#identify which entries have nsorthos for txmat pat and rsal genomes
# same as prev filtering for NC, but i use a slightly different function bc there were only 3 genomes being compared in that case (here there's many more)
nsorths<-mat_nsorths_raw %>% filter_at(vars(c("tx_mat","tx_pat","Rsal")),any_vars(grepl("\\*",.))) 
# remove any pgIDs w nsorths at any of our three 'genomes of interest'

cols<-colnames(mat_synorths_raw)
onlysyn<-mat_synorths_raw %>% filter(!(pgID%in%nsorths$pgID)) %>%
  mutate(across(all_of(cols), ~str_replace(.,"-RA$","")))

hasanortho<-onlysyn %>% 
  filter((tx_mat!=""|tx_pat!="")&Rsal!="")

# !="" means "not absent". =""means "gene absent"
# this stressed me out for a minute!
notlost<-onlysyn %>% filter((tx_mat!="" & tx_pat!="" & Rsal!=""))

## hemizygous?
tx_mat_loss<-onlysyn %>% filter(tx_mat!="" & tx_pat=="" & Rsal!="")

## gene loss by chr

summary_loss<- tx_mat_loss %>% group_by(chr) %>% summarise(n_loss = n()) #272 hemizygous
print(summary_loss)
tx_yloss<- filter(tx_mat_loss) %>% filter(chr == "X") #376
y_loss_genelist_tx_mat<-tx_yloss %>% dplyr::select(tx_mat) 

#write.table(y_loss_genelist_tx_mat,"data/genelist_tx_mat.txt",quote=FALSE,sep="\t",row.names=F,col.names = F)

###### BLAST results #######
#blast was performed - see blast_yloss_patgenome.sh
blast_yloss <- read_tsv("data/blast_pat_genome.tsv",
                        col_names = c("qseqid", "sseqid", "pident", "length", 
                                      "mismatch", "gapopen", "qstart", "qend", 
                                      "qlen", "sstart", "send", "slen", "evalue", 
                                      "bitscore", "score"))

#blast_pos<-full_join(blast_yloss,mat_synorths_raw,by=c("qseqid"="tx_mat"))  %>% add_count(qseqid) 
### which genes, out of any, are no hits in hap2 genome
#hits_x<- blast_yloss %>% filter(sseqid =="Y")

## filter by hits matching to Y only!!!!!!
hits<-blast_yloss %>% filter(sseqid == "Y") #%>% dplyr::select(qseqid) %>% distinct() # distinct is fully redundant here
nohit<- anti_join(y_loss_genelist_tx_mat,hits, by = c("tx_mat"="qseqid")) ## 222 with no hit
tophit_geneonly<- hits %>% group_by(qseqid) %>%
  mutate(the_rank  = rank(-pident, ties.method = "random")) %>%
  filter(the_rank == 1) %>%
  dplyr::select(-the_rank) %>% mutate(alignedpropn = length/qlen)


# complete and partial loss a la NC paper

tophit_partial<- tophit_geneonly %>% filter(alignedpropn <= 0.5 & alignedpropn > 0) ## 130 partially lost
tophit_notloss<-tophit_geneonly %>% filter(alignedpropn > 0.5)  ## only 24 genes "not lost"
## 376 lost before blast filtering, putatively


hemiz_blast_confirm <- y_loss_genelist_tx_mat %>% filter(!tx_mat%in%tophit_notloss$qseqid)

#write_csv(hemiz_blast_confirm, "data/tx_yloss_blastconfirmed_Aug2024.csv")

hemiz_blast_complete<- y_loss_genelist_tx_mat %>% filter(!tx_mat%in%tophit_notloss$qseqid &
                                              !tx_mat%in%tophit_partial$qseqid)
#write_csv(hemiz_blast_complete, "data/tx_yloss_blastconfirmed_complete_Aug2024.csv")

hemiz_blast_partial<-y_loss_genelist_tx_mat %>% filter(!tx_mat%in%tophit_notloss$qseqid &
                                            tx_mat%in%tophit_partial$qseqid)
#write_csv(hemiz_blast_partial, "data/tx_yloss_blastconfirmed_partial_Aug2024.csv")

heads<-c("genespace_yloss","blast_confirmed_loss","complete_loss","partial_loss","overlap_nc_completeloss")
ylo<-(filter(summary_loss,chr =="X"))$n_loss


vals<-c(376,352,222,130,125)
resTable<-data.frame(heads,vals)


## upload the completely lost genes from NC
# 167 genes lost in NC
nc_complete_loss<-read.csv("../rumex_pangenome_annotation/complete_yloss_Sept2023.csv") %>% mutate(type = "complete")
nc_partial_loss<-read.csv("../rumex_pangenome_annotation/yloss_partial_Sept2023.csv") %>% mutate(type ="partial")
nc_all_loss<-rbind(nc_complete_loss,nc_partial_loss)
#head(nc_complete_loss)
#ncGenesToKeep<-nc_complete_loss %>% filter(chr =="X") %>% select(repGene) 
ncGenesToKeep<-nc_all_loss %>% filter(chr =="X") #%>% select(repGene)

dim(ncGenesToKeep)
#167 genes completely lost in nc
#344 genes partially or completely lost nc

tx_yloss_nc_overlap<-tx_yloss %>% left_join(ncGenesToKeep, by = c("nc_hap1" = "repGene")) %>% filter(!is.na(type)) # 183 as expected
summarise()

source("scripts/generalFunctions.R")
#nc_loss_unshared <- ncGenesToKeep2 %>% filter(repGene%nin%tx_yloss$tx_mat)

## all genes
#allncGenesToKeep<-nc_complete_loss %>% select(repGene)

#write_csv(tx_yloss_nc_overlap, "tx_yloss_overlapNCcomplete_loss.csv")
#write_csv(overlap_all, "tx_")

####### further gene loss analysis ######

# load set of hemizygous genes, generated from jpop_tx_exp_geneloss.Rmd


#tx_yloss<-read_csv("data/genelist_tx_mat.txt",col_names = FALSE)


#tx_yloss_nc_overlap<-read_csv("data/tx_yloss_overlapNCcomplete_loss.csv")
#hemiz_blast_confirm<-read_csv("data/tx_yloss_blastconfirmed_Aug2024.csv")

# 376 genes lost before blast filtering
#intersect(tx_yloss_nc_overlap$repGene, hemiz_blast_confirm$tx_mat) %>% length()
#intersect(hemiz_blast_confirm$tx_mat, tx_yloss_nc_overlap$repGene) %>% length()
# all 125 nc overlap lost genes are also confirmed to be either completely or partially lost via blast
#hemiz_blast_partial<-read.csv("data/tx_yloss_blastconfirmed_partial_Aug2024.csv")
#hemiz_blast_complete<-read.csv("data/tx_yloss_blastconfirmed_complete_Aug2024.csv")

# overlaps once more
#intersect(hemiz_blast_complete$tx_mat, tx_yloss_nc_overlap$repGene) %>% length() # 112
#intersect(hemiz_blast_partial$tx_mat, tx_yloss_nc_overlap$repGene) %>% length() # 13
#setdiff(hemiz_blast_complete$tx_mat, tx_yloss_nc_overlap$repGene) %>% length() #110/222 completely lost genes are not lost in NC
#setdiff(hemiz_blast_partial$tx_mat, tx_yloss_nc_overlap$repGene) %>% length() # 117/130 partially lost genes are not lost in NC

#creating a data frame using dplyr where one column is gene name, the other is category (e.g. complete loss, partial loss), plus a third column for overlap with NC gene loss

geneloss_categorized <- hemiz_blast_confirm %>% mutate(category = case_when(
  tx_mat %in% hemiz_blast_complete$tx_mat ~ "complete loss",
  tx_mat %in% hemiz_blast_partial$tx_mat ~ "partial loss",
  TRUE ~ "no loss")) %>% 
  mutate(nc_overlap = ifelse(tx_mat %in% tx_yloss_nc_overlap$repGene, "yes", "no"))
testTable<-table(geneloss_categorized$category,geneloss_categorized$nc_overlap)
test1<-chisq.test(testTable)
#X-squared = 56.827, df = 1, p-value = 4.759e-14
geneloss_summary<- geneloss_categorized %>% group_by(category,nc_overlap) %>% summarise(n = n())              
geneloss_summary
