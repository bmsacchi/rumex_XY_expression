### gene loss
library(tidyverse, quietly = TRUE)
library(vroom, quietly = TRUE)
library(data.table, quietly = TRUE)

### gene loss and hemizyogus genes
# same approach to NC paper

mat_synorths_raw<-data.table::fread("data/pangenes_tx/tx_mat_pg_synorths.txt.gz") #%>% mutate(across(start:end, as.numeric))
mat_nsorths_raw<-data.table::fread("data/pangenes_tx/tx_mat_pg_all.txt.gz") #%>% mutate(across(start:end, as.numeric))
#identify which entries have nsorthos for txmat pat and rsal genomes
# same as prev filtering for NC, but i use a slightly different function bc there were only 3 genomes being compared in that case (here there's many more)
nsorths<-slr_rename_TX_PAR(mat_nsorths_raw) %>% filter_at(vars(c("tx_mat","tx_pat","Rsal")),any_vars(grepl("\\*",.))) 
# remove any pgIDs w nsorths at any of our three 'genomes of interest'

cols<-colnames(mat_synorths_raw)
onlysyn<-slr_rename_TX_PAR(mat_synorths_raw) %>% filter(!(pgID%in%nsorths$pgID)) %>%
  mutate(across(tx_mat:Rsal, ~str_replace(.,"-RA$","")))

hasanortho<-onlysyn %>% 
  filter((tx_mat!=""|tx_pat!="")&Rsal!="")

## n that has an ortho X and Y only
filter(hasanortho, tx_mat!="" & chr =="X") %>% nrow() # 
#filter(mat_synorths_raw, chr =="X") %>% nrow()
## 1675
# !="" means "not absent". =""means "gene absent"
# this stressed me out for a minute!
notlost<-onlysyn %>% filter((tx_mat!="" & tx_pat!="" & Rsal!=""))

## hemizygous?
tx_mat_loss<-onlysyn %>% filter(tx_mat!="" & tx_pat=="" & Rsal!="")

## gene loss by chr

summary_loss<- tx_mat_loss %>% dplyr::group_by(chr) %>% dplyr::summarise(n_loss = n()) #272 hemizygous
print(summary_loss)
tx_yloss<- filter(tx_mat_loss) %>% filter(chr == "X")  %>% #376
  mutate(region = if_else(as.numeric(start) > 220000000, "PAR", "SLR"))
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
tophit_geneonly<- hits %>% 
  dplyr::group_by(qseqid) %>%
  dplyr::mutate(the_rank = rank(-pident, ties.method = "random")) %>%
  dplyr::ungroup() %>%  # Make sure to ungroup after the operation
  filter(the_rank == 1) %>%
  dplyr::select(-the_rank) %>% mutate(alignedpropn = length/qlen)


# complete and partial loss a la NC paper

tophit_partial<- tophit_geneonly %>% filter(alignedpropn <= 0.5 & alignedpropn > 0) ## 130 partially lost
tophit_notloss<-tophit_geneonly %>% filter(alignedpropn > 0.5)  ## only 24 genes "not lost"
## 376 lost before blast filtering, putatively


hemiz_blast_confirm <- tx_yloss %>% filter(!tx_mat%in%tophit_notloss$qseqid)
#write_csv(hemiz_blast_confirm, "data/tx_yloss_blastconfirmed_Aug2024.csv")

hemiz_blast_complete<- tx_yloss %>% 
  filter(!tx_mat%in%tophit_notloss$qseqid &
        !tx_mat%in%tophit_partial$qseqid) %>% 
  mutate(txtype = "complete")

#write_csv(hemiz_blast_complete, "data/tx_yloss_blastconfirmed_complete_Aug2024.csv")

hemiz_blast_partial<-tx_yloss %>% filter(!tx_mat%in%tophit_notloss$qseqid &
          tx_mat%in%tophit_partial$qseqid) %>%  dplyr::mutate(txtype = "partial")

#write_csv(hemiz_blast_partial, "data/tx_yloss_blastconfirmed_partial_Aug2024.csv")
hemiz_blast_all<- rbind(hemiz_blast_complete,hemiz_blast_partial)

heads<-c("genespace_yloss","blast_confirmed_loss","complete_loss","partial_loss","overlap_nc_completeloss")
ylo<-(filter(summary_loss,chr =="X"))$n_loss


vals<-c(376,352,222,130,125)
resTable<-data.frame(heads,vals)
print(resTable)

## upload the completely lost genes from NC
# 167 genes lost in NC

nc_complete_loss<-read.csv("../rumex_pangenome_annotation/complete_yloss_Sept2023.csv") %>% mutate(nctype = "complete")
nc_partial_loss<-read.csv("../rumex_pangenome_annotation/yloss_partial_Sept2023.csv") %>% mutate(nctype ="partial")
nc_all_loss<-rbind(nc_complete_loss,nc_partial_loss) %>% filter(chr =="X") #%>% select(repGene)
 # 344 genes partially or completely lost nc

source("scripts/slr_rename.R")
nc_all_loss_pars<-slr_rename_PARs(nc_all_loss)
nc_all_loss <- nc_all_loss_pars %>% filter(region == "OldSLR")
#167 genes completely lost in nc
#344 genes partially or completely lost nc

tx_yloss_nc_overlap<-hemiz_blast_all %>% 
  full_join(nc_all_loss, by = c("nc_hap1" = "repGene")) %>% 
  mutate(nctype = replace_na(nctype,"no overlap")) %>%
  mutate(txtype = replace_na(txtype,"no overlap")) 

geneloss_summary<- tx_yloss_nc_overlap %>% group_by(txtype,nctype) %>% 
  dplyr::summarise(n = n())
contingency_table <- xtabs(n ~ txtype + nctype, data = geneloss_summary)
#q: how do i convert NA to something else
#a: use replace_na

#contingency_table <- xtabs(n ~ txtype + nctype, data = geneloss_summary)

print(contingency_table)

chisq_test <- chisq.test(contingency_table)
print(chisq_test) # significant relationship between the two variables
# complete and partial gene loss not independent between the cytotypes
### tukey test



#partially lost genes in tx more likely to be partially lost in nc, and vice versa
# same for completely lost genes

### chi square for overlap
tx_yloss_nc_overlap_all<-hemiz_blast_all %>% left_join(nc_all_loss, by = c("nc_hap1" = "repGene")) #%>% filter(!is.na(nctype))
## replace NA with "no overlap"
geneloss_summary_all<- tx_yloss_nc_overlap_all %>% group_by(txtype,nctype) %>% summarise(n = n())
contingency_table <- xtabs(n ~ txtype + nctype, data = geneloss_summary)
print(contingency_table)
#### what genes are "not lost" from both genomes
#### ugh this is SOOSOSOOO annoying
nc_hap1<- read_delim("data/nc_beds/hap1.bed.gz",delim = "\t",col_names = c("chr","start","end","gene")) #%>% dplyr::select(repGene)
## gotta redo the shiz from previous study UGHHHHHHHHHHHHHHHHHHHH
hap1_raw<-data.table::fread("../rumex_pangenome_annotation/pg_hap1_wide_synonly.txt")
hap1_nonsyn<-data.table::fread("../rumex_pangenome_annotation/pg_hap1_wide.txt")
nsorths_hap1<-slr_rename_PARs(hap1_nonsyn) %>% filter_all(any_vars(grepl("\\*",.)))  
onlysyn_hap1<-slr_rename_PARs(hap1_raw) %>% filter(!(pgID%in%nsorths_hap1$pgID)) # 31289 genes remain
nc_notlost1<-onlysyn_hap1 %>% filter((hap1!="" & hap2!="" & salicifolius!=""))

##nc_notlost<- nc_hap1 %>% filter(!gene%in%nc_all_loss$nc_hap1 & chr == "X") 
tx_notlost<-notlost %>% filter(chr =="X") %>% filter(!nc_hap1%in%nc_all_loss$hap1) # 1260
nc_notlost<-nc_notlost1 %>% filter(chr =="X") %>% filter(!repGene%in%hemiz_blast_all$nc_hap1) # 2610
## so many
## how many shared?
shared_notlost<-intersect(tx_notlost$nc_hap1, nc_notlost$repGene) %>% length() # 994
## add to dataframe

geneloss_summary2 <- geneloss_summary %>% ungroup() %>% add_row(txtype = "no overlap", nctype = "no overlap", n = shared_notlost)
contingency_table <- xtabs(n ~ txtype + nctype, data = geneloss_summary2)
chisq_test <- chisq.test(contingency_table)
print(chisq_test) # significant relationship between the two variables
## weeeee
contingency_table
### go back and redo add to main data for chrissakes i am soooo annoyed with this bullshit!
geneloss_summary2
write_csv(geneloss_summary2, "data/geneloss_summary_Oct2024.csv")
#####
## take into account SLRs


