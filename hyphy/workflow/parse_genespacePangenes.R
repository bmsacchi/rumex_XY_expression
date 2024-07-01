#R script

library(tidyverse)
library(vroom)

args = commandArgs(trailingOnly=TRUE)

pgMat_tx<-vroom(args[1],show_col_types = FALSE)

pgMat_PASS_ks<-filter(pgMat_tx,(flag == "PASS")) %>% 
   select(pgID,genome,id,chr,start, end, flag) %>%
   group_by(pgID,genome,flag) %>%
   filter(n() == 1) %>%
   ungroup() %>%
   distinct() %>%
   pivot_wider(names_from = genome, values_from = c(id,chr,start,end,flag)) %>%
   unite(pairID, c("id_tx_mat", "id_tx_pat"),remove = FALSE) %>%  
   distinct() %>%
   drop_na()

pgIDs_PASS<- select(pgMat_PASS_ks, pgID)

pgMat_NSOrths_ks<- pgMat_tx %>%
  select(pgID,genome,id,chr,start, end, flag) %>%
  filter(flag != "array") %>%
  group_by(pgID,genome) %>%
  filter(n() == 1) %>%
  filter(!pgID%in%pgIDs_PASS$pgID) %>%
  ungroup() %>%
  distinct() %>%
  pivot_wider(names_from = c(genome),
  	values_from = c(id,chr,start,end,flag)) %>%
  drop_na() %>%
  filter(chr_tx_mat == "X" & chr_tx_pat == "Y") %>%
  unite(pairID, c("id_tx_mat", "id_tx_pat"),remove = FALSE)

pgMatOrths_ks<-bind_rows(pgMat_PASS_ks,pgMat_NSOrths_ks) %>%
  #select(-pgID) %>%
  select(pgID,id_tx_mat:id_tx_pat)

write_csv(pgMatOrths_ks, gzfile("orthsUpdatedAnnos.csv.gz"))

