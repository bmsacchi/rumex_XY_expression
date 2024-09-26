library(tidyverse, quietly = TRUE)
library(vroom, quietly = TRUE)
library(slider)


########################### raw data ###########################################
### in X-Y ortholog ASE data

reads_pg<-read_csv("data/rna_ase_results_eqtl_sept12.csv.gz",show_col_types =FALSE) 

## read in ks dataframe 

hyphy_data<-read_csv("hyphy/all_ks_updated.csv.gz", show_col_types = FALSE) %>%  
  mutate(clean_branch = str_remove(branch, "_RA_1")) %>%
  mutate(og = str_remove(OG, "\\.json")) %>%
  select(-branch,-OG) %>% 
  dplyr::rename(.,branch = clean_branch) %>%
  dplyr::rename(.,OG = og) %>%
  mutate(taxa = ifelse(grepl("Rsag", branch), "sagittatus",
                       ifelse(grepl("TX_maternal",branch), "rhast_tx_mat",
                              ifelse(grepl("TX_paternal", branch), "rhast_tx_pat", 
                                     ifelse(grepl("buc", branch), "bucephalophorus", 
                                            ifelse(grepl("Node2", branch), "node2", 
                                                   ifelse(grepl("Node3", branch), "node3", NA))))))) 
### convert to wide format
## and remove duplicate orthogroups randomly
hyphy_wide<-hyphy_data %>% 
  pivot_wider(id_cols = OG, names_from =taxa, values_from = c("branch","dN","dS","LB","UB","MLE")) %>%
  unite("pairID", c("branch_rhast_tx_mat", "branch_rhast_tx_pat"), sep = "_", remove = FALSE) %>%
  group_by(.,across(c(branch_node2, branch_sagittatus, pairID, branch_rhast_tx_mat, branch_rhast_tx_pat, branch_bucephalophorus,branch_node3))) %>%  
  slice_sample(n = 1) %>% ungroup()


##################### Summarize dn and ds #################################
hyphy_sum <- hyphy_wide %>% 
  mutate(XYdS = dS_rhast_tx_mat+dS_rhast_tx_pat) %>% 
  mutate(dnds_calc_tx_mat = dN_rhast_tx_mat/dS_rhast_tx_mat) %>%
  mutate(dnds_calc_tx_pat = dN_rhast_tx_pat/dS_rhast_tx_pat) %>%
  mutate(YX_omegaDiff = MLE_rhast_tx_pat-MLE_rhast_tx_mat) %>%
  mutate(YX_dNdiff = dN_rhast_tx_pat-dN_rhast_tx_mat) %>%
  select(c(OG,pairID,branch_rhast_tx_mat,branch_rhast_tx_pat,dnds_calc_tx_mat,dnds_calc_tx_pat,
           dS_rhast_tx_mat,dS_rhast_tx_pat,dN_rhast_tx_mat,dN_rhast_tx_pat,
           XYdS,MLE_rhast_tx_mat,MLE_rhast_tx_pat,YX_omegaDiff,YX_dNdiff)) %>% drop_na()
########################### y/x, ds and dn analysis ###########################

###combine y/x expression data 
## sigtest variable needs fixing 
hyphy_reads_pg <- left_join(hyphy_sum, reads_pg) %>% 
  drop_na() %>% 
  mutate(sigTest = 
           case_when(res.allele.padj < 0.1 & abs(res.allele.log2FoldChange) > 1 ~ "sig",
                     TRUE ~"nonsig" )) %>% 
  mutate(yx_readRatio = male_mean.pat/male_mean.mat) %>%
  #mutate(yx_readRatio_plus = yx_readRatio+1) %>% 
  filter(yx_readRatio != 0) %>%
  mutate(log2_readRatio=log2(yx_readRatio)) 

##Ks and gametolog expression plots 

## model
ks_lm<-lm(formula = log2_readRatio ~ XYdS, data = dfKs) 
#plot(ks_lm)

summary(ks_lm)
## histogram

#ggplot(dfKs) +geom_histogram(aes(log2_readRatio))

################################################################################
################### ds expression scatterplots #################################

ggplot(hyphy_reads_pg, aes(x = XYdS, y = log2_readRatio)) + 
  geom_point(alpha = 0.5) +
  #scale_y_continuous(limits=c(-10,10)) +
  ylim(-10,10) + 
  geom_smooth(method = "lm",formula=y~x) +
  #geom_abline(intercept = coef(ks_lm)[1],slope = coef(ks_lm)[2]) +
  ylab("log2 (y/x read ratio)") +
  xlab("XYdS") +
  theme_bw(base_size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed") 
ggsave("figures/Ks_log2readRatio_scatter.png", h = 6, w = 8)

## but is it significant?
lm_yx<-lm(formula = log2_readRatio ~ XYdS, data = hyphy_reads_pg)
lm_yx
summary(lm_yx) # p value is 0.0001773
# 

###########################ks vs dn ############################################

### compare models

#lm1<-lm(formula = YX_dNdiff ~ yx_readRatio, data=df_dN)
#lm2 <- lm(formula = YX_dNdiff ~ yx_readRatio + XYdS, data=df_dN)
#lm3 <- lm(formula = YX_dNdiff ~ yx_readRatio + XYdS + XYdS:yx_readRatio, data=df_dN)
lm4<-lm(formula = YX_dNdiff ~ log2_readRatio, data=hyphy_reads_pg)
lm5<-lm(formula = YX_dNdiff ~ log2_readRatio + XYdS, data=hyphy_reads_pg)
lm6<-lm(formula = YX_dNdiff ~ log2_readRatio + XYdS + XYdS:log2_readRatio, data=hyphy_reads_pg)

AIC(lm4,lm5,lm6)
summary(lm6)
summary(lm5)
step(lm6)
#VIF(lm6) # VIF of 4 or larger means GTFO!
## yay!
AIC(lm4,lm5) # mor
###
require(ggiraph)
require(ggiraphExtra)
require(plyr)

ggPredict(lm5,se=TRUE) + theme_bw() +  
  xlab("log2 (y/x read ratio)") +
  ylab("dN Y - dN X") +
  theme_bw(base_size = 18)
ggsave("figures/Figure2_dN_log2readRatio_ggPredict.png",h=6,w=8)

summary(lm5)# geom_text(x = -7, y = 0.25, label = lm_eqn(df), parse = TRUE)

################### ks plot positions all genes ################################

## merge annotation info with Ks data

# load annotation
anno<-vroom::vroom("data/merged_TX_noMatPARlarge_txanno.gtf",
                   col_names=c("seqname","source","feature","start","end","score","strand","frame","attribute")) %>%
  filter(feature =="transcript") %>%
  separate(attribute, into =c("a","b"),sep = ";") %>% 
  separate(a, into =c("a","Geneid"),sep = " ") %>% 
  select(-c(a,b)) %>%
  mutate(Geneid = str_remove_all(Geneid, "\"")) %>%
  mutate(Geneid = str_remove_all(Geneid, "-RA"))

# merge Ks data hwere available!
hyphy_anno<- left_join(hyphy_sum, anno, by = c("branch_rhast_tx_pat" = "Geneid")) %>% 
  left_join(.,anno, by = c("branch_rhast_tx_mat" = "Geneid"), suffix = c(".pat",".mat")) %>%
  mutate(start.mat.mb = start.mat/1e6) %>%
  mutate(start.pat.mb = start.pat/1e6) 

# plot along Y position
hyphy_anno %>% filter(seqname.pat == "Y") %>%
  filter(XYdS <1) %>%
  ggplot() + 
  geom_point(aes(x=start.pat.mb, y=XYdS)) +
  xlab("Y chromosome position (Mb)") +
  geom_vline(xintercept = 45,linetype = "dashed", color = 'red') +
  theme_bw(base_size = 20)
  
ggsave("figures/ds_Ypos.png", h = 6, w = 8)
## plot along X position
hyphy_anno %>% filter(seqname.mat == "X") %>%
  filter(XYdS <1) %>%
  ggplot() + 
  geom_point(aes(x=start.mat.mb, y=XYdS)) +
  xlab("X chromosome position (Mb)") +
  theme_bw(base_size = 20) 
ggsave("figures/ds_Xpos.png", h = 6, w = 8)

## left_join with reads_pg
pg_select<-reads_pg %>% select(1:27) %>% 
  mutate(sigLFC = case_when(
    res.allele.padj < 0.1 &
      abs(res.allele.log2FoldChange) > 1 ~ "True",
    TRUE ~"False" )) %>%
  mutate(expType = case_when(sigLFC =="True" & res.allele.log2FoldChange <1 ~ "y_underexp",
                             sigLFC =="True" & res.allele.log2FoldChange >1 ~ "y_overexp",
                             sigLFC == "False" ~ NA))
 # mutate(sigLFC = as.factor(sigLFC)) %>% 
#  mutate(pairID = str_remove(pairID, "-RA"))

hyphy_anno_reads <- hyphy_anno %>% left_join(.,pg_select, by = c("pairID")) %>% 
  mutate(expType = ifelse(is.na(expType), NA, expType))

### plot points along Y and color by type of expression

hyphy_anno_reads %>% filter(seqname.pat == "Y") %>% filter(XYdS <0.5) %>%
  ggplot() + 
  geom_point(aes(x=start.pat.mb, y=XYdS, color = expType)) +
  xlab("Y chromosome position (Mb)") +
  geom_vline(xintercept = 45,linetype = "dashed", color = 'red') +
  theme_bw(base_size = 20) +
  scale_color_manual(values=alpha(c("darkorange","darkblue"),c(0.75, 0.75)), 
                     na.value = alpha("grey",0.3), labels = c("Y/X > 1", "Y/X < 1"))
ggsave("figures/ds_Ypos_exp.png", h = 6, w = 8)

write_csv(hyphy_anno_reads, "data/hyphy_anno_reads.csv.gz")

################# ks window plots ##############################################

#write_csv(ks_mb, "data/ks_mb.csv.gz")
windowmaker_mat<- function(df,win,ks){df %>% filter(XYdS<ks) %>%
    group_by(seqname.mat) %>%
    arrange(branch_rhast_tx_mat,seqname.mat, start.mat.mb) %>%
    mutate(Ks_win = slider::slide_dbl(XYdS, median, .before = win/2, .after = win/2,.step = 1,.complete =T))}

windowmaker_pat<- function(df,win,ks){df %>% filter(XYdS<ks) %>%
    group_by(seqname.pat) %>%
    arrange(branch_rhast_tx_pat,seqname.pat, start.pat.mb) %>%
    mutate(Ks_win = slider::slide_dbl(XYdS, median, .before = win/2, .after = win/2,.step = 1,.complete =T))}

ks_win100<-windowmaker_pat(hyphy_anno,win = 100,0.5) %>% drop_na(seqname.pat)
#ks_win50<-windowmaker_pat(hyphy_anno,win = 50,0.5)


ggplot(ks_win100, aes(start.pat.mb,group = seqname.pat)) +
  geom_line(aes(y=Ks_win, colour = seqname.pat)) +
  facet_grid(~seqname.pat,scales="free_x") +
  #ylim(0,0.1) + 
  scale_color_brewer(palette = "Dark2") +
  xlab("Y-bearing Haplotype\nposition(Mb)") +
  ylab("Median Ks\n(100 gene windows)") +
  theme_classic()

ggsave("figures/ds_allchr_pat_100win.png", h = 5, w = 9)

## y plot window only
ks_win100 %>% filter(seqname.pat == "Y") %>% 
  ggplot(aes(start.pat.mb,group = seqname.pat)) +
  geom_line(aes(y=Ks_win, colour = seqname.pat)) +
  facet_grid(~seqname.pat,scales="free_x") +
  #ylim(0,0.1) + 
  scale_color_brewer(palette = "Dark2") +
  xlab("Y chromosome position(Mb)") +
  ylab("Median Ks\n(100 gene windows)") +
  geom_vline(xintercept = 45,linetype = "dashed", color = 'red') +
  theme_classic()
ggsave("figures/ds_Ychr_pat_100win.png", h = 5, w = 9)


### compare with "badgenes" 
# badgenes<-read_csv("data/badgenes.csv")

rnaCountsClean<-read_csv("data/rnaCountsClean.csv.gz")
badgenes<-rnaCountsClean %>% 
  filter(female_mean>= 20 & Chr =="Y") %>% select(Geneid)
pgMatOrths<-read_csv("data/pgMatOrths.csv.gz")
# unite hyphy wide data with pangenes

hyphy_sum_filter<- hyphy_sum %>% group_by(pairID) %>% 
  slice_sample(n = 1) %>% 
  ungroup()

hyphy_sum_pos <- hyphy_sum_filter %>% 
  left_join(.,pgMatOrths, by = "pairID") %>%
  select(OG,pairID,branch_rhast_tx_mat,branch_rhast_tx_pat,dnds_calc_tx_mat,dnds_calc_tx_pat,
         dS_rhast_tx_mat,dS_rhast_tx_pat,dN_rhast_tx_mat,dN_rhast_tx_pat,
         XYdS,MLE_rhast_tx_mat,MLE_rhast_tx_pat,YX_omegaDiff,YX_dNdiff) %>% drop_na()

hyphy_badgenes<-hyphy_sum_pos %>% 
 mutate(filterType = case_when(branch_rhast_tx_pat %in% badgenes$Geneid ~ "badgene",
                               TRUE ~ "goodgene")) 
## plot distributinos of XYdS according to gene type
ggplot(hyphy_badgenes) + geom_density(aes(x=XYdS, fill = filterType), alpha = 0.5) +
  xlab("Ks") +
  ylab("Density") +
  theme_bw(base_size = 20) +
  xlim(0,0.3) #+
#  scale_fill_manual(values = c("goodgene" = "blue", "badgene" = "red"))
ggsave("figures/mismapping_kS_density.png", h = 6, w = 8)

## histogram

ggplot(hyphy_badgenes) + 
  geom_histogram(aes(x=XYdS, fill = filterType), alpha = 0.5) +
  xlab("Ks") +
  ylab("Density") +
  theme_bw(base_size = 20) +
  xlim(0,0.3) +
  scale_fill_manual(values = c("goodgene" = "blue", "badgene" = "red"))
ggsave("figures/mismapping_kS_histogram.png", h = 6, w = 8)

## no apparent difference
## odd!
rnaCountsClean<-read_csv("data/rnaCountsClean.csv.gz") %>% select(Geneid,Chr,Start,female_mean)

hyphy_allCounts<-hyphy_badgenes %>% 
  left_join(.,rnaCountsClean, by = c("branch_rhast_tx_pat" = "Geneid")) 

hyphy_allCounts %>% filter(Chr =="Y") %>%
  ggplot() + geom_point(aes(x=female_mean,y=XYdS, color = filterType)) + 
  ylim(0,1) +
  xlim(0,3000) +
  ggtitle("dS vs female read count mean on Y chromosome \nfilters: dS <1 and female mean <3000")
ggsave("figures/mismapping_kS.png", h = 6, w = 8)


