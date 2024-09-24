library(tidyverse, quietly = TRUE)
library(vroom, quietly = TRUE)
library(slider)

### read in X-Y ortholog ASE data

reads_pg<-read_csv("data/rna_ase_results_eqtl_sept12.csv.gz",show_col_types =FALSE)

### read in ks dataframe

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
hyphy_wide<-hyphy_data %>% 
  pivot_wider(id_cols = OG, names_from =taxa, values_from = c("branch","dN","dS","LB","UB","MLE")) %>%
  unite("pairID", c("branch_rhast_tx_mat", "branch_rhast_tx_pat"), sep = "_", remove = FALSE) 

### Add tx mat and tx pat Ks values
hyphy_sum <- hyphy_wide %>% 
  mutate(XYdS = dS_rhast_tx_mat+dS_rhast_tx_pat) %>% 
  mutate(dnds_calc_tx_mat = dN_rhast_tx_mat/dS_rhast_tx_mat) %>%
  mutate(dnds_calc_tx_pat = dN_rhast_tx_pat/dS_rhast_tx_pat) %>%
  mutate(YX_omegaDiff = MLE_rhast_tx_pat-MLE_rhast_tx_mat) %>%
  mutate(YX_dNdiff = dN_rhast_tx_pat-dN_rhast_tx_mat) %>%
  select(c(OG,pairID,branch_rhast_tx_mat,branch_rhast_tx_pat,dnds_calc_tx_mat,dnds_calc_tx_pat,
           dS_rhast_tx_mat,dS_rhast_tx_pat,dN_rhast_tx_mat,dN_rhast_tx_pat,
           XYdS,MLE_rhast_tx_mat,MLE_rhast_tx_pat,YX_omegaDiff,YX_dNdiff)) %>% drop_na()
#################################################################
## sigtest variable needs fixing - it's fucked up for some reason
hyphy_sum_match <- left_join(hyphy_sum, reads_pg) %>% 
  drop_na() %>% 
  mutate(sigTest = 
           case_when(res.allele.padj < 0.1 & abs(res.allele.log2FoldChange) > 1 ~ "sig",
                     TRUE ~"nonsig" )) %>% 
  mutate(yx_readRatio = male_mean.pat/male_mean.mat) %>%
  #mutate(yx_readRatio_plus = yx_readRatio+1) %>% 
  filter(yx_readRatio != 0)
### Ks plots

dfKs<- hyphy_sum_match %>% 
  select(XYdS,yx_readRatio,sigTest,res.allele.lfcSE,res.allele.log2FoldChange,res.allele.padj) %>%
  mutate(log2_readRatio=log2(yx_readRatio))
ks_lm<-lm(formula = log2_readRatio ~ XYdS, data = dfKs)
#plot(ks_lm)
summary(ks_lm)

ggplot(dfKs) +geom_histogram(aes(log2_readRatio))

ggplot(dfKs, aes(x = XYdS, y = log2_readRatio)) + 
  geom_point(alpha = 0.5) +
  scale_y_continuous(limits=c(-10,10)) +
  geom_smooth(method = "lm",formula=y~x) +
  geom_abline(intercept = coef(ks_lm)[1],slope = coef(ks_lm)[2]) +
  ylab("log2 (y/x read ratio)") +
  xlab("Ks") +
  theme_bw(base_size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed") 

df_dN<- hyphy_sum_match %>% 
  select(YX_dNdiff,yx_readRatio,sigTest,res.allele.lfcSE,res.allele.log2FoldChange,res.allele.padj,XYdS) %>%
  mutate(log2_readRatio=log2(yx_readRatio))

### compare models

#lm1<-lm(formula = YX_dNdiff ~ yx_readRatio, data=df_dN)
#lm2 <- lm(formula = YX_dNdiff ~ yx_readRatio + XYdS, data=df_dN)
#lm3 <- lm(formula = YX_dNdiff ~ yx_readRatio + XYdS + XYdS:yx_readRatio, data=df_dN)
lm4<-lm(formula = YX_dNdiff ~ log2_readRatio, data=df_dN)
lm5<-lm(formula = YX_dNdiff ~ log2_readRatio + XYdS, data=df_dN)
lm6<-lm(formula = YX_dNdiff ~ log2_readRatio + XYdS + XYdS:log2_readRatio, data=df_dN)

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
  ylab("difference in dN between Y and X gametologs") +
  theme_bw(base_size = 18)
#ggsave("figures/Figure2_dN_log2readRatio_ggPredict.png",h=6,w=8)

summary(lm5)# geom_text(x = -7, y = 0.25, label = lm_eqn(df), parse = TRUE)
################################################################################

################### ks scatterplots ############################################
ggplot(dfKs) + geom_point(aes(x=log2(yx_readRatio), y=XYdS)) +
  geom_smooth(aes(x=log2(yx_readRatio), y=XYdS), method = "lm") +
  xlab("log2 (y/x read ratio)") +
  ylab("Ks") +
  theme_bw(base_size = 20) +
  geom_vline(xintercept = 0, linetype = "dashed")
ggsave("figures/Ks_log2readRatio_scatter.png", h = 6, w = 8)

### ks across positions
ggplot(dfKs) + geom_point(aes(x=log2(yx_readRatio), y=XYdS)) +
  geom_smooth(aes(x=log2(yx_readRatio), y=XYdS), method = "lm") +
  xlab("log2 (y/x read ratio)") +
  ylab("Ks") +
  theme_bw(base_size = 20) +
  geom_vline(xintercept = 0, linetype = "dashed")

################# ks window plts ###############################################
ks_mb<-hyphy_sum_match %>% mutate_at(c('start_tx_mat', 'start_tx_pat'), as.numeric) %>% 
  mutate(start_tx_mat_mb = (start_tx_mat/1e6)) %>%
  mutate(start_tx_pat_mb = (start_tx_pat/1e6))

windowmaker_1<- function(df,win,ks){df %>% filter(XYdS<ks) %>%
    group_by(chr_tx_mat) %>%
    arrange(id_tx_mat,chr_tx_mat, start_tx_mat) %>%
    mutate(Ks_win = slider::slide_dbl(XYdS, median, .before = win/2, .after = win/2,.step = 1,.complete =T))}

# windowmaker_2<- function(df,win,ks){df %>% filter(Ks<ks) %>%
#     group_by(chr2) %>%
#     arrange(name2,chr2, start2_mb) %>%
#     mutate(Ks_win = slider::slide_dbl(Ks, median, .before = win/2, .after = win/2,.step = 1,.complete =T))}

ks_win100<-windowmaker_1(ks_mb,win = 100,1)
# ks_win100_x<-windowmaker_2(ks_mb, win = 100, 0.3)

ks_win50<-windowmaker_1(ks_mb,win = 50,1)
#ks_win20_x<-windowmaker_2(ks_mb, win = 20, 0.3)
#library(RColorBrewer)
ks_win30<-windowmaker_1(ks_mb,win = 30,1)

ggplot(ks_win100, aes(start_tx_mat_mb,group = chr_tx_mat)) +
  geom_line(aes(y=Ks_win, colour = chr_tx_mat)) +
  #facet_grid(~chr_tx_mat,scales="free_x") +
  #ylim(0,0.1) +
  scale_color_brewer(palette = "Dark2") +
  xlab("X-bearing Haplotype\nposition(Mb)") +
  ylab("Median Ks\n(100 gene windows)") +
  theme_classic()

ggsave("figures/Ks_coords_allchr_TX_mat_100win.jpg", h = 5, w = 9)


# ggplot(filter(ks_win20_x,!grepl("scaffold", chr2)), aes(start2_mb,group = chr2)) +  
#   geom_line(aes(y=Ks_win, colour = chr2)) +
#   facet_grid(~chr2,scales="free_x") +
#   ylim(0,0.05) +
#   scale_color_brewer(palette = "Dark2") +
#   xlab("X-bearing Haplotype\nposition(Mb)") +
#   ylab("Median Ks\n(100 gene windows)") 
# ggsave("Ks_coords_allchr_TX_hap1_20win.pdf", h = 5, w = 9)

#Ks of DNA badgenes ?!

# bad_ks<- filter(ks_data_separated, id_tx_mat%in%badgenes_dna$id_tx_mat & id_tx_pat%in%badgenes_dna$id_tx_pat) %>% filter(chr1  =="Y" & chr2 =="X")
# ### recall, Ks calculated b/w genes - much of these DNA reads are not genic!!!
# bad_ks_rna<- filter(ks_data_separated, id_tx_mat%in%badgenes_x$id_tx_mat & id_tx_pat%in%badgenes$Geneid) %>% filter(chr1  =="Y" & chr2 =="X")
# 
# baddf<-setdiff(badgenes_x$id_tx_mat,badgenes_dna$id_tx_mat) # 960 genes in RNA are not mismapping in the DNA data - maternal
# reversebaddf<-setdiff(badgenes_dna$id_tx_mat, badgenes_x$id_tx_mat) #404 genes in dna are not mismapping in the RNA
# 
# baddf2<-setdiff(badgenes$Geneid,badgenes_dna$id_tx_pat) #1144
# reversebaddf2<-setdiff(badgenes_dna$id_tx_pat,badgenes$Geneid) #405


# 
# 
# 

              