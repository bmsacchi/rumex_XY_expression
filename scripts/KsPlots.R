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
hyphy_wide<-hyphy_data %>% 
  pivot_wider(id_cols = OG, names_from =taxa, values_from = c("branch","dN","dS","LB","UB","MLE")) %>%
  unite("pairID", c("branch_rhast_tx_mat", "branch_rhast_tx_pat"), sep = "_", remove = FALSE) 

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

## scatterplot of ds vs read ratio
ggplot(hyphy_reads_pg, aes(x = XYdS, y = log2_readRatio)) + 
  geom_point(alpha = 0.5) +
  #scale_y_continuous(limits=c(-10,10)) +
  ylim(-10,10) + 
  geom_smooth(method = "lm",formula=y~x) +
  #geom_abline(intercept = coef(ks_lm)[1],slope = coef(ks_lm)[2]) +
  ylab("log2 (y/x read ratio)") +
  xlab("Ks") +
  theme_bw(base_size = 20) +
  geom_hline(yintercept = 0, linetype = "dashed") 


## but is it significant?
lm_yx<-lm(formula = log2_readRatio ~ XYdS, data = hyphy_reads_pg)
lm_yx
summary(lm_yx) # p value is 0.0001773
# 


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
  ylab("difference in dN between Y and X gametologs") +
  theme_bw(base_size = 18)
#ggsave("figures/Figure2_dN_log2readRatio_ggPredict.png",h=6,w=8)

summary(lm5)# geom_text(x = -7, y = 0.25, label = lm_eqn(df), parse = TRUE)

#### because I like to make life complicated: multilevel model?


################################################################################
################################################################################
################### ks scatterplots ############################################
ggplot(hyphy_reads_pg) + geom_point(aes(x=log2_readRatio, y=XYdS)) +
  geom_smooth(aes(x=log2(yx_readRatio), y=XYdS), method = "lm") +
  xlab("log2 (y/x read ratio)") +
  ylab("Ks") +
  theme_bw(base_size = 20) +
  geom_vline(xintercept = 0, linetype = "dashed")
#ggsave("figures/Ks_log2readRatio_scatter.png", h = 6, w = 8)

#### ks vs dn
ggplot(hyphy_reads_pg) + geom_point(aes(x=XYdS, y=YX_dNdiff)) + 
  theme_bw(base_size = 20) +
  xlab("dS") +
  ylab("Y - X dN") +
  geom_smooth(aes(x=XYdS, y=YX_dNdiff), method = "lm") 


### ks across positions
ggplot(hyphy_reads_pg) + geom_point(aes(x=log2_readRatio, y=XYdS)) +
  geom_smooth(aes(x=log2(yx_readRatio), y=XYdS), method = "lm") +
  xlab("log2 (y/x read ratio)") +
  ylab("Ks") +
  theme_bw(base_size = 20) +
  geom_vline(xintercept = 0, linetype = "dashed")

################# ks window plts ###############################################
ks_mb<-hyphy_sum_match %>% mutate_at(c('start_tx_mat', 'start_tx_pat'), as.numeric) %>% 
  mutate(start_tx_mat_mb = (start_tx_mat/1e6)) %>%
  mutate(start_tx_pat_mb = (start_tx_pat/1e6))
#write_csv(ks_mb, "data/ks_mb.csv.gz")
windowmaker_1<- function(df,win,ks){df %>% filter(XYdS<ks) %>%
    group_by(chr_tx_mat) %>%
    arrange(id_tx_mat,chr_tx_mat, start_tx_mat) %>%
    mutate(Ks_win = slider::slide_dbl(XYdS, median, .before = win/2, .after = win/2,.step = 1,.complete =T))}

windowmaker_1<- function(df,win,ks){df %>% filter(XYdS<ks) %>%
    group_by(chr_tx_mat) %>%
    arrange(id_tx_pat,chr_tx_pat, start_tx_mat) %>%
    mutate(Ks_win = slider::slide_dbl(XYdS, median, .before = win/2, .after = win/2,.step = 1,.complete =T))}

ks_win100<-windowmaker_1(ks_mb,win = 100,1)
# ks_win100_x<-windowmaker_2(ks_mb, win = 100, 0.3)

ks_win50<-windowmaker_1(ks_mb,win = 50,1)
#ks_win20_x<-windowmaker_2(ks_mb, win = 20, 0.3)
#library(RColorBrewer)
ks_win30<-windowmaker_1(ks_mb,win = 30,1)

ggplot(ks_win100, aes(start_tx_pat_mb,group = chr_tx_mat)) +
  geom_line(aes(y=Ks_win, colour = chr_tx_mat)) +
  #facet_grid(~chr_tx_mat,scales="free_x") +
  #ylim(0,0.1) +
  scale_color_brewer(palette = "Dark2") +
  xlab("X-bearing Haplotype\nposition(Mb)") +
  ylab("Median Ks\n(100 gene windows)") +
  theme_classic()

ggsave("figures/Ks_coords_allchr_TX_mat_100win.jpg", h = 5, w = 9)

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
  xlim(0,0.3)
  scale_fill_manual(values = c("goodgene" = "blue", "badgene" = "red"))

## boxplot
 
hyphy_badgenes%>% filter(XYdS <0.2) %>% ggplot() + geom_boxplot(aes(x=filterType, y=(XYdS), fill = filterType)) +
  xlab("Gene type") +
  ylab("Ks") +
  theme_bw(base_size = 20) +
  scale_fill_manual(values = c("goodgene" = "blue", "badgene" = "red"))
### sample part of the data

hyphy_badgenes %>% group_by(filterType) %>% slice_sample(n = 100) %>%
  filter(XYdS <0.2) %>%
  ggplot() + geom_boxplot(aes(x=filterType, y=(XYdS), fill = filterType)) +
  xlab("Gene type") +
  ylab("Ks") +
  theme_bw(base_size = 20) +
  scale_fill_manual(values = c("goodgene" = "blue", "badgene" = "red"))
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


