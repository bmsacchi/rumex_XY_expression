library(tidyverse)
source("scripts/generalFunctions.R")

reads_pg_pvals<-read_csv("data/rna_ase_results_eqtl_sept12.csv.gz")
# normcounts<-read_csv("data/norm_rna_ase_resultsOnlySept2024.csv.gz")
normcounts<-read_csv("data/norm_counts_autosomefactors_eQTLleafRhast.csv.gz")
# 

#completely_silenced<-filter(reads_pg_pvals, res.allele.padj < 0.1 & res.allele.log2FoldChange < 1 & male_mean.pat == 0) %>%
#  select(id_tx_mat)
x_overexpressed<-filter(reads_pg_pvals,
                        res.allele.padj < 0.1 & 
                          res.allele.log2FoldChange < 1) %>%
  select(id_tx_mat, id_tx_pat)

y_overexpressed<-filter(reads_pg_pvals, res.allele.padj < 0.1 &
                          res.allele.log2FoldChange > 1) %>% 
  select(id_tx_mat,id_tx_pat)

#
#
#
normcounts_clean<-#as.data.frame(normalized_counts) %>% 
  normcounts %>%
  dplyr::rename("Geneid" ="rowname") %>%
  rowwise %>% 
  mutate(maleMean = mean(c_across(contains("M",ignore.case = FALSE)))) %>%
  mutate(maleVar = sd(c_across(contains("M",ignore.case = FALSE)))) %>%
  mutate(maleSE = plotrix::std.error(c_across(contains("M",ignore.case = FALSE)))) %>%
  mutate(femaleMean = mean(c_across(contains("F",ignore.case = FALSE)))) %>%
  mutate(femaleVar = sd(c_across(contains("F",ignore.case = FALSE)))) %>%
  mutate(femaleSE = plotrix::std.error(c_across(contains("F",ignore.case = FALSE)))) %>%
  mutate(totalCount =sum(c_across("39dF":"43bM"))) %>%
  mutate(totalVar = sd(c_across("39dF":"43bM"))) %>%
  ungroup() %>%
  filter(totalCount >20) #%>%
#
#
#

hemiz_blast_confirm<-read_csv("data/tx_yloss_blastconfirmed_Aug2024.csv.gz") %>% 
  select(tx_mat)
dc_combo <- normcounts_clean %>% mutate(genetype = case_when(
  #Geneid%in%completely_silenced$id_tx_mat ~ "completeYsilenced",
  Geneid%in%x_overexpressed$id_tx_mat ~ "xOverexp",
  Geneid%in%y_overexpressed$id_tx_mat ~ "yOverexp",
  Geneid%in%hemiz_blast_confirm$tx_mat ~ "hemizygous", TRUE ~ "other")) %>% 
  filter(genetype != "other")
#
#
# ggplot(dc_combo, aes(x = femaleMean, y = maleMean, color = genetype)) + 
#   geom_point(alpha = 0.5, size = 3) +
#   ylim(0,500)+
#   xlim(0,500)+
#   geom_errorbar(aes(ymin = maleMean-maleSE, ymax = maleMean+maleSE, color = genetype)) +
#   geom_errorbarh(aes(xmin = femaleMean-femaleSE, xmax = femaleMean+femaleSE, color = genetype)) +
#   geom_smooth(method="lm", fullrange = TRUE) +
#   geom_abline(slope = 1, intercept = 0, linetype=1) + 
#   geom_abline(slope = 0.5, intercept = 0, linetype=2) +
#   theme_bw() +
#   xlab("Average female expression") +
#   ylab("Average male expression") +
#   labs(color= NULL) +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 16),
#         legend.text= element_text(size =12)) +
#   scale_color_discrete(labels = c("hemizygous"="Hemizygous", "xOverexp" = "Partial or complete\n Y silencing"))
# 
# #ggsave("plots/combined_dc_plot_hemi_ysilence_updated.png", h = 6, w =9)
# #ggsave("plots/combined_dc_plot_hemi_ysilence_limit500.png", h = 6, w =9)

dc_combo %>% filter(genetype == "yOverexp") %>% 
  ggplot(aes(x = femaleMean, y = maleMean)) + 
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = maleMean-maleSE, ymax = maleMean+maleSE)) +
  geom_errorbarh(aes(xmin = femaleMean-femaleSE, xmax = femaleMean+femaleSE)) +
  geom_smooth(method="lm", fullrange = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype=1) + 
  geom_abline(slope = 0.5, intercept = 0, linetype=2) +
  theme_bw() +
  xlab("Average female expression") +
  ylab("Average male expression") +
  labs(color= NULL) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text= element_text(size =12))

#ggsave("plots/yOverexp_dc_plot.png", h = 5, w = 7)

dc_combo %>% filter(genetype != "yOverexp") %>% ggplot(aes(x = femaleMean, y = maleMean, color = genetype)) + 
  geom_point(alpha = 0.5, size = 3) +
  geom_errorbar(aes(ymin = maleMean-maleSE, ymax = maleMean+maleSE, color = genetype)) +
  geom_errorbarh(aes(xmin = femaleMean-femaleSE, xmax = femaleMean+femaleSE, color = genetype)) +
  geom_smooth(method="lm", fullrange = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype=1) + 
  geom_abline(slope = 0.5, intercept = 0, linetype=2) +
  theme_bw() +
  xlab("Average female expression") +
  ylab("Average male expression") +
  labs(color= NULL) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text= element_text(size =12)) +
  scale_color_discrete(labels = c("hemizygous"="Hemizygous", "xOverexp" = "Partial or complete\n Y silencing"))
ggsave("figures/Fig3_hemiz_xOverexp_dc_plot.png", h = 5, w = 7)



ggplot(dc_combo, aes(x = femaleMean, y = maleMean, color = genetype)) + 
  geom_point(alpha = 0.5, size = 3) +
  ylim(0,500)+
  xlim(0,500)+
  geom_errorbar(aes(ymin = maleMean-maleSE, ymax = maleMean+maleSE, color = genetype)) +
  geom_errorbarh(aes(xmin = femaleMean-femaleSE, xmax = femaleMean+femaleSE, color = genetype)) +
  geom_smooth(method="lm", fullrange = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype=1) + 
  geom_abline(slope = 0.5, intercept = 0, linetype=2) +
  theme_bw() +
  xlab("Average female expression") +
  ylab("Average male expression") +
  labs(color= NULL) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text= element_text(size =12)) +
  scale_color_discrete(labels = c("hemizygous"="Hemizygous", "xOverexp" = "Partial or complete\n Y silencing"))


####

dc_combo %>% filter(genetype != "yOverexp") %>% ggplot(aes(x = femaleMean, y = maleMean, color = genetype)) + 
  geom_point(alpha = 0.5, size = 3) +
  geom_errorbar(aes(ymin = maleMean-maleSE, ymax = maleMean+maleSE, color = genetype)) +
  geom_errorbarh(aes(xmin = femaleMean-femaleSE, xmax = femaleMean+femaleSE, color = genetype)) +
  geom_smooth(method="lm", fullrange = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype=1) + 
  geom_abline(slope = 0.5, intercept = 0, linetype=2) +
  theme_bw() +
  xlab("Average female expression") +
  ylab("Average male expression") +
  labs(color= NULL) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text= element_text(size =12)) +
  scale_color_discrete(labels = c("hemizygous"="Hemizygous", "xOverexp" = "Partial or complete\n Y silencing"))
ggsave("figures/combined_dc_plot_hemi_ysilence_updated.png", h = 6, w =9)
