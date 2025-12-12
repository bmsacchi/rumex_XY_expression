library(ggplot2)
library(tidyverse)
library(cowplot)
#library(RColorBrewer)
library(khroma)

## use palette
vibrant<- color("vibrant")
#vibrant
#vibrant(2)
vibrant()
reads_pg_pvals<-read_csv("data/rna_ase_results_summary_pollen_April2025.csv")

reads_pg_pvals2 <- reads_pg_pvals %>% mutate(sigLFC = case_when(
    res.allele.padj < 0.1 &
      abs(res.allele.log2FoldChange) > 1 ~ "True",
    TRUE ~ "False")) %>%
  mutate(log2male_mean.mat = log2(pollen_mean.mat)) %>%
  mutate(log2male_mean.pat = log2(pollen_mean.pat)) %>%
  #filter(male_mean.mat > 10 & male_mean.pat >10) #%>% 
  filter(log2male_mean.pat != "-Inf" & log2male_mean.mat != "-Inf") 

  
  
reads_pg_pvals2$sigLFC<-factor(reads_pg_pvals2$sigLFC, levels = c("True","False"))
### retransforming the data bc I hate myself!
colnames(reads_pg_pvals)

df<-reads_pg_pvals  %>% mutate(sigLFC = case_when(
  res.allele.padj < 0.1 &
    abs(res.allele.log2FoldChange) > 1 ~ "True",
  TRUE ~ "False"))

# readcounts_long<- df %>% pivot_longer(cols = 14:311, names_to = c("sample",".value"),
#                     names_sep = "\\.") #%>% 
#   mutate(sex = case_when(grepl("M",sample) ~ "Male", 
#                        grepl("F", sample) ~ "Female", TRUE ~ NA_character_)) 


# sum_filt_readcounts_long <- readcounts_long %>% 
#   filter(sex == "Male") %>%
#   #filter(pat >=0 & mat >=0) %>%
#   group_by(pairID,sigTest,sigLFC) %>% 
#   summarise(n = n(), totalReadCount = sum(pat+mat), male_mean.pat = mean(pat),male_mean.mat = mean(mat), 
#             male_SE.pat = plotrix::std.error(pat),male_SE.mat = plotrix::std.error(mat))# %>%
#   filter(male_mean.pat !="-Inf")
# readcounts_long$sigLFC<-factor(readcounts_long$sigLFC, levels = c("True","False"))
  

colors <- vibrant(7)[c(4,7)]

figure1a <- ggplot(reads_pg_pvals2,aes(x=pollen_mean.mat,y=pollen_mean.pat)) + 
  #geom_point(aes(size = res.allele.lfcSE),alpha = 0.60) +
  geom_point(alpha = 0.5, size =1,aes(color = sigLFC)) +
  #geom_errorbar(aes(color = sigLFC,ymin = ((male_mean.pat)-(male_SE.pat)), ymax = ((male_mean.pat)+(male_SE.pat)))) +
  #geom_errorbarh(aes(color = sigLFC,xmin = male_mean.mat-male_SE.mat, xmax = male_mean.mat+male_SE.mat)) +
  geom_abline(slope = 1, intercept = 0, linetype=1) +
  theme_classic()+
  geom_smooth(method ="lm", color =  vibrant(7)[7]) + # all points regression line
  geom_smooth(method = "lm", data = function(x) { filter(x, sigLFC == "True") }, se = T, fullrange = TRUE, color =  vibrant(7)[4]) +
  # regression sig only
  xlab("mean expression of X\n gametolog in males") +
  ylab("mean expression of Y\n gametolog in males") +
  scale_x_continuous(transform = "log2",breaks = c(0.0625,2,64,2048,32768)) +
  scale_y_continuous(transform = "log2",breaks = c(0.0625,2,64,2048,32768)) +
  #labs(color= "p-val < 0.1 & |log2FC| >1 ") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.position = "none") +
  scale_color_manual(values= colors) # from the vibrant palette, see vignette("tol")
## not the long version
## no filter
figure1a
ggsave("figures/pollen_xy_expression.png", h =6, w =10)
##

colnames(reads_pg_pvals2) # rm NaN

lm_all<-lm(formula = log2male_mean.pat ~ log2male_mean.mat, data = reads_pg_pvals2) 
#plot(ks_lm)
lm_sig<-lm(formula = log2male_mean.pat ~ log2male_mean.mat, data = reads_pg_pvals2 %>% filter(sigLFC == "True"))
summary(lm_all)
summary(lm_sig)


