library(tidyverse)
library(mcp)
library(rjags)
########################### Ks step change ##################################
hyphy_anno_reads <- read_csv("data/hyphy_anno_reads.csv.gz")

hyphy_anno_reads %>% filter(seqname.pat =="Y" & XYdS <0.5) %>% 
  ggplot() + geom_point(aes(x = start.pat.mb, y = XYdS)) +
  geom_step(aes(x = start.pat.mb, y = XYdS)) +
#  geom_smooth(aes(x = start.pat.mb, y = XYdS), method = "loess", se = FALSE) +
  xlab("X-bearing Haplotype\nposition(Mb)") +
  ylab("Ks") +
  theme_classic()

## define the model ## exampl3 ##
model = list(
  response ~ 1,  # plateau (int_1)
  ~ 0 + time,    # joined slope (time_2) at cp_1
  ~ 1 + time     # disjoined slope (int_3, time_3) at cp_2
)
#get data and fit it
ex = mcp_example("demo")
fit = mcp(model, ex$data)
## run this on a server my god!
plot(fit)

### ok my real model
model1 <- model2 = list(
  XYdS ~ 1,  
  ~ 0 + start.pat.mb,   
  ~ 1 + start.pat.mb
)
model2 <- list(
  XYdS ~ 1,  
  ~ 0 + start.pat.mb,   
  ~ 1 + start.pat.mb, 
  ~ 0 + start.pat.mb,
  ~ 1 + start.pat.mb 
)
df<-hyphy_anno_reads %>% filter(seqname.pat =="Y" & XYdS <0.5) %>% select(start.pat.mb, XYdS) %>% drop_na()
fit1 <- mcp(model1, df)
fit2 <- mcp (model2, df)
s1<-summary(fit1)
s2<-summary(fit2)
p1<-plot(fit1)
p2<-plot(fit2)
p2
cowplot::plot_grid(p1,p2)

## George suggests 1) unsupervised clustering positions by Ks, 
# 2) compare clusterings vs position - strata fall out? e.g. does the ks cluster match the positions? idk
#confusing
#k means clustering?

library(factoextra)
library(cluster)

km.3<-kmeans(df$XYdS, centers = 3, nstart = 20)
km.3
fviz_cluster(km.3, data = df,ellipse.type = "norm")
km.2<-kmeans(df$XYdS, 2, 20)
km.2
fviz_cluster(km.2, data = df,ellipse.type = "norm")
km.4<-kmeans(df$XYdS, 4, 20)
km.4
fviz_cluster(km.4, data = df, ellipse.type = "norm")

## example from george
library(MASS)
kmod<- kmeans(df$XYdS, centers =2)
df$cluster<-kmod$cluster
mylogit<-glm(as.factor(cluster)~start.pat.mb, data =df, family = "binomial")
summary(mylogit)
a<-ggplot(df) + geom_point(aes(x=start.pat.mb, y = cluster,color = XYdS)) 

kmod<- kmeans(df$XYdS, centers =3)
df$cluster<-kmod$cluster
mylogit<-glm(as.factor(cluster)~start.pat.mb, data =df, family = "binomial")
summary(mylogit)
b<-ggplot(df) + geom_point(aes(x=start.pat.mb, y = cluster,color = XYdS)) 

kmod<- kmeans(df$XYdS, centers =4)
df$cluster<-kmod$cluster
mylogit<-glm(as.factor(cluster)~start.pat.mb, data =df, family = "binomial")
summary(mylogit)
c<-ggplot(df) + geom_point(aes(x=start.pat.mb, y = cluster,color = XYdS)) 
cowplot::plot_grid(a,b,c)
ggsave("figures/ks_kmeans_genes.png", width = 8, height = 6, units = "in")
### do it on windowed data
## ks_win100_filt<-read_csv("ks_win100_windowNumsY.csv.gz")
kmod2<-kmeans(ks_win100_filt$Ks_win, centers =2)
ks_win100_filt$cluster2<-kmod2$cluster
mylogit2<-glm(as.factor(cluster2)~start.pat.mb, data =ks_win100_filt, family = "binomial")
summary(mylogit2)
a<-ggplot(ks_win100_filt) + geom_point(aes(x=start.pat.mb, y = cluster2,color = Ks_win)) 
## awesome!!!
kmod3<-kmeans(ks_win100_filt$Ks_win, centers =3)
ks_win100_filt$cluster3<-kmod3$cluster
mylogit3<-glm(as.factor(cluster3)~start.pat.mb, data =ks_win100_filt, family = "binomial")
summary(mylogit3)
b<-ggplot(ks_win100_filt) + geom_point(aes(x=start.pat.mb, y = cluster3,color = Ks_win)) 
## next = 4
kmod4<-kmeans(ks_win100_filt$Ks_win, centers =4)
ks_win100_filt$cluster4<-kmod4$cluster
mylogit4<-glm(as.factor(cluster4)~start.pat.mb, data =ks_win100_filt, family = "binomial")
summary(mylogit4)
c<-ggplot(ks_win100_filt) + geom_point(aes(x=start.pat.mb, y = cluster4, color = Ks_win)) 


cowplot::plot_grid(a,b,c)
ggsave("figures/ks_win100_kmeans_clusters.png", width = 8, height = 6, units = "in")
## sliding windows, makes sense that it would be like this.
## what about non-sliding windows
ks_win100_noslide<-hyphy_anno_reads %>% filter(seqname.pat =="Y" & XYdS <0.5) %>% #1250 obs - so divide by 12
  select(start.pat.mb, XYdS) %>% group_by(group = gl(length(XYdS)/10, 10)) %>%
    summarise(med_val = median(XYdS))

kmod2<-kmeans(ks_win100_noslide$med_val, centers =2)
ks_win100_noslide$cluster2<-kmod2$cluster
mylogit2<-glm(as.factor(cluster2)~group, data =ks_win100_noslide, family = "binomial")
summary(mylogit2)
a<-ggplot(ks_win100_noslide) + geom_point(aes(x=group, y = cluster2,color = med_val)) 


kmod3<-kmeans(ks_win100_noslide$med_val, centers =3)
ks_win100_noslide$cluster3<-kmod3$cluster
mylogit3<-glm(as.factor(cluster3)~group, data =ks_win100_noslide, family = "binomial")
summary(mylogit3)
b<-ggplot(ks_win100_noslide) + geom_point(aes(x=group, y = cluster3,color = med_val))

kmod4<-kmeans(ks_win100_noslide$med_val, centers =4)
ks_win100_noslide$cluster4<-kmod4$cluster
mylogit4<-glm(as.factor(cluster4)~group, data =ks_win100_noslide, family = "binomial")
summary(mylogit4)
c<-ggplot(ks_win100_noslide) + geom_point(aes(x=group, y = cluster4,color = med_val))
cowplot::plot_grid(a,b,c)
ggsave("figures/ks_winNoSlide_kmeans.png", width = 8, height = 6, units = "in")

#### time to give up!!!

# okey :)
#for now
