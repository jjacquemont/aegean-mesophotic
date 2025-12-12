library(tidyverse)
library(ggplot2)
#remotes::install_github("ropensci/rfishbase")
library(dplyr)
library(FD)
library(mgcv) #gam model
library(modelr) #data manipulation such as data_grid
library(ggpubr) #combining ggplot figures
library(viridis)
library(vegan) # ecological and environmental data analysis
library(clustsig) 
library(ggdendro)
library(labdsv) # indval - indicator species
#install.packages("mvpart", repos = NULL, type = "source")
  ## had to manually replace Sint by int
library(mvpart) # multivariate regression tree
library(ggrepel)
library(car) # to perform homogeneity of variance tests 

fish.data <- read.csv("fish-transects.csv")

#### 1. Site analysis ####
#### 1.1 Trikeri ####
Trikeri.fish <- read.csv("fish-transects.csv") %>%
  filter(area=="Trikeri" & species !="fishing line")%>%
  mutate(depth=case_when(depth==41~40,
                         depth==19~20,
                         TRUE~depth)) 


#### species richness ####
Trikeri.fish %>%
  group_by(species)%>%
  summarize(abundance=sum(number))

no.fish <- data.frame(depth=c(50,60),forest=c("out","out"),sp.mean=c(0,0),
                      sp.se=c(0,0),replicates=c(3,3))

Trikeri.sp <-  Trikeri.fish %>%
  group_by(depth,forest,transect.nb,species) %>%
  summarize(abundance=sum(number)) %>%
  group_by(depth,forest,transect.nb) %>%
  summarize(sp.nb=n())%>%
  group_by(depth,forest) %>%
  summarize(sp.mean=mean(sp.nb),sp.se=sd(sp.nb)*1.96/n(),replicates=n())%>%
  bind_rows(no.fish) %>%
  mutate(depth=factor(depth,levels=c("10","20","30","40","50","55","60")),
         forest=factor(forest,levels=c("out","in")))%>%
  arrange(depth,forest)


ggplot(Trikeri.sp,aes(x=depth,y=sp.mean,fill=forest))+
  geom_bar(width=0.6,stat="identity", position=position_dodge())+
  #position=position_dodge2(preserve = "single"))+
  theme_classic()+
  geom_errorbar(width=0.1,mapping=aes(x=c(1,2,3,4,4.85,5.15,5.85,6.15,6.85,7.15),
                                      ymin = sp.mean-sp.se, 
                                      ymax = sp.mean+sp.se), 
                position=position_dodge())+
  #position=position_dodge2(width=2,preserve = "single"))+
  scale_fill_manual(values=c("#FEB72DFF","#D14E72FF"))+
  ylab("species richness (ind. transect-1)")


#### abundance ####
no.fish <- data.frame(depth=c(50,60),forest=c("out","out"),abundance.mean=c(5,5),
                      abundance.se=c(0,0),replicates=c(3,3))

Trikeri.abu <- Trikeri.fish %>%
  group_by(depth,forest,transect.nb) %>%
  summarize(abundance=sum(number)) %>%
  group_by(depth,forest) %>%
  summarize(abundance.mean=mean(abundance),abundance.se=sd(abundance)*1.96/n(),replicates=n())%>%
  bind_rows(no.fish) %>%
  mutate(depth=factor(depth,levels=c("10","20","30","40","50","55","60")),
         forest=factor(forest,levels=c("out","in")))%>%
  arrange(depth,forest)


ggplot(Trikeri.abu,aes(x=depth,y=abundance.mean,fill=forest))+
  geom_bar(width=0.6,stat="identity", position=position_dodge())+
  #position=position_dodge2(preserve = "single"))+
  theme_classic()+
  geom_errorbar(width=0.1,mapping=aes(x=c(1,2,3,4,4.85,5.15,5.85,6.15,6.85,7.15),ymin = abundance.mean-abundance.se, 
                                      ymax = abundance.mean+abundance.se), 
                position=position_dodge())+
  #position=position_dodge2(width=2,preserve = "single"))+
  scale_fill_manual(values=c("#FEB72DFF","#D14E72FF"))+
  ylab("abundance (ind. transect-1)")

#### total length ####
no.fish <- data.frame(depth=c(50,60),forest=c("out","out"),length.mean=c(0,0),
                      length.se=c(0,0),replicates=c(3,3))%>%
  mutate(depth=factor(depth,levels=c("10","20","30","40","50","55","60")),
         forest=factor(forest,levels=c("out","in")))

Trikeri.length <- Trikeri.fish %>%
  group_by(depth,forest,transect.nb) %>%
  summarize(length=sum(size)) %>%
  mutate(depth=factor(depth,levels=c("10","20","30","40","50","55","60")),
         forest=factor(forest,levels=c("out","in")))

Trikeri.length.mean <- Trikeri.length %>%
  group_by(depth,forest) %>%
  summarize(length.mean=mean(length),length.se=sd(length)*1.96/n(),replicates=n()) %>%
  bind_rows(no.fish)%>%
  arrange(depth,forest)


ggplot()+
  geom_bar(Trikeri.length.mean,mapping=aes(x=depth,y=length.mean,fill=forest),
           width=0.6,stat="identity", position=position_dodge())+
  geom_point(Trikeri.length,mapping=aes(x=depth,y=length,fill=forest), shape=21)+
  theme_classic()+
  geom_errorbar(Trikeri.length.mean,width=0.1,mapping=aes(x=c(1,2,3,4,4.85,5.15,5.85,6.15,6.85,7.15),ymin = length.mean-length.se, 
                                                         ymax = length.mean+length.se), 
                position=position_dodge())+
  #position=position_dodge2(width=2,preserve = "single"))+
  scale_fill_manual(values=c("#FEB72DFF","#D14E72FF"))+
  ylab("length (cm. transect-1)")

#### 1.2. Fourni ####
Fourni.fish <- read.csv("fish-transects.csv") %>%
  filter(area=="Fourni" & species !="fishing line")

#### species richness ####
Fourni.fish %>%
  group_by(species)%>%
  summarize(abundance=sum(number))

no.fish <- data.frame(depth=c(85,85),forest=c("out","out"),
                      transect_nb=c("85_2","85_3"),sp.nb=c(0,0))

Fourni.sp <-  Fourni.fish %>%
  group_by(depth,forest,transect.nb,species) %>%
  summarize(abundance=sum(number)) %>%
  group_by(depth,forest,transect.nb) %>%
  summarize(sp.nb=n())%>%
  bind_rows(no.fish) %>%
  group_by(depth,forest) %>%
  summarize(sp.mean=mean(sp.nb),sp.se=sd(sp.nb)*1.96/n(),replicates=n())%>%
  mutate(depth=factor(depth,levels=c("10","20","30","40","80","85","90")),
         forest=factor(forest,levels=c("out","in")))%>%
  arrange(depth,forest)


ggplot(Fourni.sp,aes(x=depth,y=sp.mean,fill=forest))+
  geom_bar(width=0.6,stat="identity", position=position_dodge())+
  #position=position_dodge2(preserve = "single"))+
  theme_classic()+
  geom_errorbar(width=0.1,mapping=aes(x=c(1,2,3,4,4.85,5.15,5.85,6.15,6.85,7.15),
                                      ymin = sp.mean-sp.se, 
                                      ymax = sp.mean+sp.se), 
                position=position_dodge())+
  #position=position_dodge2(width=2,preserve = "single"))+
  scale_fill_manual(values=c("#FEB72DFF","#D14E72FF"))+
  ylab("species richness (ind. transect-1)")

#### abundance ####
no.fish <- data.frame(depth=c(85,85),forest=c("out","out"),
                      transect.nb=c("85_2","85_3"),abundance=c(0,0))

Fourni.abu <- Fourni.fish %>%
  group_by(depth,forest,transect.nb) %>%
  summarize(abundance=sum(number)) %>%
  bind_rows(no.fish)%>%
  mutate(depth=factor(depth,levels=c("10","20","30","40","80","85","90")),
         forest=factor(forest,levels=c("out","in")))%>%
  arrange(depth,forest)
  
Fourni.abu.mean <- Fourni.abu %>%
  group_by(depth,forest) %>%
  summarize(abundance.mean=mean(abundance),abundance.se=sd(abundance)*1.96/n(),replicates=n())


ggplot()+
  geom_bar(Fourni.abu.mean,mapping=aes(x=depth,y=abundance.mean,fill=forest),
           width=0.6,stat="identity", position=position_dodge())+
  geom_point(Fourni.abu,mapping=aes(x=depth,y=abundance,fill=forest), shape=21)+
  theme_classic()+
  geom_errorbar(Fourni.abu.mean,width=0.1,mapping=aes(x=c(1,2,3,4,4.85,5.15,5.85,6.15,6.85,7.15),ymin = abundance.mean-abundance.se, 
                                      ymax = abundance.mean+abundance.se), 
                position=position_dodge())+
                #position=position_dodge2(width=2,preserve = "single"))+
  scale_fill_manual(values=c("#FEB72DFF","#D14E72FF"))+
  ylab("abundance (ind. transect-1)")

#### total length ####
no.fish <- data.frame(depth=c(85,85),forest=c("out","out"),
                      transect.nb=c("85_2","85_3"),length=c(0,0))

Fourni.length <- Fourni.fish %>%
  group_by(depth,forest,transect.nb) %>%
  summarize(length=sum(size)) %>%
  bind_rows(no.fish)%>%
  mutate(depth=factor(depth,levels=c("10","20","30","40","80","85","90")),
         forest=factor(forest,levels=c("out","in")))%>%
  arrange(depth,forest)

Fourni.length.mean <- Fourni.length %>%
  group_by(depth,forest) %>%
  summarize(length.mean=mean(length),length.se=sd(length)*1.96/n(),replicates=n())


ggplot()+
  geom_bar(Fourni.length.mean,mapping=aes(x=depth,y=length.mean,fill=forest),
           width=0.6,stat="identity", position=position_dodge())+
  geom_point(Fourni.length,mapping=aes(x=depth,y=length,fill=forest), shape=21)+
  theme_classic()+
  geom_errorbar(Fourni.length.mean,width=0.1,mapping=aes(x=c(1,2,3,4,4.85,5.15,5.85,6.15,6.85,7.15),ymin = length.mean-length.se, 
                                                      ymax = length.mean+length.se), 
                position=position_dodge())+
  #position=position_dodge2(width=2,preserve = "single"))+
  scale_fill_manual(values=c("#FEB72DFF","#D14E72FF"))+
  ylab("length (cm. transect-1)")


#### 1.3. Fourni + Trikeri ####
#### species richness ####
no.fish <- data.frame(area=c(rep("Trikeri",6),rep("Fourni",2)),habitat="mesophotic barren",
                      transect_nb=1:8,sp.nb=0)

greece.sp1 <- read.csv("fish-transects.csv") %>%
  filter(species !="fishing line")%>%
  mutate(depth=case_when(depth==41~40,
                         depth==19~20,
                         TRUE~depth)) %>%
  mutate(habitat=case_when(depth<=40 ~"shallow",
                           depth>40 & forest=="in" ~ "MAF",
                           depth>40 & forest=="out"~ "mesophotic barren"))%>%
  group_by(area,habitat,transect.nb,species) %>%
  summarize(abundance=sum(number)) %>%
  group_by(area,habitat,transect.nb) %>%
  summarize(sp.nb=n())%>%
  bind_rows(no.fish)%>%
  mutate(habitat=factor(habitat,levels=c("shallow","MAF","mesophotic barren")))

# means per habitat and site
greece.sp2 <- greece.sp1 %>%
  group_by(habitat,area) %>%
  summarize(sp.mean=mean(sp.nb),sp.se=sd(sp.nb)/sqrt(n()),replicates=n())

# means per habitat 
greece.sp3 <- greece.sp1 %>%
  group_by(habitat) %>%
  summarize(sp.mean=mean(sp.nb),sp.se=sd(sp.nb)/sqrt(n()),replicates=n())

## grouping by site
sp.richness <- ggplot(greece.sp1,aes(x=area,y=sp.nb,fill=habitat))+
  #geom_violin(position="dodge",alpha=0.8)+
  geom_boxplot(alpha=0.8, width = 0.5,position=position_dodge((width = 0.7)))+
  #geom_violin( position="dodge",outlier.shape = NA)+
  #geom_violin(position = position_dodge(width = 0.8), alpha = 0.5, trim = FALSE) +
  #geom_bar(width=0.6,stat="identity", position=position_dodge((width = 0.8)))+
  geom_point(greece.sp1,mapping=aes(x=area,y=sp.nb,fill=habitat), shape=21,
             position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7))+
  geom_text(greece.sp2, mapping=aes(x=area,y=-1,label=paste("n=",replicates)),
            position=position_dodge((width = 0.7)),size=3)+
  scale_x_discrete(expand = c(0, 0))+
  theme(axis.text.x = element_text(color = "black"))+
  #position=position_dodge2(preserve = "single"))+
  theme_classic()+
  #geom_errorbar(width=0.1,mapping=aes(ymin = sp.mean-sp.se, ymax = sp.mean+sp.se), 
   #             position=position_dodge(width=0.8))+
  #position=position_dodge2(width=2,preserve = "single"))+
  scale_fill_manual(values=c("#1FBFC0","#B589D6","#D6B588"))+
  labs(y=expression(paste("species. 120 m"^-2)),x="",title="species richness")
sp.richness

## grouping by habitat
 ggplot(greece.sp1,aes(x=habitat,y=sp.nb,fill=habitat,color=area))+
  geom_boxplot(alpha=0.8, width = 0.5,position=position_dodge((width = 0.7)))+
  geom_point(greece.sp1,mapping=aes(x=habitat,y=sp.nb,fill=habitat), shape=21,
             position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7))+
  geom_text(greece.sp2, mapping=aes(x=habitat,y=-1,label=paste("n=",replicates)),
            position=position_dodge((width = 0.7)),size=3)+
  scale_x_discrete(expand = c(0, 0))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"))+
   #     axis.text.x = element_text(color = "black"))+
  #position=position_dodge2(preserve = "single"))+
  scale_fill_manual(values=c("#1FBFC0","#B589D6","#D6B588"))+
  scale_color_manual(values=c("black","black"))+
  labs(y=expression(paste("species. 120 m"^-2)),x="",title="species richness")

summary(aov(data=greece.sp1, sp.nb~habitat)) # signif param
summary(aov(data=greece.sp1, sp.nb~as.factor(area))) # not signif
aov_sp <- aov(data=greece.sp1, sp.nb~as.factor(area)*habitat) # signif, only habitat
summary(aov_sp)

shapiro.test(residuals(aov_sp)) # p =0.27, normality of residuals OK
leveneTest(sp.nb~area*habitat,greece.sp1) #p=0.07, homogeneity of variance OK 



kruskal.test(data=greece.sp1, sp.nb~habitat) #signif non param
TukeyHSD(aov(data=greece.sp1, sp.nb~habitat)) # pairwise comparisons

kruskal.test(data=greece.sp1 %>% filter(area=="Trikeri"), sp.nb~habitat) # signif Trikeri
TukeyHSD(aov(data=greece.sp1 %>% filter(area=="Trikeri"), sp.nb~habitat)) # signif Trikeri

kruskal.test(data=greece.sp1 %>% filter(area=="Fourni"), sp.nb~habitat) # signif Fourni
TukeyHSD(aov(data=greece.sp1 %>% filter(area=="Fourni"), sp.nb~habitat)) # signif Trikeri

# non parametric tests
pairwise.wilcox.test(greece.sp1$sp.nb, greece.sp1$habitat,
                     p.adjust.method = "BH")

trikeri.sp <- greece.sp1 %>% filter(area=="Trikeri")
pairwise.wilcox.test(trikeri.sp$sp.nb, trikeri.sp$habitat,
                     p.adjust.method = "BH")

fourni.sp <- greece.sp1 %>% filter(area=="Fourni")
pairwise.wilcox.test(fourni.sp$sp.nb, fourni.sp$habitat,
                     p.adjust.method = "BH")

#### abundance ####
no.fish <- data.frame(area=c(rep("Trikeri",6),rep("Fourni",2)),habitat="deep bare",
                      transect.nb=as.character(1:8),abundance=0)

greece.abu1 <- read.csv("fish-transects.csv") %>%
  filter(species !="fishing line")%>%
  mutate(depth=case_when(depth==41~40,
                         depth==19~20,
                         TRUE~depth)) %>%
  mutate(habitat=case_when(depth<=40 ~"shallow",
                           depth>40 & forest=="in" ~ "deep forest",
                           depth>40 & forest=="out"~ "deep bare"))%>%
  group_by(area,habitat,transect.nb) %>%
  summarize(abundance=sum(number))%>% 
  bind_rows(no.fish)%>%
  mutate(habitat=factor(habitat,levels=c("shallow","deep forest","deep bare")))

# per habitat and location
greece.abu2 <- greece.abu1 %>%
  mutate(abu.log=log(abundance+1))%>%
  group_by(habitat,area) %>%
  summarize(abu.mean=mean(abundance),abu.se=sd(abundance)/n(),replicates=n(),
            abu.log.mean=mean(abu.log),abu.log.se=sd(abu.log)/n())

# per habitat 
greece.abu3 <- greece.abu1 %>%
  mutate(abu.log=log(abundance+1))%>%
  group_by(habitat) %>%
  summarize(abu.mean=mean(abundance),abu.se=sd(abundance)/n(),replicates=n(),
            abu.log.mean=mean(abu.log),abu.log.se=sd(abu.log)/n())

# grouped by site
abundance <- ggplot(greece.abu1,mapping=aes(x=area,y=log(abundance+1),fill=habitat))+
  geom_boxplot(alpha=0.8, width = 0.5,position=position_dodge((width = 0.7)))+
  geom_point(shape=21,position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7))+
  geom_text(greece.abu2, mapping=aes(x=area,y=-0.5,label=paste("n=",replicates)),
           position=position_dodge(width = 0.7),size=3)+
  scale_x_discrete(expand = c(0, 0))+
  theme_classic()+
 # scale_y_continuous(limits = c(-100, 2000)) +
  theme(axis.text.x = element_text(color = "black"))+
  scale_fill_manual(values=c("#1FBFC0","#B589D6","#D6B588"))+
  labs(y=expression(paste("individuals (log). 120 m"^-2)),x="",title="abundance")
abundance

# grouped by habitat
ggplot(greece.abu1,mapping=aes(x=habitat,y=log(abundance+1),fill=habitat,color=area))+
  geom_boxplot(alpha=0.8, width = 0.5,position=position_dodge((width = 0.7)),outliers = F)+
  geom_point(shape=21,position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7))+
  geom_text(greece.abu2, mapping=aes(x=habitat,y=-1,label=paste("n=",replicates)),
            position=position_dodge(width = 0.7),size=3)+
  scale_x_discrete(expand = c(0, 0))+
  theme_classic()+
  scale_y_continuous(limits = c(-1, 8)) +
  theme(axis.text = element_text(color = "black"))+
  scale_fill_manual(values=c("#6CBDE9","#D6B588","#FD7C6E"))+
  scale_color_manual(values=c("black","black"))+
  labs(y="log individuals transect-1",x="",title="abundance")

##
summary(aov(data=greece.abu1, log(abundance+1)~habitat)) # signif param
summary(aov(data=greece.abu1, abundance~habitat)) # signif param
summary(aov(data=greece.abu1, abundance~as.factor(area))) # not signif
abu.aov <- aov(data=greece.abu1, log(abundance+1)~as.factor(area)*habitat) # signif, only area and area*habitat
summary(abu.aov)
# pairwise comparisons: signif between deep bare and the others
TukeyHSD(aov(data=greece.abu1, abundance~habitat))
shapiro.test(residuals(abu.aov)) # p =0.27, normality of residuals OK
leveneTest(log(abundance+1)~area*habitat,greece.abu1) #p=0.07, homogeneity of variance OK 

kruskal.test(data=greece.abu1, log(abundance+1)~habitat) #signif non param
kruskal.test(data=greece.abu1, abundance~habitat) #signif non param
TukeyHSD(aov(data=greece.abu1, log(abundance+1)~habitat))

kruskal.test(data=greece.abu1 %>% filter(area=="Trikeri"), log(abundance+1)~habitat) # signif Trikeri
TukeyHSD(aov(data=greece.abu1 %>% filter(area=="Trikeri"), log(abundance+1)~habitat)) # signif Trikeri


kruskal.test(data=greece.abu1 %>% filter(area=="Fourni"), log(abundance+1)~habitat) # non-signif Fourni
TukeyHSD(aov(data=greece.abu1 %>% filter(area=="Fourni"), log(abundance+1)~habitat)) # non-signif Fourni

# non parametric tests
pairwise.wilcox.test(log(greece.abu1$abundance+1), greece.abu1$habitat,
                     p.adjust.method = "BH")

trikeri.abu <- greece.abu1 %>% filter(area=="Trikeri")
pairwise.wilcox.test(log(trikeri.abu$abundance+1), trikeri.abu$habitat,
                     p.adjust.method = "BH")

fourni.abu <- greece.abu1 %>% filter(area=="Fourni")
pairwise.wilcox.test(log(fourni.abu$abundance+1), fourni.abu$habitat,
                     p.adjust.method = "BH")

## comparing data one site at a time
TukeyHSD(aov(data=(greece.abu1 %>% filter (area=="Fourni")) , abundance~habitat)) # no signif diff
TukeyHSD(aov(data=(greece.abu1 %>% filter (area=="Trikeri")) , abundance~habitat)) # diff between shallow and deep
TukeyHSD(aov(data=(greece.abu1 %>% filter (area=="Trikeri")) , log(abundance+1)~habitat)) # 

#### total length ####
no.fish <- data.frame(area=c(rep("Trikeri",6),rep("Fourni",2)),habitat="deep bare",
                      transect.nb=as.character(1:8),total.length=0)

greece.size1 <- read.csv("fish-transects.csv") %>%
  filter(species !="fishing line")%>%
  mutate(depth=case_when(depth==41~40,
                         depth==19~20,
                         TRUE~depth)) %>%
  mutate(habitat=case_when(depth<=40 ~"shallow",
                           depth>40 & forest=="in" ~ "deep forest",
                           depth>40 & forest=="out"~ "deep bare"))%>%
  group_by(area,habitat,transect.nb) %>%
  summarize(total.length=sum(size))%>% 
  bind_rows(no.fish)%>%
  mutate(habitat=factor(habitat,levels=c("shallow","deep forest","deep bare")))

greece.size2 <- greece.size1 %>%
  group_by(habitat,area) %>%
  summarize(size.mean=mean(total.length),size.se=sd(total.length)*1.96/n(),replicates=n())

# by site
length <-ggplot(greece.size1,aes(x=area,y=total.length,fill=habitat))+
  geom_boxplot(alpha=0.8, width = 0.5,position=position_dodge((width = 0.7)),
               outlier.shape = NA)+
  geom_point(shape=21,position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7))+
  geom_text(greece.size2, mapping=aes(x=area,y=-20,label=paste("n=",replicates)),
        position=position_dodge((width = 0.7)),size=3)+
  theme_classic()+
  scale_x_discrete(expand = c(0, 0))+
  theme(axis.text.x = element_text(color = "black"))+
  scale_fill_manual(values=c("#1FBFC0","#B589D6","#D6B588"))+
  labs(y=expression(paste("total length (cm). 120 m"^-2)),x="",title="total length")
length

# by habitat
 ggplot(greece.size1,aes(x=habitat,y=total.length,fill=habitat,color=area))+
  geom_boxplot(alpha=0.8, width = 0.5,position=position_dodge((width = 0.7)),
               outlier.shape = NA)+
  geom_point(shape=21,position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7))+
  geom_text(greece.size2, mapping=aes(x=habitat,y=-20,label=paste("n=",replicates)),
            position=position_dodge((width = 0.7)),size=3)+
  theme_classic()+
  scale_x_discrete(expand = c(0, 0))+
  theme(axis.text.x = element_text(color = "black"))+
  scale_fill_manual(values=c("#6CBDE9","#D6B588","#FD7C6E"))+
  scale_color_manual(values=c("black","black"))+
  labs(y="total length. transect-1",x="",title="total length")


summary(aov(data=greece.size1, total.length~habitat)) # signif param
summary(aov(data=greece.size1, total.length~as.factor(area))) # not signif
aov_size <- aov(data=greece.size1, total.length~as.factor(area)*habitat) # signif, only habitat
summary(aov_size)

kruskal.test(data=greece.size1, log(total.length+1)~habitat) #signif non param

shapiro.test(residuals(abu.aov)) # p =0.02, not normal
leveneTest(total.length~area*habitat,greece.size1) #p=0.01, no homogeneity of variance 


# pairwise comparisons: signif between deep bare and the others
TukeyHSD(aov(data=greece.size1, log(total.length+1)~habitat)) 
TukeyHSD(aov(data=greece.size1 %>% filter(area=="Trikeri"), log(total.length+1)~habitat)) 
TukeyHSD(aov(data=greece.size1 %>% filter(area=="Fourni"), log(total.length+1)~habitat)) 

# non parametric tests
pairwise.wilcox.test(log(greece.size1$total.length+1), greece.size1$habitat,
                     p.adjust.method = "BH")

trikeri.size <- greece.size1 %>% filter(area=="Trikeri")
pairwise.wilcox.test(log(trikeri.size$total.length+1), trikeri.size$habitat,
                     p.adjust.method = "BH")

fourni.size <- greece.size1 %>% filter(area=="Fourni")
pairwise.wilcox.test(log(fourni.size$total.length+1), fourni.size$habitat,
                     p.adjust.method = "BH")


#### three plots together - FIGURE 3 ####
ggarrange(abundance, sp.richness, length, nrow=3,ncol=1,legend="right",
          common.legend=T)


#### 1.4. Species list ####
species.list <- fish.data %>%
  distinct(species)

#### 2. Trans-site analysis ####
#### 2.1. Depth x site ####
fish.com <- read.csv("fish-transects.csv") %>%
  filter(species !="fishing line") %>%
  mutate(depth=case_when(depth==41~40,
                         depth==19~20,
                         TRUE~depth)) %>%
  select(species, depth, area, number) %>%
  group_by(species, depth, area) %>%
  summarize(abundance = sqrt(sum(number))) %>%
  spread(species, abundance, fill =  0) %>%
  as.data.frame()

rownames(fish.com) <- paste(fish.com$area,fish.com$depth)
fish.com.mds <- fish.com[-c(1:2)]

## SIMPROF
greece.clust <- simprof(data = fish.com.mds, method.cluster = "complete", method.distance="braycurtis",
                     method.transform="identity", alpha=0.0000001,
                     sample.orientation="row", const=0,
                     silent=TRUE, increment=100,
                     undef.zero=TRUE, warn.braycurtis=TRUE)

#### dendrogram ####
# 2 groups from shallow to deep
greece.res.tib <- tibble(cl1 = list(greece.clust$significantclusters[[1]]),
                      cl2 = list(greece.clust$significantclusters[[2]]),
                      cl3 = list(greece.clust$significantclusters[[3]])) %>%
  unnest_longer(cl1) %>%
  unnest_longer(cl2) %>%
  unnest_longer(cl3) %>%
  gather(key = "cluster", value = "label") %>% # transform wide to long format 
  distinct()  # table with cluster and associated dbands
  #mutate(label = as.factor(site)) %>%
  #arrange(dband) %>%
  #mutate(run_num = 1:27) %>%
  #group_by(cluster) %>%
  #mutate(ord_num = min(run_num))

# vegdist() : distance between each depth bin given relative abundance of species
# hlcust(): hierarchical clustering. 
#takes a dissimilarity table as produced by vegdist()
#returns the order of observations by similarity ($order), and the order of merges ($merge)
greece.clust.veg <- hclust(vegdist(fish.com.mds, method = "bray"), method = "complete")

# changes class of object from hclust to dendogram 
greece.dend <- as.dendrogram(greece.clust.veg) 

# length and coordinates of each branch in the dendogram ($segments)
greece.dend.data <- dendro_data(greece.dend, type = "rectangle")

# joining with values of dbands
greece.dend.helper <- greece.dend.data$labels %>%
  inner_join(greece.res.tib) 
greece.dend.data$labels <- greece.dend.helper

# get colors for each level using fishualize and tinter package
carib.pal <- fish(4, option = "Bodianus_rufus")
cur.pal <- tinter(carib.pal[1], steps = 5)
blue_palette <- c("#BAE0F3","#87CEEB","#6CBDE9","#50ABE7", "#4895EF", "#4361EE", "#2835AF", "#12086F")
mesoP_palette <- c("#DCEEF3","#C2E2EA","#A7D5E1","#8DC8D8","#72BBCE")
rariP_palette <- c("#DDC5EB","#BC90DB","#831EB6")
pal <- c(blue_palette[3],blue_palette[3],rariP_palette[1],rariP_palette[2])

ggplot(greece.dend.data$segments) + 
  #geom_rect(xmin=22.7,xmax=27.3,ymin=-0.13, ymax=0.95,fill=blue_palette[3])+
  #geom_rect(xmin=18.7,xmax=22.3,ymin=-0.13, ymax=0.95,fill=blue_palette[1])+
  #geom_rect(xmin=12.7,xmax=18.3,ymin=-0.13, ymax=0.95,fill=rariP_palette[1])+
  #geom_rect(xmin=0.7,xmax=12.3,ymin=-0.13, ymax=0.95,fill=rariP_palette[2])+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), color = "black")+
  geom_text(data = greece.dend.data$labels, aes(x, y, label = label), color = "black", #below for yellow gradient
            #color = as.factor(ord_num)), 
            hjust = 1, angle = 90, size = 3, fontface = "bold") +
  #scale_color_manual(values = cur.pal[-1]) +
  #theme of plot
  theme(legend.position = "none",
        line = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  # clusters highlight
  #geom_segment(aes(x=0.8,xend=3.2,y=-0.13,yend=-0.13),size=1)+
  #geom_segment(aes(x=3.8,xend=8.2,y=-0.13,yend=-0.13),size=1)+
  #geom_segment(aes(x=8.8,xend=12.2,y=-0.13,yend=-0.13),size=1)+
  #geom_segment(aes(x=12.8,xend=18.2,y=-0.13,yend=-0.13),size=1)+
  #geom_segment(aes(x=18.8,xend=22.2,y=-0.13,yend=-0.13),size=1)+
  #geom_segment(aes(x=22.8,xend=27.2,y=-0.13,yend=-0.13),size=1)+
  ylab("")+
  ylim(-0.2,1)

  

#### 2.2. Depth x site x MAF - FIG 2 LEFT ####
fish.com <- read.csv("fish-transects.csv") %>%
  filter(species !="fishing line") %>%
  mutate(depth=case_when(depth==41~40,
                         depth==19~20,
                         TRUE~depth)) %>%
  group_by(species, depth, area, forest) %>%
  summarize(abundance = sqrt(sum(number))) %>%
  spread(species, abundance, fill =  0) %>%
  as.data.frame()

rownames(fish.com) <- paste(fish.com$area,fish.com$depth,fish.com$forest)
fish.com.mds <- fish.com[-c(1:3)]

greece.clust <- simprof(data = fish.com.mds, method.cluster = "complete", method.distance="braycurtis",
                        method.transform="identity", alpha=0.0000001,
                        sample.orientation="row", const=0,
                        silent=TRUE, increment=100,
                        undef.zero=TRUE, warn.braycurtis=TRUE)

# 3 groups from shallow to deep
greece.res.tib <- tibble(cl1 = list(greece.clust$significantclusters[[1]]),
                         cl2 = list(greece.clust$significantclusters[[2]]),
                         cl3 = list(greece.clust$significantclusters[[3]])) %>%
  unnest_longer(cl1) %>%
  unnest_longer(cl2) %>%
  unnest_longer(cl3) %>%
  gather(key = "cluster", value = "label") %>% # transform wide to long format 
  distinct() 


  # vegdist() : distance between each depth bin given relative abundance of species
  # hlcust(): hierarchical clustering. 
  #takes a dissimilarity table as produced by vegdist()
  #returns the order of observations by similarity ($order), and the order of merges ($merge)
greece.clust.veg <- hclust(vegdist(fish.com.mds, method = "bray"), method = "complete")

# changes class of object from hclust to dendogram 
greece.dend <- as.dendrogram(greece.clust.veg) 

# length and coordinates of each branch in the dendogram ($segments)
greece.dend.data <- dendro_data(greece.dend, type = "rectangle")

# joining with values of dbands
greece.dend.helper <- greece.dend.data$labels %>%
  inner_join(greece.res.tib) 
greece.dend.data$labels <- greece.dend.helper

ggplot(greece.dend.data$segments) + 
  #geom_rect(xmin=22.7,xmax=27.3,ymin=-0.13, ymax=0.95,fill=blue_palette[3])+
  geom_rect(xmin=9.7,xmax=17.3,ymin=-0.8, ymax=0.4,fill="#1FBFC0",alpha=0.7)+
  geom_rect(xmin=6.7,xmax=9.3,ymin=-0.8, ymax=0.4,fill="#D6B588",alpha=0.7)+
  geom_rect(xmin=0.7,xmax=6.3,ymin=-0.8, ymax=0.4,fill="#B589D6",alpha=0.7)+
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), color = "black")+
  geom_text(data = greece.dend.data$labels, aes(x, y, label = label), color = "black", #below for yellow gradient
            #color = as.factor(ord_num)), 
            hjust = 0, angle = 0, size = 3, fontface = "bold") +
  geom_point(data=greece.dend.helper,aes(x=1:17,y=-0.35,shape=word(greece.dend.helper$label, 1)))+
  #scale_color_manual(values = cur.pal[-1]) +
  #theme of plot
  theme(legend.position = "none",
        line = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  # clusters highlight
  #geom_segment(aes(x=0.8,xend=6.2,y=-0.38,yend=-0.38),size=1)+
  #geom_segment(aes(x=6.8,xend=9.2,y=-0.38,yend=-0.38),size=1)+
  #geom_segment(aes(x=9.8,xend=17.2,y=-0.38,yend=-0.38),size=1)+

  ylab("")+
  coord_flip()+
  scale_y_reverse(limits=c(1,-1))

#### MDS joining pseudoreplicates #### 
mds_greece <- cmdscale(vegdist(fish.com.mds, method = "bray")) %>%
  as_tibble() %>%
  mutate(habitat=row.names(cmdscale(vegdist(fish.com.mds, method = "bray")))) %>%
  left_join(greece.res.tib,by=c("habitat"="label")) %>%
  left_join(fish.com%>%
              select(1:3)%>%
              mutate(habitat=row.names(.)))%>%
  mutate(habitat=case_when(depth<=40 ~"shallow",
                           depth>40 & forest=="in" ~ "deep forest",
                           depth>40 & forest=="out"~ "deep bare"))%>%
  mutate(habitat=factor(habitat,levels=c("deep forest","deep bare","shallow")))
  

pal=c("#B589D6","#D6B588","#1FBFC0")

greece.mds.plot <- ggscatter(mds_greece, x = "V1", y = "V2", 
                          palette = pal,
                          size = 1.5,
                          fill="habitat",
                          shape="area",
                          color = "black",
                          ellipse.alpha = 0.7,
                          label=paste(mds_greece$area,mds_greece$depth,"m"),
                          ellipse = TRUE,
                          ellipse.type = "convex",
                          repel = TRUE)
greece.mds.plot

## with ggplot 
greece.hull.mds <- mds_greece %>% 
  group_by(cluster) %>% # using cluster from hierarchial clustering hclust()
  # chull: coordinates of points which constitute the outside border (=convex hull) of a cluster of points
  slice(chull(V1, V2)) 

ggplot(mds_greece) +
  # polygon shape of each 7 cluster
  geom_polygon(data = greece.hull.mds,
               aes(x = V1, y = V2, fill = habitat), alpha = 0.9) +
  geom_jitter(data = mds_greece, aes(x = V1, y = V2,fill = habitat,shape=area),
              color="black",size=2) +
  geom_text_repel(aes(x = V1, y = V2, label = paste(depth,"m")), 
                  size = 3, segment.color = "grey69")+
  theme_classic() +
  xlim(-0.5,0.5)+
  ylim(-0.65,0.4)+
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) 



#### MDS - FIGURE 2 - RIGHT ####
transect.com <- read.csv("fish-transects.csv") %>%
  filter(species !="fishing line") %>%
  mutate(depth=case_when(depth==41~40,
                         depth==19~20,
                         TRUE~depth)) %>%
  group_by(species, depth, transect.nb,area, forest) %>%
  summarize(abundance = sqrt(sum(number))) %>%
  spread(species, abundance, fill =  0) %>%
  as.data.frame()

rownames(transect.com) <- paste(transect.com$area,transect.com$depth,transect.com$forest,transect.com$transect.nb)
transect.com.mds <- transect.com[-c(1:4)]

## SIMPROF
transect.clust <- simprof(data = transect.com.mds, method.cluster = "complete", method.distance="braycurtis",
                        method.transform="identity", alpha=0.0000001,
                        sample.orientation="row", const=0,
                        silent=TRUE, increment=100,
                        undef.zero=TRUE, warn.braycurtis=TRUE)

# 2 groups from shallow to deep
transect.res.tib <- tibble(cl1 = list(transect.clust$significantclusters[[1]]),
                         cl2 = list(transect.clust$significantclusters[[2]]),
                         cl3 = list(transect.clust$significantclusters[[3]]),
                         cl4 = list(transect.clust$significantclusters[[4]]),
                         cl5 = list(transect.clust$significantclusters[[5]]),) %>%
  unnest_longer(cl1) %>%
  unnest_longer(cl2) %>%
  unnest_longer(cl3) %>%
  unnest_longer(cl4) %>%
  unnest_longer(cl5) %>%
  gather(key = "cluster", value = "label") %>% # transform wide to long format 
  distinct() 


# vegdist() : distance between each depth bin given relative abundance of species
# hlcust(): hierarchical clustering. 
#takes a dissimilarity table as produced by vegdist()
#returns the order of observations by similarity ($order), and the order of merges ($merge)
transect.clust.veg <- hclust(vegdist(transect.com.mds, method = "bray"), method = "complete")

mds_stat <- cmdscale(vegdist(transect.com.mds, method = "bray"), eig=TRUE) 

mds_transect <- cmdscale(vegdist(transect.com.mds, method = "bray")) %>%
  as_tibble() %>%
  mutate(habitat=row.names(cmdscale(vegdist(transect.com.mds, method = "bray")))) %>%
  left_join(transect.res.tib,by=c("habitat"="label")) %>%
  left_join(transect.com%>%
              select(1:4)%>%
              mutate(habitat=row.names(.)))%>%
  mutate(habitat=case_when(depth<=40 ~"shallow",
                           depth>40 & forest=="in" ~ "deep forest",
                           depth>40 & forest=="out"~ "deep bare"))%>%
  mutate(habitat=factor(habitat,levels=c("deep forest","deep bare","shallow")))

transect.hull.mds <- mds_transect %>% 
  filter(!(forest=="out"&depth=="80"))%>% # two outlier points
  group_by(habitat) %>% # using cluster from hierarchial clustering hclust()
  # chull: coordinates of points which constitute the outside border (=convex hull) of a cluster of points
  slice(chull(V1, V2)) 

ggplot() +
  # polygon shape of each 7 cluster
  geom_polygon(data = transect.hull.mds,
               mapping=aes(x = V1, y = V2, fill = habitat,color=habitat), alpha = 0.7) +
  geom_jitter(data = mds_transect, mapping=aes(x = V1, y = V2,shape=area,fill=habitat),
              color="black",size=2) +
  geom_text_repel(mds_transect,mapping=aes(x = V1, y = V2, label = paste(depth,"m")), 
                  size = 3, segment.color = "grey69")+
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 8, color = "black"),
    axis.title.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.title.y = element_text(size = 8, color = "black"),
  )+
  labs(x="nMDS axis 1", y="nMDS axis 2")+
  scale_fill_manual(values = pal)+
  scale_color_manual(values = pal)+
  scale_shape_manual(values=c(21,24))

#### 3.1. Indval ####
## per transect - 5 clusters
clusters <- transect.com.mds %>%
  mutate(label=row.names(.)) %>%
  left_join(transect.res.tib)

indval <- indval(x=transect.com.mds,clustering=clusters$cluster)
summary(indval)

## joining pseudoreplicates - 3 clusters
clusters2 <- fish.com.mds %>%
  mutate(label=row.names(.)) %>%
  left_join(greece.res.tib)

indval2 <- indval(x=fish.com.mds,clustering=clusters2$cluster)
summary(indval2, type="short")
summary(indval2, type="long")

#### multivariate regression tree #### 
mvpart(data.matrix(fish.com[,4:ncol(fish.com)])~depth+area+forest,fish.com,
       size =3)
  # R2= 1-Error = 1 - 0.423 = 0.567
mvpart(data.matrix(fish.com[,4:ncol(fish.com)])~depth+area+forest,fish.com,
       size =4)
  # R2= 0.665
greece.mrt <- mvpart(data.matrix(fish.com[,4:ncol(fish.com)])~depth+area+forest,fish.com,
       size =5)
  # R2=0.773

indval.mrt <- indval(fish.com.mds, greece.mrt$where)
summary(indval.mrt)
#### 3. Statistical tests ####
#### abundance ####
## total abundance
read.csv("fish-transects.csv") %>%
  filter(species!="fishing line") %>%
  summarize(tot.abu=sum(number))

## effect of site, depth, forest
no.fish <- data.frame(area=c(rep(c("Trikeri"),6),rep(c("Fourni"),2)),
                      depth=c(rep(c(50,60),each=3),85,85),
                              forest=rep(c("out"),8),
                      transect.nb=c("50_1","50_2","50_3","60_1","60_2","60_3","85_2","85_3"),
                      abundance=rep(c(1),8))

stat.abu <- read.csv("fish-transects.csv") %>%
  filter(species!="fishing line") %>%
  group_by(area,depth,forest,transect.nb) %>%
  summarize(abundance=sum(sqrt(number))) %>%
  bind_rows(no.fish)

abu.aov <- aov(abundance~area*depth+area+forest,stat.abu)
summary(abu.aov)
shapiro.test(residuals(abu.aov)) # p =0.27, normality of residuals OK
leveneTest(abundance~area*forest,stat.abu) #p=0.07, homogeneity of variance OK 

#### species richness ####
## total species richness
tot.richness <- read.csv("fish-transects.csv") %>%
  filter(species!="fishing line") %>%
  group_by(species)%>%
  summarize(abundance=sum(number))





#### 4. Traits per site ####
#### 4.1. Diets ####
fish.diets <- read.csv("fish-transects.csv") %>%
  filter(species!="fishing line") %>%
  mutate(depth=case_when(depth==41~40,
                         depth==19~20,
                         TRUE~depth),
         species=case_when(species=="gobius fallax/geniporus"~"gobius fallax",
                           TRUE~species)) %>%
  mutate(habitat=case_when(depth<=40 ~"shallow reef",
                           depth>40 & forest=="in" ~ "mesophotic forest",
                           depth>40 & forest=="out"~ "mesophotic barren")) %>%
  group_by(area,habitat,species) %>%
  summarize(abundance=sum(number))%>%
  left_join(read.csv("Quimbayo_trait_database.csv") %>% 
              select (species,Diet,Diel_activity,Size_group,Level_water,Size_class,Spawning,Body_size_max))%>%
  mutate(habitat=factor(habitat,levels=c("shallow reef","mesophotic forest","mesophotic barren")),
         Body_size_max=as.numeric(Body_size_max),
         Size_class2=case_when( Body_size_max>=100 ~ "xl",
                                Body_size_max>=50 & Body_size_max < 100 ~ "l",
                                Body_size_max>=30 & Body_size_max < 50 ~ "m",
                                Body_size_max>=10 & Body_size_max < 30 ~ "s",
                                Body_size_max<10~"xs"),
         Diet=factor(Diet, levels=rev(c("hd","pk","om","im","fc"))),
         Level_water=factor(Level_water,levels=rev(c("bottom","low","high"))),
         Spawning=factor(Spawning,levels=rev(c("demersal","oral","pelagic"))),
         Size_group=factor(Size_group,levels=rev(c("sol","pair","smallg","medg","largeg"))),
         Size_class2=factor(Size_class2,levels=rev(c("xs","s","m","l","xl"))))%>%
  arrange(desc(abundance))

diets <- ggplot(fish.diets, aes(x=habitat,y=sqrt(abundance),fill=Diet,color=area))+
  geom_bar(stat="identity",position="fill",color="black",size=0.04,width=0.8) +
  facet_wrap(~area)+
  scale_fill_manual(labels = c("fc" = "high predator", "pk" = "planktivore",
                               "om"="omnivore","im"="mobile invertivore","hd"="herbivore-detritivore"),
                    values=c("#d48c84","#EBC662","#143896","#ebdbfd", "#9bddb1"))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(labels=c(0,25,50,75,100))+
  labs(x="",y="% total abundance",fill="",title="diet")+
  #scale_fill_manual(values=c("#c03a76","grey","#fde5a7"))+
  #labs(y="number of species",x="",fill="IUCN status")+
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black",size=12),
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 9,color="black"),
        axis.text.y = element_text(size = 10,color="black"),
        legend.position = "bottom")

#### position in water column ####

position <- ggplot(fish.diets, aes(x=habitat,y=sqrt(abundance),fill=Level_water,color=area))+
  geom_bar(stat="identity",position="fill",color="black",size=0.04,width=0.8) +
  facet_wrap(~area)+
  scale_fill_manual(labels = c("fc" = "high predator", "pk" = "planktivore",
                               "om"="omnivore","im"="mobile invertivore","hd"="herbivore-detritivore"),
                    values=rev(c("#EBC662","#d48c84","#143896")))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(labels=c(0,25,50,75,100))+
  labs(x="",y="% total abundance",fill="",title="position in water column")+
  #scale_fill_manual(values=c("#c03a76","grey","#fde5a7"))+
  #labs(y="number of species",x="",fill="IUCN status")+
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black",size=12),
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 9,color="black"),
        axis.text.y = element_text(size = 10,color="black"),
        legend.position = "bottom")

#### Reproduction ####
reproduction <- ggplot(fish.diets, aes(x=habitat,y=sqrt(abundance),fill=Spawning,color=area))+
  geom_bar(stat="identity",position="fill",color="black",size=0.04,width=0.8) +
  facet_wrap(~area)+
  scale_fill_manual(labels = c("fc" = "high predator", "pk" = "planktivore",
                               "om"="omnivore","im"="mobile invertivore","hd"="herbivore-detritivore"),
                    values=rev(c("#EBC662","#d48c84","#143896")))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(labels=c(0,25,50,75,100))+
  labs(x="",y="% total abundance",fill="",title="reproduction mode")+
  #scale_fill_manual(values=c("#c03a76","grey","#fde5a7"))+
  #labs(y="number of species",x="",fill="IUCN status")+
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black",size=12),
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 9,color="black"),
        axis.text.y = element_text(size = 10,color="black"),
        legend.position = "bottom")


#### Gregariousness ####
gregariousness <- ggplot(fish.diets, aes(x=habitat,y=sqrt(abundance),fill=Size_group,color=area))+
  geom_bar(stat="identity",position="fill",color="black",size=0.04,width=0.8) +
  facet_wrap(~area)+
  scale_fill_manual(labels = c("fc" = "high predator", "pk" = "planktivore",
                               "om"="omnivore","im"="mobile invertivore","hd"="herbivore-detritivore"),
                    values=rev(c("#d48c84","#143896","#EBC662","#ebdbfd", "#9bddb1")))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(labels=c(0,25,50,75,100))+
  labs(x="",y="% total abundance",fill="",title="gregariousness")+
  #scale_fill_manual(values=c("#c03a76","grey","#fde5a7"))+
  #labs(y="number of species",x="",fill="IUCN status")+
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black",size=12),
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 9,color="black"),
        axis.text.y = element_text(size = 10,color="black"),
        legend.position = "bottom")

#### Size ####
size <- ggplot(fish.diets, aes(x=habitat,y=sqrt(abundance),fill=Size_class2,color=area))+
  geom_bar(stat="identity",position="fill",color="black",size=0.04,width=0.8) +
  facet_wrap(~area)+
  scale_fill_manual(labels = c("fc" = "high predator", "pk" = "planktivore",
                               "om"="omnivore","im"="mobile invertivore","hd"="herbivore-detritivore"),
                    values=c("#d48c84","#EBC662", "#9bddb1","#ebdbfd","#143896"))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(labels=c(0,25,50,75,100))+
  labs(x="",y="% total abundance",fill="",title="body size")+
  #scale_fill_manual(values=c("#c03a76","grey","#fde5a7"))+
  #labs(y="number of species",x="",fill="IUCN status")+
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black",size=12),
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 9,color="black"),
        axis.text.y = element_text(size = 10,color="black"),
        legend.position = "bottom")
size

ggarrange(size, diets, position,gregariousness,reproduction, nrow=3,ncol=2)

#### 4. Traits per habitat ####
#### 4.1. Diets ####
fish.diets <- read.csv("fish-transects.csv") %>%
  filter(species!="fishing line") %>%
  mutate(depth=case_when(depth==41~40,
                         depth==19~20,
                         TRUE~depth),
         species=case_when(species=="gobius fallax/geniporus"~"gobius fallax",
                           TRUE~species)) %>%
  mutate(habitat=case_when(depth<=40 ~"shallow reef",
                           depth>40 & forest=="in" ~ "mesophotic forest",
                           depth>40 & forest=="out"~ "mesophotic barren")) %>%
  group_by(habitat,species) %>%
  summarize(abundance=sum(number))%>%
  left_join(read.csv("Quimbayo_trait_database.csv") %>% 
              select (species,Diet,Diel_activity,Size_group,Level_water,Size_class,Spawning,Body_size_max))%>%
  mutate(habitat=factor(habitat,levels=c("shallow reef","mesophotic forest","mesophotic barren")),
         Body_size_max=as.numeric(Body_size_max),
         Size_class2=case_when( Body_size_max>=100 ~ "xl",
                                Body_size_max>=50 & Body_size_max < 100 ~ "l",
                                Body_size_max>=30 & Body_size_max < 50 ~ "m",
                                Body_size_max>=10 & Body_size_max < 30 ~ "s",
                                Body_size_max<10~"xs"),
         Diet=factor(Diet, levels=rev(c("hd","pk","om","im","fc"))),
         Level_water=factor(Level_water,levels=rev(c("bottom","low","high"))),
         Spawning=factor(Spawning,levels=rev(c("demersal","oral","pelagic"))),
         Size_group=factor(Size_group,levels=rev(c("sol","pair","smallg","medg","largeg"))),
         Size_class2=factor(Size_class2,levels=rev(c("xs","s","m","l","xl"))))%>%
  arrange(desc(abundance))

diets <- ggplot(fish.diets, aes(x=habitat,y=sqrt(abundance),fill=Diet))+
  geom_bar(stat="identity",position="fill",color="black",size=0.04,width=0.8) +
  scale_fill_manual(labels = c("fc" = "high predator", "pk" = "planktivore",
                               "om"="omnivore","im"="mobile invertivore","hd"="herbivore-detritivore"),
                    values=c("#d48c84","#EBC662","#143896","#ebdbfd", "#9bddb1"))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(labels=c(0,25,50,75,100))+
  labs(x="",y="% total abundance",fill="",title="diet")+
  #scale_fill_manual(values=c("#c03a76","grey","#fde5a7"))+
  #labs(y="number of species",x="",fill="IUCN status")+
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black",size=12),
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 9,color="black"),
        axis.text.y = element_text(size = 10,color="black"),
        legend.position = "bottom")

#### position in water column ####

position <- ggplot(fish.diets, aes(x=habitat,y=sqrt(abundance),fill=Level_water))+
  geom_bar(stat="identity",position="fill",color="black",size=0.04,width=0.8) +
  scale_fill_manual(labels = c("fc" = "high predator", "pk" = "planktivore",
                               "om"="omnivore","im"="mobile invertivore","hd"="herbivore-detritivore"),
                    values=rev(c("#EBC662","#d48c84","#143896")))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(labels=c(0,25,50,75,100))+
  labs(x="",y="% total abundance",fill="",title="position in water column")+
  #scale_fill_manual(values=c("#c03a76","grey","#fde5a7"))+
  #labs(y="number of species",x="",fill="IUCN status")+
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black",size=12),
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 9,color="black"),
        axis.text.y = element_text(size = 10,color="black"),
        legend.position = "bottom")

position <- ggplot(fish.diets, aes(x=habitat,y=sqrt(abundance),fill=Level_water))+
  geom_bar(stat="identity",position="stack",color="black",size=0.04,width=0.8) +
  scale_fill_manual(labels = c("fc" = "high predator", "pk" = "planktivore",
                               "om"="omnivore","im"="mobile invertivore","hd"="herbivore-detritivore"),
                    values=rev(c("#EBC662","#d48c84","#143896")))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  #scale_y_continuous(labels=c(0,25,50,75,100))+
  labs(x="",y="total abundance",fill="",title="position in water column")+
  #scale_fill_manual(values=c("#c03a76","grey","#fde5a7"))+
  #labs(y="number of species",x="",fill="IUCN status")+
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black",size=12),
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 9,color="black"),
        axis.text.y = element_text(size = 10,color="black"),
        legend.position = "bottom")

#### Reproduction ####
reproduction <- ggplot(fish.diets, aes(x=habitat,y=sqrt(abundance),fill=Spawning))+
  geom_bar(stat="identity",position="fill",color="black",size=0.04,width=0.8) +
  scale_fill_manual(labels = c("fc" = "high predator", "pk" = "planktivore",
                               "om"="omnivore","im"="mobile invertivore","hd"="herbivore-detritivore"),
                    values=rev(c("#EBC662","#d48c84","#143896")))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(labels=c(0,25,50,75,100))+
  labs(x="",y="% total abundance",fill="",title="reproduction mode")+
  #scale_fill_manual(values=c("#c03a76","grey","#fde5a7"))+
  #labs(y="number of species",x="",fill="IUCN status")+
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black",size=12),
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 9,color="black"),
        axis.text.y = element_text(size = 10,color="black"),
        legend.position = "bottom")


#### Gregariousness ####
gregariousness <- ggplot(fish.diets, aes(x=habitat,y=sqrt(abundance),fill=Size_group))+
  geom_bar(stat="identity",position="fill",color="black",size=0.04,width=0.8) +
  scale_fill_manual(labels = c("fc" = "high predator", "pk" = "planktivore",
                               "om"="omnivore","im"="mobile invertivore","hd"="herbivore-detritivore"),
                    values=rev(c("#d48c84","#143896","#EBC662","#ebdbfd", "#9bddb1")))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(labels=c(0,25,50,75,100))+
  labs(x="",y="% total abundance",fill="",title="gregariousness")+
  #scale_fill_manual(values=c("#c03a76","grey","#fde5a7"))+
  #labs(y="number of species",x="",fill="IUCN status")+
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black",size=12),
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 9,color="black"),
        axis.text.y = element_text(size = 10,color="black"),
        legend.position = "bottom")

#### Size ####
size <- ggplot(fish.diets, aes(x=habitat,y=sqrt(abundance),fill=Size_class2))+
  geom_bar(stat="identity",position="fill",color="black",size=0.04,width=0.8) +
  scale_fill_manual(labels = c("fc" = "high predator", "pk" = "planktivore",
                               "om"="omnivore","im"="mobile invertivore","hd"="herbivore-detritivore"),
                    values=c("#d48c84","#EBC662", "#9bddb1","#ebdbfd","#143896"))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(labels=c(0,25,50,75,100))+
  labs(x="",y="% total abundance",fill="",title="body size")+
  #scale_fill_manual(values=c("#c03a76","grey","#fde5a7"))+
  #labs(y="number of species",x="",fill="IUCN status")+
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black",size=12),
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        legend.title = element_text(size = 10),
        title= element_text(size = 10),
        axis.text.x = element_text(size = 9,color="black"),
        axis.text.y = element_text(size = 10,color="black"),
        legend.position = "bottom")
size

ggarrange(size, diets + rremove("ylab"), position+ rremove("ylab"),gregariousness,
          reproduction+ rremove("ylab"), nrow=2,ncol=3,align="h")
