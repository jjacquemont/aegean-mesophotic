library(tidyverse)
library(ggplot2)
#remotes::install_github("ropensci/rfishbase")
library(rfishbase)
library(dplyr)
library(FD)
library(fishualize)
library(mgcv) #gam model
library(modelr) #data manipulation such as data_grid
library(ggpubr) #combining ggplot figures
library(viridis)
library(vegan) # ecological and environmental data analysis
library(clustsig) 
library(ggdendro)
library(labdsv) # indval - indicator species
library(mvpart) # multivariate regression tree

cryptos <- read.csv("crypto.captured.csv") %>%
  mutate(species=paste(Genus,Species))

## total nb of species 
cryptos %>%
  distinct(species) #20

# total species per site and habitat
cryptos %>%
  filter(site!="ciotat") %>%
  mutate(habitat=case_when(Depth<=40 ~"shallow",
                           Depth>40 & forest=="in" ~ "deep forest",
                           Depth>40 & forest=="out"~ "deep bare")) %>%
  mutate(habitat=factor(habitat,levels=c("shallow","deep bare","deep forest"))) %>%
  group_by(habitat,site,species) %>%
  summarize(abundance=n()) %>%
  group_by(habitat,site) %>%
  summarize(species=n())

## most abundant species
abundant.cryptos <- cryptos %>%
  #filter(site!="ciotat") %>%
  mutate(habitat=case_when(Depth<=40 ~"shallow",
                           Depth>40 & forest=="in" ~ "deep forest",
                           Depth>40 & forest=="out"~ "deep barren")) %>%
  mutate(habitat=factor(habitat,levels=c("shallow","deep forest","deep barren"))) %>%
  mutate(species=case_when(str_detect(species, "Scorpaena") ~ "Scorpaena spp.",
                            TRUE~ species)) %>%
  group_by(species, habitat, site) %>%
  summarize(abundance=n())

ggplot(abundant.cryptos,aes(x=habitat,y=abundance,fill=species))+
  geom_bar(stat="identity",position="fill",width=0.8) +
  facet_wrap(~site)+
  #scale_fill_viridis_d(option="magma",begin=0.2,end=0.9)+
  #scale_fill_manual(values=rev(c("#c53c74","#fc9065","#EBC662","#99c5c4","grey")))+
  labs(y="proportion of total abundance",x="")+
  theme(strip.placement = "outside",
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 9,color="black"),
        axis.text.y = element_text(size = 10,color="black"))

#### station characteristics ####
#### total abundance ####
stations <- read.csv("crypto-station.csv")%>%
  filter(site!="La Ciotat") %>%
  mutate(habitat=case_when(depth<=40 ~"shallow",
                           depth>40 & forest=="in" ~ "deep forest",
                           depth>40 & forest=="out"~ "deep barren")) %>%
  mutate(habitat=factor(habitat,levels=c("shallow","deep forest","deep barren"))) 

mean.abu <- stations %>%
  group_by(habitat,site) %>%
  summarize(mean.abu=mean(abundance),se=sd(abundance)/n(),replicate=n())

mean.abu2 <- stations %>%
  group_by(habitat) %>%
  summarize(mean.abu=mean(abundance),se=sd(abundance)/n(),replicate=n())

abundance <- ggplot(stations,aes(x=site,y=abundance,fill=habitat))+
  geom_boxplot(alpha=0.8, width = 0.55,position=position_dodge((width = 0.7)),
               outlier.shape = NA)+
  geom_point(shape=21,position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7))+
  geom_text(mean.abu, mapping=aes(x=site,y=-0.5,label=paste("n=",replicate)),
            position=position_dodge((width = 0.7)),size=3)+
  scale_x_discrete(expand = c(0, 0))+
  theme_classic()+
  theme(axis.text.x = element_text(color = "black"))+
  scale_fill_manual(values=c("#1FBFC0","#B589D6","#D6B588"))+
  labs(y=expression(paste("individuals. m"^-2)),x="",title="abundance")
abundance

ggplot(stations,aes(x=habitat,y=abundance,fill=habitat,color=site))+
  geom_boxplot(alpha=0.8, width = 0.55,position=position_dodge((width = 0.7)),
               outlier.shape = NA)+
  geom_point(shape=21,position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7))+
  geom_text(mean.abu, mapping=aes(x=habitat,y=-0.5,label=paste("n=",replicate)),
            position=position_dodge((width = 0.7)),size=3)+
  scale_x_discrete(expand = c(0, 0))+
  theme_classic()+
  theme(axis.text.x = element_text(color = "black"))+
  scale_fill_manual(values=c("#6CBDE9","#D6B588","#FD7C6E"))+
  scale_color_manual(values=c("black","black"))+
  labs(y="individuals. station-1",x="",title="cryptic abundance")

summary(aov(data=stations, abundance~habitat)) # not signif - param
summary(aov(data=stations, abundance~as.factor(site))) # signif
summary(aov(data=stations, abundance~as.factor(site)*habitat)) # signif, only site
TukeyHSD(aov(data=stations, abundance~habitat)) # pairwise comparisons

kruskal.test(data=stations, abundance~habitat) #not signif - non param

## non-parametric tests
pairwise.wilcox.test(stations$abundance, stations$habitat,
                     p.adjust.method = "BH")

trikeri.abu <- stations %>% filter(site=="Trikeri")
pairwise.wilcox.test(trikeri.abu$abundance, trikeri.abu$habitat,
                     p.adjust.method = "BH")

fourni.abu <- stations %>% filter(site=="Fourni")
pairwise.wilcox.test(fourni.abu$abundance, fourni.abu$habitat,
                     p.adjust.method = "BH")

#### richness ####
mean.richness <- stations %>%
  group_by(habitat,site) %>%
  summarize(mean.richness=mean(richness),se=sd(richness)/n(),replicate=n())

richness <- ggplot(stations,aes(x=site,y=richness,fill=habitat))+
  geom_boxplot(alpha=0.8, width = 0.55,position=position_dodge((width = 0.7)),
               outlier.shape = NA)+
  geom_point(shape=21,position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7))+
  geom_text(mean.abu, mapping=aes(x=site,y=-0.5,label=paste("n=",replicate)),
            position=position_dodge((width = 0.7)),size=3)+
  scale_x_discrete(expand = c(0, 0))+
  theme_classic()+
  theme(axis.text.x = element_text(color = "black"))+
  scale_fill_manual(values=c("#1FBFC0","#B589D6","#D6B588"))+
  labs(y=expression(paste("species. m"^-2)),x="",title="species richness")
richness

ggplot(stations,aes(x=habitat,y=richness,fill=habitat,color=site))+
  geom_boxplot(alpha=0.8, width = 0.55,position=position_dodge((width = 0.7)),
               outlier.shape = NA)+
  geom_point(shape=21,position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7))+
  geom_text(mean.abu, mapping=aes(y=-0.5,label=paste("n=",replicate)),
            position=position_dodge((width = 0.7)),size=3)+
  scale_x_discrete(expand = c(0, 0))+
  theme_classic()+
  theme(axis.text.x = element_text(color = "black"))+
  scale_fill_manual(values=c("#1FBFC0","#B589D6","#D6B588"))+
  scale_color_manual(values=c("black","black"))+
  labs(y="species station-1",x="",title="cryptic richness")

summary(aov(data=stations, richness~habitat)) # not signif habitat - param
summary(aov(data=stations, richness~as.factor(site))) # signif site
summary(aov(data=stations, richness~as.factor(site)*habitat)) # signif, only site

kruskal.test(data=stations, richness~habitat) #not signif habitat - non param

## non-parametric tests
pairwise.wilcox.test(stations$richness, stations$habitat,
                     p.adjust.method = "BH")

trikeri.sp <- stations %>% filter(site=="Trikeri")
pairwise.wilcox.test(trikeri.sp$richness, trikeri.sp$habitat,
                     p.adjust.method = "BH")

fourni.abu <- stations %>% filter(site=="Fourni")
pairwise.wilcox.test(fourni.abu$richness, fourni.sp$habitat,
                     p.adjust.method = "BH")

#### length ####
length <- read.csv("crypto.captured.csv") %>%
  mutate(species=paste(Genus,Species))%>%
  filter(site!="ciotat") %>%
  mutate(habitat=case_when(Depth<=40 ~"shallow",
                           Depth>40 & forest=="in" ~ "mesophotic forest",
                           Depth>40 & forest=="out"~ "mesophotic barren")) %>%
  mutate(habitat=factor(habitat,levels=c("shallow","mesophotic forest","mesophotic barren"))) 

mean.length <- length %>%
  filter(!is.na(length)) %>%
  group_by(habitat,site) %>%
  summarize(mean.length=mean(length),se=sd(length)/n(),replicate=n())


length.plot <-ggplot(length,aes(x=site,y=length/10,fill=habitat))+
  geom_boxplot(alpha=0.8, width = 0.55,position=position_dodge(width = 0.7),
               outlier.shape = NA)+
  geom_point(shape=21,position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7))+
  theme_classic()+
  scale_x_discrete(expand=expansion(mult = c(0.05, 0.05)))+
  theme(axis.text.x = element_text(color = "black"))+
  geom_text(mean.length, mapping=aes(x=site,y=-0.5,label=paste("n=",replicate)),
            position=position_dodge((width = 0.7)),size=3)+
  scale_fill_manual(values=c("#1FBFC0","#B589D6","#D6B588"))+
  labs(y="cm",x="",title="individual length")
length.plot

ggplot(length,aes(x=habitat,y=length,fill=habitat,color=site))+
  geom_boxplot(alpha=0.8, width = 0.55,position=position_dodge(width = 0.7),
               outlier.shape = NA)+
  geom_point(shape=21,position=position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7))+
  theme_classic()+
  scale_x_discrete(expand=expansion(mult = c(0.05, 0.05)))+
  theme(axis.text.x = element_text(color = "black"))+
  geom_text(mean.length, mapping=aes(y=-0.5,label=paste("n=",replicate)),
            position=position_dodge((width = 0.7)),size=3)+
  scale_fill_manual(values=c("#6CBDE9","#D6B588","#FD7C6E"))+
  scale_color_manual(values=c("black","black"))+
  labs(y="individual length (mm)",x="",title="cryptic length")

summary(aov(data=length, length~habitat)) # not signif habitat - param
summary(aov(data=length, length~as.factor(site))) # pvalue for site: 0.06
summary(aov(data=length, length~as.factor(site)*habitat)) #  pvalue for site: 0.06
TukeyHSD(aov(data=length, length~habitat)) # pairwise comparisons - not signif habitat

kruskal.test(data=length, length~habitat) #not signif habitat - non param

## non-parametric tests
pairwise.wilcox.test(length$length, length$habitat,
                     p.adjust.method = "BH")

trikeri.length <- length %>% filter(site=="Trikeri")
pairwise.wilcox.test(trikeri.length$length, trikeri.length$habitat,
                     p.adjust.method = "BH")

fourni.length <- length %>% filter(site=="Fourni")
pairwise.wilcox.test(fourni.length$length, fourni.length$habitat,
                     p.adjust.method = "BH")

#### three plots together ####
ggarrange(abundance, richness, length.plot, nrow=3,ncol=1,legend="right",
          common.legend=T)