library(tidyverse)
library(ggplot2)
library(dplyr)
library(fishualize)
library(mgcv) #gam model
library(modelr) #data manipulation such as data_grid
library(ggpubr) #combining ggplot figures
library(viridis)

fish.data <- read.csv("fish-transects.csv") %>%
  filter(species !="fishing line")%>%
  mutate(depth=case_when(depth==41~40,
                         depth==19~20,
                         TRUE~depth)) %>%
  mutate(habitat=case_when(depth<=40 ~"shallow",
                           depth>40 & forest=="in" ~ "mesophotic forest",
                           depth>40 & forest=="out"~ "mesophotic barren"),
         habitat=factor(habitat,levels=c("shallow","mesophotic forest","mesophotic barren")))

#### 1. Fishing lines ####
fishing.lines <- read.csv("fish-transects.csv") %>%
  filter(species=="fishing line") %>%
  mutate(depth=case_when(depth==41~40,
                         depth==19~20,
                         TRUE~depth)) %>%
  group_by(area, transect.nb,depth,forest) %>%
  summarize (fishing.lines=sum(number))%>%
  bind_rows(data_frame(area=rep("Trikeri",3),depth=rep(55,3),forest=rep("out",3),
                       fishing.lines=rep(0,3))) %>%
  mutate(habitat=case_when(depth<=40 ~"shallow",
                           depth>40 & forest=="in" ~ "mesophotic forest",
                           depth>40 & forest=="out"~ "mesophotic barren"),
         habitat=factor(habitat,levels=c("shallow","mesophotic forest","mesophotic barren")),
         depth=factor(depth)) %>%
  group_by(area,depth,habitat) %>%
  summarize (replicates=n(), fishing.lines.mean=mean(fishing.lines))

fishing.line.means <- read.csv("fish-transects.csv") %>%
  filter(species=="fishing line") %>%
  mutate(depth=case_when(depth==41~40,
                         depth==19~20,
                         TRUE~depth)) %>%
  group_by(area, transect.nb,depth,forest) %>%
  summarize (fishing.lines=sum(number))%>%
  bind_rows(data_frame(area=rep("Trikeri",3),depth=rep(55,3),forest=rep("out",3),
                       fishing.lines=rep(0,3))) %>%
  mutate(habitat=case_when(depth<=40 ~"shallow",
                           depth>40 & forest=="in" ~ "mesophotic forest",
                           depth>40 & forest=="out"~ "mesophotic barren"),
         habitat=factor(habitat,levels=c("shallow","mesophotic forest","mesophotic barren"))) %>%
  group_by(area, habitat) %>%
  summarize (replicates=n(), fishing.lines.mean=mean(fishing.lines)) %>%
  mutate(depth="overall") %>%
  bind_rows(fishing.lines)%>%
  mutate(depth=factor(depth,levels=rev(c("10","20","30","40","50","55","60","80","85","90","overall"))))
  

fishing.line.means2 <- read.csv("fish-transects.csv") %>%
  filter(species=="fishing line") %>%
  mutate(depth=case_when(depth==41~40,
                         depth==19~20,
                         TRUE~depth)) %>%
  group_by(area, transect.nb,depth,forest) %>%
  summarize (fishing.lines=sum(number))%>%
  bind_rows(data_frame(area=rep("Trikeri",3),depth=rep(55,3),forest=rep("out",3),
                       fishing.lines=rep(0,3))) %>%
  mutate(habitat=case_when(depth<=40 ~"shallow",
                           depth>40 & forest=="in" ~ "deep forest",
                           depth>40 & forest=="out"~ "deep baren")) %>%
  group_by(area, habitat) %>%
  #group_by(area) %>%
  summarize (replicates=n(), fishing.lines.mean=mean(fishing.lines), se=sd(fishing.lines)/sqrt(n()))

## stat test
fishing.line.stat <- read.csv("fish-transects.csv") %>%
  filter(species=="fishing line") %>%
  mutate(depth=case_when(depth==41~40,
                         depth==19~20,
                         TRUE~depth)) %>%
  group_by(area, transect.nb,depth,forest) %>%
  summarize (fishing.lines=sum(number))%>%
  mutate(habitat=case_when(depth<=40 ~"shallow",
                           depth>40 & forest=="in" ~ "deep forest",
                           depth>40 & forest=="out"~ "deep barren"),
         habitat=factor(habitat))

lines.aov <- aov(fishing.lines~habitat+area,data=fishing.line.stat)
summary(lines.aov)

## in Trikeri
fishing.line.stat.pigadi <- fishing.line.stat%>% filter(area=="Trikeri")
pairwise.t.test(fishing.line.stat.pigadi$fishing.lines, 
                                    fishing.line.stat.pigadi$habitat, paired = F,
                                    p.adjust.method = "bonferroni")

## in Fourni
fishing.line.stat.fourni <- fishing.line.stat%>% filter(area=="Fourni")
pairwise.t.test(fishing.line.stat.fourni$fishing.lines, 
                fishing.line.stat.fourni$habitat, paired = F,
                p.adjust.method = "bonferroni")


## heatmap
color_scale<- c("white", viridis(6, option = "magma", direction = -1,begin=0.4)[2:5])

ggplot(data=fishing.line.means, aes(x=depth, y=habitat,fill=fishing.lines.mean))+
  geom_tile(data = fishing.line.means %>% filter(depth != "overall"),
            color="black",size=0.2,height=0.95)+
  geom_tile(data = fishing.line.means %>% filter(depth == "overall"), 
            aes(x=depth, y=habitat), 
            width = 0.9,height=0.95,size=0.2,color="black") +
  #scale_fill_viridis_c(option="magma", begin=1, end=0.4)+
  scale_fill_gradientn(colors = color_scale, values = scales::rescale(c(0, 1, 2, 3, 4, 5)))+
  facet_wrap(~area, strip.position = "top")+
  theme(strip.placement = "outside",
    panel.background = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(color = "black",size=10),
        axis.line.y = element_line(),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 10,color="black"),
        axis.text.y = element_text(size = 10,color="black"),
       legend.position = "bottom")+
        #panel.grid.major = element_blank())+
        #panel.grid.minor = element_blank())+
   scale_y_discrete(position="right",labels=str_wrap(levels(fishing.line.means$habitat),width=8))+
  labs(y="",x="depth",fill=expression(paste("fishing lines. 120 m"^{-2})))+
  coord_flip()
  

## avg in forest
read.csv("fish-transects.csv") %>%
  filter(species=="fishing line" & area=="pigadi" & forest=="in") %>%
  group_by(transect.nb) %>%
  summarize (number=sum(number)) %>%
  summarize (fishing.lines=mean(number), se=sd(number)*1.96/n(),replicates=n())

read.csv("fish-transects.csv") %>%
  filter(species=="fishing line" & area=="fourni" & forest=="in") %>%
  group_by(transect.nb) %>%
  summarize (number=sum(number)) %>%
  summarize (fishing.lines=mean(number), se=sd(number)*1.96/n(),replicates=n())




#### 2. vulnerability & fishing data ####

#### IUCN status ####
fishing.targets <- fish.data %>%
  filter(species!="ind") %>%
  group_by(area,habitat,species) %>%
  summarize(abundance=sum(number))%>%
  left_join(read.csv("IUCN_fishMed.csv"))%>%
  group_by(area,habitat,IUCN)%>%
  summarize(count=n()) %>%
  mutate(IUCN=factor(IUCN, levels=c("VU","DD","LC")))

ggplot(fishing.targets, aes(x=habitat,y=count,fill=IUCN))+
  geom_bar(stat="identity",position="stack") +
  facet_wrap(~area)+
  scale_fill_manual(values=c("#c03a76","grey","#fde5a7"))+
  labs(y="number of species",x="",fill="IUCN status")+
  theme(strip.placement = "outside",
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 9,color="black"),
        axis.text.y = element_text(size = 10,color="black"))

#### targeted by fishery - IUCN ####
fishing.targets <- fish.data %>%
  filter(species!="ind") %>%
  group_by(area,habitat,species) %>%
  summarize(abundance=sum(number))%>%
  left_join(read.csv("IUCN_fishMed.csv"))%>%
  group_by(area,habitat,fishery)%>%
  summarize(count=n()) %>%
  mutate(fishery=factor(fishery, levels=c("primary","secondary","no")))

ggplot(fishing.targets, aes(x=habitat,y=count,fill=fishery))+
  geom_bar(stat="identity",position="stack") +
  facet_wrap(~area)+
  scale_fill_manual(values=c("#c03a76","#fecd90","#BAE0F3"))+
  labs(y="number of species",x="",fill="fishery target")+
  theme(strip.placement = "outside",
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 9,color="black"),
        axis.text.y = element_text(size = 10,color="black"))

## target by aquarium trade
fishing.targets <- fish.data %>%
  filter(species!="ind") %>%
  group_by(area,habitat,species) %>%
  summarize(abundance=sum(number))%>%
  left_join(read.csv("IUCN_fishMed.csv"))%>%
  group_by(area,habitat,aquarium)%>%
  summarize(count=n()) %>%
  mutate(aquarium=factor(aquarium, levels=c("yes","rare","no")))

ggplot(fishing.targets, aes(x=habitat,y=count,fill=aquarium))+
  geom_bar(stat="identity",position="stack") +
  facet_wrap(~area)+
  scale_fill_manual(values=c("#c03a76","#fecd90","#BAE0F3"))+
  labs(y="number of species",x="",fill="aquarium trade")+
  theme(strip.placement = "outside",
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 9,color="black"),
        axis.text.y = element_text(size = 10,color="black"))

#### 3. targeted by fishery - Sini et al., 2019 ####
#### 3.1. By species richness ####

fishing.targets <- fish.data %>%
  filter(species!="ind") %>%
  group_by(area,habitat,species) %>%
  summarize(abundance=sum(number))%>%
  left_join(read.csv("IUCN_fishMed.csv")) %>%
  group_by(area,habitat,commercial_status,species)%>%
  summarize(count=n()) %>%
  mutate(commercial_status=factor(commercial_status, levels=c("NC","aquarium-trade","LC","C-sp","C")))

fish.data %>%
  filter(species!="ind") %>%
  group_by(species) %>%
  summarize(abundance=sum(number))%>%
  left_join(read.csv("IUCN_fishMed.csv")) %>%
  group_by(commercial_status) %>%
  summarize(nb_species=n())

ggplot(fishing.targets, aes(x=habitat,y=count,fill=commercial_status))+
  geom_bar(stat="identity",position="stack",width=0.8) +
  facet_wrap(~area)+
  #scale_fill_viridis_d(option="magma",begin=0.2,end=0.9)+
  scale_fill_manual(values=rev(c("#c53c74","#fc9065","#EBC662","#99c5c4","grey")))+
  labs(y="number of species",x="",fill="commercial status")+
  scale_x_discrete(labels=str_wrap(levels(fishing.targets$habitat),width=8))+
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black",size=12),
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 9,color="black"),
        axis.text.y = element_text(size = 10,color="black"))

#### 3.2. By abundance ####
fishing.targets <- fish.data %>%
  filter(species!="ind") %>%
  group_by(area,habitat,species) %>%
  summarize(abundance=sum(number))%>%
  left_join(read.csv("IUCN_fishMed.csv")) %>%
  group_by(area,habitat,commercial_status,species)%>%
  summarize(abundance=sum(abundance)) %>%
  mutate(commercial_status=factor(commercial_status, levels=c("NC","aquarium-trade","LC","C-sp","C")))


ggplot(fishing.targets, aes(x=habitat,y=abundance,fill=commercial_status))+
  geom_bar(stat="identity",position="stack",width=0.8) +
  facet_wrap(~area)+
  #scale_fill_viridis_d(option="magma",begin=0.2,end=0.9)+
  scale_fill_manual(values=rev(c("#c53c74","#fc9065","#EBC662","#99c5c4","grey")))+
  labs(y="abundance",x="",fill="commercial status")+
  scale_x_discrete(labels=str_wrap(levels(fishing.targets$habitat),width=8))+
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(color = "black",size=12),
        panel.background = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(),
        legend.title = element_text(size = 10),
        axis.text.x = element_text(size = 9,color="black"),
        axis.text.y = element_text(size = 10,color="black"))
