setwd("~/Desktop/R/RNAseq_MDII/July25-reanalysis/")

library('ggplot2')
library('ggrepel')
library('ggbreak')
library('dplyr')
library('tidyr')
library('stringr')
library('data.table')
library('ggprism')
library('ggh4x')
library('dplyr') 
library('Rmisc')

spu <- read.table("~/Desktop/TAMU/Spurp_reference_files/Spurp_genome/conversion_references/refseqLocus_spu_IDmapping.txt", 
                  sep = "\t", header = TRUE, quote = "")
colnames(spu) <- c("names", "SPU", "sp_name")

grn <- read.table("~/Desktop/Maroon_Debut/MD-ref-files/Davidson2022_supp1_dGRNlist.txt", sep = "\t") #192 genes
colnames(grn) = c("SPU", "sp_name_GW", "short_GW", "sp_description_GW")

LOCgrn <- merge(grn, spu, by = "SPU") 
#note: increases in length due to multiple LOC IDs fitting some spu IDs
#same thing happens the other way, where some of the same LOC IDs have multiple spuIDs

RW <- read.csv("HMRvsLMR_DEG_725.csv") 
RW <- subset(RW, select = -c(baseMean, lfcSE, stat, pvalue)) 
colnames(RW) <- c("names", "log2FoldChange", "padj")
RW$DE <- "No change"
RW$DE[RW$log2FoldChange > 0.5 & RW$padj < 0.05] <- "Up"
RW$DE[RW$log2FoldChange < -0.5 & RW$padj < 0.05] <- "Down"
head(RW)

RT <- read.csv("18vs14_DEG_725.csv")
RT <- subset(RT, select = -c(baseMean, lfcSE, stat, pvalue)) 
colnames(RT) <- c("names", "log2FoldChange", "padj")
RT$DE <- "No change"
RT$DE[RT$log2FoldChange > 0.5 & RT$padj < 0.05] <- "Up"
RT$DE[RT$log2FoldChange < -0.5 & RT$padj < 0.05] <- "Down"
head(RT)

all <- merge(RW, RT, by = "names", all = TRUE)
all <- unite(all, DE_WT, "DE.x", "DE.y", sep = "-")

dgrn <- merge(all, LOCgrn, by = "names")
colnames(dgrn) = c("names", "lfcW", "padjW", "DE_WT", 
                   "lfcTemp", "padjTemp", 
                   "SPU", "sp_name_GW", "shortGW", 
                   "sp_description_GW", "sp_name")
length(unique(dgrn$names)) #177
length(unique(dgrn$SPU)) #179
#dgrn dataframe has 212 rows 

dgrnPlot <- dgrn %>% 
  group_by(names, lfcW, DE_WT, lfcTemp) %>% 
  dplyr::summarise(spu_count = n())
#177 total on list 

dewtColors <- c("#CC0033", "#33CCCC", "#CC0033", 
                "azure4", "#FFCC00", "azure4", "#FFCC00",
                "#CC0033", "#33CCCC", "#CC0033")
names(dewtColors) <- c("Down-Down", "Down-No change", "Down-Up", 
                       "NA-No change", "No change-Down", 
                       "No change-No change", "No change-Up", 
                       "Up-Down", "Up-No change", "Up-Up")

devPlot <- ggplot() + 
  theme_prism() + 
  #ylim(-6, 4) +
  #xlim(-4, 4) +
  geom_vline(xintercept = c(0, 0), col = "red") +
  geom_hline(yintercept = c(0, 0), col = "red") +
  geom_point(data = dgrnPlot,
             aes(x = lfcW, y = lfcTemp, color = DE_WT)) +
  scale_color_manual(values = dewtColors) + 
  theme(legend.position = "none") +
  geom_rect(data = dgrnPlot, 
            aes(xmin = 0 - 0.5, 
                xmax = 0 + 0.5, 
                ymin = 0 - 0.5, 
                ymax = 0 + 0.5), 
            fill = "transparent", color = "purple", size = 0.2) +
  geom_text_repel(data = subset(dgrnPlot, !(DE_WT == "No change-No change" | 
                                           names %like% "LOC")), 
                  aes(x = lfcW, y = lfcTemp, label = names), 
                  size = 4, color = "black", max.overlaps = Inf) +  
  xlab(expression(bold("Microbial Richness"~"-"~"log"[2]~Fold~Change~"\n"))) +
  ylab(expression(bold(Temp~"-"~"log"[2]~Fold~Change~"\n"))) +
  theme(panel.grid.major = element_line(color = "#666666",
                                        size = 0.3,
                                        linetype = 2)) + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) + 
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))
print(devPlot)

ggsave(plot = devPlot, 
       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/devRnaPlot725.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

#lists of DEGs
#177 genes present in RNAseq dataset that are dGRN
degGrn <- subset(dgrnPlot, !(DE_WT == "No change-No change")) 
#70 of these are DEG

#write.csv(degGrn, "degGrn_list_725.csv")

ncDown <- subset(dgrnPlot, DE_WT == "No change-Down") 
length(unique(ncDown$names)) #9 LocIDs down only from temp
#12.9% (9/70)
ncDownSpu <- merge(ncDown, LOCgrn, by = "names")

ncUp <- subset(dgrnPlot, DE_WT == "No change-Up") 
length(unique(ncUp$names)) #36 LocIDs Up only from temp 
#51.4%

UpNC <- subset(dgrnPlot, DE_WT == "Up-No change") 
length(unique(UpNC$names)) #6 LocIDs Up only from MBE 
#8.6%
 
DownNC <- subset(dgrnPlot, DE_WT == "Down-No change") 
length(unique(DownNC$names)) #6 LocIDs Down only from MBE 
6/70 #8.6%

UpDown <- subset(dgrnPlot, DE_WT == "Up-Down")
length(unique(UpDown$names)) #9 LocIDs Up in MBE, Down in temp 
9/111 #8.1%
UpDownSpu <- merge(UpDown, LOCgrn, by = "names")

DownUp <- subset(dgrnPlot, DE_WT == "Down-Up")
length(unique(DownUp$names)) #8 LocIDs Down in MBE, Up in temp 
8/70 #11.4% 

# x-plots of specific dGRNs ----
vsdpvals <- read.csv("MDII_rldpvals_July25.csv")

colData <- data.frame(row.names = colnames(vsdpvals),
                      sample = c(colnames(vsdpvals)))
colData <- colData %>% 
  mutate(water = case_when(
    str_detect(sample, "St") ~ "low", #low richness
    str_detect(sample, "F") ~ "mbe"  #high richness
  )) %>% 
  mutate(temp = case_when(
    str_detect(sample, "14") ~ "14", 
    str_detect(sample, "18") ~ "18"
  )) 
colData$treat = paste(colData$water, colData$temp, sep = "") 

vsdpvals <- as.data.frame(vsdpvals)
head(vsdpvals)

vsdpvals[grepl("blimp1", vsdpvals$X), ]
cand1= vsdpvals[grepl("hh", vsdpvals$X), ]
cand2= vsdpvals[grepl("LOC576047", vsdpvals$X), ] #smo
cand3= vsdpvals[grepl("LOC593924", vsdpvals$X), ] #ptch1
cand4= vsdpvals[grepl("FGFR", vsdpvals$X), ]
skelCands = rbind(cand1, cand2, cand3, cand4)

geneL= skelCands %>%
  gather("sample", "expr", -X)

geneL <- geneL %>% 
  mutate(water = case_when(
    str_detect(sample, "St") ~ "low", #low richness
    str_detect(sample, "F") ~ "mbe"  #high richness
  )) %>% 
  mutate(temp = case_when(
    str_detect(sample, "14") ~ "14", 
    str_detect(sample, "18") ~ "18"
  )) 
geneL <- subset(geneL, !(is.na(water)))
#geneL$treat = paste(geneL$water, geneL$temp, sep = "")
#geneL$time=str_sub(geneL$sample, 6,6)
#geneL$time<-sub("2","24", geneL$time)
#geneL$time=as.factor(geneL$time)

summ=summarySE(data=geneL,measurevar="expr",groupvars=c("X","temp","water"))

pd <- position_dodge(0.1)
skelGenes=ggplot(summ,aes(x= temp,y=expr, color= water))+
  geom_point(aes(shape= water),size=3,position=pd)+
  geom_line(aes(group= water,linetype= water),position=pd)+
  geom_errorbar(aes(ymin=expr-se,ymax=expr+se),lwd=0.4,width=0.3,position=pd)+
  scale_color_manual(values=c("darkorange","darkred", "darkgreen","darkblue"))+
  scale_shape_manual(values=c(15,16,17,18))+
  theme_minimal()+
  facet_wrap(~ X,scales="free_y", ncol=5)+
  theme(legend.text=element_text(size=10)) +
  theme(legend.key = element_blank())+
  theme(legend.direction = 'horizontal', legend.position = 'top')
print(skelGenes)

# other named genes x-plot ---
vsdpvals[grepl("delta", vsdpvals$X), ]
cand1= vsdpvals[grepl("LOC764088", vsdpvals$X), ] #numb
cand2= vsdpvals[grepl("blimp1", vsdpvals$X), ]
cand3= vsdpvals[grepl("GATAe", vsdpvals$X), ]
cand4= vsdpvals[grepl("Gsc", vsdpvals$X), ]
cand5= vsdpvals[grepl("LOC585551", vsdpvals$X), ] #gemin2
cand6= vsdpvals[grepl("LOC592057", vsdpvals$X), ] #hesC
cand7= vsdpvals[grepl("Six1", vsdpvals$X), ]
cand8= vsdpvals[grepl("snail", vsdpvals$X), ]
cand9= vsdpvals[grepl("Nanos2", vsdpvals$X), ]
cand10= vsdpvals[grepl("Nanos2", vsdpvals$X), ]
cand11= vsdpvals[grepl("delta", vsdpvals$X), ]

Cands = rbind(cand1, cand2, cand3, cand4, cand5, 
              cand6, cand7, cand8, cand9, cand10, cand11)

geneL= Cands %>%
  gather("sample", "expr", -X)

geneL <- geneL %>% 
  mutate(water = case_when(
    str_detect(sample, "St") ~ "low", #low richness
    str_detect(sample, "F") ~ "mbe"  #high richness
  )) %>% 
  mutate(temp = case_when(
    str_detect(sample, "14") ~ "14", 
    str_detect(sample, "18") ~ "18"
  )) 
geneL <- subset(geneL, !(is.na(water)))

summ=summarySE(data=geneL,measurevar="expr",groupvars=c("X","temp","water"))

pd <- position_dodge(0.1)
Genes=ggplot(summ,aes(x= temp,y=expr, color= water))+
  geom_point(aes(shape= water),size=3,position=pd)+
  geom_line(aes(group= water,linetype= water),position=pd)+
  geom_errorbar(aes(ymin=expr-se,ymax=expr+se),lwd=0.4,width=0.3,position=pd)+
  scale_color_manual(values=c("darkorange","darkred", "darkgreen","darkblue"))+
  scale_shape_manual(values=c(15,16,17,18))+
  theme_minimal()+
  facet_wrap(~ X,scales="free_y", ncol=5)+
  theme(legend.text=element_text(size=10)) +
  theme(legend.key = element_blank())+
  theme(legend.direction = 'horizontal', legend.position = 'top')
print(Genes)
