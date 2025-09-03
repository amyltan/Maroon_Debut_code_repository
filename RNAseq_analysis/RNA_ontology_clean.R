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

spu <- read.table("~/Desktop/TAMU/Spurp_reference_files/Spurp_genome/conversion_references/refseqLocus_spu_IDmapping.txt", 
                  sep = "\t", header = TRUE, quote = "")
colnames(spu) <- c("names", "SPU", "sp_name")

onto <- read.csv("~/Desktop/TAMU/Spurp_reference_files/Spurp_genome/Tu_Ontology_SuppTableS2.csv")

LOConto <- merge(onto, spu, by = "SPU") 

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

unchar <- read.csv("~/Desktop/Maroon_Debut/phmmer-out/phmmer-proteinList-all.csv")
unchar <- subset(unchar, select = -c(X))
colnames(unchar)[colnames(unchar) == "locID"] <- "names"
colnames(unchar)[colnames(unchar) == "protein"] <- "description_unchar"

all <- merge(RW, RT, by = "names", all = TRUE)
all <- unite(all, DE_WT, "DE.x", "DE.y", sep = "-")

conv <- read.csv("~/Desktop/TAMU/Spurp_reference_files/Spurp_genome/AmyGeneNameConversions.csv")
colnames(conv) <- c("names", "description") 
all <- merge(all, conv, by = "names", all = TRUE)
all <- all[!is.na(all$DE_WT),]

allOnto <- merge(all, LOConto, by = "names")
table(allOnto$Class.L1) 
#8650 ontology entries representing 5530 genes 

colnames(allOnto) <- c("names", "lfc_W", "padj_W", "DE_WT", "lfc_temp", 
                   "padj_temp", "description","SPU", "Class.L1", 
                   "Class.L2", "Class.L3", "sp_name")

#allOnto <- allOnto %>% 
#  group_by(names, Class.L1, DE_WT, lfc_MBE, padj_MBE,
#           lfc_temp, padj_temp, description_temp) %>% 
#  summarise(count = n()) 

RWonto <- merge(RW, LOConto, by = "names")
table(RWonto$Class.L1)

RTonto <- merge(RT, LOConto, by = "names")
table(RTonto$Class.L1)

dewtColors <- c("#CC0033", "#33CCCC", "#CC0033", 
                "azure4", "#FFCC00", "azure4", "#FFCC00",
                "#CC0033", "#33CCCC", "#CC0033")
names(dewtColors) <- c("Down-Down", "Down-No change", "Down-Up", 
                       "NA-No change", "No change-Down", 
                       "No change-No change", "No change-Up", 
                       "Up-Down", "Up-No change", "Up-Up")

# Level 1 onto numbers plot - these are not in publication #####
plotOntoW <- RWonto %>% 
  group_by(Class.L1, DE) %>% 
  summarise(count = n())

plotOntoT <- RTonto %>% 
  group_by(Class.L1, DE) %>% 
  summarise(count = n())

ggplot() +
  geom_col(data = plotOntoW, 
           aes(x = Class.L1, y = count, fill = DE)) + 
  scale_fill_manual(values = c("blue", "azure4", "coral2")) +
  coord_flip() +
  xlab("Ontology class (level 1)") + 
  ylab ("# of genes") 

ggplot() +
  geom_col(data = plotOntoT, 
           aes(x = Class.L1, y = count, fill = DE)) + 
  scale_fill_manual(values = c("blue", "azure4", "coral2")) +
  coord_flip() +
  xlab("Ontology class (level 1)") + 
  ylab ("# of genes")

deg <- merge(RW, RT, by = "names", all = TRUE)
colnames(deg) <- c("names", "lfc_W", "padj_W", "DE_W", 
                   "lfc_temp", "padj_temp", "DE_temp")
deg <- merge(deg, LOConto, by = "names")

plotRW <- RWonto %>% 
  group_by(Class.L1, DE) %>% 
  summarise(count = n()) %>% 
  mutate(count = if_else(DE == "Down", -count, count))
plotRW$cond <- "Microbial Richness"

plotRT <- RTonto %>% 
  group_by(Class.L1, DE) %>% 
  summarise(count = n()) %>% 
  mutate(count = if_else(DE == "Down", -count, count))
plotRT$cond <- "Temperature"

ggplot() + 
  theme_bw() + 
  ylim(-400, 650) +
  scale_y_break(c(250,500), scale = 0.1) +
  geom_col(data = subset(plotRW, !(DE == "No change")), 
           aes(x = Class.L1, y = count, fill = cond),
           width = 0.4,
           position = ggplot2::position_nudge(x = -0.2)) + 
  geom_col(data = subset(plotRT, !(DE == "No change")),
           aes(x = Class.L1, y = count, fill = cond), 
           width = 0.4, 
           position = ggplot2::position_nudge(x = 0.2)) + 
  scale_fill_manual("Condition", 
                    values = c("#336666", "#CC9933")) +
  geom_hline(yintercept = c(0), col = "black") + 
  coord_flip() + 
  xlab(NULL) + 
  ylab("# of genes")


# NOTE TO SELF: Removed the defunct IMMUNITY plot code ----

# NOTE TO SELF: Removed the defunct METABOLISM plot code ----

# Transcription factors = TF ##### 
tf <- subset(allOnto, Class.L1 == "TF")

#without znf
ggplot() + 
  theme_prism() + 
  #ylim(-3, 3) +
  #xlim(-3, 3) +
  geom_vline(xintercept = c(0, 0), col = "red") +
  geom_hline(yintercept = c(0, 0), col = "red") +
  geom_point(data = subset(tf, !(Class.L2 %like% "ZNF")),
             aes(x = lfc_W, y = lfc_temp, color = DE_WT)) + 
  geom_rect(data = tf, 
            aes(xmin = 0 - 0.5, 
                xmax = 0 + 0.5, 
                ymin = 0 - 0.5, 
                ymax = 0 + 0.5), 
            fill = "transparent", color = "purple", size = 0.2) +
  geom_text_repel(data = subset(tf, (DE_WT == "Up-Up" | 
                                        DE_WT == "Down-Down" | 
                                        DE_WT == "Up-Down" | 
                                        DE_WT == "Down-Up")), 
                  aes(x = lfc_W, y = lfc_temp, label = names), 
                  size = 2, color = "black", max.overlaps = Inf) + 
  xlab(expression(MBE~"-"~"log"[2]~Fold~Change)) +
  ylab(expression(Temp~"-"~"log"[2]~Fold~Change)) +
  scale_color_manual("TF gene\nexpression", values = dewtColors, 
                     labels = 
                       c("MBE Down, Temp Down", 
                         "MBE Down", 
                         "MBE Down, Temp Up", 
                         "Temp NC", 
                         "Temp Down", 
                         "No change", 
                         "Temp Up", 
                         "MBE Up, Temp Down", 
                         "MBE Up", 
                         "MBE Up, Temp Up")) +
  theme(legend.position = "none")

#finding other TFs not included in Tu ontology - this is a list of the names searched
#bHLH
all[all$description %like% "helix-loop-helix",] 
unchar[unchar$protein %like% "bHLH",] #"Helix-loop-helix"

#bZIP - leucine zipper
all[all$description.y %like% "leucine zipper",] 
unchar[unchar$protein %like% "leucine zipper",] 


#ETS
all[all$description.y %like% "Ets",] #ETS, ets, Ets
unchar[unchar$protein %like% "Ets",] #none

#Forkhead
all[all$description.y %like% "fox",] #Forkhead, forkhead, fox
unchar[unchar$protein %like% "Forkhead",] #none

#GATA
all[all$description.y %like% "GATA",]
#GATAc, GATAe, LOC105436739, LOC115926363, LOC591657
unchar[unchar$protein %like% "GATA",] #none

#Homeodomain
all[all$description.y %like% "homeobox",] #homeodomain, homeobox
unchar[unchar$protein %like% "Hox",] 

#Nuclear receptor
all[all$description.y %like% "nuclear receptor",]
unchar[unchar$protein %like% "nuclear receptor",]

#Sox
all[all$description.y %like% "Sox",]
unchar[unchar$protein %like% "Sox",] #none 

#T-box
all[all$description.y %like% "T-box",]
#LOC576774, LOC592389, Tbx2/3
unchar[unchar$protein %like% "T-box",] #none

#creating larger list - 9 families
bhlh <- subset(all, description %like% "helix-loop-helix")
bhlh2 <- subset(unchar, (description_unchar %like% "bHLH" |
                              description_unchar %like% "Helix-loop-helix")) 
bhlh$TFfam <- "bHLH"
bhlh <- select(bhlh, c(names, TFfam))
bhlh2$TFfam <- "bHLH"
bhlh2 <- select(bhlh2, c(names, TFfam))
bhlhR <- rbind(bhlh, bhlh2)

bzip <- subset(all, description %like% "leucine zipper")
bzip2 <- subset(unchar, description_unchar %like% "leucine zipper" )
bzip$TFfam <- "bZIP"
bzip <- select(bzip, c(names, TFfam))
bzip2$TFfam <- "bZIP"
bzip2 <- select(bzip2, c(names, TFfam))
bzipR <- rbind(bzip, bzip2)

ets <- subset(all,description %like% "ETS" | 
                description %like% "ets" | 
                description %like% "Ets")
ets$TFfam <- "Ets"
etsR <- select(ets, c(names, TFfam))

fox <- subset(all, description %like% "Forkhead" |
                description %like% "forkhead" |
                description %like% "fox")
fox$TFfam <- "Fox"
foxR <- select(fox, c(names, TFfam))

gata <- subset(all, description %like% "GATA")
gata$TFfam <- "GATA"
gataR <- select(gata, c(names, TFfam))


hox <- subset(all, description %like% "homeobox" |
                description %like% "homeodomain")
hox2 <- subset(unchar, description_unchar %like% "Hox")
hox$TFfam <- "Hox"
hox <- select(hox, c(names, TFfam))
hox2$TFfam <- "Hox"
hox2 <- select(hox2, c(names, TFfam))
hoxR <- rbind(hox, hox2)

nr <- subset(all, description %like% "nuclear receptor")
nr2 <- subset(unchar, description_unchar %like% "nuclear receptor")
nr$TFfam <- "NR"
nr <- select(nr, c(names, TFfam))
nr2$TFfam <- "NR"
nr2 <- select(nr2, c(names, TFfam))
nrR <- rbind(nr, nr2)

sox <- subset(all, description %like% "Sox")
sox$TFfam <- "Sox"
soxR <- select(sox, c(names, TFfam))

tbox <- subset(all, description %like% "T-box")
tbox$TFfam <- "Tbox"
tboxR <- select(tbox, c(names, TFfam))

tfRef <- rbind(bhlhR, bzipR, etsR, foxR, gataR, hoxR, 
               nrR, soxR, tboxR)

tfAll <- subset(all, (description %like% "helix-loop-helix" | 
                        description %like% "leucine zipper" | 
                        description %like% "ETS" | 
                        description %like% "ets" | 
                        description %like% "Ets" |
                        description %like% "Forkhead" |
                        description %like% "forkhead" |
                        description %like% "fox" |
                        description %like% "GATA" |
                        description %like% "homeobox" |
                        description %like% "homeodomain" |
                        description %like% "nuclear receptor" |
                        description %like% "Sox" |
                        description %like% "T-box"))

tfcomp <- merge(tf, tfAll, by = "names", all = TRUE) 

tfUnchar <- subset(unchar, (description_unchar %like% "bHLH" |
                              description_unchar %like% "Helix-loop-helix" |
                              description_unchar %like% "leucine zipper" |
                              description_unchar %like% "Hox" |
                              description_unchar %like% "nuclear receptor"))
tfU <- all[all$names %in% tfUnchar$names,]

tfcomp <- merge(tfcomp, tfU, by = "names", all = TRUE)
tfcomp <- merge(tfcomp, tfRef, by = "names", all = TRUE) #450
length(unique(tfcomp$names)) #394

tfcomp <- tfcomp %>% 
  group_by(names, Class.L1, Class.L2, TFfam, 
           DE_WT.x, lfc_W, padj_W,
           lfc_temp, padj_temp, description.x) %>% 
  summarise(count = n()) 

tfcomp$Class.L2 <- gsub("TF_", "", tfcomp$Class.L2) 
tfcomp$Class.L2 <- gsub("Homeo", "Hox", tfcomp$Class.L2) 
tfcomp$Class.L2 <- gsub("Forkhead", "Fox", tfcomp$Class.L2) 
tfcomp$Class.L2 <- gsub("bzip", "bZIP", tfcomp$Class.L2) 
tfcomp$Class.L2 <- gsub("Sox/HMG", "Sox", tfcomp$Class.L2) 

tfcomp$TFfam <- ifelse(is.na(tfcomp$TFfam), tfcomp$Class.L2, tfcomp$TFfam)
tfnoznf <- subset(tfcomp, !(TFfam == "ZNF")) #note: this set already ignored ZNFs, but this is pulling out all the ZNFs added via Tu's sheet
#write.csv(tfnoznf, "~/Desktop/R/RNAseq_MDII/July25-reanalysis/supp-info725/allTF_DEG_Fig4L_725.csv")

tfplot <- ggplot() + 
  theme_prism() + 
  #ylim(-5, 3) +
  #xlim(-3, 3) +
  geom_vline(xintercept = c(0, 0), col = "red") +
  geom_hline(yintercept = c(0, 0), col = "red") +
  geom_point(data = tfnoznf,
             aes(x = lfc_W, y = lfc_temp, color = DE_WT.x)) + 
  geom_rect(data = tf, 
            aes(xmin = 0 - 0.5, 
                xmax = 0 + 0.5, 
                ymin = 0 - 0.5, 
                ymax = 0 + 0.5), 
            fill = "transparent", color = "purple", size = 0.2) +
  geom_text_repel(data = subset(tfnoznf, ((DE_WT.x == "Up-Up" | 
                                           DE_WT.x == "Down-Down" | 
                                           DE_WT.x == "Up-Down" |
                                           DE_WT.x == "Down-Up" |
                                            abs(lfc_temp > 1.25) | 
                                           abs(lfc_W > 1)) &
                                            !(TFfam == "Other"))), 
                  aes(x = lfc_W, y = lfc_temp, label = TFfam), 
                  size = 5, color = "black", max.overlaps = Inf) + 
  xlab(expression(bold("Microbial Richness"~"-"~"log"[2]~Fold~Change))) +
  ylab(expression(bold(Temp~"-"~"log"[2]~Fold~Change))) +
  scale_color_manual("TF gene\nexpression", values = dewtColors, 
                     labels = 
                       c("MBE Down, Temp Down", 
                         "MBE Down", 
                         "MBE Down, Temp Up", 
                         "Temp NC", 
                         "Temp Down", 
                         "No change", 
                         "Temp Up", 
                         "MBE Up, Temp Down", 
                         "MBE Up", 
                         "MBE Up, Temp Up")) +
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
plot(tfplot)
ggsave(plot = tfplot, 
       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/tf-expression-Plot725.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

table(tfnoznf$TFfam)
#bHLH  bZIP   Ets   Fox  GATA   Hox    NR Other  Smad   Sox  Tbox 
#  47    31    14    26     5   103    44    56     8    10     5 

# TF with matching motifs ####
#note: only Sox appears in Tu ontology list
#Nfia (nucelar factor IA) - yes "nuclear factor", NFkB !(inhibitor)
#Rfx - yes "regulatory factor X"
#Mtf (metal-responsive) - yes "metal regulatory", "metal-response"
#MSANTD(Myb/sant dna bd containing) - yes "myb", "Myb" (also unchar-Myb)
#FIGLA (factor in germline alpha) - no 
#TCF - yes "Tcf"
#Klf - yes "Kruppel"
#Sox - yes "Sox" 
#Znf - purposefully not including: too many with unknown functions
all[all$description %like% "Tcf",]
unchar[unchar$description_unchar %like% "Myb",]

nfia <- subset(all, description %like% "nuclear factor")
nfkb <- subset(all, description %like% "NFkB") 
nfia <- rbind(nfia, nfkb)
nfia$TFfam <- "Nfia" 
nfia <- subset(nfia, !(description %like% "inhibitor"))

rfx <- subset(all, description %like% "regulatory factor X")
rfx$TFfam <- "Rfx"

mtf <- subset(all, description %like% "metal regulatory")
mtf2 <- subset(all, description %like% "metal-response")
mtf <- rbind(mtf, mtf2)
mtf$TFfam <- "Mtf"

msantd <- subset(all, description %like% "myb")
msantd2 <- subset(all, description %like% "Myb")
unchar[unchar$description_unchar %like% "Myb",]
msantd3 <- subset(all, names %like% "LOC100889035")
msantd <- rbind(msantd, msantd2, msantd3)
msantd$TFfam <- "MSANTD"

tcf <- subset(all, description %like% "Tcf")
tcf$TFfam <- "Tcf"

klf <- subset(all, description %like% "Kruppel")
klf$TFfam <- "Klf"

sox2 <- subset(all, description %like% "Sox")
sox2$TFfam <- "Sox"

tfplot2 <- rbind(nfia, rfx, mtf, msantd, tcf, klf, sox2)

colnames(tfplot2) <- c("names", "lfc_W", "padj_W", "DE_WT", 
                       "lfc_temp", 
                       "padj_temp", "description", "TFfam")

table(tfplot2$DE_WT)
#write.csv(tfplot2, "~/Desktop/R/RNAseq_MDII/July25-reanalysis/supp-info725/TFwithmotifs_DEG_Fig5D725.csv")

tfMotifPlot <- ggplot() + 
  theme_prism() + 
  #ylim(-3, 3) +
  #xlim(-3, 3) +
  geom_vline(xintercept = c(0, 0), col = "red") +
  geom_hline(yintercept = c(0, 0), col = "red") +
  geom_point(data = tfplot2,
             aes(x = lfc_W, y = lfc_temp, color = DE_WT)) + 
  geom_rect(data = tfplot2, 
            aes(xmin = 0 - 0.5, 
                xmax = 0 + 0.5, 
                ymin = 0 - 0.5, 
                ymax = 0 + 0.5), 
            fill = "transparent", color = "purple", size = 0.2) +
  geom_text_repel(data = subset(tfplot2, 
                                !(DE_WT == "No change-No change" | 
                                    TFfam == "MSANTD")), 
                  aes(x = lfc_W, y = lfc_temp, 
                      label = TFfam), 
                  size = 4, color = "black", max.overlaps = Inf) + 
  xlab(expression(bold("Microbial Richness"~"-"~"log"[2]~Fold~Change))) +
  ylab(expression(bold(Temp~"-"~"log"[2]~Fold~Change))) +
  scale_color_manual("TF gene\nexpression", values = dewtColors, 
                     labels = 
                       c("MBE Down, Temp Down", 
                         "MBE Down", 
                         "MBE Down, Temp Up", 
                         "Temp NC", 
                         "Temp Down", 
                         "No change", 
                         "Temp Up", 
                         "MBE Up, Temp Down", 
                         "MBE Up", 
                         "MBE Up, Temp Up")) +
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
plot(tfMotifPlot)

ggsave(plot = tfMotifPlot, 
       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/tf-Motif-Plot725.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

#znf quick check 
znf <- subset(all, description %like% "zinc finger")
table(znf$DE_WT) #223 no change; 112 DE
ggplot() + 
  theme_prism() + 
  ylim(-3, 3) +
  xlim(-3, 3) +
  geom_vline(xintercept = c(0, 0), col = "red") +
  geom_hline(yintercept = c(0, 0), col = "red") +
  geom_point(data = znf,
             aes(x = log2FoldChange.x, 
                 y = log2FoldChange.y, color = DE_WT)) + 
  geom_rect(data = tf, 
            aes(xmin = 0 - 0.5, 
                xmax = 0 + 0.5, 
                ymin = 0 - 0.5, 
                ymax = 0 + 0.5), 
            fill = "transparent", color = "purple", size = 0.2) +
  #geom_text_repel(data = znf, 
  #                aes(x = lfc_MBE, y = lfc_temp, label = names), 
  #                size = 2, color = "black", max.overlaps = Inf) + 
  xlab(expression(bold(MBE~"-"~"log"[2]~Fold~Change))) +
  ylab(expression(bold(Temp~"-"~"log"[2]~Fold~Change))) +
  scale_color_manual("TF gene\nexpression", values = dewtColors, 
                     labels = 
                       c("MBE Down, Temp Down", 
                         "MBE Down", 
                         "MBE Down, Temp Up", 
                         "Temp NC", 
                         "Temp Down", 
                         "No change", 
                         "Temp Up", 
                         "MBE Up, Temp Down", 
                         "MBE Up", 
                         "MBE Up, Temp Up")) +
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) +  
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))

# Metabolism - clean #####
met <- subset(allOnto, Class.L1 == "Metabolism")
#write.csv(met, "~/Desktop/R/RNAseq_MDII/July25-reanalysis/supp-info725/metabolismDEG_Fig5B725.csv")

met <- met %>% 
  group_by(names, Class.L1, DE_WT, lfc_W, padj_W,
           lfc_temp, padj_temp, description) %>% 
  summarise(count = n()) 
table(met$DE_WT) 

#Down-Down      Down-No change             Down-Up      No change-Down 
#       20                 141                  34                 177 
#No change-No change        No change-Up             Up-Down        Up-No change 
#               1025                 284                  34                  68 
#Up-Up 
#    9 
#1792
#DEG: 767

metPlot <- ggplot() + 
  theme_prism() + 
  #ylim(-6, 4) +
  #xlim(-4, 4) +
  geom_vline(xintercept = c(0, 0), col = "red") +
  geom_hline(yintercept = c(0, 0), col = "red") +
  geom_point(data = met,
             aes(x = lfc_W, y = lfc_temp, color = DE_WT)) + 
  geom_rect(data = met, 
            aes(xmin = 0 - 0.5, 
                xmax = 0 + 0.5, 
                ymin = 0 - 0.5, 
                ymax = 0 + 0.5), 
            fill = "transparent", color = "purple", size = 0.2) +
  geom_text_repel(data = subset(met, (DE_WT == "Up-Up" | 
                                           DE_WT == "Down-Down" | 
                                           DE_WT == "Up-Down" |
                                           DE_WT == "Down-Up") & 
                                  !(names %like% "LOC")), 
                  aes(x = lfc_W, 
                      y = lfc_temp, label = names), 
                  size = 4, color = "black", max.overlaps = Inf) + 
  xlab(expression(bold("Microbial Richness"~"-"~"log"[2]~Fold~Change))) +
  ylab(expression(bold(Temp~"-"~"log"[2]~Fold~Change))) +
  scale_color_manual("TF gene\nexpression", values = dewtColors, 
                     labels = 
                       c("MBE Down, Temp Down", 
                         "MBE Down", 
                         "MBE Down, Temp Up", 
                         "Temp NC", 
                         "Temp Down", 
                         "No change", 
                         "Temp Up", 
                         "MBE Up, Temp Down", 
                         "MBE Up", 
                         "MBE Up, Temp Up")) +
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
plot(metPlot)
ggsave(plot = metPlot, 
       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/metabolism-Plot725.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

metlookup <- subset(met, !(names %like% "LOC" | DE_WT == "No change-No change"))
metlookup2 <- subset(met, !(DE_WT %like% "No change")) #97
table(metlookup2$DE_WT)
#Down-Down   Down-Up   Up-Down     Up-Up 
#       20        34        34         9 
metDU <- subset(metlookup, DE_WT == "Down-Up") #3 with names
metUD <- subset(metlookup, DE_WT == "Up-Down") #1
metUU <- subset(metlookup, DE_WT == "Up-Up") #none
metDD <- subset(metlookup, DE_WT == "Down-Down") #1

#met GO output - needs stat added back 
head(met)
RW <- read.csv("HMRvsLMR_DEG_725.csv") 
RW <- subset(RW, select = -c(baseMean, log2FoldChange, lfcSE, pvalue, padj)) 
colnames(RW) <- c("names", "stat")

RT <- read.csv("18vs14_DEG_725.csv") 
RT <- subset(RT, select = -c(baseMean, log2FoldChange, lfcSE, pvalue, padj)) 
colnames(RT) <- c("names", "stat")

metRW <- merge(met, RW, by = "names") 
metRW <- subset(metRW, select = c("names", "stat"))
colnames(metRW) <- c("geneId", "stat")
metRT <- merge(met, RT, by = "names")
metRT <- subset(metRT, select = c("names", "stat"))
colnames(metRT) <- c("geneId", "stat")

write.csv(metRW, "~/Desktop/R/MDII_GO_MWU_725/metRW-waldStat-725.csv", 
          row.names = FALSE, quote = F)
write.csv(metRT, "~/Desktop/R/MDII_GO_MWU_725/metRT-waldStat-725.csv", 
          row.names = FALSE, quote = F)

# Immunity - clean #### 
immune <- subset(allOnto, Class.L1 == "Immunity")
#write.csv(immune, "~/Desktop/R/RNAseq_MDII/July25-reanalysis/supp-info725/immuneAll_DEG_Fig6F725.csv")

immune <- immune %>% 
  group_by(names, Class.L1, DE_WT, lfc_W, padj_W,
           lfc_temp, padj_temp, description) %>% 
  summarise(count = n()) 

table(immune$DE_WT)
#Down-Down      Down-No change             Down-Up      No change-Down 
#        3                  21                   7                  52 
#No change-No change        No change-Up             Up-Down        Up-No change 
#                368                 115                  14                  21 
#Up-Up 
#    1 
#total: 602
#DEG: 234
#temp only: 167
#microbes only: 42 
#both: 25

sigimmune <- subset(immune, !(DE_WT == "No change-No change"))
table(sigimmune$DE_WT)

immunePlot <- ggplot() + 
  theme_prism() + 
  #ylim(-6, 4) +
  #xlim(-4, 4) +
  geom_vline(xintercept = c(0, 0), col = "red") +
  geom_hline(yintercept = c(0, 0), col = "red") +
  geom_point(data = immune,
             aes(x = lfc_W, 
                 y = lfc_temp, color = DE_WT)) + 
  geom_rect(data = immune, 
            aes(xmin = 0 - 0.5, 
                xmax = 0 + 0.5, 
                ymin = 0 - 0.5, 
                ymax = 0 + 0.5), 
            fill = "transparent", color = "purple", size = 0.2) +
  geom_rect(data = immune, 
            aes(xmin = 0 - 2.5, 
                xmax = 0 + 1, 
                ymin = 0 - 2.5, 
                ymax = 0 + 1), 
            fill = "transparent", color = "#333333", size = 0.5, linetype = 2) +
  geom_text_repel(data = subset(immune, ((abs(lfc_W > 1)) |
                                  (abs(lfc_temp > 1))) & 
                                  !(names %like% "LOC")), 
                  aes(x = lfc_W, 
                      y = lfc_temp, label = names), 
                  size = 4, color = "black", max.overlaps = Inf) + 
  xlab(expression(bold("Microbial Richness"~"-"~"log"[2]~Fold~Change))) +
  ylab(expression(bold(Temp~"-"~"log"[2]~Fold~Change))) +
  scale_color_manual("TF gene\nexpression", values = dewtColors, 
                     labels = 
                       c("MBE Down, Temp Down", 
                         "MBE Down", 
                         "MBE Down, Temp Up", 
                         "Temp NC", 
                         "Temp Down", 
                         "No change", 
                         "Temp Up", 
                         "MBE Up, Temp Down", 
                         "MBE Up", 
                         "MBE Up, Temp Up")) +
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
plot(immunePlot)
ggsave(plot = immunePlot, 
       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/immune-Plot725.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

# pigment cells - clean #####

# note: gcm, pks1, sult1c2, abcg11, mif5, oprk are all markers of differentiated pigment cells 
cand1= all[grepl("gcm", all$names), ] #gcm 
cand2= all[grepl("LOC588806", all$names), ] #pks1 
cand3= all[grepl("LOC582581", all$names), ] #fmo2  #no change
cand4= all[grepl("LOC586240", all$names), ] #fmo3 #essential for pigmentation - Wessel 2020
cand5= all[grepl("LOC589643", all$names), ] #fmo2-2 
cand6= all[grepl("LOC588371", all$names), ] #fmo5-1 
cand7= all[grepl("LOC592500", all$names), ] #sult1c2 
cand8= all[grepl("LOC100891288", all$names), ] #oprk1c/L 
cand9= all[grepl("LOC754218", all$names), ] #mif5 #In already
cand10= all[grepl("LOC579779", all$names), ] #glur6    #up in temp
cand11= all[grepl("LOC590101", all$names), ] #abcg11 

Cands=rbind(cand1, cand2, cand3, cand4, cand5, cand6, cand7, 
            cand8, cand9, cand10, cand11)

pigment <- subset(all, names %in% Cands$names)
colnames(pigment) <- c("names", "lfc_W", "padj_W", "DE_WT", 
                       "lfc_temp", 
                       "padj_temp", "description")

#adding common name - ignore these 2 lines
#immuneLookUp <- read.csv("immune_genes_presentMDRNA.csv")
#pigment <- merge(pigment, immuneLookUp, by = "names")

#write.csv(pigment, "~/Desktop/R/RNAseq_MDII/July25-reanalysis/supp-info725/pigmentCell_DEG_Fig6G725.csv")

pigmentPlot <- ggplot() + 
  theme_prism() + 
  ylim(-2.5, 1) +
  xlim(-2.5, 1) +
  geom_vline(xintercept = c(0, 0), col = "red") +
  geom_hline(yintercept = c(0, 0), col = "red") +
  geom_point(data = pigment,
             aes(x = lfc_W, y = lfc_temp, color = DE_WT)) + 
  geom_rect(data = pigment, 
            aes(xmin = 0 - 0.5, 
                xmax = 0 + 0.5, 
                ymin = 0 - 0.5, 
                ymax = 0 + 0.5), 
            fill = "transparent", color = "purple", size = 0.2) +
  geom_text_repel(data = pigment, 
                  aes(x = lfc_W, y = lfc_temp, label = names), 
                  size = 4, color = "black", max.overlaps = Inf) + 
  xlab(expression(bold("Microbial Richness"~"-"~"log"[2]~Fold~Change))) +
  ylab(expression(bold(Temp~"-"~"log"[2]~Fold~Change))) +
  scale_color_manual("Immune gene\nexpression", values = dewtColors, 
                     labels = 
                       c("MBE Down, Temp Down", 
                         "MBE Down", 
                         "MBE Down, Temp Up", 
                         "Temp NC", 
                         "Temp Down", 
                         "No change", 
                         "Temp Up", 
                         "MBE Up, Temp Down", 
                         "MBE Up", 
                         "MBE Up, Temp Up")) +
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
plot(pigmentPlot)
ggsave(plot = pigmentPlot, 
       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/pigmentCell-Plot725.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)
# Removed dGRN #####

# Defensome - ABC transporters #### 
abc <- subset(allOnto, Class.L2 == "Defensome_TransporterABC")
unique(abc$Class.L3)

abc[abc$sp_name %like% "Abc", ]

abc <- abc %>% 
  group_by(names, Class.L1, DE_WT, lfc_W, padj_W,
           lfc_temp, padj_temp, description) %>% 
  summarise(count = n()) 

abcPlot <- ggplot() + 
  theme_prism() + 
  #ylim(-2.5, 1) +
  #xlim(-2, 1) +
  geom_vline(xintercept = c(0, 0), col = "red") +
  geom_hline(yintercept = c(0, 0), col = "red") +
  geom_point(data = abc,
             aes(x = lfc_W, y = lfc_temp, color = DE_WT)) + 
  geom_rect(data = abc, 
            aes(xmin = 0 - 0.5, 
                xmax = 0 + 0.5, 
                ymin = 0 - 0.5, 
                ymax = 0 + 0.5), 
            fill = "transparent", color = "purple", size = 0.2) +
  geom_text_repel(data = subset(abc, !(DE_WT %like% "No change")), 
                  aes(x = lfc_W, y = lfc_temp, label = names), 
                  size = 4, color = "black", max.overlaps = Inf) + 
  xlab(expression(bold(MBE~"-"~"log"[2]~Fold~Change))) +
  ylab(expression(bold(Temp~"-"~"log"[2]~Fold~Change))) +
  scale_color_manual("Immune gene\nexpression", values = dewtColors, 
                     labels = 
                       c("MBE Down, Temp Down", 
                         "MBE Down", 
                         "MBE Down, Temp Up", 
                         "Temp NC", 
                         "Temp Down", 
                         "No change", 
                         "Temp Up", 
                         "MBE Up, Temp Down", 
                         "MBE Up", 
                         "MBE Up, Temp Up")) +
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
print(abcPlot)
ggsave(plot = abcPlot, 
       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/abc-transporter-Plot725.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

#note: haven't written this file
#write.csv(abc, "~/Desktop/abc-transporter-expression.csv")

# TESTING INTERACTION INFORMATION ----
# Note: this hasn't been changed past the TF stuff because I'm mentally stuck on whether it's useful
I <- read.csv("interaction_DEG_725.csv") 
I <- subset(I, select = -c(baseMean, lfcSE, stat, pvalue)) 
colnames(I) <- c("names", "log2FoldChange", "padj")
I$DE <- "No change"
I$DE[I$log2FoldChange > 0.5 & I$padj < 0.05] <- "Up"
I$DE[I$log2FoldChange < -0.5 & I$padj < 0.05] <- "Down"
head(I)

WI <- read.csv("WI_DEG_725.csv") 
WI <- subset(WI, select = -c(baseMean, lfcSE, stat, pvalue)) 
colnames(WI) <- c("names", "log2FoldChange", "padj")
WI$DE <- "No change"
WI$DE[WI$log2FoldChange > 0.5 & WI$padj < 0.05] <- "Up"
WI$DE[WI$log2FoldChange < -0.5 & WI$padj < 0.05] <- "Down"
head(WI)

TI <- read.csv("TI_DEG_725.csv")
TI <- subset(TI, select = -c(baseMean, lfcSE, stat, pvalue)) 
colnames(TI) <- c("names", "log2FoldChange", "padj")
TI$DE <- "No change"
TI$DE[TI$log2FoldChange > 0.5 & TI$padj < 0.05] <- "Up"
TI$DE[TI$log2FoldChange < -0.5 & TI$padj < 0.05] <- "Down"
head(TI)

allI <- merge(I, WI, by = "names", all = TRUE) 
allI <- merge(allI, TI, by = "names", all = TRUE)
#all <- unite(all, DE_WT, "DE.x", "DE.y", sep = "-")

#conv <- read.csv("~/Desktop/TAMU/Spurp_reference_files/Spurp_genome/AmyGeneNameConversions.csv")
#colnames(conv) <- c("names", "description") 
allI <- merge(allI, conv, by = "names", all = TRUE)
allI <- allI[!is.na(allI$log2FoldChange.x),]

allOntoI <- merge(allI, LOConto, by = "names")
table(allOntoI$Class.L1) 
#8650 ontology entries representing 5530 genes 

colnames(allOntoI) <- c("names", "lfc_I", "padj_I", "DE_I", 
                        "lfc_WI", "padj_WI", "DE_WI",
                        "lfc_TI", "padj_TI", "DE_TI", "description",
                        "SPU", "Class.L1", 
                       "Class.L2", "Class.L3", "sp_name")

#allOnto <- allOnto %>% 
#  group_by(names, Class.L1, DE_WT, lfc_MBE, padj_MBE,
#           lfc_temp, padj_temp, description_temp) %>% 
#  summarise(count = n()) 

degIColors <- c("#CC0033", "#FFCC00", "#33CCCC")
names(degIColors) <- c("F18", "F18vsF14", "F18vsS18")

# Transcription factors = TF ##### 
tfI <- subset(allOntoI, Class.L1 == "TF")

#without znf
ggplot() + 
  theme_prism() + 
  #ylim(-3, 3) +
  #xlim(-3, 3) +
  geom_vline(xintercept = c(0, 0), col = "red") +
  geom_hline(yintercept = c(0, 0), col = "red") +
  geom_point(data = subset(tfI, !(Class.L2 %like% "ZNF")),
             aes(x = lfc_W, y = lfc_temp, color = DE_WT)) + 
  geom_rect(data = tf, 
            aes(xmin = 0 - 0.5, 
                xmax = 0 + 0.5, 
                ymin = 0 - 0.5, 
                ymax = 0 + 0.5), 
            fill = "transparent", color = "purple", size = 0.2) +
  geom_text_repel(data = subset(tf, (DE_WT == "Up-Up" | 
                                       DE_WT == "Down-Down" | 
                                       DE_WT == "Up-Down" | 
                                       DE_WT == "Down-Up")), 
                  aes(x = lfc_W, y = lfc_temp, label = names), 
                  size = 2, color = "black", max.overlaps = Inf) + 
  xlab(expression(MBE~"-"~"log"[2]~Fold~Change)) +
  ylab(expression(Temp~"-"~"log"[2]~Fold~Change)) +
  scale_color_manual("TF gene\nexpression", values = dewtColors, 
                     labels = 
                       c("MBE Down, Temp Down", 
                         "MBE Down", 
                         "MBE Down, Temp Up", 
                         "Temp NC", 
                         "Temp Down", 
                         "No change", 
                         "Temp Up", 
                         "MBE Up, Temp Down", 
                         "MBE Up", 
                         "MBE Up, Temp Up")) +
  theme(legend.position = "none")

#finding other TFs not included in Tu ontology - this is a list of the names searched
#bHLH
all[all$description %like% "helix-loop-helix",] 
unchar[unchar$protein %like% "bHLH",] #"Helix-loop-helix"

#bZIP - leucine zipper
all[all$description.y %like% "leucine zipper",] 
unchar[unchar$protein %like% "leucine zipper",] 


#ETS
all[all$description.y %like% "Ets",] #ETS, ets, Ets
unchar[unchar$protein %like% "Ets",] #none

#Forkhead
all[all$description.y %like% "fox",] #Forkhead, forkhead, fox
unchar[unchar$protein %like% "Forkhead",] #none

#GATA
all[all$description.y %like% "GATA",]
#GATAc, GATAe, LOC105436739, LOC115926363, LOC591657
unchar[unchar$protein %like% "GATA",] #none

#Homeodomain
all[all$description.y %like% "homeobox",] #homeodomain, homeobox
unchar[unchar$protein %like% "Hox",] 

#Nuclear receptor
all[all$description.y %like% "nuclear receptor",]
unchar[unchar$protein %like% "nuclear receptor",]

#Sox
all[all$description.y %like% "Sox",]
unchar[unchar$protein %like% "Sox",] #none 

#T-box
all[all$description.y %like% "T-box",]
#LOC576774, LOC592389, Tbx2/3
unchar[unchar$protein %like% "T-box",] #none

#creating larger list - 9 families
bhlh <- subset(all, description %like% "helix-loop-helix")
bhlh2 <- subset(unchar, (description_unchar %like% "bHLH" |
                           description_unchar %like% "Helix-loop-helix")) 
bhlh$TFfam <- "bHLH"
bhlh <- select(bhlh, c(names, TFfam))
bhlh2$TFfam <- "bHLH"
bhlh2 <- select(bhlh2, c(names, TFfam))
bhlhR <- rbind(bhlh, bhlh2)

bzip <- subset(all, description %like% "leucine zipper")
bzip2 <- subset(unchar, description_unchar %like% "leucine zipper" )
bzip$TFfam <- "bZIP"
bzip <- select(bzip, c(names, TFfam))
bzip2$TFfam <- "bZIP"
bzip2 <- select(bzip2, c(names, TFfam))
bzipR <- rbind(bzip, bzip2)

ets <- subset(all,description %like% "ETS" | 
                description %like% "ets" | 
                description %like% "Ets")
ets$TFfam <- "Ets"
etsR <- select(ets, c(names, TFfam))

fox <- subset(all, description %like% "Forkhead" |
                description %like% "forkhead" |
                description %like% "fox")
fox$TFfam <- "Fox"
foxR <- select(fox, c(names, TFfam))

gata <- subset(all, description %like% "GATA")
gata$TFfam <- "GATA"
gataR <- select(gata, c(names, TFfam))


hox <- subset(all, description %like% "homeobox" |
                description %like% "homeodomain")
hox2 <- subset(unchar, description_unchar %like% "Hox")
hox$TFfam <- "Hox"
hox <- select(hox, c(names, TFfam))
hox2$TFfam <- "Hox"
hox2 <- select(hox2, c(names, TFfam))
hoxR <- rbind(hox, hox2)

nr <- subset(all, description %like% "nuclear receptor")
nr2 <- subset(unchar, description_unchar %like% "nuclear receptor")
nr$TFfam <- "NR"
nr <- select(nr, c(names, TFfam))
nr2$TFfam <- "NR"
nr2 <- select(nr2, c(names, TFfam))
nrR <- rbind(nr, nr2)

sox <- subset(all, description %like% "Sox")
sox$TFfam <- "Sox"
soxR <- select(sox, c(names, TFfam))

tbox <- subset(all, description %like% "T-box")
tbox$TFfam <- "Tbox"
tboxR <- select(tbox, c(names, TFfam))

tfRef <- rbind(bhlhR, bzipR, etsR, foxR, gataR, hoxR, 
               nrR, soxR, tboxR)

tfAllI <- subset(allI, (description %like% "helix-loop-helix" | 
                        description %like% "leucine zipper" | 
                        description %like% "ETS" | 
                        description %like% "ets" | 
                        description %like% "Ets" |
                        description %like% "Forkhead" |
                        description %like% "forkhead" |
                        description %like% "fox" |
                        description %like% "GATA" |
                        description %like% "homeobox" |
                        description %like% "homeodomain" |
                        description %like% "nuclear receptor" |
                        description %like% "Sox" |
                        description %like% "T-box"))

tfcompI <- merge(tfI, tfAllI, by = "names", all = TRUE) 

tfUnchar <- subset(unchar, (description_unchar %like% "bHLH" |
                              description_unchar %like% "Helix-loop-helix" |
                              description_unchar %like% "leucine zipper" |
                              description_unchar %like% "Hox" |
                              description_unchar %like% "nuclear receptor"))
tfUI <- allI[allI$names %in% tfUnchar$names,]

tfcompI <- merge(tfcompI, tfUI, by = "names", all = TRUE)
tfcompI <- merge(tfcompI, tfRef, by = "names", all = TRUE) #450
length(unique(tfcompI$names)) #394

tfcomp <- tfcomp %>% 
  group_by(names, Class.L1, Class.L2, TFfam, 
           DE_WT.x, lfc_W, padj_W,
           lfc_temp, padj_temp, description.x) %>% 
  summarise(count = n()) 

tfcomp$Class.L2 <- gsub("TF_", "", tfcomp$Class.L2) 
tfcomp$Class.L2 <- gsub("Homeo", "Hox", tfcomp$Class.L2) 
tfcomp$Class.L2 <- gsub("Forkhead", "Fox", tfcomp$Class.L2) 
tfcomp$Class.L2 <- gsub("bzip", "bZIP", tfcomp$Class.L2) 
tfcomp$Class.L2 <- gsub("Sox/HMG", "Sox", tfcomp$Class.L2) 

tfcompI$TFfam <- ifelse(is.na(tfcompI$TFfam), tfcompI$Class.L2, tfcompI$TFfam)
tfnoznfI <- subset(tfcompI, !(TFfam == "ZNF")) #note: this set already ignored ZNFs, but this is pulling out all the ZNFs added via Tu's sheet
#write.csv(tfnoznf, "~/Desktop/R/RNAseq_MDII/July25-reanalysis/supp-info725/allTF_DEG_Fig4L_725.csv")

tfplot <- ggplot() + 
  theme_prism() + 
  #ylim(-5, 3) +
  #xlim(-3, 3) +
  geom_vline(xintercept = c(0, 0), col = "red") +
  geom_hline(yintercept = c(0, 0), col = "red") +
  geom_point(data = tfnoznfI,
             aes(x = lfc_I, y = lfc_WI, color = DE_I)) + 
  geom_rect(data = tfI, 
            aes(xmin = 0 - 0.5, 
                xmax = 0 + 0.5, 
                ymin = 0 - 0.5, 
                ymax = 0 + 0.5), 
            fill = "transparent", color = "purple", size = 0.2) +
  #geom_text_repel(data = subset(tfnoznfI, ((DE_WT.x == "Up-Up" | 
  #                                           DE_WT.x == "Down-Down" | 
  #                                           DE_WT.x == "Up-Down" |
  #                                           DE_WT.x == "Down-Up" |
  ##                                           abs(lfc_temp > 1.25) | 
  #                                           abs(lfc_W > 1)) &
  #                                          !(TFfam == "Other"))), 
  #                aes(x = lfc_W, y = lfc_temp, label = TFfam), 
  #                size = 4, color = "black", max.overlaps = Inf) + 
  xlab(expression(bold("F18"~"-"~"log"[2]~Fold~Change))) +
  ylab(expression(bold(WI~"-"~"log"[2]~Fold~Change))) +
  scale_color_manual("TF gene\nexpression", values = dewtColors, 
                     labels = 
                       c("MBE Down, Temp Down", 
                         "MBE Down", 
                         "MBE Down, Temp Up", 
                         "Temp NC", 
                         "Temp Down", 
                         "No change", 
                         "Temp Up", 
                         "MBE Up, Temp Down", 
                         "MBE Up", 
                         "MBE Up, Temp Up")) +
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
plot(tfplot)
ggsave(plot = tfplot, 
       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/tf-expression-Plot725.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

table(tfnoznf$TFfam)
#bHLH  bZIP   Ets   Fox  GATA   Hox    NR Other  Smad   Sox  Tbox 
#  47    31    14    26     5   103    44    56     8    10     5 

# TF with matching motifs ####
#note: only Sox appears in Tu ontology list
#Nfia (nucelar factor IA) - yes "nuclear factor", NFkB !(inhibitor)
#Rfx - yes "regulatory factor X"
#Mtf (metal-responsive) - yes "metal regulatory", "metal-response"
#MSANTD(Myb/sant dna bd containing) - yes "myb", "Myb" (also unchar-Myb)
#FIGLA (factor in germline alpha) - no 
#TCF - yes "Tcf"
#Klf - yes "Kruppel"
#Sox - yes "Sox" 
#Znf - purposefully not including: too many with unknown functions
all[all$description %like% "Tcf",]
unchar[unchar$description_unchar %like% "Myb",]

nfia <- subset(all, description %like% "nuclear factor")
nfkb <- subset(all, description %like% "NFkB") 
nfia <- rbind(nfia, nfkb)
nfia$TFfam <- "Nfia" 
nfia <- subset(nfia, !(description %like% "inhibitor"))

rfx <- subset(all, description %like% "regulatory factor X")
rfx$TFfam <- "Rfx"

mtf <- subset(all, description %like% "metal regulatory")
mtf2 <- subset(all, description %like% "metal-response")
mtf <- rbind(mtf, mtf2)
mtf$TFfam <- "Mtf"

msantd <- subset(all, description %like% "myb")
msantd2 <- subset(all, description %like% "Myb")
unchar[unchar$description_unchar %like% "Myb",]
msantd3 <- subset(all, names %like% "LOC100889035")
msantd <- rbind(msantd, msantd2, msantd3)
msantd$TFfam <- "MSANTD"

tcf <- subset(all, description %like% "Tcf")
tcf$TFfam <- "Tcf"

klf <- subset(all, description %like% "Kruppel")
klf$TFfam <- "Klf"

sox2 <- subset(all, description %like% "Sox")
sox2$TFfam <- "Sox"

tfplot2 <- rbind(nfia, rfx, mtf, msantd, tcf, klf, sox2)

colnames(tfplot2) <- c("names", "lfc_W", "padj_W", "DE_WT", 
                       "lfc_temp", 
                       "padj_temp", "description", "TFfam")

table(tfplot2$DE_WT)
#write.csv(tfplot2, "~/Desktop/R/RNAseq_MDII/July25-reanalysis/supp-info725/TFwithmotifs_DEG_Fig5D725.csv")

tfMotifPlot <- ggplot() + 
  theme_prism() + 
  #ylim(-3, 3) +
  #xlim(-3, 3) +
  geom_vline(xintercept = c(0, 0), col = "red") +
  geom_hline(yintercept = c(0, 0), col = "red") +
  geom_point(data = tfplot2,
             aes(x = lfc_W, y = lfc_temp, color = DE_WT)) + 
  geom_rect(data = tfplot2, 
            aes(xmin = 0 - 0.5, 
                xmax = 0 + 0.5, 
                ymin = 0 - 0.5, 
                ymax = 0 + 0.5), 
            fill = "transparent", color = "purple", size = 0.2) +
  geom_text_repel(data = subset(tfplot2, 
                                !(DE_WT == "No change-No change")), 
                  aes(x = lfc_W, y = lfc_temp, 
                      label = TFfam), 
                  size = 4, color = "black", max.overlaps = Inf) + 
  xlab(expression(bold("Microbial Richness"~"-"~"log"[2]~Fold~Change))) +
  ylab(expression(bold(Temp~"-"~"log"[2]~Fold~Change))) +
  scale_color_manual("TF gene\nexpression", values = dewtColors, 
                     labels = 
                       c("MBE Down, Temp Down", 
                         "MBE Down", 
                         "MBE Down, Temp Up", 
                         "Temp NC", 
                         "Temp Down", 
                         "No change", 
                         "Temp Up", 
                         "MBE Up, Temp Down", 
                         "MBE Up", 
                         "MBE Up, Temp Up")) +
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
plot(tfMotifPlot)

ggsave(plot = tfMotifPlot, 
       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/tf-Motif-Plot725.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

#znf quick check 
znf <- subset(all, description %like% "zinc finger")
table(znf$DE_WT) #223 no change; 112 DE
ggplot() + 
  theme_prism() + 
  ylim(-3, 3) +
  xlim(-3, 3) +
  geom_vline(xintercept = c(0, 0), col = "red") +
  geom_hline(yintercept = c(0, 0), col = "red") +
  geom_point(data = znf,
             aes(x = log2FoldChange.x, 
                 y = log2FoldChange.y, color = DE_WT)) + 
  geom_rect(data = tf, 
            aes(xmin = 0 - 0.5, 
                xmax = 0 + 0.5, 
                ymin = 0 - 0.5, 
                ymax = 0 + 0.5), 
            fill = "transparent", color = "purple", size = 0.2) +
  #geom_text_repel(data = znf, 
  #                aes(x = lfc_MBE, y = lfc_temp, label = names), 
  #                size = 2, color = "black", max.overlaps = Inf) + 
  xlab(expression(bold(MBE~"-"~"log"[2]~Fold~Change))) +
  ylab(expression(bold(Temp~"-"~"log"[2]~Fold~Change))) +
  scale_color_manual("TF gene\nexpression", values = dewtColors, 
                     labels = 
                       c("MBE Down, Temp Down", 
                         "MBE Down", 
                         "MBE Down, Temp Up", 
                         "Temp NC", 
                         "Temp Down", 
                         "No change", 
                         "Temp Up", 
                         "MBE Up, Temp Down", 
                         "MBE Up", 
                         "MBE Up, Temp Up")) +
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) +  
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))

# Metabolism - clean #####
met <- subset(allOnto, Class.L1 == "Metabolism")
#write.csv(met, "~/Desktop/R/RNAseq_MDII/July25-reanalysis/supp-info725/metabolismDEG_Fig5B725.csv")
table(met$DE_WT) 
# Down-Down      Down-No change             Down-Up 
#        21                 171                  40 
# No change-Down No change-No change        No change-Up 
#            199                1163                 325 
# Up-Down        Up-No change               Up-Up 
#      39                  79                   9 
#total: 2046 
#deg: 883 (-no change:no change)

metPlot <- ggplot() + 
  theme_prism() + 
  #ylim(-6, 4) +
  #xlim(-4, 4) +
  geom_vline(xintercept = c(0, 0), col = "red") +
  geom_hline(yintercept = c(0, 0), col = "red") +
  geom_point(data = met,
             aes(x = log2FoldChange.x, y = log2FoldChange.y, color = DE_WT)) + 
  geom_rect(data = met, 
            aes(xmin = 0 - 0.5, 
                xmax = 0 + 0.5, 
                ymin = 0 - 0.5, 
                ymax = 0 + 0.5), 
            fill = "transparent", color = "purple", size = 0.2) +
  geom_text_repel(data = subset(met, (DE_WT == "Up-Up" | 
                                        DE_WT == "Down-Down" | 
                                        DE_WT == "Up-Down" |
                                        DE_WT == "Down-Up") & 
                                  !(names %like% "LOC")), 
                  aes(x = log2FoldChange.x, 
                      y = log2FoldChange.y, label = names), 
                  size = 4, color = "black", max.overlaps = Inf) + 
  xlab(expression(bold("Microbial Richness"~"-"~"log"[2]~Fold~Change))) +
  ylab(expression(bold(Temp~"-"~"log"[2]~Fold~Change))) +
  scale_color_manual("TF gene\nexpression", values = dewtColors, 
                     labels = 
                       c("MBE Down, Temp Down", 
                         "MBE Down", 
                         "MBE Down, Temp Up", 
                         "Temp NC", 
                         "Temp Down", 
                         "No change", 
                         "Temp Up", 
                         "MBE Up, Temp Down", 
                         "MBE Up", 
                         "MBE Up, Temp Up")) +
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
plot(metPlot)
ggsave(plot = metPlot, 
       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/metabolism-Plot725.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

metlookup <- subset(met, !(names %like% "LOC" | DE_WT == "No change-No change"))
metlookup2 <- subset(met, !(DE_WT %like% "No change"))
metDU <- subset(metlookup, DE_WT == "Down-Up")
metUD <- subset(metlookup, DE_WT == "Up-Down")
metUU <- subset(metlookup, DE_WT == "Up-Up") #none
metDD <- subset(metlookup, DE_WT == "Down-Down") 


# Immunity - clean #### 
immune <- subset(allOnto, Class.L1 == "Immunity")
#write.csv(immune, "~/Desktop/R/RNAseq_MDII/July25-reanalysis/supp-info725/immuneAll_DEG_Fig6F725.csv")
table(immune$DE_WT)
# Down-Down      Down-No change             Down-Up 
#         3                  24                   8 
# No change-Down No change-No change        No change-Up 
#             64                 475                 142 
# Up-Down        Up-No change               Up-Up 
#      15                  36                   1 

sigimmune <- subset(immune, !(DE_WT == "No change-No change"))
table(sigimmune$DE_WT)
#293 DEG 
#Temp only = 64+142 = 206
#MBE only = 24+36 = 60
#both = 3+8+15+1 = 27


immunePlot <- ggplot() + 
  theme_prism() + 
  #ylim(-6, 4) +
  #xlim(-4, 4) +
  geom_vline(xintercept = c(0, 0), col = "red") +
  geom_hline(yintercept = c(0, 0), col = "red") +
  geom_point(data = immune,
             aes(x = log2FoldChange.x, 
                 y = log2FoldChange.y, color = DE_WT)) + 
  geom_rect(data = immune, 
            aes(xmin = 0 - 0.5, 
                xmax = 0 + 0.5, 
                ymin = 0 - 0.5, 
                ymax = 0 + 0.5), 
            fill = "transparent", color = "purple", size = 0.2) +
  geom_rect(data = immune, 
            aes(xmin = 0 - 2.5, 
                xmax = 0 + 1, 
                ymin = 0 - 2.5, 
                ymax = 0 + 1), 
            fill = "transparent", color = "#333333", size = 0.5, linetype = 2) +
  geom_text_repel(data = subset(immune, ((abs(log2FoldChange.x > 1)) |
                                           (abs(log2FoldChange.y > 1))) & 
                                  !(names %like% "LOC")), 
                  aes(x = log2FoldChange.x, 
                      y = log2FoldChange.y, label = names), 
                  size = 4, color = "black", max.overlaps = Inf) + 
  xlab(expression(bold("Microbial Richness"~"-"~"log"[2]~Fold~Change))) +
  ylab(expression(bold(Temp~"-"~"log"[2]~Fold~Change))) +
  scale_color_manual("TF gene\nexpression", values = dewtColors, 
                     labels = 
                       c("MBE Down, Temp Down", 
                         "MBE Down", 
                         "MBE Down, Temp Up", 
                         "Temp NC", 
                         "Temp Down", 
                         "No change", 
                         "Temp Up", 
                         "MBE Up, Temp Down", 
                         "MBE Up", 
                         "MBE Up, Temp Up")) +
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
plot(immunePlot)
ggsave(plot = immunePlot, 
       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/immune-Plot725.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

# pigment cells - clean #####

cand1= all[grepl("gcm", all$names), ] #gcm 
cand2= all[grepl("LOC588806", all$names), ] #pks1 
cand3= all[grepl("LOC582581", all$names), ] #fmo2 
cand4= all[grepl("LOC586240", all$names), ] #fmo3 
cand5= all[grepl("LOC589643", all$names), ] #fmo2-2 
cand6= all[grepl("LOC588371", all$names), ] #fmo5-1 
cand7= all[grepl("LOC592500", all$names), ] #sult1c2 
cand8= all[grepl("LOC100891288", all$names), ] #oprk1c/L 
cand9= all[grepl("LOC754218", all$names), ] #mif5 #In already
cand10= all[grepl("LOC579779", all$names), ] #glur6
cand11= all[grepl("LOC590101", all$names), ] #abcg11 

Cands=rbind(cand1, cand2, cand3, cand4, cand5, cand6, cand7, 
            cand8, cand9, cand10, cand11)

pigment <- subset(all, names %in% Cands$names)
colnames(pigment) <- c("names", "lfc_W", "padj_W", "DE_WT", 
                       "lfc_temp", 
                       "padj_temp", "description")

#adding common name - ignore these 2 lines
#immuneLookUp <- read.csv("immune_genes_presentMDRNA.csv")
#pigment <- merge(pigment, immuneLookUp, by = "names")

#write.csv(pigment, "~/Desktop/R/RNAseq_MDII/July25-reanalysis/supp-info725/pigmentCell_DEG_Fig6G725.csv")

pigmentPlot <- ggplot() + 
  theme_prism() + 
  ylim(-2.5, 1) +
  xlim(-2.5, 1) +
  geom_vline(xintercept = c(0, 0), col = "red") +
  geom_hline(yintercept = c(0, 0), col = "red") +
  geom_point(data = pigment,
             aes(x = lfc_W, y = lfc_temp, color = DE_WT)) + 
  geom_rect(data = pigment, 
            aes(xmin = 0 - 0.5, 
                xmax = 0 + 0.5, 
                ymin = 0 - 0.5, 
                ymax = 0 + 0.5), 
            fill = "transparent", color = "purple", size = 0.2) +
  geom_text_repel(data = pigment, 
                  aes(x = lfc_W, y = lfc_temp, label = names), 
                  size = 4, color = "black", max.overlaps = Inf) + 
  xlab(expression(bold("Microbial Richness"~"-"~"log"[2]~Fold~Change))) +
  ylab(expression(bold(Temp~"-"~"log"[2]~Fold~Change))) +
  scale_color_manual("Immune gene\nexpression", values = dewtColors, 
                     labels = 
                       c("MBE Down, Temp Down", 
                         "MBE Down", 
                         "MBE Down, Temp Up", 
                         "Temp NC", 
                         "Temp Down", 
                         "No change", 
                         "Temp Up", 
                         "MBE Up, Temp Down", 
                         "MBE Up", 
                         "MBE Up, Temp Up")) +
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
plot(pigmentPlot)
ggsave(plot = pigmentPlot, 
       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/supp-info725//pigmentCell-Plot725.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)
# Removed dGRN #####

# Defensome - ABC transporters #### 
abc <- subset(allOnto, Class.L2 == "Defensome_TransporterABC")
unique(abc$Class.L3)

abc[abc$sp_name %like% "Abc", ]

abcPlot <- ggplot() + 
  theme_prism() + 
  #ylim(-2.5, 1) +
  #xlim(-2, 1) +
  geom_vline(xintercept = c(0, 0), col = "red") +
  geom_hline(yintercept = c(0, 0), col = "red") +
  geom_point(data = abc,
             aes(x = lfc_W, y = lfc_temp, color = DE_WT)) + 
  geom_rect(data = abc, 
            aes(xmin = 0 - 0.5, 
                xmax = 0 + 0.5, 
                ymin = 0 - 0.5, 
                ymax = 0 + 0.5), 
            fill = "transparent", color = "purple", size = 0.2) +
  geom_text_repel(data = subset(abc, !(DE_WT %like% "No change")), 
                  aes(x = lfc_W, y = lfc_temp, label = sp_name), 
                  size = 4, color = "black", max.overlaps = Inf) + 
  xlab(expression(bold(MBE~"-"~"log"[2]~Fold~Change))) +
  ylab(expression(bold(Temp~"-"~"log"[2]~Fold~Change))) +
  scale_color_manual("Immune gene\nexpression", values = dewtColors, 
                     labels = 
                       c("MBE Down, Temp Down", 
                         "MBE Down", 
                         "MBE Down, Temp Up", 
                         "Temp NC", 
                         "Temp Down", 
                         "No change", 
                         "Temp Up", 
                         "MBE Up, Temp Down", 
                         "MBE Up", 
                         "MBE Up, Temp Up")) +
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
print(abcPlot)
ggsave(plot = abcPlot, 
       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/abc-transporter-Plot725.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

#note: haven't written this file
#write.csv(abc, "~/Desktop/abc-transporter-expression.csv")

