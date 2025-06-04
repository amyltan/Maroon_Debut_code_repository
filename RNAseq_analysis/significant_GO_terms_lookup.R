setwd("~/Desktop/R/MDII_GO_MWU/")

library('ggprism')
library('ggplot2')
library('ggrepel')
library('VennDiagram')
library('tidyr')
library('dplyr')
library('data.table')
library('stringr')
library('eulerr')

go <- read.csv("significant_MWU_GO_MDRNA.csv")
table(go$name)

dup <- go %>% 
  group_by(name) %>% 
  summarize(comparison = paste(comparison, collapse = ";"))
table(dup$comparison)

tonly <- subset(dup, comparison == "temp")
dim(tonly) #29
wonly <- subset(dup, comparison == "water")
dim(wonly) #26
tw <- subset(dup, comparison == "temp;water") 
dim(tw) #27

twfull <- subset(go, name %in% tw$name)
twfull[twfull$name %like% "translation",]
count(twfull[twfull$GO_category %like% "MF",]) #6
count(twfull[twfull$GO_category %like% "CC",]) #28
count(twfull[twfull$GO_category %like% "BP",]) #20

tfull <- subset(go, name %in% tonly$name)
count(tfull[tfull$GO_category %like% "MF",]) #17
count(tfull[tfull$GO_category %like% "CC",]) #3
count(tfull[tfull$GO_category %like% "BP",]) #9

wfull <- subset(go, name %in% wonly$name)
count(wfull[wfull$GO_category %like% "MF",]) #0
count(wfull[wfull$GO_category %like% "CC",]) #20
count(wfull[wfull$GO_category %like% "BP",]) #6

# SHARED ONLY - Graphing direction of GO regulation #####
library('ggh4x')

twcolor <- c("#FFCC00", "#33CCCC")
names(twcolor) <- c("temp", "water")

BPgo <- ggplot() + 
  theme_prism() + 
  geom_col(data = subset(twfull, GO_category == "BP"), 
           aes(x = name, y = delta.rank, fill = comparison)) +
  scale_fill_manual(values = twcolor, 
                    labels = c("Temp", "MBE")) +
  geom_hline(yintercept = c(0), col = "black") + 
  coord_flip() + 
  ylim(c(-475, 350)) +
  xlab(NULL) + 
  ylab("Delta rank") + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) +
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(2, "in")) #note, these are not flipped. Col is still horizontal axis.
print(BPgo) 

ggsave(plot = BPgo, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/BPgo.jpg", 
       scale = 1, width = 7.5, height = 6, units = c("in"), 
       dpi = 300)

CCgo <- ggplot() + 
  theme_prism() + 
  geom_col(data = subset(twfull, GO_category == "CC"), 
           aes(x = name, y = delta.rank, fill = comparison)) +
  scale_fill_manual("Condition", 
                    values = c("#FFCC00", "#33CCCC"), 
                    labels = c("Temp", "MBE")) +
  geom_hline(yintercept = c(0), col = "black") + 
  coord_flip() + 
  ylim(c(-475, 350)) +
  xlab(NULL) + 
  ylab("Delta rank") + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) +
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(2, "in")) #note, these are not flipped. Col is still horizontal axis.
print(CCgo)
ggsave(plot = CCgo, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/CCgo.jpg", 
       scale = 1, width = 7, height = 6, units = c("in"), 
       dpi = 300)

MFgo <- ggplot() + 
  theme_prism() + 
  geom_col(data = subset(twfull, GO_category == "MF"), 
           aes(x = name, y = delta.rank, fill = comparison)) +
  scale_fill_manual("Condition", 
                    values = c("#FFCC00", "#33CCCC"), 
                    labels = c("Temp", "MBE")) +
  geom_hline(yintercept = c(0), col = "black") + 
  coord_flip() + 
  ylim(c(-475, 350)) +
  xlab(NULL) + 
  ylab("Delta rank") + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) +
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(2, "in")) #note, these are not flipped. Col is still horizontal axis.
print(MFgo)
ggsave(plot = MFgo, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/MFgo.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

# SINGULAR ONLY - Graphing direction of GO regulation #####
tempBP <-
  ggplot() + 
  theme_prism() + 
  geom_col(data = subset(tfull, GO_category == "BP"), 
           aes(x = name, y = delta.rank, fill = comparison)) +
  scale_fill_manual(values = twcolor, 
                    labels = c("Temp")) +
  geom_hline(yintercept = c(0), col = "black") + 
  coord_flip() + 
  ylim(c(-475, 350)) +
  xlab(NULL) + 
  ylab("Delta rank") + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) +
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(2, "in")) #note, these are not flipped. Col is still horizontal axis.
print(tempBP) 

ggsave(plot = tempBP, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/tempBPgo-FigS3.jpg", 
       scale = 1, width = 8, height = 6, units = c("in"), 
       dpi = 300)

mbeBP <- ggplot() + 
  theme_prism() + 
  geom_col(data = subset(wfull, GO_category == "BP"), 
           aes(x = name, y = delta.rank, fill = comparison)) +
  scale_fill_manual(values = twcolor, 
                    labels = c("MBE")) +
  geom_hline(yintercept = c(0), col = "black") + 
  coord_flip() + 
  ylim(c(-475, 350)) +
  xlab(NULL) + 
  ylab("Delta rank") + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) +
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(2, "in")) #note, these are not flipped. Col is still horizontal axis.
print(mbeBP) 

ggsave(plot = mbeBP, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/mbeBPgo-FigS3.jpg", 
       scale = 1, width = 8, height = 6, units = c("in"), 
       dpi = 300)

tempCC <- ggplot() + 
  theme_prism() + 
  geom_col(data = subset(tfull, GO_category == "CC"), 
           aes(x = name, y = delta.rank, fill = comparison)) +
  scale_fill_manual(values = twcolor, 
                    labels = c("Temp")) +
  geom_hline(yintercept = c(0), col = "black") + 
  coord_flip() + 
  ylim(c(-475, 350)) +
  xlab(NULL) + 
  ylab("Delta rank") + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) +
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(2, "in")) #note, these are not flipped. Col is still horizontal axis.
print(tempCC) 

ggsave(plot = tempCC, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/tempCCgo-FigS3.jpg", 
       scale = 1, width = 8, height = 6, units = c("in"), 
       dpi = 300)


mbeCC <- ggplot() + 
  theme_prism() + 
  geom_col(data = subset(wfull, GO_category == "CC"), 
           aes(x = name, y = delta.rank, fill = comparison)) +
  scale_fill_manual(values = twcolor, 
                    labels = c("MBE")) +
  geom_hline(yintercept = c(0), col = "black") + 
  coord_flip() + 
  ylim(c(-475, 350)) +
  xlab(NULL) + 
  ylab("Delta rank") + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) +
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(2, "in")) #note, these are not flipped. Col is still horizontal axis.
print(mbeCC) 

ggsave(plot = mbeCC, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/mbeCCgo-FigS3.jpg", 
       scale = 1, width = 8, height = 6, units = c("in"), 
       dpi = 300)

tempMF <- ggplot() + 
  theme_prism() + 
  geom_col(data = subset(tfull, GO_category == "MF"), 
           aes(x = name, y = delta.rank, fill = comparison)) +
  scale_fill_manual(values = twcolor, 
                    labels = c("Temp")) +
  geom_hline(yintercept = c(0), col = "black") + 
  coord_flip() + 
  ylim(c(-475, 400)) +
  xlab(NULL) + 
  ylab("Delta rank") + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) +
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(2, "in")) #note, these are not flipped. Col is still horizontal axis.
print(tempMF) 

ggsave(plot = tempMF, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/tempMFgo-FigS3.jpg", 
       scale = 1, width = 8, height = 6, units = c("in"), 
       dpi = 300)

#reminder: MF full doesn't have it

# Venn ##### 
temp <- subset(go, comparison == "temp")
wat <- subset(go, comparison == "water")

allGO = list("Microbial\nexposure" = (wat$name), 
                     "Temperature" = (temp$name))

quartz()
vennGO = venn.diagram(
  x = allGO, 
  filename = NULL,
  col = "transparent", 
  fill = c("#00CCCC", "#FFCC33"),
  cex = 3, 
  fontfamily = "sans", 
  cat.default.pos = "text", 
  cat.col = c("#996633", "#003333"), 
  cat.cex = 0, #size to 0 to remove labels, also hashed out irrelevant lines
  #cat.fontfamily = "sans", 
  #cat.dist = c(0.15, 0.15), 
  #cat.pos = 1, 
  rotation.degree = 180
); 
grid.draw(vennGO) 

tBP <- subset(temp, GO_category == "BP")
wBP <- subset(wat, GO_category == "BP")
BP = list("Microbial\nexposure" = (wBP$name), 
             "Temperature" = (tBP$name))
quartz()
vennBP = venn.diagram(
  x = BP, 
  filename = NULL, 
  col = "transparent", 
  fill = c("#33CCCC","#FFCC00"),
  cex = 3, 
  fontfamily = "sans", 
  cat.default.pos = "text", 
  cat.col = c("#996633", "#003333"), 
  cat.cex = 0, #size to 0 to remove labels, also hashed out irrelevant lines
  #cat.fontfamily = "sans", 
  #cat.dist = c(0.15, 0.15), 
  #cat.pos = 1, 
); 
grid.draw(vennBP) 

tcc <- subset(temp, GO_category == "CC")
wcc <- subset(wat, GO_category == "CC")
CC = list("Microbial\nexposure" = (wcc$name), 
          "Temperature" = (tcc$name))
quartz()
vennCC = venn.diagram(
  x = CC, 
  filename = NULL, 
  col = "transparent", 
  fill = c("#33CCCC","#FFCC00"),
  cex = 3, 
  fontfamily = "sans", 
  cat.default.pos = "text", 
  cat.col = c("#996633", "#003333"), 
  cat.cex = 0, #size to 0 to remove labels, also hashed out irrelevant lines
  #cat.fontfamily = "sans", 
  #cat.dist = c(0.15, 0.15), 
  #cat.pos = 1,
  inverted = TRUE
); 
grid.draw(vennCC) 

tmf <- subset(temp, GO_category == "MF")
wmf <- subset(wat, GO_category == "MF")
MF = list("Microbial\nexposure" = (wmf$name), 
          "Temperature" = (tmf$name))
quartz()
vennMF = venn.diagram(
  x = MF, 
  filename = NULL, 
  col = "transparent", 
  fill = c("#33CCCC","#FFCC00"),
  cex = 3, 
  fontfamily = "sans", 
  cat.default.pos = "text", 
  cat.col = c("#996633", "#003333"), 
  cat.cex = 0, #size to 0 to remove labels, also hashed out irrelevant lines
  #cat.fontfamily = "sans", 
  #cat.dist = c(0.15, 0.15), 
  #cat.pos = 1,
); 
grid.draw(vennMF) 

gridExtra::grid.arrange(vennGO, vennBP, vennCC, vennMF)
#this prints all four venns but not scaled to each other.


# List of GO terms to look at genes #####      
go[go$name %like% "name",] #

#translation
go[go$name %like% "translation",] #GO:0006412
go[go$name %like% "protein-DNA complex assembly",] #GO:0006334;GO:0034728;GO:0065004;GO:0006338
go[go$name %like% "ribosome",] #GO:0005840
go[go$name %like% "structural constituent of ribosome",] #GO:0003735
go[go$name %like% "ribosomal subunit",] #GO:0044391 
go[go$name %like% "ribonucleoprotein complex",] #GO:1990904
go[go$name %like% "protein-DNA complex",] #GO:0000786;GO:0032993
go[go$name %like% "RNA binding",] #GO:0003723 

#nuc acid met process 
go[go$name %like% "nucleic acid metabolic process",] #GO:0090304
go[go$name %like% "DNA metabolic process",] #GO:0006259 
go[go$name %like% "RNA metabolic process",] #GO:0016070

#chromosomes
go[go$name %like% "chromosome organization",] #GO:0051276
go[go$name %like% "chromosomal region",] #GO:0098687
go[go$name %like% "structural constituent of chromatin",] #GO:0030527;GO:0046982


# Combining lists of GO for easy searching ##### 
tempBP <- read.csv("BP_temp_waldStat_MDRNA.csv", sep = "\t")
tempBP$GO_category <- "BP"
tempCC <- read.csv("CC_temp_waldStat_MDRNA.csv", sep = "\t")
tempCC$GO_category <- "CC"
tempMF <- read.csv("MF_temp_waldStat_MDRNA.csv", sep = "\t")
tempMF$GO_category <- "MF"

tempGO <- rbind(tempBP, tempCC, tempMF)
head(tempGO)
tempGO$comparison <- "temp"

watBP <- read.csv("BP_water_waldStat_MDRNA.csv", sep = "\t")
watBP$GO_category <- "BP"
watCC <- read.csv("CC_water_waldStat_MDRNA.csv", sep = "\t")
watCC$GO_category <- "CC"
watMF <- read.csv("MF_water_waldStat_MDRNA.csv", sep = "\t")
watMF$GO_category <- "MF"

waterGO <- rbind(watBP, watCC, watMF)
waterGO$comparison <- "water"

#creating list of the higher level GO terms that match the gene names ##### 
combo <- rbind(tempGO, waterGO)
combo <- subset(combo, select = -c(name, lev, value, GO_category, comparison))
single <- combo[!duplicated(combo),]

ref <- single %>% 
  group_by(seq) %>% 
  summarize(term = paste(term, collapse = ";"))

#double checking the frame condensed properly
#combo[combo$seq %like% "LOC100888908",]
#single[single$seq %like% "LOC100888908",]
#ref[ref$seq %like% "LOC100888908",]

#write.csv(ref, "GO_highLevel_MDII.csv")

# Pulling genes matching GO terms #####
RW <- read.csv("~/Desktop/R/RNAseq_MDII/microbe-exposureDEgenes-all.csv") 
RW <- subset(RW, select = -c(X, baseMean, lfcSE, stat, pvalue)) 
colnames(RW)[colnames(RW) == "names"] <- "v5name" #same as "names" in other files, but called v5name here to distinguish from GO column called "name"
RT <- read.csv("~/Desktop/R/RNAseq_MDII/temp-DEgenes-all.csv")
RT <- subset(RT, select = -c(X, baseMean, lfcSE, stat, pvalue)) 
colnames(RT)[colnames(RT) == "names"] <- "v5name"
#conv <- read.csv("~/Desktop/TAMU/Spurp_reference_files/Spurp_genome/v5conversions-LocSpuDescripGO.csv")

colnames(ref)[colnames(ref) == "seq"] <- "v5name" 
rna <- merge(RW, RT, by = "v5name", all = TRUE) ###NOTE need to re run these now that "all = TRUE is added"
colnames(rna) <- c("v5name", "lfcMBE", "padjMBE", "DE-MBE", 
                   "descriptionMBE", "lfcTemp", "padjTemp", 
                   "DE-Temp", "descriptionTemp") 
rna <- unite(rna, DE_WT, "DE-MBE", "DE-Temp", sep = "-")

# lists of GO terms to match #####

#transcription GO terms
glist1 <- ref %>% filter(
  str_detect(term, "GO:0006412|GO:0006334|GO:0034728|GO:0065004|GO:0006338|GO:0005840|GO:0003735|GO:0044391|GO:1990904|GO:0000786|GO:0032993|GO:0003723"))

# verification 
#ref[ref$term %like% "0006338",]
#glist1[glist1$term %like% "0006338",]

#nucleic acid metabolic processes
glist2 <- ref %>% filter(
  str_detect(term, "GO:0090304|GO:0006259|GO:0016070"))

#chromosomes
glist3 <- ref %>% filter(
  str_detect(term, "GO:0051276|GO:0098687|GO:0030527|GO:0046982"))

# scatterplots #### 
dewtColors <- c("#CC0033", "#33CCCC", "#CC0033", 
                "azure4", "#FFCC00", "azure4", "#FFCC00",
                "#CC0033", "#33CCCC", "#CC0033")
names(dewtColors) <- c("Down-Down", "Down-No change", "Down-Up", 
                       "NA-No change", "No change-Down", 
                       "No change-No change", "No change-Up", 
                       "Up-Down", "Up-No change", "Up-Up")

abs(min(rna$lfcMBE)) #3.464577 
abs(max(rna$lfcMBE)) #5.008481
abs(min(rna$lfcTemp)) #7.331406 
abs(max(rna$lfcTemp)) #4.133199

#translation related
trGO <- subset(rna, v5name %in% glist1$v5name)
#write.csv(trGO, "~/Desktop/Maroon_Debut/supp_info/transcriptranslaGO_DEG_Fig4K.csv")
translationPlot <- ggplot() + 
  theme_prism() + 
  ylim(-8, 5) +
  xlim(-8, 5) +
  geom_vline(xintercept = c(0, 0), col = "red") +
  geom_hline(yintercept = c(0, 0), col = "red") +
  geom_point(data = trGO,
             aes(x = lfcMBE, y = lfcTemp, color = DE_WT)) + 
  geom_text_repel(data = subset(trGO, !(v5name %like% "LOC" | 
                                          DE_WT == "No change-No change")), 
                  aes(x = lfcMBE, y = lfcTemp, label = v5name), 
                  size = 4, color = "black", max.overlaps = Inf) +
  geom_rect(data = rna, 
            aes(xmin = 0 - 0.5, 
                xmax = 0 + 0.5, 
                ymin = 0 - 0.5, 
                ymax = 0 + 0.5), 
            fill = "transparent", color = "purple", size = 0.2) +
  xlab(expression(bold(MBE~"-"~"log"[2]~Fold~Change))) +
  ylab(expression(bold(Temp~"-"~"log"[2]~Fold~Change))) +
  scale_color_manual("transcription gene\nexpression", values = dewtColors, 
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
plot(translationPlot)
ggsave(plot = translationPlot, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/GO-translation-expression-Plot.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)


namGO <- subset(rna, v5name %in% glist2$v5name)
#write.csv(namGO, "~/Desktop/Maroon_Debut/supp_info/nucacidmetGO_DEG_Fig4J.csv")
nametPlot <- ggplot() + 
  theme_prism() + 
  ylim(-8, 5) +
  xlim(-8, 5) +
  geom_vline(xintercept = c(0, 0), col = "red") +
  geom_hline(yintercept = c(0, 0), col = "red") +
  geom_point(data = namGO,
             aes(x = lfcMBE, y = lfcTemp, color = DE_WT)) +
  geom_text_repel(data = subset(namGO, !(v5name %like% "LOC")), 
                  aes(x = lfcMBE, y = lfcTemp, label = v5name), 
                  size = 4, color = "black", max.overlaps = Inf) +
  geom_rect(data = rna, 
            aes(xmin = 0 - 0.5, 
                xmax = 0 + 0.5, 
                ymin = 0 - 0.5, 
                ymax = 0 + 0.5), 
            fill = "transparent", color = "purple", size = 0.2) +
  xlab(expression(bold(MBE~"-"~"log"[2]~Fold~Change))) +
  ylab(expression(bold(Temp~"-"~"log"[2]~Fold~Change))) +
  scale_color_manual("transcription gene\nexpression", values = dewtColors, 
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
plot(nametPlot)
ggsave(plot = nametPlot, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/GO-nucacidmet-expression-Plot.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

chrGO <- subset(rna, v5name %in% glist3$v5name)
#write.csv(chrGO, "~/Desktop/Maroon_Debut/supp_info/chromatinGO_DEG_Fig4I.csv")
chrPlot <- ggplot() + 
  theme_prism() + 
  ylim(-8, 5) +
  xlim(-8, 5) +
  geom_vline(xintercept = c(0, 0), col = "red") +
  geom_hline(yintercept = c(0, 0), col = "red") +
  geom_point(data = chrGO,
             aes(x = lfcMBE, y = lfcTemp, color = DE_WT)) + 
  geom_text_repel(data = subset(chrGO, !(v5name %like% "LOC")), 
                  aes(x = lfcMBE, y = lfcTemp, label = v5name), 
                  size = 4, color = "black", max.overlaps = Inf) +
  geom_rect(data = rna, 
            aes(xmin = 0 - 0.5, 
                xmax = 0 + 0.5, 
                ymin = 0 - 0.5, 
                ymax = 0 + 0.5), 
            fill = "transparent", color = "purple", size = 0.2) +
  xlab(expression(bold(MBE~"-"~"log"[2]~Fold~Change))) +
  ylab(expression(bold(Temp~"-"~"log"[2]~Fold~Change))) +
  scale_color_manual("transcription gene\nexpression", values = dewtColors, 
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
plot(chrPlot)
ggsave(plot = chrPlot, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/GO-chromosome-expression-Plot.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)
