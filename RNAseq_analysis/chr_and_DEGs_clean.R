setwd("~/Desktop/R/RNAseq_MDII/July25-reanalysis/")

library('dplyr')
library('data.table')
library('ggplot2')
library('ggprism')
library('geomtextpath')
library('ggh4x')


atac <- read.csv("~/Desktop/R/atac_seq_MDII/DAinDB-ChIPseeker/normReads_withAnnotations_MDATAC.csv")
#69932
atac <- subset(atac, select = c(annotation, geneId, distanceToTSS,
                                feature)) 
colnames(atac) <- c("full_annotation", "names", "distanceToTSS", "feature_code")

RW <- read.csv("HMRvsLMR_DEG_725.csv") 
colnames(RW)[colnames(RW) == "X"] <- "names"
RW$DE <- "No change"
RW$DE[RW$log2FoldChange > 0.5 & RW$padj < 0.05] <- "Up"
RW$DE[RW$log2FoldChange < -0.5 & RW$padj < 0.05] <- "Down"
head(RW)

RT <- read.csv("18vs14_DEG_725.csv")
colnames(RT)[colnames(RT) == "X"] <- "names"
RT$DE <- "No change"
RT$DE[RT$log2FoldChange > 0.5 & RT$padj < 0.05] <- "Up"
RT$DE[RT$log2FoldChange < -0.5 & RT$padj < 0.05] <- "Down"
head(RT)

rldpvals <- read.csv("MDII_rldpvals_July25.csv") #18115
colnames(rldpvals)[colnames(rldpvals) == "X"] <- "names"

# number of genes with multiple peaks <1kb to TSS ----
P <- subset(atac, full_annotation == "Promoter(<=1kb)")
dim(P)
countP <- P %>% 
  group_by(names) %>%
  summarize(n = n())
countP <- subset(countP, n > 1)
max(countP$n) #6 
mean(countP$n) #2.2
median(countP$n) #2
length(unique(P$names)) #11826
(1714/11826) *100 #14.5%

## condensing list of ATAC regions 
patac <- P %>% 
  group_by(names, full_annotation) %>% 
  slice(which.min(abs(distanceToTSS)))

#sanity check 
names(which(table(P$names) > 4))
P[P$names %like% "LOC588892",] 
patac[patac$names %like% "LOC588892",]

# Open vs Closed Promoter by condition####
#MBE genes 
length(intersect(RW$names, patac$names)) #9599

mbeP <- merge(RW, patac, by = "names", all = TRUE)
dim(mbeP) #20342
mbeC <- subset(mbeP, is.na(full_annotation))
dim(mbeC) #8516
mbeC$chr_state <- "Closed"
length(intersect(mbeC$names, patac$names)) #0 - as expected
min(mbeC$log2FoldChange) #-4.314407 
max(mbeC$log2FoldChange) #6.659325
mbeO <- subset(mbeP, !(is.na(full_annotation) | is.na(log2FoldChange))) 
dim(mbeO) #9599 
mbeO$chr_state <- "Open"

de_colors <- c("blue", "coral2", "azure4")
names(de_colors) <- c("Down", "Up", "No change")

lfcTSSmbe <- ggplot() + 
  theme_prism() + 
  geom_point(data = subset(mbeO, !(DE == "No change")), 
             aes(x = distanceToTSS, y = abs(log2FoldChange), 
                 color = DE), 
             position = "jitter") + 
  scale_color_manual(values = de_colors) + 
  xlab(expression(bold("Distance from ATAC peak to TSS"))) + 
  ylab(expression(bold("Absolute log"[2]~fold~change~of~DEG))) +
  ylim(0, 8) + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) +  
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))
print(lfcTSSmbe)

ggsave(plot = lfcTSSmbe, 
       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/lfc-TSS-MBE725.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)


levels = c("Open", "Closed")

lfcPromMBE <- ggplot() + 
  theme_prism() + 
  geom_point(data = subset(mbeO, !(DE == "No change")), 
             aes(x = chr_state, y = abs(log2FoldChange),
                 color = DE), 
             position = "jitter") + 
  geom_boxplot(data = subset(mbeO, !(DE == "No change")), 
               aes(x = chr_state, y = abs(log2FoldChange)), 
               show.legend = FALSE, width = 0.5, 
               position = position_dodge(width = 1))+ 
  geom_point(data = subset(mbeC, !(DE == "No change")), 
             aes(x = chr_state, y = abs(log2FoldChange),
                 color = DE), 
             position = "jitter") + 
  geom_boxplot(data = subset(mbeC, !(DE == "No change")), 
               aes(x = chr_state, y = abs(log2FoldChange)), 
               show.legend = FALSE, width = 0.5, 
               position = position_dodge(width = 1))+ 
  scale_color_manual(values = de_colors) + 
  xlab(expression(bold("Promoter accessibility"))) + 
  scale_x_discrete(limits = levels) +
  ylab(expression(bold("Absolute log"[2]~fold~change~of~DEG))) +
  ylim(0, 8) + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) +  
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))
print(lfcPromMBE)

ggsave(plot = lfcPromMBE, 
       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/lfc-Promoters-MBE725.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

#Temp genes
length(intersect(RT$names, patac$names)) #9599

tP <- merge(RT, patac, by = "names", all = TRUE)
dim(tP) #20342
tC <- subset(tP, is.na(full_annotation))
dim(tC) #8516
tC$chr_state <- "Closed"
length(intersect(tC$names, patac$names)) #0 - as expected
min(tC$log2FoldChange) #-7.021326 
max(tC$log2FoldChange) #4.496332
tO <- subset(tP, !(is.na(full_annotation) | is.na(log2FoldChange))) 
dim(tO) #9599 
tO$chr_state <- "Open"

#de_colors <- c("blue", "coral2", "azure4")
#names(de_colors) <- c("Down", "Up", "No change")

lfcTSStemp <- ggplot() + 
  theme_prism() + 
  geom_point(data = subset(tO, !(DE == "No change")), 
             aes(x = distanceToTSS, y = abs(log2FoldChange), 
                 color = DE), 
             position = "jitter") + 
  scale_color_manual(values = de_colors) + 
  xlab(expression(bold("Distance from ATAC peak to TSS"))) + 
  ylab(expression(bold("Absolute log"[2]~fold~change~of~DEG))) +
  ylim(0, 8) + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) +  
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))
print(lfcTSStemp)

ggsave(plot = lfcTSStemp, 
       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/lfc-TSS-Temp725.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

levels = c("Open", "Closed")

lfcPromTemp <- ggplot() + 
  theme_prism() + 
  geom_point(data = subset(tO, !(DE == "No change")), 
             aes(x = chr_state, y = abs(log2FoldChange),
                 color = DE), 
             position = "jitter") + 
  geom_boxplot(data = subset(tO, !(DE == "No change")), 
               aes(x = chr_state, y = abs(log2FoldChange)), 
               show.legend = FALSE, width = 0.5, 
               position = position_dodge(width = 1))+ 
  geom_point(data = subset(tC, !(DE == "No change")), 
             aes(x = chr_state, y = abs(log2FoldChange),
                 color = DE), 
             position = "jitter") + 
  geom_boxplot(data = subset(tC, !(DE == "No change")), 
               aes(x = chr_state, y = abs(log2FoldChange)), 
               show.legend = FALSE, width = 0.5, 
               position = position_dodge(width = 1))+ 
  scale_color_manual(values = de_colors) + 
  xlab(expression(bold("Promoter accessibility"))) + 
  scale_x_discrete(limits = levels) +
  ylab(expression(bold("Absolute log"[2]~fold~change~of~DEG))) +
  ylim(0, 8) + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) +  
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))
print(lfcPromTemp)

ggsave(plot = lfcPromTemp, 
       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/lfc-Promoter-Temp725.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)


# Base mean of all DEGs ####
# Base mean from DESeq calculation
rna <- merge(RW, RT, by = "names", all = TRUE)
length(intersect(rldpvals$names, rna$names)) 
missing <- subset(rldpvals, !(names %in% rna$names)) #NONE
colnames(rna) <- c("names", "baseMean.MBE", "log2FoldChange.MBE", 
                   "lfcSE.MBE", "stat.MBE", "pvalue.MBE", "padj.MBE", 
                   "DE.MBE", "baseMean.Temp", 
                   "log2FoldChange.Temp", "lfcSE.Temp", "stat.Temp", 
                   "pvalue.Temp", "padj.Temp", "DE.Temp")

# Note: there were steps checking on different length dataframes but that's irrelevant now

rna$log10baseMean <- log10(rna$baseMean.Temp)

rna$DEsum <- "NC"
rna$DEsum[rna$DE.MBE == "No change" & !(rna$DE.Temp == "No change")] <- "Temp-only"
rna$DEsum[!(rna$DE.MBE == "No change") & rna$DE.Temp == "No change"] <- "MBE-only"
rna$DEsum[!(rna$DE.MBE == "No change") & !(rna$DE.Temp == "No change")] <- "DE-both"

dewtColors <- c("#CC0033", "#33CCCC", "azure4", "#FFCC00")
names(dewtColors) <- c("DE-both", "Temp-only", "NC", "MBE-only")

length(intersect(rna$names, patac$names)) #9598

rnaP <- merge(rna, patac, by = "names", all = TRUE)
dim(rnaP) #20342
rnaC <- subset(rnaP, is.na(full_annotation))
dim(rnaC) #8516
rnaC$chr_state <- "Closed"
length(intersect(rnaC$names, patac$names)) #0 - as expected
rnaO <- subset(rnaP, !(is.na(full_annotation) | is.na(baseMean.Temp))) 
dim(rnaO) #9599 
rnaO$chr_state <- "Open"

basemeanTSS <- ggplot() + 
  theme_prism() + 
  geom_point(data = subset(rnaO, !(DEsum == "NC")), 
             aes(x = distanceToTSS, y = log10baseMean), 
             position = "jitter", color = "#666666") + 
  xlab(expression(bold("Distance from ATAC peak to TSS"))) + 
  ylab(expression(bold("log"[10]~base~mean~expression))) + 
  ylim(0, 8) + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) +  
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in")) 
print(basemeanTSS)

ggsave(plot = basemeanTSS, 
       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/baseMeanTSS725.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

levels = c("Open", "Closed")

basemeanProm <- ggplot() + 
  theme_prism() + 
  geom_point(data = subset(rnaO, !(DEsum == "NC")), 
             aes(x = chr_state, y = log10baseMean), 
             position = "jitter", color = "#666666") + 
  geom_boxplot(data = subset(rnaO, !(DEsum == "NC")), 
               aes(x = chr_state, y = log10baseMean), 
               show.legend = FALSE)+ 
  geom_point(data = subset(rnaC, !(DEsum == "NC")), 
             aes(x = chr_state, y = log10baseMean), 
             position = "jitter", color = "#666666") + 
  geom_boxplot(data = subset(rnaC, !(DEsum == "NC")), 
               aes(x = chr_state, y = log10baseMean), 
               show.legend = FALSE)+ 
  #scale_color_manual(values = dewtColors) + 
  xlab(expression(bold("Promoter accessibility"))) + 
  scale_x_discrete(limits = levels) +
  ylab(expression(bold("log"[10]~base~mean~expression))) +
  ylim(0, 8) + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) +  
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))
print(basemeanProm)

ggsave(plot = basemeanProm, 
       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/baseMeanPromoter725.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

#are there any RNA genes NOT in ATAC set? ####
RNAonly <- subset(rna, !(names %in% atac$names))
dim(RNAonly)
dim(rna)
(5552/18115) *100 #30.6%

closedDEG <- ggplot() + 
  theme_prism() + 
  geom_point(data = subset(RNAonly, !(is.na(DE.MBE) | DE.MBE == "No change")), 
             aes(x = "HMR vs. LMR", y = abs(log2FoldChange.MBE), color = DE.MBE), 
             position = "jitter") + 
  geom_point(data = subset(RNAonly, !(DE.MBE == "Up" | DE.MBE == "Down")), 
             aes(x = "HMR vs. LMR ", y = abs(log2FoldChange.MBE), color = DE.MBE), 
             position = "jitter") + 
  geom_point(data = subset(RNAonly, !(is.na(DE.Temp) | DE.Temp == "No change")), 
             aes(x = "18ºC vs. 14ºC", y = abs(log2FoldChange.Temp), color = DE.Temp), 
             position = "jitter") + 
  geom_point(data = subset(RNAonly, !(DE.Temp == "Up" | DE.Temp == "Down")), 
             aes(x = "18ºC vs. 14ºC ", y = abs(log2FoldChange.Temp), color = DE.Temp), 
             position = "jitter") + 
  scale_color_manual(values = de_colors) + 
  xlab(expression("DEGs without ATAC peaks")) + 
  ylab(expression("Abs. log"[2]~Fold~Change)) +
  ylim(0, 8) + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) +  
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))
print(closedDEG)

ggsave(plot = closedDEG, 
       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/closedDEG725.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)


#Note: oldest chr file also has code for calculating nearest open chr region (for closed chr) against LFC

# Trying stats ##### 
#with LFC and all values
t.test(mbeO$log2FoldChange, mbeC$log2FoldChange) # sig - but has "no change"
#  mean of x   mean of y 
#-0.22004809  0.05702788 
#change: 3.86x (based on abs. LFC)

t.test(tO$log2FoldChange, tC$log2FoldChange) # sig - but has "no change"
# mean of x  mean of y 
#0.15841823 0.03015542 
#change: 5x

mbeO$absLFC <- abs(mbeO$log2FoldChange)
mbeODE <- subset(mbeO, !(DE == "No change"))
mbeC$absLFC <- abs(mbeC$log2FoldChange)
mbeCDE <- subset(mbeC, !(DE == "No change"))

t.test(mbeODE$absLFC, mbeCDE$absLFC) #p-value < 2.2e-16
#mean of x mean of y 
#1.557115  1.944611 
#change: 1.24x
t.test(mbeODE$log2FoldChange, mbeCDE$log2FoldChange) #p-value < 2.2e-16
#  mean of x   mean of y 
#-0.04222982  0.79881076 

tO$absLFC <- abs(tO$log2FoldChange)
tODE <- subset(tO, !(DE == "No change"))
tC$absLFC <- abs(tC$log2FoldChange)
tCDE <- subset(tC, !(DE == "No change"))

t.test(tODE$absLFC, tCDE$absLFC) #p-value < 2.2e-16
#mean of x mean of y 
#1.177719  1.387962 
#change: 1.2x
t.test(tODE$log2FoldChange, tCDE$log2FoldChange)
#  mean of x   mean of y 
#0.26529639 -0.02401869 

rnaO$absBM <- abs(rnaO$log10baseMean)
rnaODE <- subset(rnaO, !(DEsum == "NC"))
rnaC$absBM <- abs(rnaC$log10baseMean)
rnaCDE <- subset(rnaC, !(DEsum == "NC"))

t.test(rnaODE$absBM, rnaCDE$absBM) #p-value < 2.2e-16
#mean of x mean of y 
#1.941663  1.624440 
#change: 1.2x

RWonly$absLFC <- abs(RWonly$log2FoldChange)
RWonlyDE <- subset(RWonly, !(DE == "No change"))
RTonly$absLFC <- abs(RTonly$log2FoldChange)
RTonlyDE <- subset(RTonly, !(DE == "No change"))

t.test(RWonlyDE$absLFC, RTonlyDE$absLFC) #p-value < 2.2e-16
#mean of x mean of y 
#1.967584  1.400852  
#change: 1.4x 

#creating df for anova PROBABLY REMOVE everything below - I don't think I used this ----
RWonly$chr_state <- "Closed"
RTonly$chr_state <- "Closed"

RWdf <- subset(RWonly, select = c("names", "baseMean", "log2FoldChange", 
                                "DE", "absLFC", "chr_state"))
RWdf$condition <- "rnaOnly"
RTdf <- subset(RTonly, select = c("names", "baseMean", "log2FoldChange", 
                                  "DE", "absLFC", "chr_state"))
RTdf$condition <- "rnaOnly"
mbeOdf <- subset(mbeO, select = c("names", "baseMean", "log2FoldChange", 
                                  "DE", "absLFC", "chr_state"))
mbeOdf$condition <- "MBE_DEG"
mbeCdf <- subset(mbeC, select = c("names", "baseMean", "log2FoldChange", 
                                  "DE", "absLFC", "chr_state"))
mbeCdf$condition <- "MBE_DEG"
tOdf <- subset(tO, select = c("names", "baseMean", "log2FoldChange", 
                              "DE", "absLFC", "chr_state"))
tOdf$condition <- "Temp_DEG"
tCdf <- subset(tC, select = c("names", "baseMean", "log2FoldChange", 
                              "DE", "absLFC", "chr_state"))
tCdf$condition <- "Temp_DEG"

df <- rbind(RWdf, RTdf, mbeOdf, mbeCdf, tOdf, tCdf) 

aov <- aov(data = df, absLFC ~ chr_state * condition * DE)
summary(aov)
TukeyHSD(aov)
#                                  p adj
#Open:MBE_DEG-Closed:MBE_DEG    0.0000000
#Closed:rnaOnly-Closed:MBE_DEG  0.0000023
#Closed:Temp_DEG-Closed:MBE_DEG 0.0000000
#Open:Temp_DEG-Closed:MBE_DEG   0.9837123
#Closed:rnaOnly-Open:MBE_DEG    0.0000000
#Closed:Temp_DEG-Open:MBE_DEG   0.0000000
#Open:Temp_DEG-Open:MBE_DEG     0.0000000
#Closed:Temp_DEG-Closed:rnaOnly 0.1598535
#Open:Temp_DEG-Closed:rnaOnly   0.0000000
#Open:Temp_DEG-Closed:Temp_DEG  0.0000000

#Open:rnaOnly-Closed:MBE_DEG           NA
#Open:rnaOnly-Open:MBE_DEG             NA
#Open:rnaOnly-Closed:rnaOnly           NA
#Closed:Temp_DEG-Open:rnaOnly          NA
#Open:Temp_DEG-Open:rnaOnly            NA

#chr_state:condition : all sig
#chr_state:DE only ns Closed:Up-Open:Down             0.9978853
#3 way: ns are:
#Open:MBE_DEG:Up-Closed:MBE_DEG:Down                1.0000000 
#Open:Temp_DEG:Up-Closed:MBE_DEG:Down               0.4243459
#Closed:MBE_DEG:Up-Open:Temp_DEG:Down               1.0000000
#Closed:rnaOnly:No change-Closed:MBE_DEG:No change  0.9450885
#Open:Temp_DEG:No change-Open:MBE_DEG:No change     0.9999130
#Closed:Temp_DEG:No change-Closed:rnaOnly:No change 0.1982241
#Open:Temp_DEG:Up-Open:MBE_DEG:Up                   0.0731161 
#Closed:Temp_DEG:Up-Closed:rnaOnly:Up               1.0000000

df2 <- subset(df, !(DE == "No change"))
aov2 <- aov(data = df2, absLFC ~ chr_state * condition * DE)
summary(aov2)
TukeyHSD(aov2)

# trying plotting everything together 
df2$conchr <- paste(df2$condition, df2$chr_state, sep = "-")
df2$conDE <- paste(df2$condition, df2$DE, sep = "-")
df2$all <- paste(df2$conchr, df2$DE, sep = "-")

chrde <- c("#33CCCC", "azure4", "#FFCC00")
names(chrde) <- c("MBE_DEG", "rnaOnly", "Temp_DEG")

ggplot() + 
  theme_prism() + 
  geom_point(data = df2, 
             aes(x = all, y = abs(log2FoldChange), 
                 color = condition), 
             position = "jitter") + 
  geom_boxplot(data = df2, 
               aes(x = all, y = abs(log2FoldChange))) + 
  scale_color_manual(values = chrde) + 
  theme(axis.text.x = element_text(angle = 25))

ggplot() + 
  theme_prism() + 
  geom_point(data = df2, 
             aes(x = all, y = abs(log2FoldChange)), 
             position = "jitter", color = "#666666") + 
  geom_boxplot(data = df2, 
               aes(x = all, y = abs(log2FoldChange), 
                   fill = DE), 
               width = 0.5) + 
  scale_fill_manual(values = de_colors) + 
  theme(axis.text.x = element_text(angle = 25)) + 
  scale_x_discrete(labels = c("Closed", 
                              "Closed", 
                              "Open", 
                              "Open", 
                              "No open chr", 
                              "No open chr", 
                              "Closed", 
                              "Closed", 
                              "Open", 
                              "Open")) + 
  xlab(NULL)

#                                                 upr     p adj
#Open:MBE_DEG:Down-Closed:MBE_DEG:Down    -0.021776938 0.0067453
#Closed:rnaOnly:Down-Closed:MBE_DEG:Down   0.540454592 0.0000000
#Closed:Temp_DEG:Down-Closed:MBE_DEG:Down  0.626723445 0.0000000
#Open:Temp_DEG:Down-Closed:MBE_DEG:Down    0.381348807 0.0000000
#Closed:MBE_DEG:Up-Closed:MBE_DEG:Down     0.395636737 0.0000000
#Closed:rnaOnly:Up-Closed:MBE_DEG:Down     0.240983079 0.0132396
#Closed:Temp_DEG:Up-Closed:MBE_DEG:Down    0.234524897 0.0130813
#Closed:rnaOnly:Down-Open:MBE_DEG:Down     0.659355464 0.0000000
#Closed:Temp_DEG:Down-Open:MBE_DEG:Down    0.744769902 0.0000000
#Open:Temp_DEG:Down-Open:MBE_DEG:Down      0.499402523 0.0000000
#Closed:MBE_DEG:Up-Open:MBE_DEG:Down       0.519410735 0.0000000
#Open:MBE_DEG:Up-Open:MBE_DEG:Down         0.249332386 0.0006743
#Closed:rnaOnly:Up-Open:MBE_DEG:Down       0.359698548 0.0000000
#Closed:Temp_DEG:Up-Open:MBE_DEG:Down      0.352264964 0.0000000
#Open:Temp_DEG:Up-Open:MBE_DEG:Down        0.281486851 0.0000000
#Closed:Temp_DEG:Down-Closed:rnaOnly:Down  0.159578959 0.0021951
#Open:Temp_DEG:Down-Closed:rnaOnly:Down   -0.085782523 0.0000000
#Closed:MBE_DEG:Up-Closed:rnaOnly:Down    -0.061791035 0.0000121
#Open:MBE_DEG:Up-Closed:rnaOnly:Down      -0.332627619 0.0000000
#Closed:rnaOnly:Up-Closed:rnaOnly:Down    -0.224959978 0.0000000
#Closed:Temp_DEG:Up-Closed:rnaOnly:Down   -0.233177385 0.0000000
#Open:Temp_DEG:Up-Closed:rnaOnly:Down     -0.305092360 0.0000000
#Open:Temp_DEG:Down-Closed:Temp_DEG:Down  -0.179402217 0.0000000
#Closed:MBE_DEG:Up-Closed:Temp_DEG:Down   -0.153890684 0.0000000
#Open:MBE_DEG:Up-Closed:Temp_DEG:Down     -0.424995624 0.0000000
#Closed:rnaOnly:Up-Closed:Temp_DEG:Down   -0.318360803 0.0000000
#Closed:Temp_DEG:Up-Closed:Temp_DEG:Down  -0.326906494 0.0000000
#Open:Temp_DEG:Up-Closed:Temp_DEG:Down    -0.399326323 0.0000000
#Open:MBE_DEG:Up-Open:Temp_DEG:Down       -0.179570898 0.0000000
#Closed:rnaOnly:Up-Open:Temp_DEG:Down     -0.072927108 0.0000000
#Closed:Temp_DEG:Up-Open:Temp_DEG:Down    -0.081469870 0.0000000
#Open:Temp_DEG:Up-Open:Temp_DEG:Down      -0.153885109 0.0000000
#Open:MBE_DEG:Up-Closed:MBE_DEG:Up        -0.143117283 0.0000000
#Closed:rnaOnly:Up-Closed:MBE_DEG:Up      -0.030260841 0.0014225
#Closed:Temp_DEG:Up-Closed:MBE_DEG:Up     -0.037046224 0.0004575
#Open:Temp_DEG:Up-Closed:MBE_DEG:Up       -0.106942014 0.0000000
#Closed:rnaOnly:Up-Open:MBE_DEG:Up         0.224300890 0.0005163
#Closed:Temp_DEG:Up-Open:MBE_DEG:Up        0.217209231 0.0003886
#Open:Temp_DEG:Up-Closed:rnaOnly:Up       -0.007372743 0.0135165
#Open:Temp_DEG:Up-Closed:Temp_DEG:Up      -0.009923043 0.0066346

#n.s. 
#Open:MBE_DEG:Up-Closed:MBE_DEG:Down       0.126778999 1.0000000
#Open:Temp_DEG:Up-Closed:MBE_DEG:Down      0.165058861 0.8256646
#Closed:MBE_DEG:Up-Open:Temp_DEG:Down      0.091531764 1.0000000
#Open:Temp_DEG:Up-Open:MBE_DEG:Up          0.146902122 0.4910951
#Closed:Temp_DEG:Up-Closed:rnaOnly:Up      0.064646099 1.0000000
