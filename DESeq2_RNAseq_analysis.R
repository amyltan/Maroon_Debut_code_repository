setwd("~/Desktop/R/RNAseq_MDII/")

library('DESeq2')
library('arrayQualityMetrics')
library('vegan')
library('rgl')
library('ape')
library('pheatmap')
library('dplyr')
library('adegenet')
library('VennDiagram')
library('tidyverse')
library('ggplot2')
library('data.table')
library('ggprism')
library('ggh4x')

######
# Load saved DESeq matrix 
ll = load("RNAseq_ddsDESeq2_MDII.Rdata")

gcountsAll = read.table("v5_counts_featureCounts/counts_v5gtf.txt", header = T)
row.names(gcountsAll) <- gcountsAll$Geneid 

gcounts = as.data.frame(gcountsAll[-c(1:6)]) #removing columns for Geneid, Chr, Start, End, Strand, Length

length(gcounts[,1]) #27447 
dim(gcounts) #27447 12 

#rename columns 
colnames(gcounts) #checking order
colnames(gcounts) = sub("X.scratch.user.atan.atan_Terra_transfer.MDII.RNAseq.RNAseq.processing.MDII.12OCT23.Aligned_Reads.",
                        "", colnames(gcounts))
colnames(gcounts) = sub("_R.C.*.final.bam", "", colnames(gcounts))
colnames(gcounts) = sub("_R2.C.*.final.bam", "", colnames(gcounts))

#note: running the line starting with sub() will report what the transformation will do; good way to check without making alterations to the dataframe in case file names don't edit as expected

summary(colSums(gcounts))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#327483 1939502 5248144 4334613 6419272 7314617 

totalCounts = colSums(gcounts)
totalCounts
#C8St14  C9St14 C10St14  C8St18  C9St18 C10St18  C11F14  C12F14  C13F14  C11F18 
#5264557 3660229 2137777 7039467 5946565 7314617 1344675  516047  327483 7009760 
#C12F18  C13F18 
#5231731 6222443 

min(totalCounts) #327483
mean(totalCounts) #4334613
max(totalCounts) #7314617 


###Remove genes with low mean counts 

mns = apply(gcounts, 1, mean)

gcounts = gcounts[mns>3,] #removing genes with little or no expression

table(mns>3)
#FALSE  TRUE 
#9332 18115 

dim(gcounts)
#18115    12

##Building a dataframe to associate sample name with treatment conditions
colData <- data.frame(row.names = colnames(gcounts), 
                      chamber = c("C8St14","C9St14","C10St14","C8St18","C9St18",
                                  "C10St18","C11F14","C12F14",
                                  "C13F14","C11F18","C12F18","C13F18" ))
colData <- colData %>%
  mutate(water = case_when(
    str_detect(chamber, "St") ~ "sterile",
    str_detect(chamber, "F") ~ "filtered" 
  )) %>%
  mutate(temp = case_when(
    str_detect(chamber, "14") ~ "14",
    str_detect(chamber, "18") ~ "18"
  ))

colData$treat = paste(colData$water, colData$temp, sep = "")


###DESeq dataset and outliers 
dds <- DESeqDataSetFromMatrix(countData = gcounts, 
                              colData = colData, 
                              design = ~ treat)
vsd = varianceStabilizingTransformation(dds, blind = T)
e = ExpressionSet(assay(vsd), AnnotatedDataFrame(as.data.frame(colData(vsd))))
#arrayQualityMetrics(e, intgroup = c("treat"), force = T,
#                    outdir = "outliers_report_for_genes") #note: creates a new folder with name. DO give it an output dir--creates several reports and figures.
#C13F14 identified as an outlier by both boxplot and distances between arrays heatmap
#C12F14 identified as an outlier by boxplot 
#Because neither of these is identified as a major outlier by all three tests, both are being left.



### Wald Test - full model ### 
dds <- DESeqDataSetFromMatrix(gcounts, colData = colData, 
                              design = formula(~ water + temp))
rld = rlog(dds)
rld.df = assay(rld) #extracts matrix of normalized values

#save(rld.df, colData, file="RNAseq_ddsDESeq2_MDII.Rdata")
#write.csv(rld.df, "normalized_counts_allGenes_MDRNA.csv") #for baseline gene names

## Wald Test for WATER: filtered vs. sterile
ddsW <- DESeq(dds, minReplicatesForReplace = Inf)
resW <- results(ddsW, contrast = c('water', 'filtered', 'sterile')) #this defines contrasting conditions
#note that order for contrasts is contrast = c("cond", "B", "A") and return is log2(B/A)
#second listed should be control condition
mcols(resW, use.names = TRUE)
summary(resW)

#out of 18115 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1679, 9.3%
#LFC < 0 (down)     : 1912, 11%
#outliers [1]       : 4, 0.022%
#low counts [2]     : 703, 3.9%
#(mean count < 2)

sum(is.na(resW)) #711 

#output CSV 
resWdf <- as.data.frame(resW)
row.names(resWdf) <- gsub("gene-", "", row.names(resWdf))
#write.csv(resWdf, "~/Desktop/R/RNAseq_MDII/MBE_DEG.csv")

## Wald Test for TEMP: 18 vs. 14 
ddsT <- DESeq(dds, minReplicatesForReplace = Inf)
resT <- results(ddsT, contrast = c('temp', '18', '14'))
mcols(resT, use.names = TRUE)
summary(resT)

#out of 18115 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 6149, 34%
#LFC < 0 (down)     : 4141, 23%
#outliers [1]       : 4, 0.022%
#low counts [2]     : 0, 0%
#(mean count < 1) 

sum(is.na(resT)) #8

#output CSV
resTdf <- as.data.frame(resT)
row.names(resTdf) <- gsub("gene-", "", row.names(resTdf))
#write.csv(resTdf, "~/Desktop/R/RNAseq_MDII/Temp_DEG.csv")

#save(rld.df, colData, resW, resT, file = "MDRNA_rld&generesults_July2024.Rdata")

# Running with grouping 
de_colors <- c("blue", "coral2", "azure4")
names(de_colors) <- c("Down", "Up", "No change")

dds$group <- factor(paste0(dds$water, dds$temp))
design(dds) <- ~group
dds <- DESeq(dds)
resultsNames(dds)
ddsFS14 <- results(dds, contrast = c("group", "filtered14", "sterile14"))
summary(ddsFS14)
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1845, 10%
#LFC < 0 (down)     : 1830, 10%
FS14 <- subset(ddsFS14, ddsFS14$padj <= 0.05)
FS14 <- as.data.frame(FS14)
FS14$comparison <- "FS14"
FS14$names <- sub("gene-", "", row.names(FS14)) 
FS14$DE <- "No change"
FS14$DE[FS14$log2FoldChange > 0.5 & FS14$padj < 0.05] <- "Up"
FS14$DE[FS14$log2FoldChange < -0.5 & FS14$padj < 0.05] <- "Down" 
FS14deg <- subset(FS14, !(DE == "No change")) 
table(FS14deg$DE)
#Down   Up 
#1200 1478
#write.csv(FS14deg, "~/Desktop/Maroon_Debut/supp_info/FS14deg-FigS2Aa.csv")

FS14plot <- ggplot(data = FS14, 
                        aes(x = log2FoldChange, y = -log10(padj), col = DE)) + 
  geom_point() + scale_color_manual(values = de_colors) + 
  geom_vline(xintercept = c(-0.5, 0.5), col = "brown1") +
  geom_hline(yintercept = -log10(0.05), col = "brown1") + 
  xlab(expression(bold("log"[2]~Fold~Change))) +
  ylab(expression(bold("-log"[10]~"Adjusted p-value"))) +
  theme_prism() + 
  xlim(-9, 9) +
  #geom_text_repel(data = subset(FS14, 
  #                              DE == "Up" & (log2FoldChange > 2.5 | -log10(padj) > 2.5) &
  #                                !(names %like% "LOC")), 
  #                aes(x = log2FoldChange, y = -log10(padj), label = names), 
  #                size = 3, color = "black", max.overlaps = Inf) + 
  #geom_text_repel(data = subset(FS14, 
  #                              DE == "Down" & (log2FoldChange < -2.5 | -log10(padj) > 2.5) &
  #                                !(names %like% "LOC")), 
  #                aes(x = log2FoldChange, y = -log10(padj), label = names), 
  #                size = 3, color = "black", max.overlaps = Inf) +
  #theme(panel.grid.major = element_line(color = "#666666",
  #                                      size = 0.3,
  #                                      linetype = 2)) + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) + 
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))
print(FS14plot)
ggsave(plot = FS14plot, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/FS14plot-volcano.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)


ddsFS18 <- results(dds, contrast = c("group", "filtered18", "sterile18"))
summary(ddsFS18)
#adjusted p-value < 0.1
#LFC > 0 (up)       : 458, 2.5%
#LFC < 0 (down)     : 598, 3.3%
FS18 <- subset(ddsFS18, ddsFS18$padj <= 0.05)
FS18 <- as.data.frame(FS18)
FS18$comparison <- "FS18"
FS18$names <- sub("gene-", "", row.names(FS18))
FS18$DE <- "No change"
FS18$DE[FS18$log2FoldChange > 0.5 & FS18$padj < 0.05] <- "Up"
FS18$DE[FS18$log2FoldChange < -0.5 & FS18$padj < 0.05] <- "Down"
FS18deg <- subset(FS18, !(DE == "No change")) 
table(FS18deg$DE)
#Down   Up 
#405  311
#write.csv(FS18deg, "~/Desktop/Maroon_Debut/supp_info/FS18deg-FigS2Ba.csv")
FS18plot <- ggplot(data = FS18, 
                   aes(x = log2FoldChange, y = -log10(padj), col = DE)) + 
  geom_point() + scale_color_manual(values = de_colors) + 
  geom_vline(xintercept = c(-0.5, 0.5), col = "brown1") +
  geom_hline(yintercept = -log10(0.05), col = "brown1") + 
  xlab(expression(bold("log"[2]~Fold~Change))) +
  ylab(expression(bold("-log"[10]~"Adjusted p-value"))) +
  theme_prism() + 
  xlim(-9, 9) +
  #geom_text_repel(data = subset(FS18, 
  #                              DE == "Up" & (log2FoldChange > 2.5 | -log10(padj) > 2.5) &
  #                                !(names %like% "LOC")), #(log2FoldChange > 3.5 | -log10(padj) > 12)
  #                aes(x = log2FoldChange, y = -log10(padj), label = names), 
  #                size = 3, color = "black", max.overlaps = Inf) + 
  #geom_text_repel(data = subset(FS18, 
  #                              DE == "Down" & (log2FoldChange < -2.5 | -log10(padj) > 2.5) &
  #                                !(names %like% "LOC")), #(log2FoldChange < -2.5 | -log10(padj) > 10)
  #                aes(x = log2FoldChange, y = -log10(padj), label = names), 
  #                size = 3, color = "black", max.overlaps = Inf) +
  #theme(panel.grid.major = element_line(color = "#666666",
  #                                      size = 0.3,
  #                                      linetype = 2)) + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) + 
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))
print(FS18plot)
ggsave(plot = FS18plot, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/FS18plot-volcano.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

ddsF18S14 <- results(dds, contrast = c("group", "filtered18", "sterile14"))
summary(ddsF18S14)
#adjusted p-value < 0.1
#LFC > 0 (up)       : 3373, 19%
#LFC < 0 (down)     : 2644, 15%
F18S14 <- subset(ddsF18S14, ddsF18S14$padj <= 0.05)
F18S14 <- as.data.frame(F18S14)
F18S14$comparison <- "F18S14"
F18S14$names <- sub("gene-", "", row.names(F18S14))
F18S14$DE <- "No change"
F18S14$DE[F18S14$log2FoldChange > 0.5 & F18S14$padj < 0.05] <- "Up"
F18S14$DE[F18S14$log2FoldChange < -0.5 & F18S14$padj < 0.05] <- "Down"
F18S14deg <- subset(F18S14, !(DE == "No change")) 
table(F18S14deg$DE)
#Down   Up 
#1930 2395 
#write.csv(F18S14deg, "~/Desktop/Maroon_Debut/supp_info/F18S14deg-FigS2Da.csv")
F18S14plot <- ggplot(data = F18S14, 
                   aes(x = log2FoldChange, y = -log10(padj), col = DE)) + 
  geom_point() + scale_color_manual(values = de_colors) + 
  geom_vline(xintercept = c(-0.5, 0.5), col = "brown1") +
  geom_hline(yintercept = -log10(0.05), col = "brown1") + 
  xlab(expression(bold("log"[2]~Fold~Change))) +
  ylab(expression(bold("-log"[10]~"Adjusted p-value"))) +
  theme_prism() + 
  xlim(-9, 9) +
  #geom_text_repel(data = subset(F18S14, 
  #                              DE == "Up" & (log2FoldChange > 2.5 | -log10(padj) > 2.5) &
  #                                !(names %like% "LOC")), #(log2FoldChange > 3.5 | -log10(padj) > 12)
  #                aes(x = log2FoldChange, y = -log10(padj), label = names), 
  #                size = 3, color = "black", max.overlaps = Inf) + 
  #geom_text_repel(data = subset(F18S14, 
  #                              DE == "Down" & (log2FoldChange < -2.5 | -log10(padj) > 2.5) &
  #                                !(names %like% "LOC")), #(log2FoldChange < -2.5 | -log10(padj) > 10)
  #                aes(x = log2FoldChange, y = -log10(padj), label = names), 
  #                size = 3, color = "black", max.overlaps = Inf) +
  #theme(panel.grid.major = element_line(color = "#666666",
  #                                      size = 0.3,
  #                                      linetype = 2)) + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) + 
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))
print(F18S14plot)
ggsave(plot = F18S14plot, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/F18S14plot-volcano.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

ddsF14S18 <- results(dds, contrast = c("group", "filtered14", "sterile18"))
summary(ddsF14S18)
#adjusted p-value < 0.1
#LFC > 0 (up)       : 3731, 21%
#LFC < 0 (down)     : 5064, 28%
F14S18 <- subset(ddsF14S18, ddsF14S18$padj <= 0.05)
F14S18 <- as.data.frame(F14S18)
F14S18$comparison <- "F14S18"
F14S18$names <- sub("gene-", "", row.names(F14S18))
F14S18$DE <- "No change"
F14S18$DE[F14S18$log2FoldChange > 0.5 & F14S18$padj < 0.05] <- "Up"
F14S18$DE[F14S18$log2FoldChange < -0.5 & F14S18$padj < 0.05] <- "Down"
F14S18deg <- subset(F14S18, !(DE == "No change")) 
table(F14S18deg$DE)
#Down   Up 
#4184 3293
#write.csv(F14S18deg, "~/Desktop/Maroon_Debut/supp_info/F14S18deg-FigS2Ca.csv")

F14S18plot <- ggplot(data = F14S18, 
                     aes(x = log2FoldChange, y = -log10(padj), col = DE)) + 
  geom_point() + scale_color_manual(values = de_colors) + 
  geom_vline(xintercept = c(-0.5, 0.5), col = "brown1") +
  geom_hline(yintercept = -log10(0.05), col = "brown1") + 
  xlab(expression(bold("log"[2]~Fold~Change))) +
  ylab(expression(bold("-log"[10]~"Adjusted p-value"))) +
  theme_prism() + 
  xlim(-9, 9) +
  #geom_text_repel(data = subset(F14S18, 
  #                              DE == "Up" & (log2FoldChange > 2.5 | -log10(padj) > 2.5) &
  #                                !(names %like% "LOC")), #(log2FoldChange > 3.5 | -log10(padj) > 12)
  #                aes(x = log2FoldChange, y = -log10(padj), label = names), 
  #                size = 3, color = "black", max.overlaps = Inf) + 
  #geom_text_repel(data = subset(F14S18, 
  #                              DE == "Down" & (log2FoldChange < -2.5 | -log10(padj) > 2.5) &
  #                                !(names %like% "LOC")), #(log2FoldChange < -2.5 | -log10(padj) > 10)
  #                aes(x = log2FoldChange, y = -log10(padj), label = names), 
  #                size = 3, color = "black", max.overlaps = Inf) +
  #theme(panel.grid.major = element_line(color = "#666666",
  #                                     size = 0.3,
  #                                      linetype = 2)) + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) + 
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))
print(F14S18plot)
ggsave(plot = F14S18plot, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/F14S18plot-volcano.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

comp <- rbind(FS14deg, FS18deg, F18S14deg, F14S18deg) #15196 
length(unique(comp$names)) #9572 

ref <- comp %>% 
  group_by(names) %>% 
  summarise(comparison = toString(comparison))

refFS14 <- subset(ref, comparison == "FS14")
refFS18 <- subset(ref, comparison == "FS18")
refF18S14 <- subset(ref, comparison == "F18S14")
refF14S18 <- subset(ref, comparison == "F14S18")

FS14only <- subset(FS14deg, names %in% refFS14$names)
table(FS14only$DE)
#Down   Up 
#32   44 
#write.csv(FS14only, "~/Desktop/Maroon_Debut/supp_info/FS14only-FigS2Ab.csv")

FS14onlyplot <- ggplot(data = FS14only, 
                     aes(x = log2FoldChange, y = -log10(padj), col = DE)) + 
  geom_point() + scale_color_manual(values = de_colors) + 
  geom_vline(xintercept = c(-0.5, 0.5), col = "brown1") +
  geom_hline(yintercept = -log10(0.05), col = "brown1") + 
  xlab(expression(bold("log"[2]~Fold~Change))) +
  ylab(expression(bold("-log"[10]~"Adjusted p-value"))) +
  theme_prism() + 
  xlim(-9, 9) +
  #geom_text_repel(data = subset(FS14only, 
  #                              DE == "Up" & (log2FoldChange > 2.5 | -log10(padj) > 2.5) &
  #                                !(names %like% "LOC")),
  #                aes(x = log2FoldChange, y = -log10(padj), label = names), 
  #                size = 3, color = "black", max.overlaps = Inf) + 
  #geom_text_repel(data = subset(FS14only, 
  #                              DE == "Down" & (log2FoldChange < -2.5 | -log10(padj) > 2.5) &
  #                                !(names %like% "LOC")), 
  #                aes(x = log2FoldChange, y = -log10(padj), label = names), 
  #                size = 3, color = "black", max.overlaps = Inf) +
  #theme(panel.grid.major = element_line(color = "#666666",
  #                                      size = 0.3,
  #                                      linetype = 2)) + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) + 
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))
print(FS14onlyplot)
ggsave(plot = FS14onlyplot, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/FS14onlyplot-volcano.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

FS18only <- subset(FS18deg, names %in% refFS18$names)
table(FS18only$DE)
#Down   Up 
#77   37
#write.csv(FS18only, "~/Desktop/Maroon_Debut/supp_info/FS18only-FigS2Bb.csv")

FS18onlyplot <- ggplot(data = FS18only, 
                       aes(x = log2FoldChange, y = -log10(padj), col = DE)) + 
  geom_point() + scale_color_manual(values = de_colors) + 
  geom_vline(xintercept = c(-0.5, 0.5), col = "brown1") +
  geom_hline(yintercept = -log10(0.05), col = "brown1") + 
  xlab(expression(bold("log"[2]~Fold~Change))) +
  ylab(expression(bold("-log"[10]~"Adjusted p-value"))) +
  theme_prism() + 
  xlim(-9, 9) +
  #geom_text_repel(data = subset(FS18only, 
  #                              DE == "Up" & (log2FoldChange > 2.5 | -log10(padj) > 2.5) &
  #                                !(names %like% "LOC")),
  #                aes(x = log2FoldChange, y = -log10(padj), label = names), 
  #                size = 3, color = "black", max.overlaps = Inf) + 
  #geom_text_repel(data = subset(FS18only, 
  #                              DE == "Down" & (log2FoldChange < -2.5 | -log10(padj) > 2.5) &
  #                                !(names %like% "LOC")), 
  #                aes(x = log2FoldChange, y = -log10(padj), label = names), 
  #                size = 3, color = "black", max.overlaps = Inf) +
  #theme(panel.grid.major = element_line(color = "#666666",
  #                                      size = 0.3,
  #                                      linetype = 2)) + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) + 
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))
print(FS18onlyplot)
ggsave(plot = FS18onlyplot, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/FS18onlyplot-volcano.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

F18S14only <- subset(F18S14deg, names %in% refF18S14$names)
table(F18S14only$DE)
#Down   Up 
#626  905 
#write.csv(F18S14only, "~/Desktop/Maroon_Debut/supp_info/F18S14only-FigS2Db.csv")

F18S14onlyplot <- ggplot(data = F18S14only, 
                       aes(x = log2FoldChange, y = -log10(padj), col = DE)) + 
  geom_point() + scale_color_manual(values = de_colors) + 
  geom_vline(xintercept = c(-0.5, 0.5), col = "brown1") +
  geom_hline(yintercept = -log10(0.05), col = "brown1") + 
  xlab(expression(bold("log"[2]~Fold~Change))) +
  ylab(expression(bold("-log"[10]~"Adjusted p-value"))) +
  theme_prism() + 
  xlim(-9, 9) +
  #geom_text_repel(data = subset(F18S14only, 
  #                              DE == "Up" & (log2FoldChange > 2.5 | -log10(padj) > 2.5) &
  #                                !(names %like% "LOC")),
  #                aes(x = log2FoldChange, y = -log10(padj), label = names), 
  #                size = 3, color = "black", max.overlaps = Inf) + 
  #geom_text_repel(data = subset(F18S14only, 
  #                              DE == "Down" & (log2FoldChange < -2.5 | -log10(padj) > 2.5) &
  #                                !(names %like% "LOC")), 
  #                aes(x = log2FoldChange, y = -log10(padj), label = names), 
  #                size = 3, color = "black", max.overlaps = Inf) +
  #theme(panel.grid.major = element_line(color = "#666666",
  #                                      size = 0.3,
  #                                      linetype = 2)) + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) + 
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))
print(F18S14onlyplot)
ggsave(plot = F18S14onlyplot, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/F18S14onlyplot-volcano.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

F14S18only <- subset(F14S18deg, names %in% refF14S18$names)
table(F14S18only$DE)
#Down   Up 
#1908 1153 
#write.csv(F14S18only, "~/Desktop/Maroon_Debut/supp_info/F14S18only-FigS2Cb.csv")

F14S18onlyplot <- ggplot(data = F14S18only, 
                         aes(x = log2FoldChange, y = -log10(padj), col = DE)) + 
  geom_point() + scale_color_manual(values = de_colors) + 
  geom_vline(xintercept = c(-0.5, 0.5), col = "brown1") +
  geom_hline(yintercept = -log10(0.05), col = "brown1") + 
  xlab(expression(bold("log"[2]~Fold~Change))) +
  ylab(expression(bold("-log"[10]~"Adjusted p-value"))) +
  theme_prism() + 
  xlim(-9, 9) +
  #geom_text_repel(data = subset(F14S18only, 
  #                              DE == "Up" & (log2FoldChange > 2.5 | -log10(padj) > 2.5) &
  #                                !(names %like% "LOC")),
  #                aes(x = log2FoldChange, y = -log10(padj), label = names), 
  #                size = 3, color = "black", max.overlaps = Inf) + 
  #geom_text_repel(data = subset(F14S18only, 
  #                              DE == "Down" & (log2FoldChange < -2.5 | -log10(padj) > 2.5) &
  #                                !(names %like% "LOC")), 
  #                aes(x = log2FoldChange, y = -log10(padj), label = names), 
  #                size = 3, color = "black", max.overlaps = Inf) +
  #theme(panel.grid.major = element_line(color = "#666666",
  #                                      size = 0.3,
  #                                      linetype = 2)) + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) + 
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))
print(F14S18onlyplot)
ggsave(plot = F14S18onlyplot, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/F14S18onlyplot-volcano.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

#
#######
ll2 = load("MDRNA_rld&generesults_July2024.Rdata")

head(rld.df)
head(resW)
head(resT)

vals = cbind(resW$stat, resW$pvalue, resW$padj, 
             resT$stat, resT$pvalue, resT$padj)
#collect pvalues for each gene 
colnames(vals) = c("stat_resW", "pval_resW", "padj_resW", 
                   "stat_resT", "pval_resT", "padj_resT")
rldpvals = as.data.frame(cbind(rld.df, vals)) #combine RLD with pvalues 
row.names(rldpvals) <- gsub("gene-", "", row.names(rldpvals))
#write.csv(rldpvals, "MDII_rldpvals.csv")

#subset rld and pvals table to only those significant for each comparison 
sigpvals_W = rldpvals[rldpvals$padj_resW <= 0.05 & !is.na (rldpvals$padj_resW),]
dim(sigpvals_W) #2408   18 #significant only based on p-value!!! Volcano plot further narrows based on lfc.


sigpvals_T = rldpvals[rldpvals$padj_resT <= 0.05 & !is.na (rldpvals$padj_resT),]
dim(sigpvals_T) #9117   18

#save pvalue data 
#save(rldpvals, sigpvals_W, sigpvals_T, file = "rldpvals_MDIIRNA_July2024.Rdata")

#Venn Diagram - only p-values (does not match numbers in Volcano)#####
#see below Volcano plots for Venn diagram that matches
### Venn Diagram ### ####Note: not excluding low fold change values
library('VennDiagram')
rldpvals = as.data.frame(rldpvals)
wat = row.names(rldpvals[rldpvals$padj_resW<0.05 & !is.na(rldpvals$padj_resW),])
temp = row.names(rldpvals[rldpvals$padj_resT < 0.05 & !is.na (rldpvals$padj_resT),])
degs05 = union(wat,temp)
length(degs05) #9979
candidates = list("Microbial\nexposure" = wat, "Temperature" = temp)

#"#FFCC00", "#CC9933", "#336666", "#33CCCC"

quartz()
venn = venn.diagram(
  x = candidates, 
  filename = NULL, 
  col = "transparent", 
  fill = c("#00CCCC","#FFCC33"),
  cex = 3, 
  fontfamily = "sans", 
  cat.default.pos = "text", 
  #cat.col = c("#336666", "#996600"), 
  cat.cex = 0, #setting font size to 0 to remove labels, also hashed out other irrelevant lines
  #cat.fontfamily = "sans", 
  #cat.dist = c(0.20, 0.13), #dist from edge of circle 
  #cat.pos = c(20, -40), #degrees along circle, 0 = 12 o'clock
  rotation.degree = 180
); 
grid.draw(venn) 

#Principal Coordinate Calculation#####
### Principal coordinate calculation
##PCoA works off dissimilarity matrix (calculated by vegdist) vs. PCA that calculates similarity

dd.veg = vegdist(t(rld.df), "manhattan") #t() transposes frame so that chambers are the rows and genes are columns. Columns are used to compare different samples.
div.dd.veg = dd.veg/1000
head(div.dd.veg)

dd.pcoa = pcoa(div.dd.veg)
head(dd.pcoa)
scores = dd.pcoa$vectors

conditions = colData
conditions$water = as.factor(conditions$water)
conditions$temp = as.factor(conditions$temp)
conditions$treat = as.factor(conditions$treat)

dd.pcoa
#$values
#Eigenvalues Relative_eig Broken_stick Cumul_eig Cumul_br_stick
#1   528.496463  0.634660302  0.274534304 0.6346603      0.2745343
#2    92.200933  0.110722164  0.183625213 0.7453825      0.4581595
#3    62.288934  0.074801473  0.138170668 0.8201839      0.5963302

### plotting PCoA 

#plotting first and second principal axes 
colorsW <- c("#336666","#CC9933","#33CCCC","#FFCC00") 
conditions$treat_c <- colorsW[as.factor(conditions$treat)]
legendNames <- c("14ºC microbiome-exposed", "18ºC microbiome-exposed", 
                 "14ºC sterile ASW", "18ºC sterile ASW")

par(mai = c(0.8,0.9,0.1,0.1))
plot(scores[,1], scores[,2], col = conditions$treat_c, pch = 19, cex = 3,
     xlab = "PCo1", ylab = "PCo2", cex.axis = 1.2, cex.lab = 1.9, 
     xlim = c(-12, 17), ylim = c(-8, 12))
ordispider(scores, conditions$treat, label = F)
#ordiellipse(ord = scores, groups = conditions$treat, 
#            conf = 0.95)
legend(legend = c("14ºC MBE", "18ºC MBE", 
                             "14ºC SASW", "18ºC SASW"), 
       col = colorsW, pch = 19, cex = 1.3, bty = "n",
       #x = c(-7.8, -4), y = c(6.7, 3)) 
       x = c(-11, -11), y = c(12, 12))

#ggplot PCoA 
colorsPC <- c("#336666","#336666","#336666",
              "#CC9933","#CC9933","#CC9933",
              "#33CCCC","#33CCCC","#33CCCC",
              "#FFCC00","#FFCC00","#FFCC00") 
names(colorsPC) <- c("C11F14", "C12F14", "C13F14", 
                     "C11F18", "C12F18", "C13F18",
                     "C8St14", "C9St14", "C10St14", 
                     "C8St18", "C9St18", "C10St18")
scoresplot <- as.data.frame(scores)

#add ordispider but ggplot equivalent
#calculate centroids for each treatment group 
scoresplot$SampleID <- rownames(scoresplot)
scoresplot$treat <- colData$treat 
centroids <- aggregate(cbind(Axis.1, Axis.2) ~ treat, 
                       data = scoresplot, FUN = mean)

treatShapes = c(16, 17, 16, 17)
names(treatShapes) = c("sterile14", "sterile18", 
                       "filtered14", "filtered18")

pcoa <- ggplot() + 
  theme_prism() +
  geom_point(data = scoresplot, 
             aes(x = scoresplot[,1], y = scoresplot[,2], 
                 color = conditions$chamber, shape = treat), 
             size = 4) + 
  scale_color_manual(values = colorsPC) + 
  scale_shape_manual(values = treatShapes) +
  xlab("PCoA1 (63%)") + 
  ylab("PCoA2 (11%)") + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) +
  geom_segment(data = scoresplot, aes(x = Axis.1, y = Axis.2,
                                      xend = centroids$Axis.1[match(treat, centroids$treat)],
                                      yend = centroids$Axis.2[match(treat, centroids$treat)],
                                      color = conditions$chamber)) +
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))
print(pcoa)
ggsave(plot = pcoa, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/rna-pcoa.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

#plotting only pcoa legend
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', 
     xlim=0:1, ylim=0:1) #creates a blank plot
plotSym <- c(16, 16, 16, 16) #shapes list
par(mai = c(0.5,0.5,0.5,0.5))
legend("center", 
       legend = c("14ºC SASW", 
                  "18ºC SASW", 
                  "14ºC MBE", 
                  "18ºC MBE"), 
       pch = plotSym, pt.cex = 5.1,
       col = c("#33CCCC", "#FFCC00", "#336666", "#CC9933"), 
       cex = 2.5, bty = 'n', par(mai = c(0.1,0.1,0.1,0.1))
) 

#plotting second and third principal axes 
plot(scores[,2], scores[,3], col = conditions$treat_c, pch = 19, 
     xlab = "PCo2", ylab = "PCo3")
ordispider(scores[,2:3], conditions$treat, label = F)
legend("topleft", legend = levels(conditions$treat), col = colorsW, pch = 19)

ad = adonis2(t(rld.df)~water*temp, data = conditions, method = "manhattan")
ad 
#note: p-values vary slightly between runs due to permutations. 
#Copied below are the results in the paper (as of 12Feb2025):
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#adonis2(formula = t(rld.df) ~ water * temp, data = conditions, method = "manhattan")
#            Df  SumOfSqs      R2       F Pr(>F)    
#water       1 114318905 0.13728  3.3978  0.045 *  
#temp        1 391233700 0.46982 11.6282  0.001 ***
#water:temp  1  58007650 0.06966  1.7241  0.180    
#Residual    8 269163111 0.32323                   
#Total      11 832723365 1.00000                   
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

labels = c("Microbe\ncondition", "Temperature\ncondition", "Interaction", "Residuals")
cols = c("#669966","#FF9900","gray30","grey80") 
labels2 = paste(labels, round(ad$R2[1:4]*100, digits = 1))
pie(ad$R2[1:4], labels = labels2, col = cols, main = "Gene Expression")

#Volcano plots#####
library('ggrepel')
library('ggprism')
library('data.table')

#checking min and max fold changes in significant regions 
resWsig <- subset(resW, resW$padj <= 0.05)
min(abs(resWsig$log2FoldChange)) #0.2436684
max(abs(resWsig$log2FoldChange)) #5.008481 
sum(resWsig$padj <= 0.5 & (abs(resWsig$log2FoldChange) > 0.5)) #2185
sum(resWsig$log2FoldChange > 0.5) #0.5 - 1062, 1 - 565, 2 - 79, 3 - 17
sum(resWsig$log2FoldChange < -0.5) #1123
plotW <- na.omit(resW) #removes if NA value; NAs occur in pval and padj if B/A didn't have a divisible number

resTsig <- subset(resT, resT$padj <= 0.05)
min(abs(resTsig$log2FoldChange)) #0.2032912
max(abs(resTsig$log2FoldChange)) #7.331406
sum(resTsig$padj <= 0.5 & (abs(resTsig$log2FoldChange) > 0.5)) #8630
sum(resTsig$log2FoldChange > 0.5) #0.5 - 5186, #1 - 2903, 2 - 193, 3 - 14
sum(resTsig$log2FoldChange < -0.5) #3444
plotT <- na.omit(resT)

# WATER
plotW$DE <- "No change"
plotW$DE[plotW$log2FoldChange > 0.5 & plotW$padj < 0.05] <- "Up"
plotW$DE[plotW$log2FoldChange < -0.5 & plotW$padj < 0.05] <- "Down"
plotW$geneId <- gsub("gene-", "", row.names(plotW)) 
#write.csv(plotW, "water_ATAC_DAgenes.csv")

de_colors <- c("blue", "coral2", "azure4")
names(de_colors) <- c("Down", "Up", "No change")

water_volcano <- ggplot(data = plotW, 
       aes(x = log2FoldChange, y = -log10(padj), col = DE)) + 
  geom_point() + scale_color_manual(values = de_colors) + 
  geom_vline(xintercept = c(-0.5, 0.5), col = "brown1") +
  geom_hline(yintercept = -log10(0.05), col = "brown1") + 
  xlab(expression(bold("log"[2]~Fold~Change))) +
  ylab(expression(bold("-log"[10]~"Adjusted p-value"))) +
  theme_prism() + 
  xlim(-9, 9) +
  #geom_text_repel(data = subset(plotW, 
  #                              DE == "Up" & (log2FoldChange > 2.5 | -log10(padj) > 2.5) &
  #                                !(geneId %like% "LOC")), #(log2FoldChange > 3.5 | -log10(padj) > 12)
  #                aes(x = log2FoldChange, y = -log10(padj), label = geneId), 
  #                size = 5, color = "black", max.overlaps = Inf) + 
  #geom_text_repel(data = subset(plotW, 
  #                              DE == "Down" & (log2FoldChange < -2.5 | -log10(padj) > 2.5) &
  #                                !(geneId %like% "LOC")), #(log2FoldChange < -2.5 | -log10(padj) > 10)
  #                aes(x = log2FoldChange, y = -log10(padj), label = geneId), 
  #                size = 5, color = "black", max.overlaps = Inf) +
  #theme(panel.grid.major = element_line(color = "#666666",
  #                                      size = 0.3,
  #                                      linetype = 2)) + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) + 
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))
print(water_volcano)

ggsave(plot = water_volcano, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/rna-water-volcano.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

# TEMPERATURE
plotT$DE <- "No change"
plotT$DE[plotT$log2FoldChange > 0.5 & plotT$padj < 0.05] <- "Up"
plotT$DE[plotT$log2FoldChange < -0.5 & plotT$padj < 0.05] <- "Down"
plotT$names <- gsub("gene-", "", row.names(plotT))
#write.csv(plotT, "temp_ATAC_DAgenes.csv")

temp_volcano <- ggplot(data = plotT, 
       aes(x = log2FoldChange, y = -log10(padj), col = DE)) + 
  geom_point() + scale_color_manual(values = de_colors) + 
  geom_vline(xintercept = c(-0.5, 0.5), col = "brown1") +
  geom_hline(yintercept = -log10(0.05), col = "brown1") + 
  xlab(expression(bold("log"[2]~Fold~Change))) +
  ylab(expression(bold("-log"[10]~"Adjusted p-value"))) +
  theme_prism() + 
  xlim(-9, 9) +
  #geom_text_repel(data = subset(plotT, 
  #                              DE == "Up" & (log2FoldChange > 3 | -log10(padj) > 7) & 
  #                                !(names %like% "LOC")), 
  #                aes(x = log2FoldChange, y = -log10(padj), label = names), 
  #                size = 3, color = "black", max.overlaps = Inf) + 
  #geom_text_repel(data = subset(plotT, 
  #                              DE == "Down" & (log2FoldChange < -4.5 | -log10(padj) > 9) &
  #                                !(names %like% "LOC")), 
  #                aes(x = log2FoldChange, y = -log10(padj), label = names), 
  #                size = 3, color = "black", max.overlaps = Inf) + 
  #theme(panel.grid.major = element_line(color = "#666666",
  #                                      size = 0.3,
  #                                      linetype = 2)) +
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) +
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))
print(temp_volcano)

ggsave(plot = temp_volcano, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/rna-temp-volcano.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

#plotting only legend
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', 
     xlim=0:1, ylim=0:1) #creates a blank plot
plotSym <- c(16, 16) #shapes list
par(mai = c(0.5,0.2,0.5,0.5)) #bottom, left, top, right
legend("center", 
       legend = c("Up-regulated", 
                  "Down-regulated", 
                  "No change"), 
       pch = plotSym, pt.cex = 5.1,
       col = c("coral2", "blue", "azure4"), 
       cex = 2.5, bty = 'n', par(mai = c(0.1,0.1,0.1,0.1))
)

# Plotting Venn based on both lfc and padj 
library('VennDiagram')
plotW = as.data.frame(plotW)
plotT = as.data.frame(plotT)
wat = row.names(plotW[plotW$padj<0.05 & (abs(plotW$log2FoldChange) > 0.5),])
temp = row.names(plotT[plotT $padj < 0.05 & (abs(plotT$log2FoldChange) > 0.5),])
degs05 = union(wat,temp)
length(degs05) #9424
candidates = list("Microbial\nexposure" = wat, "Temperature" = temp)

#"#FFCC00", "#CC9933", "#336666", "#33CCCC"
quartz()
venn = venn.diagram(
  x = candidates, 
  filename = NULL, 
  col = "transparent", 
  fill = c("#00CCCC","#FFCC33"),
  cex = 3, 
  fontfamily = "sans", 
  cat.default.pos = "text", 
  #cat.col = c("#336666", "#996600"), 
  cat.cex = 0, #setting font size to 0 to remove labels, also hashed out other irrelevant lines
  #cat.fontfamily = "sans", 
  #cat.dist = c(0.20, 0.13), #dist from edge of circle 
  #cat.pos = c(20, -40), #degrees along circle, 0 = 12 o'clock
  rotation.degree = 180
); 
grid.draw(venn) 

### Adding gene name conversion information #####
conv <- read.csv("~/Desktop/TAMU/Spurp_reference_files/Spurp_genome/AmyGeneNameConversions.csv")
colnames(conv) <- c("names", "description") 
#resWdf <- read.csv("water_ATAC_DAgenes.csv")
#resTdf <- read.csv("temp_ATAC_DAgenes.csv")

resWdf <- as.data.frame(plotW) 
colnames(resWdf)[colnames(resWdf) == "geneId"] <- "names"
resWdes <- merge(resWdf, conv, by = c("names"), all = TRUE)
resWdes <- resWdes[!is.na(resWdes$log2FoldChange),]
#write.csv(resWdes, "microbe-exposureDEgenes-all.csv")
W <- read.csv("microbe-exposureDEgenes-all.csv") 
W <- subset(W, select = -c(X)) 
Wup <- subset(W, DE == "Up") 
write.csv(Wup, "microbe-exposureDEgenes-UP.csv")
Wdown <- subset(W, DE == "Down") 
write.csv(Wdown, "microbe-exposureDEgenes-DOWN.csv")

resTdf <- as.data.frame(plotT)
colnames(resTdf)[colnames(resTdf) == "geneId"] <- "names"
resTdes <- merge(resTdf, conv, by = c("names"), all = TRUE)
resTdes <- resTdes[!is.na(resTdes$log2FoldChange),]
#write.csv(resTdes, "temp-DEgenes-all.csv") 
T <- read.csv("temp-DEgenes-all.csv")
T <- subset(T, select = -c(X)) 
Tup <- subset(T, DE == "Up") 
write.csv(Tup, "temp-DEgenes-UP.csv")
Tdown <- subset(T, DE == "Down") 
write.csv(Tdown, "temp-DEgenes-DOWN.csv")

# subsetting genes that were labeled in volcano plot 
WdesLabeled = subset(resWdes, ((DE == "Up" & (log2FoldChange > 3 | -log10(padj) > 11)) | 
                       (DE == "Down" & (log2FoldChange < -2.5 | -log10(padj) > 10))))
WdesLabeled <- subset(WdesLabeled, select = -c(baseMean, lfcSE, stat, pvalue))                  
WdesID <- WdesLabeled %>% filter(!str_detect(description, "uncharacterized protein"))
WdesID <- subset(WdesID, select = -c(log2FoldChange, padj)) 
WdesID
#write.csv(WdesID, "water-DEgeneConvertedNames.csv")

TdesLabeled = subset(resWdes, ((DE == "Up" & (log2FoldChange > 3 | -log10(padj) > 11)) | 
                                 (DE == "Down" & (log2FoldChange < -4.5 | -log10(padj) > 15))))
TdesLabeled <- subset(TdesLabeled, select = -c(baseMean, lfcSE, stat, pvalue))                  
TdesID <- TdesLabeled %>% filter(!str_detect(description, "uncharacterized protein"))
TdesID <- subset(TdesID, select = -c(log2FoldChange, padj)) 
TdesID
#write.csv(TdesID, "temp-DEgeneConvertedNames.csv")

## Creating list of genes altered by both conditions 
lap <- merge(resWdf, resTdf, by = c("names"), all = TRUE)
colnames(lap) <- c("names", "baseMean_W", "log2FoldChange_W", "lfcSE_W", "stat_W", "pavlue_W", "padj_W", "DE_W", 
                   "baseMean_T", "log2FoldChange_T", "lfcSE_T", "stat_T", "pavlue_T", "padj_T", "DE_T")
lap <- subset(lap, select = -c(baseMean_W, lfcSE_W, stat_W, pavlue_W, baseMean_T, lfcSE_T, stat_T, pavlue_T))
head(lap)
lapDes <- merge(lap, conv, by = c("names"), all = TRUE)

sigBoth <- subset(lapDes, lapDes$padj_W <= 0.05 & lapDes$padj_T <= 0.05)
sum(str_detect(sigBoth$description, "uncharacterized protein")) #357
sum(str_detect(sigBoth$DE_W, "Up") & str_detect(sigBoth$DE_T, "Up")) #22 
sum(str_detect(sigBoth$DE_W, "Down") & str_detect(sigBoth$DE_T, "Down")) #53 
sum(str_detect(sigBoth$DE_W, "Up") & str_detect(sigBoth$DE_T, "Down")) #681 
sum(str_detect(sigBoth$DE_W, "Down") & str_detect(sigBoth$DE_T, "Up")) #635 
sum(str_detect(sigBoth$DE_W, "No change") & str_detect(sigBoth$DE_T, "Up")) #24
sum(str_detect(sigBoth$DE_W, "Up") & str_detect(sigBoth$DE_T, "No change")) #19
sum(str_detect(sigBoth$DE_W, "No change") & str_detect(sigBoth$DE_T, "Down")) #13
sum(str_detect(sigBoth$DE_W, "Down") & str_detect(sigBoth$DE_T, "No change")) #42
sum(str_detect(sigBoth$DE_W, "Down") & 
      (str_detect(sigBoth$DE_T, "No change") | str_detect(sigBoth$DE_T, "Up"))) #677
sum(str_detect(sigBoth$DE_W, "Up") & 
      (str_detect(sigBoth$DE_T, "No change") | str_detect(sigBoth$DE_T, "Down"))) #700
sum(str_detect(sigBoth$DE_T, "Up") & 
      (str_detect(sigBoth$DE_W, "No change") | str_detect(sigBoth$DE_W, "Down"))) #659
sum(str_detect(sigBoth$DE_T, "Down") & 
      (str_detect(sigBoth$DE_W, "No change") | str_detect(sigBoth$DE_W, "Up"))) #694

sigBothFiltered <- sigBoth %>% 
  filter(!str_detect(DE_W, "No change"
                     )) %>% 
  filter(!str_detect(DE_T, "No change"
                     )) %>% 
  filter(!str_detect(description, "uncharacterized protein"
  ))

write.csv(sigBothFiltered, "overlapGenesFiltered-RNA.csv", row.names = FALSE)

# Venn diagrams of overlaps between conditions#####
## Venn diagram of overlaps - remove for publish
####Note: only contains genes that are significant in BOTH conditions (lower total count than earlier Venn that plotted # sig in EITHER)
row.names(sigBoth) <- sigBoth[,1]
UW = row.names(sigBoth[(str_detect(sigBoth$DE_W, "Up")),])
DW = row.names(sigBoth[(str_detect(sigBoth$DE_W, "Down")),])
UT = row.names(sigBoth[(str_detect(sigBoth$DE_T, "Up")),])
DT = row.names(sigBoth[(str_detect(sigBoth$DE_T, "Down")),])

candidates3 = list("Up-regulated\nin 18ºC" = UT, 
                   "Up-regulated\nby microbial\n-exposure" = UW, 
                   "Down-regulated\nin 18ºC" = DT, 
                   "Down-regulated\n by microbial\n -exposure" = DW)

quartz()
venn4 = venn.diagram(
  x = candidates3, 
  filename = NULL, 
  col = "transparent", 
  fill = c("orange", "darkred","#3333CC", "#6633CC"),
  cex = 1.5, 
  fontfamily = "sans", 
  cat.default.pos = "text", 
  cat.col = c( "#FF6633", "darkred", "#000066", "#330066"), 
  cat.cex = 1.5, 
  cat.fontfamily = "sans", 
  cat.dist = c(0.18, 0.2, 0.12, 0.15), 
  cat.pos = 1, 
); 
grid.draw(venn4)

