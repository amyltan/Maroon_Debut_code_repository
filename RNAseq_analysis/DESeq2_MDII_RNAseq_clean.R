#DESeq2 script for MDII 

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

# Saved DESeq matrices, etc ---- 
#ll = load("RNAseq_ddsDESeq2_MDII_July25.Rdata")

# Organizing count data ---- 
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

# Removing genes with low mean counts 
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

# DESeq dataset and checking outliers ---- 
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

# DESeq matrix and model #### 
dds <- DESeqDataSetFromMatrix(gcounts, colData = colData, 
                              design = formula(~ water * temp)) 

#checking reference levels. Ref. condition is listed first. 
dds$water #Levels: filtered sterile 
dds$water = relevel(dds$water, "sterile") #list the condition to use as ref. 
dds$water #Levels: sterile filtered
dds$temp #Levels: 14 18 #good to go

rld = rlog(dds)
rld.df = assay(rld) #extracts matrix of normalized values

#save(rld.df, colData, file="RNAseq_ddsDESeq2_MDII_July25.Rdata")
#write.csv(rld.df, "normalized_counts_allGenes_MDRNA_July25.csv") #for baseline gene names

# PERMANOVA - Adonis2 test and results ----
conditions = colData
conditions$water = as.factor(conditions$water)
conditions$temp = as.factor(conditions$temp)
conditions$treat = as.factor(conditions$treat)
ad <- adonis2(t(rld.df) ~ water * temp, data = conditions, 
              method = "manhattan")
ad #note p-values will change slightly due to permutations 
#           Df  SumOfSqs      R2       F Pr(>F)    
#water       1 114318905 0.13728  3.3978  0.039 *  
#temp        1 391233700 0.46982 11.6282  0.001 ***
#water:temp  1  58007650 0.06966  1.7241  0.175    
#Residual    8 269163111 0.32323 

# Wald Test Results ---- 
ddsR <- DESeq(dds, minReplicatesForReplace = Inf) #note: model was set up earlier when creating the DESeqDataSet
resultsNames(ddsR) #second name listed is the Ref. condition
#[1] "Intercept"                 "water_filtered_vs_sterile"
#[3] "temp_18_vs_14"             "waterfiltered.temp18"  

# to check model matrix (either of these lines should work) 
attr(ddsR, "modelMatrix") 
model.matrix(~ water * temp, colData(ddsR)) 
#note: for more complicated models, might need to double check betaPrior (FALSE)

#        Intercept water_filtered_vs_sterile temp_18_vs_14 waterfiltered.temp18
#C8St14          1                         0             0                    0
#C9St14          1                         0             0                    0
#C10St14         1                         0             0                    0
#C8St18          1                         0             1                    0
#C9St18          1                         0             1                    0
#C10St18         1                         0             1                    0
#C11F14          1                         1             0                    0
#C12F14          1                         1             0                    0
#C13F14          1                         1             0                    0
#C11F18          1                         1             1                    1
#C12F18          1                         1             1                    1
#C13F18          1                         1             1                    1


resW <- results(ddsR, name = c("water_filtered_vs_sterile"))
#results(ddsR, contrast = c('water', 'filtered', 'sterile')) <- note: functionally the same as looking at the results with contrast specified here
mcols(resW, use.names = TRUE)
summary(resW)
sum(is.na(resW)) #2112
#out of 18115 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 1845, 10%
#LFC < 0 (down)     : 1830, 10%
#outliers [1]       : 2, 0.011%
#low counts [2]     : 2108, 12%
#(mean count < 5)

resT <- results(ddsR, name = c("temp_18_vs_14"))
mcols(resT, use.names = TRUE)
summary(resT)
sum(is.na(resT)) #4
#out of 18115 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 4522, 25%
#LFC < 0 (down)     : 2996, 17%
#outliers [1]       : 2, 0.011%
#low counts [2]     : 0, 0%
#(mean count < 1)

#I (interaction) is effect in F18 compared against all the other conditions
resI <- results(ddsR, name = c("waterfiltered.temp18")) 
mcols(resI, use.names = TRUE)
summary(resI)
sum(is.na(resI)) #9486
#out of 18115 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 375, 2.1%
#LFC < 0 (down)     : 313, 1.7%
#outliers [1]       : 2, 0.011%
#low counts [2]     : 9482, 52%
#(mean count < 35)

#WI contrasts F14 and F18
resWI <- results(ddsR, list(c("water_filtered_vs_sterile", "waterfiltered.temp18")))
mcols(resWI, use.names = TRUE)
summary(resWI)
sum(is.na(resWI)) #4
#out of 18115 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 458, 2.5%
#LFC < 0 (down)     : 598, 3.3%
#outliers [1]       : 2, 0.011%
#low counts [2]     : 0, 0%
#(mean count < 1)

#TI contrasts S18 and F18 
resTI <- results(ddsR, list(c("temp_18_vs_14", "waterfiltered.temp18")))
mcols(resTI, use.names = TRUE)
summary(resTI)
sum(is.na(resTI)) #1058
#out of 18115 with nonzero total read count
#adjusted p-value < 0.1
#LFC > 0 (up)       : 4262, 24%
#LFC < 0 (down)     : 3321, 18%
#outliers [1]       : 2, 0.011%
#low counts [2]     : 1054, 5.8%
#(mean count < 3)

#save(rld.df, colData, resW, resT, resI, resWI, resTI, 
#     file = "MDRNA_rld&geneResults_July2025.Rdata")

# output CSV lines for Wald Test Results ---- 
resWdf <- as.data.frame(resW)
row.names(resWdf) <- gsub("gene-", "", row.names(resWdf))
#write.csv(resWdf, "~/Desktop/R/RNAseq_MDII/July25-reanalysis/HMRvsLMR_DEG_725.csv")

resTdf <- as.data.frame(resT)
row.names(resTdf) <- gsub("gene-", "", row.names(resTdf))
#write.csv(resTdf, "~/Desktop/R/RNAseq_MDII/July25-reanalysis/18vs14_DEG_725.csv")

resIdf <- as.data.frame(resI)
row.names(resIdf) <- gsub("gene-", "", row.names(resIdf))
#write.csv(resIdf, "~/Desktop/R/RNAseq_MDII/July25-reanalysis/interaction_DEG_725.csv")

resWIdf <- as.data.frame(resWI)
row.names(resWIdf) <- gsub("gene-", "", row.names(resWIdf))
#write.csv(resWIdf, "~/Desktop/R/RNAseq_MDII/July25-reanalysis/WI_DEG_725.csv")

resTIdf <- as.data.frame(resTI)
row.names(resTIdf) <- gsub("gene-", "", row.names(resTIdf))
#write.csv(resTIdf, "~/Desktop/R/RNAseq_MDII/July25-reanalysis/TI_DEG_725.csv")

# NOTE TO SELF: CHANGE HERE. "Running with grouping" was here. Move it to later if we're going to still do something with that ----

# Formatting rld and results files for downstream work ---- 
ll2 = load("MDRNA_rld&geneResults_July2025.Rdata") 

head(rld.df)
head(resW)
head(resT)
head(resI)
head(resWI)
head(resTI)

#collect pvalues for each gene 
vals = cbind(resW$stat, resW$pvalue, resW$padj, 
             resT$stat, resT$pvalue, resT$padj, 
             resI$stat, resI$pvalue, resI$padj, 
             resWI$stat, resWI$pvalue, resWI$padj,
             resTI$stat, resTI$pvalue, resTI$padj)

colnames(vals) = c("stat_resW", "pval_resW", "padj_resW", 
                   "stat_resT", "pval_resT", "padj_resT",
                   "stat_resI", "pval_resI", "padj_resI", 
                   "stat_resWI", "pval_resWI", "padj_resWI", 
                   "stat_resTI", "pval_resTI", "padj_resTI")
rldpvals = as.data.frame(cbind(rld.df, vals)) #combine RLD with pvalues 
row.names(rldpvals) <- gsub("gene-", "", row.names(rldpvals))
#write.csv(rldpvals, "July25-reanalysis/MDII_rldpvals_July25.csv")

# Collecting sig. pvalues from each results set 
sigpvals_W = rldpvals[rldpvals$padj_resW <= 0.05 & !is.na (rldpvals$padj_resW),]
dim(sigpvals_W) #2681   27 (note: this number differs from Volcano plots because it is only based on padj, not LFC)

sigpvals_T = rldpvals[rldpvals$padj_resT <= 0.05 & !is.na (rldpvals$padj_resT),]
dim(sigpvals_T) #5779   27

sigpvals_I = rldpvals[rldpvals$padj_resI <= 0.05 & !is.na (rldpvals$padj_resI),]
dim(sigpvals_I) #360  27 

sigpvals_WI = rldpvals[rldpvals$padj_resWI <= 0.05 & !is.na (rldpvals$padj_resWI),]
dim(sigpvals_WI) #747  27 

sigpvals_TI = rldpvals[rldpvals$padj_resTI <= 0.05 & !is.na (rldpvals$padj_resTI),]
dim(sigpvals_TI) #6183   27

# NOTE TO SELF: Venn diagrams was here. Move that elsewhere so that it can be plotted as significant in both padj and LFC---- 

# Principal Coordinates Calculation ---- 
##PCoA works off dissimilarity matrix (calculated by vegdist)

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
#      Eigenvalues Relative_eig Broken_stick Cumul_eig Cumul_br_stick
#1 63%  528.496463  0.634660302  0.274534304 0.6346603      0.2745343
#2 11%   92.200933  0.110722164  0.183625213 0.7453825      0.4581595
#3       62.288934  0.074801473  0.138170668 0.8201839      0.5963302

### Plotting PCoA 
colorsPC <- c("#336666","#336666","#336666",
              "#CC9933","#CC9933","#CC9933",
              "#33CCCC","#33CCCC","#33CCCC",
              "#FFCC00","#FFCC00","#FFCC00") 
names(colorsPC) <- c("C11F14", "C12F14", "C13F14", 
                     "C11F18", "C12F18", "C13F18",
                     "C8St14", "C9St14", "C10St14", 
                     "C8St18", "C9St18", "C10St18")
scoresplot <- as.data.frame(scores)

## plotting first and second principal axes - ggplot 

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
  xlab("PCo1 (63%)") + 
  ylab("PCo2 (11%)") + 
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

#ggsave(plot = pcoa, 
#       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/rna-pcoa-725.jpg", 
#       scale = 1, width = 6, height = 6, units = c("in"), 
#       dpi = 300)

#plotting only pcoa legend
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', 
     xlim=0:1, ylim=0:1) #creates a blank plot
plotSym <- c(16, 16, 16, 16) #shapes list
par(mai = c(0.5,0.5,0.5,0.5))
legend("center", 
       legend = c("14ºC LMR", 
                  "18ºC LMR", 
                  "14ºC HMR", 
                  "18ºC HMR"), 
       pch = plotSym, pt.cex = 5.1,
       col = c("#33CCCC", "#FFCC00", "#336666", "#CC9933"), 
       cex = 2.5, bty = 'n', par(mai = c(0.1,0.1,0.1,0.1))
) 

# plotting axes 2 and 3
centroids23 <- aggregate(cbind(Axis.2, Axis.3) ~ treat, 
                       data = scoresplot, FUN = mean)

pcoa23 <- ggplot() + 
  theme_prism() +
  geom_point(data = scoresplot, 
             aes(x = scoresplot[,2], y = scoresplot[,3], 
                 color = conditions$chamber, shape = treat), 
             size = 4) + 
  scale_color_manual(values = colorsPC) + 
  scale_shape_manual(values = treatShapes) +
  xlab("PCo2") + 
  ylab("PCo3") + 
  theme(axis.text.x = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) +
  geom_segment(data = scoresplot, aes(x = Axis.2, y = Axis.3,
                                      xend = centroids23$Axis.2[match(treat, centroids23$treat)],
                                      yend = centroids23$Axis.3[match(treat, centroids23$treat)],
                                      color = conditions$chamber)) +
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))
print(pcoa23)


# VOLCANO PLOTS ---- 
library('ggrepel')
library('ggprism')
library('data.table')

# setting up dataframes for plotting 
resWsig <- subset(resW, resW$padj <= 0.05)
min(abs(resWsig$log2FoldChange)) #0.4126948
max(abs(resWsig$log2FoldChange)) #6.659325
sum(resWsig$padj <= 0.5 & (abs(resWsig$log2FoldChange) > 0.5)) #2678
sum(resWsig$log2FoldChange > 0.5) #1478 
sum(resWsig$log2FoldChange < -0.5) #1200
plotW <- na.omit(resW) #removes if NA value; NAs occur in pval and padj if B/A didn't have a divisible number

resTsig <- subset(resT, resT$padj <= 0.05)
min(abs(resTsig$log2FoldChange)) #0.2855876
max(abs(resTsig$log2FoldChange)) #7.250701
sum(resTsig$padj <= 0.5 & (abs(resTsig$log2FoldChange) > 0.5)) #5607
sum(resTsig$log2FoldChange > 0.5) #3351 
sum(resTsig$log2FoldChange < -0.5) #2256
plotT <- na.omit(resT) 

resIsig <- subset(resI, resI$padj <= 0.05)
min(abs(resIsig$log2FoldChange)) #0.59097
max(abs(resIsig$log2FoldChange)) #4.410615
sum(resIsig$padj <= 0.5 & (abs(resIsig$log2FoldChange) > 0.5)) #360
sum(resIsig$log2FoldChange > 0.5) #173 
sum(resIsig$log2FoldChange < -0.5) #187
plotI <- na.omit(resI)

resWIsig <- subset(resWI, resWI$padj <= 0.05)
min(abs(resWIsig$log2FoldChange)) #0.3665892
max(abs(resWIsig$log2FoldChange)) #6.023607
sum(resWIsig$padj <= 0.5 & (abs(resWIsig$log2FoldChange) > 0.5)) #716
sum(resWIsig$log2FoldChange > 0.5) #311 
sum(resWIsig$log2FoldChange < -0.5) #405
plotWI <- na.omit(resWI)

resTIsig <- subset(resTI, resTI$padj <= 0.05)
min(abs(resTIsig$log2FoldChange)) #0.3964666
max(abs(resTIsig$log2FoldChange)) #7.350015
sum(resTIsig$padj <= 0.5 & (abs(resTIsig$log2FoldChange) > 0.5)) #6170
sum(resTIsig$log2FoldChange > 0.5) #3294 
sum(resTIsig$log2FoldChange < -0.5) #2876
plotTI <- na.omit(resTI)

# Colors for volcano plots
de_colors <- c("blue", "coral2", "azure4")
names(de_colors) <- c("Down", "Up", "No change")

# WATER
plotW$DE <- "No change"
plotW$DE[plotW$log2FoldChange > 0.5 & plotW$padj < 0.05] <- "Up"
plotW$DE[plotW$log2FoldChange < -0.5 & plotW$padj < 0.05] <- "Down"
plotW$geneId <- gsub("gene-", "", row.names(plotW)) 

water_volcano <- ggplot(data = plotW, 
                        aes(x = log2FoldChange, y = -log10(padj), col = DE)) + 
  geom_point() + scale_color_manual(values = de_colors) + 
  geom_vline(xintercept = c(-0.5, 0.5), col = "brown1") +
  geom_hline(yintercept = -log10(0.05), col = "brown1") + 
  xlab(expression(bold("log"[2]~Fold~Change))) +
  ylab(expression(bold("-log"[10]~"Adjusted p-value"))) +
  theme_prism() + 
  xlim(-9, 9) + 
  ylim(0, 40) +
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

#ggsave(plot = water_volcano, 
#       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/rna-water-volcano725.jpg", 
#       scale = 1, width = 6, height = 6, units = c("in"), 
#       dpi = 300)

# TEMPERATURE
plotT$DE <- "No change"
plotT$DE[plotT$log2FoldChange > 0.5 & plotT$padj < 0.05] <- "Up"
plotT$DE[plotT$log2FoldChange < -0.5 & plotT$padj < 0.05] <- "Down"
plotT$names <- gsub("gene-", "", row.names(plotT))

temp_volcano <- ggplot(data = plotT, 
                       aes(x = log2FoldChange, y = -log10(padj), col = DE)) + 
  geom_point() + scale_color_manual(values = de_colors) + 
  geom_vline(xintercept = c(-0.5, 0.5), col = "brown1") +
  geom_hline(yintercept = -log10(0.05), col = "brown1") + 
  xlab(expression(bold("log"[2]~Fold~Change))) +
  ylab(expression(bold("-log"[10]~"Adjusted p-value"))) +
  theme_prism() + 
  xlim(-9, 9) +
  ylim(0, 40) +
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

#ggsave(plot = temp_volcano, 
#       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/rna-temp-volcano725.jpg", 
#       scale = 1, width = 6, height = 6, units = c("in"), 
#       dpi = 300)

# INTERACTION
plotI$DE <- "No change"
plotI$DE[plotI$log2FoldChange > 0.5 & plotI$padj < 0.05] <- "Up"
plotI$DE[plotI$log2FoldChange < -0.5 & plotI$padj < 0.05] <- "Down"
plotI$names <- gsub("gene-", "", row.names(plotI))

interaction_volcano <- ggplot(data = plotI, 
                       aes(x = log2FoldChange, y = -log10(padj), col = DE)) + 
  geom_point() + scale_color_manual(values = de_colors) + 
  geom_vline(xintercept = c(-0.5, 0.5), col = "brown1") +
  geom_hline(yintercept = -log10(0.05), col = "brown1") + 
  xlab(expression(bold("log"[2]~Fold~Change))) +
  ylab(expression(bold("-log"[10]~"Adjusted p-value"))) +
  theme_prism() + 
  xlim(-9, 9) +
  ylim(0, 40) +
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
print(interaction_volcano)

#ggsave(plot = interaction_volcano, 
#       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/rna-interaction-volcano725.jpg", 
#       scale = 1, width = 6, height = 6, units = c("in"), 
#       dpi = 300)

# WI plot 
plotWI$DE <- "No change"
plotWI$DE[plotWI$log2FoldChange > 0.5 & plotWI$padj < 0.05] <- "Up"
plotWI$DE[plotWI$log2FoldChange < -0.5 & plotWI$padj < 0.05] <- "Down"
plotWI$names <- gsub("gene-", "", row.names(plotWI))

WI_volcano <- ggplot(data = plotWI, 
                              aes(x = log2FoldChange, y = -log10(padj), col = DE)) + 
  geom_point() + scale_color_manual(values = de_colors) + 
  geom_vline(xintercept = c(-0.5, 0.5), col = "brown1") +
  geom_hline(yintercept = -log10(0.05), col = "brown1") + 
  xlab(expression(bold("log"[2]~Fold~Change))) +
  ylab(expression(bold("-log"[10]~"Adjusted p-value"))) +
  theme_prism() + 
  xlim(-9, 9) +
  ylim(0, 40) +
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
print(WI_volcano) 

#ggsave(plot = WI_volcano, 
#       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/rna-WI-volcano725.jpg", 
#       scale = 1, width = 6, height = 6, units = c("in"), 
#       dpi = 300)

# TI plot 
plotTI$DE <- "No change"
plotTI$DE[plotTI$log2FoldChange > 0.5 & plotTI$padj < 0.05] <- "Up"
plotTI$DE[plotTI$log2FoldChange < -0.5 & plotTI$padj < 0.05] <- "Down"
plotTI$names <- gsub("gene-", "", row.names(plotTI))

TI_volcano <- ggplot(data = plotTI, 
                     aes(x = log2FoldChange, y = -log10(padj), col = DE)) + 
  geom_point() + scale_color_manual(values = de_colors) + 
  geom_vline(xintercept = c(-0.5, 0.5), col = "brown1") +
  geom_hline(yintercept = -log10(0.05), col = "brown1") + 
  xlab(expression(bold("log"[2]~Fold~Change))) +
  ylab(expression(bold("-log"[10]~"Adjusted p-value"))) +
  theme_prism() + 
  xlim(-9, 9) +
  ylim(0, 40) +
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
print(TI_volcano)

#ggsave(plot = TI_volcano, 
#       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/rna-TI-volcano725.jpg", 
#       scale = 1, width = 6, height = 6, units = c("in"), 
#       dpi = 300)

# Venn diagrams (significance based on both LFC & Padj) ----
library('VennDiagram')
plotW = as.data.frame(plotW)
plotT = as.data.frame(plotT)
wat = row.names(plotW[plotW$padj<0.05 & (abs(plotW$log2FoldChange) > 0.5),])
temp = row.names(plotT[plotT $padj < 0.05 & (abs(plotT$log2FoldChange) > 0.5),])
degs05 = union(wat,temp)
length(degs05) #7259
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

#NOTE: Leaving off at adding all the gene name conversion information  ---- 

# Fold change and Fisher's test data for GO-MWU ---- 

### WATER
row.names(resW) <- sub("gene-", "", row.names(resW))
water_stat_RNA = data.frame(cbind("geneId" = row.names(resW), "stat" = resW$stat))
write.csv(water_stat_RNA, 
          file = "~/Desktop/R/MDII_GO_MWU/July25-reanalysis/water_waldStat_MDRNA_725.csv", 
          quote = F, row.names = F)

### TEMP
row.names(resT) <- sub("gene-", "", row.names(resT))
temp_stat_RNA = data.frame(cbind("geneId" = row.names(resT), "stat" = resT$stat))
write.csv(temp_stat_RNA, 
          file = "~/Desktop/R/MDII_GO_MWU/July25-reanalysis/temp_waldStat_MDRNA_725.csv", 
          quote = F, row.names = F)

### Interaction
row.names(resI) <- sub("gene-", "", row.names(resI))
interaction_stat_RNA = data.frame(cbind("geneId" = row.names(resI), "stat" = resI$stat))
write.csv(interaction_stat_RNA, 
          file = "~/Desktop/R/MDII_GO_MWU/July25-reanalysis/interaction_waldStat_MDRNA_725.csv", 
          quote = F, row.names = F)

### WI
row.names(resWI) <- sub("gene-", "", row.names(resWI))
WI_stat_RNA = data.frame(cbind("geneId" = row.names(resWI), "stat" = resWI$stat))
write.csv(WI_stat_RNA, 
          file = "~/Desktop/R/MDII_GO_MWU/July25-reanalysis/WI_waldStat_MDRNA_725.csv", 
          quote = F, row.names = F)

### TI
row.names(resTI) <- sub("gene-", "", row.names(resTI))
TI_stat_RNA = data.frame(cbind("geneId" = row.names(resTI), "stat" = resTI$stat))
write.csv(TI_stat_RNA, 
          file = "~/Desktop/R/MDII_GO_MWU/July25-reanalysis/TI_waldStat_MDRNA_725.csv", 
          quote = F, row.names = F)

#Plotting LFC of all DEG (scatterplot) ---- 
RW <- read.csv("July25-reanalysis/HMRvsLMR_DEG_725.csv") 
RW <- subset(RW, select = -c(baseMean, lfcSE, stat, pvalue)) 
RW$DE <- "No change"
RW$DE[RW$log2FoldChange > 0.5 & RW$padj < 0.05] <- "Up"
RW$DE[RW$log2FoldChange < -0.5 & RW$padj < 0.05] <- "Down"

RT <- read.csv("July25-reanalysis/18vs14_DEG_725.csv")
RT <- subset(RT, select = -c(baseMean, lfcSE, stat, pvalue)) 
RT$DE <- "No change"
RT$DE[RT$log2FoldChange > 0.5 & RT$padj < 0.05] <- "Up"
RT$DE[RT$log2FoldChange < -0.5 & RT$padj < 0.05] <- "Down"

DEG <- merge(RW, RT, by = "X", all = TRUE)
colnames(DEG) <- c("names", "lfc_W", "padj_W", 
                   "DE_W", "lfc_temp", 
                   "padj_temp", "DE_temp")
DEG$lfc_W[is.na(DEG$lfc_W)] <- 0
DEG$lfc_W[is.na(DEG$lfc_W)] <- "No change"
DEG$DE_W[is.na(DEG$DE_W)] <- 0
DEG$DE_temp[is.na(DEG$DE_temp)] <- "No change"
DEG$DE_W <- sub("No change", "Nc", DEG$DE_W)
DEG$DE_temp <- sub("No change", "Nc", DEG$DE_temp)
DEG <- unite(DEG, DE_WT, "DE_W", "DE_temp", sep = "-")

dewtColors <- c("#CC0033", "#33CCCC", "#CC0033", 
                "azure4", "#FFCC00", "azure4", "#FFCC00",
                "#CC0033", "#33CCCC", "#CC0033")
names(dewtColors) <- c("Down-Down", "Down-Nc", "Down-Up", 
                       "NA-Nc", "Nc-Down", 
                       "Nc-Nc", "Nc-Up", 
                       "Up-Down", "Up-Nc", "Up-Up")

library('ggrepel')
min(DEG$lfc_temp) #-7.250701
max(DEG$lfc_temp) #4.496332
min(DEG$lfc_W) #-0.000159411167426084
max(DEG$lfc_W) #6.65932475617044
DEG$lfc_W <- as.numeric(DEG$lfc_W)

degPlot <- ggplot() + 
  theme_prism() + 
  ylim(-8, 5) +
  xlim(-6, 6) +
  geom_vline(xintercept = c(0, 0), col = "red") +
  geom_hline(yintercept = c(0, 0), col = "red") +
  geom_point(data = subset(DEG, DE_WT == "Nc-Nc"),
             aes(x = lfc_W, y = lfc_temp, color = DE_WT)) + 
  geom_point(data = subset(DEG, !(DE_WT == "Nc-Nc")),
             aes(x = lfc_W, y = lfc_temp, color = DE_WT)) + 
  geom_rect(data = DEG, 
            aes(xmin = 0 - 0.5, 
                xmax = 0 + 0.5, 
                ymin = 0 - 0.5, 
                ymax = 0 + 0.5), 
            fill = "transparent", color = "purple", size = 0.2) +
  #geom_text_repel(data = subset(DEG, ((DE_WT == "Down-Down" |
  #                                      DE_WT == "Up-Up"| 
  #                                      DE_WT == "Up-Down"| 
  #                                      DE_WT == "Down-Up"|
  #                                      abs(lfc_MBE) > 1 | 
  #                                  abs(lfc_temp) > 2) & 
  #                                !(names %like% "LOC"))), 
  #                aes(x = lfc_MBE, y = lfc_temp, label = names), 
  #                size = 3.5, color = "black", max.overlaps = Inf) + 
  xlab(expression(bold("Microbial Richness"~"-"~"log"[2]~Fold~Change))) +
  ylab(expression(bold(Temp~"-"~"log"[2]~Fold~Change))) +
  scale_color_manual("DEG expression", values = dewtColors) +
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
plot(degPlot)

ggsave(plot = degPlot, 
       filename = "~/Desktop/R/RNAseq_MDII/July25-reanalysis/deg-scatter-725.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

table(DEG$DE_WT)
#Down-Down   Down-Nc   Down-Up   Nc-Down     Nc-Nc     Nc-Up   Up-Down     Up-Nc 
#      110       783       307      1586     10856      2995       560       869 
#Up-Up 
#   49
