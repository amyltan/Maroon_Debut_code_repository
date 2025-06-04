setwd("~/Desktop/R/atac_seq_MDII/DAinDB-ChIPseeker/")

library("DiffBind")
library("csaw")

### Load saved analysis of normalized count data: 
multiDBA <- dba.load(file = 'ATACseq-MDII-analysis-wContrasts', 
                     dir = '.', pre = 'dba_', ext = 'Rdata')

### To save changes to analysis step: 
#dba.save(multiDBA, file = 'ATACseq-MDII-analysis-wContrasts', 
#         dir = '.', pre = 'dba_', ext = 'Rdata', bRemoveAnalysis = FALSE)

### Read in sample list containing file names and locations
samples <- read.csv("v5-int25bampe-sampleList.csv")
MDII <- dba(sampleSheet = samples)

counts <- dba.count(MDII)  #note: this step takes a long time
counts
#12 Samples, 69932 sites in matrix: 

countinfo <- dba.show(counts)
libsizes <- cbind(LibReads=countinfo$Reads, FRiP=countinfo$FRiP,
                  PeakReads=round(countinfo$Reads * countinfo$FRiP))
rownames(libsizes) <- countinfo$ID
countinfo

#ID Condition Treatment Caller Intervals    Reads FRiP
#1   C8St14        14   Sterile counts     69932  8671800 0.25
#2   C9St14        14   Sterile counts     69932  7453326 0.26
#3  C10St14        14   Sterile counts     69932 10325132 0.25
#4   C11F14        14  Filtered counts     69932  1499280 0.29
#5   C12F14        14  Filtered counts     69932  4826209 0.31
#6   C13F14        14  Filtered counts     69932 10607227 0.27
#7   C8St18        18   Sterile counts     69932  5180118 0.25
#8   C9St18        18   Sterile counts     69932  8319984 0.26
#9  C10St18        18   Sterile counts     69932  7423490 0.24
#10  C11F18        18  Filtered counts     69932  6174982 0.27
#11  C12F18        18  Filtered counts     69932  7126832 0.28
#12  C13F18        18  Filtered counts     69932 10033066 0.25
 
### Normalize by RLE on background windows 
normcounts <- dba.normalize(counts, method=DBA_DESEQ2, normalize=DBA_NORM_NATIVE, 
                            library=DBA_LIBSIZE_BACKGROUND, background=TRUE)
norm <- dba.normalize(normcounts, bRetrieve = TRUE)
norm 

#$norm.factors
#C8St14    C9St14   C10St14    C11F14    C12F14    C13F14    C8St18    C9St18   C10St18 
#1.4147111 1.1475333 1.6808231 0.1953563 0.6500328 1.6955338 0.7964571 1.2828065 1.1486505 
#C11F18    C12F18    C13F18 
#0.9394497 1.0593054 1.6157561 

#creating table of normalized reads for PCoA (note: table created another way at the end, but results are the same. That one has been exported)
norm_ret <- dba.count(normcounts, peaks=NULL, score=DBA_SCORE_NORMALIZED)
norm.table <- dba.peakset(norm_ret, bRetrieve = TRUE) 
norm.table <- data.frame(norm.table)
library('tidyr')
library('vegan')
library('ape')
norm.table <- unite(norm.table, id_start_end, "seqnames", "start", "end", sep = "-")
norm.table <- subset(norm.table, select = -c(width, strand))
tmp <- as.data.frame(t(norm.table[,-1]))
colnames(tmp) <- norm.table$id_start_end

veg.tmp <- vegdist(tmp, "manhattan")
veg.tmp.div = veg.tmp/1000
norm.pcoa = pcoa(veg.tmp.div) 
head(norm.pcoa)
scores = norm.pcoa$vectors

#$values
#   Eigenvalues Relative_eig Broken_stick Cumul_eig Cumul_br_stick
#1   1225013.43   0.47690197  0.274534304 0.4769020      0.2745343
#2    639282.51   0.24887489  0.183625213 0.7257769      0.4581595
#3    140707.73   0.05477801  0.138170668 0.7805549      0.5963302

colData <- data.frame(row.names = rownames(scores), 
                      chamber = c("C8St14","C9St14","C10St14","C8St18","C9St18",
                                  "C10St18","C11F14","C12F14",
                                  "C13F14","C11F18","C12F18","C13F18" ))
library('dplyr')
library('stringr')
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

conditions = colData
conditions$water = as.factor(conditions$water) 
conditions$temp = as.factor(conditions$temp) 
conditions$treat = as.factor(conditions$treat)

colorsW <- c("#336666","#CC9933","#99FFFF","#FFCC00")
conditions$treat_c <- colorsW[as.factor(conditions$treat)]
legendNames <- c("14ºC microbiome-exposed", "18ºC microbiome-exposed", 
                 "14ºC sterile ASW", "18ºC sterile ASW")

par(mai = c(0.8,0.9,0.1,0.1))
plot(scores[,1], scores[,2], col = conditions$treat_c, pch = 19, cex = 3,
     xlab = "PCo1", ylab = "PCo2", cex.axis = 1.2, cex.lab = 1.9, 
     xlim = c(-400,900), ylim = c(-500, 700))
ordispider(scores, conditions$treat, label = F) 
#ordiellipse(ord = scores, groups = conditions$treat, 
#            conf = 0.95)
legend(legend = c("14ºC MBE", "18ºC MBE", 
                  "14ºC SASW", "18ºC SASW"), 
       col = colorsW, pch = 19, cex = 1.3, bty = "n",
       xlim = c(-400,900), ylim = c(-500, 700)) 

# trying to plot PCo with ggplot ####
library('ggplot2')
library('ggprism')
library('ggh4x')

#to redo without only this plot
#scoresplot <- read.csv("pcoa-atac-scores.csv")


colorsPC <- c("#336666","#336666","#336666",
             "#CC9933","#CC9933","#CC9933",
             "#33CCCC","#33CCCC","#33CCCC",
             "#FFCC00","#FFCC00","#FFCC00") 
names(colorsPC) <- c("C11F14", "C12F14", "C13F14", 
                   "C11F18", "C12F18", "C13F18",
                   "C8St14", "C9St14", "C10St14", 
                   "C8St18", "C9St18", "C10St18")

#scoresplot <- as.data.frame(scores)
#calculate centroids for each treatment group 
scoresplot$SampleID <- rownames(scoresplot)
scoresplot$treat <- colData$treat 
centroids <- aggregate(cbind(Axis.1, Axis.2) ~ treat, 
                       data = scoresplot, FUN = mean)
#write.csv(scoresplot, "pcoa-atac-scores.csv")

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
  xlab("PCoA1 (48%)") + 
  ylab("PCoA2 (25%)") + 
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
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/atac-pcoa.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

#plotting only pcoa legend
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', 
     xlim=0:1, ylim=0:1) #creates a blank plot
plotSym <- c(16, 17, 16, 17) #shapes list
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


#####
# to make a nice table 
norm_table <- cbind(FullLibSize = norm$lib.sizes, 
              NormFacts = norm$norm.factors, 
              NormLibSize = round(norm$lib.sizes/norm$norm.factors))

#       FullLibSize NormFacts NormLibSize
#C8St14      8639019 1.4147111     6106561
#C9St14      7125890 1.1475333     6209746
#C10St14    10116668 1.6808231     6018877
#C11F14      1321487 0.1953563     6764498
#C12F14      4382175 0.6500328     6741467
#C13F14     10551786 1.6955338     6223283
#C8St18      5029423 0.7964571     6314745
#C9St18      7949110 1.2828065     6196656
#C10St18     6968301 1.1486505     6066511
#C11F18      5948116 0.9394497     6331490
#C12F18      6694258 1.0593054     6319479
#C13F18      9954157 1.6157561     6160681

### Multifactorial differential accessibility analysis 
normContrasts <- dba.contrast(normcounts, design = "~ Condition + Treatment")
multiDBA <- dba.analyze(normContrasts) 

dba.show(multiDBA, bContrasts = TRUE)
#     Factor    Group Samples  Group2 Samples2 DB.DESeq2
#1 Condition       14       6      18        6       134
#2 Treatment Filtered       6 Sterile        6      4816

### To save any changes in the analysis: 
#dba.save(multiDBA, file = 'ATACseq-MDII-analysis-wContrasts', 
#         dir = '.', pre = 'dba_', ext = 'Rdata', bRemoveAnalysis = FALSE)

### Load saved analysis of normalized count data: 
#multiDBA <- dba.load(file = 'ATACseq-MDII-analysis-wContrasts', 
#                     dir = '.', pre = 'dba_', ext = 'Rdata')

dba.plotMA(multiDBA, contrast = 1, bFlip = TRUE)
dba.plotMA(multiDBA, contrast = 2)

## Retrieving reports for differentially accessible sites based on condition

# Temperature: 18 vs. 14
temp_report <- dba.report(multiDBA, contrast = 1, bFlip = TRUE, th = 1) 
#th = 1 reports all sites 

sum(temp_report$FDR<0.05) #significant: 134 
sum(temp_report$FDR<0.05 & temp_report$Fold>0) #3 sites more accessible
sum(temp_report$FDR<0.05 & temp_report$Fold<0) #131 sites less accessible 

temp_sig = temp_report[temp_report$FDR<0.05] #significant results only
#checking min and max significant fold changes
min(abs(temp_sig$Fold)) #0.1950321 
max(abs(temp_sig$Fold)) #0.4724931 

# Microbial exposure: Filtered vs. Sterile
water_report <- dba.report(multiDBA, contrast = 2, bFlip = FALSE, th = 1)
sum(water_report$FDR<0.05) #significant: 4816 
sum(water_report$FDR<0.05 & water_report$Fold>0) #4668 sites more accessible
sum(water_report$FDR<0.05 & water_report$Fold<0) #148 sites less accessible 

water_sig = water_report[water_report$FDR<0.05] #significant results only
#checking min and max significant fold changes
min(abs(water_sig$Fold)) #0.2546238 
max(abs(water_sig$Fold)) #0.8442567

#####
## Exporting for ChIPseeker 
temp_DF = data.frame(temp_report)
write.table(temp_DF, "tempReport_ATAC_DBanalysis.txt", quote = FALSE)
water_DF = data.frame(water_report)
write.table(water_DF, "waterReport_ATAC_DBanalysis.txt", quote = FALSE)

## Exporting BEDs 
library("rtracklayer")
export.bed(temp_report, "temp_ATAC_allSites.bed") 
export.bed(temp_sig, "temp_ATAC_sig.bed")
export.bed(water_report, "water_ATAC_allSites.bed")
export.bed(water_sig, "water_ATAC_sig.bed")

temp.UP = temp_sig[temp_sig$Fold>0]
temp.DOWN = temp_sig[temp_sig$Fold<0]
water.UP = water_sig[water_sig$Fold>0]
water.DOWN = water_sig[water_sig$Fold<0]

export.bed(temp.UP, "temp_ATAC_UP.bed")
export.bed(temp.DOWN, "temp_ATAC_DOWN.bed")
export.bed(water.UP, "water_ATAC_UP.bed")
export.bed(water.DOWN, "water_ATAC_DOWN.bed")

##### 
#Creating dataframe of normalized reads by each chamber and matching up to gene names 
# for output to mixOmics
 
normOut <- dba.peakset(normcounts, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
head(normOut)
library('rtracklayer')
export.bed(normOut, "normReads_MDATAC.bed") #bed has CHR-START-END locations but no read info

library('ChIPseeker')
library('GenomicFeatures')
library('GenomicRanges')
library('tidyr')
spurp_txdb = loadDb("~/Desktop/R/atac_seq_MDII/DAinDB-ChIPseeker/Spurp5_TxDb")
ll = load("~/Desktop/R/atac_seq_MDII/DAinDB-ChIPseeker/Spurp5_txdb_Annotations.Rdata")

#read in bed file to ChIPseeker and add geneId annotations to peak locations
# hashmarked lines add single letter feature code. Unhash these lines to use them.
normIn <- readPeakFile("normReads_MDATAC.bed")
normAnno <- annotatePeak(normIn, tssRegion = c(-3000, 3000), TxDb = spurp_txdb)
annoDF <- data.frame(normAnno)
#annoDF$feature <- as.factor(str_sub(annoDF$annotation, 1, 1))
head(annoDF)
annoDF <- subset(annoDF, select = -c(width, strand, V4, V5, V6, annotation, 
                                     geneChr, geneStart, geneEnd, geneLength, 
                                     geneStrand, transcriptId, distanceToTSS))
annoDF <- unite(annoDF, "id-start-end", seqnames, start, end)
normOut <- unite(normOut, "id-start-end", CHR, START, END)

#Combine normalized read count dataframe with annotated peak locations
normReadsAnno <- merge(annoDF, normOut, by = "id-start-end", all = TRUE)
normReadsAnno <- subset(normReadsAnno, select = -c(`id-start-end`))
#normReadsAnno <- unite(normReadsAnno, "geneId", geneId, feature)
write.csv(normReadsAnno, "~/Desktop/R/MDII_mixOmics/MDATAC_normalizedReadsPerGene_wFeature.csv", 
          row.names = FALSE)
