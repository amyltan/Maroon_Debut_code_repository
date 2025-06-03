setwd("~/Desktop/R/atac_seq_MDII/DAinDB-ChIPseeker/")

library('stringr')
library('plyr')
library('ChIPseeker')
library('GenomicFeatures')
library('GenomicRanges')
library('Repitools')

#load Spurp GFF
#spurp_txdb <- makeTxDbFromGFF(file = "~/Desktop/R/atac_seq/peak-counts/sp5_0_GCF_new.gff3", 
#                              organism = "Strongylocentrotus purpuratus")


#partition Spurp genomic features 
#promoter <- getPromoters(TxDb = spurp_txdb,
#                         upstream = 3000, downstream = 3000)
#introns <- getBioRegion(TxDb = spurp_txdb, by = "intron")
#exons <- getBioRegion(TxDb = spurp_txdb, by = "exon")
#genes <- getBioRegion(TxDb = spurp_txdb, by = "gene")

#save 
#save(promoter, introns, exons, genes, file = "Spurp5_txdb_Annotations.Rdata")
#saveDb(spurp_txdb, file = "Spurp5_TxDb")

#load 
spurp_txdb = loadDb("Spurp5_TxDb")
ll = load("Spurp5_txdb_Annotations.Rdata")
ll2 = load("temp_AllPeakAnnoGR_ATAC-MDII.Rdata")
ll3 = load("TempATAC_peakAnnoAll_df.Rdata") 
ll4 = load("temp_SigPeakAnno_ATAC-MDII.Rdata")
ll5 = load("TempATAC_peakAnno_Sig_df.Rdata")

#Read in BED files from DiffBind results 
peaks <- readPeakFile("bed-track-files/temp_ATAC_allSites.bed") 
peaks_df <- data.frame(peaks)
head(peaks_df)

peakSig <- readPeakFile("bed-track-files/temp_ATAC_sig.bed")
peakSig_df <- data.frame(peakSig)
head(peakSig_df) 

#Color list for figures 
fullregionColList = c("#FF0000", #promoter <1kb 
                      "#CC0000", #promoter 1-2kb 
                      "#993300", #promoter 2-3kb 
                      "#660033", #5' UTR 
                      "#3399FF", #3' UTR 
                      "#FF9933", #1st exon
                      "#FFCC33", #other exon
                      "#006600", #1st intron
                      "#003300", #other intron 
                      "#CC00CC", #Downstream 
                      "#99FFCC" #distal intergenic
)

#Annotate peaks 
peakAnno_All <- annotatePeak(peaks, tssRegion = c(-3000,3000), TxDb = spurp_txdb) #note: gave extensive errors
plotAnnoPie(peakAnno_All, cex = 0.6, radius = 0.8, col = fullregionColList)
peakAnno_GR <- as.GRanges(peakAnno_All)
peakAnno_All_df <- data.frame(peakAnno_GR)
#save(peakAnno_GR, file = "temp_AllPeakAnnoGR_ATAC-MDII.Rdata")
#write.csv(peakAnno_All_df, "temp_peakAnno_AllSites.csv", quote = TRUE)

peakAnno_Sig <- annotatePeak(peakSig, tssRegion = c(-3000,3000), TxDb = spurp_txdb)
peakAnnoSig_GR <- as.GRanges(peakAnno_Sig)
peakAnno_Sig_df <- data.frame(peakAnnoSig_GR)
#save(peakAnnoSig_GR, file = "temp_SigPeakAnnoGR_ATAC-MDII.Rdata")
#write.csv(peakAnno_Sig_df, "temp_peakAnno_SigSites.csv", quote = TRUE)

## Combining DA from DiffBind with annotations from ChIPseeker
library('tidyr')

anno <- read.csv("temp_peakAnno_AllSites.csv")
anno = subset(anno, select = -c(X, width, strand, V4, V5, V6))
anno <- unite(anno, id_start_end, "seqnames", "start", "end", sep = "-")

DAtable <-read.table("tempReport_ATAC_DBanalysis.txt")
DAtable <- unite(DAtable, id_start_end, "seqnames", "start", "end", sep = "-")

DAanno <- merge(anno, DAtable, by = c("id_start_end"), all = TRUE)

#write.csv(DAanno, "~/Desktop/R/atac_seq_MDII/Temp_DAR_ATACseq.csv")

## Marking up Sig peaks
### Peaks +/- 500 bp from TSS 
tss_Sig <-  peakAnnoSig_GR[abs(peakAnno_Sig_df$distanceToTSS) <500]
tss_Sig_df = data.frame(tss_Sig)
tss_Sig_df$feature <- 'P' #add feature line to TSS peaks data frame

length(unique(tss_Sig$transcriptId)) #Filter for unique: temp - 10 genes 

peakCov_tss <- count(tss_Sig_df$geneId) #Sum counts per TSS by geneId 
names(peakCov_tss) = c("common_name", "peak_density_TSS")
head(peakCov_tss)

### Wrangling rest of peak information
#Clean up annotation column so that single letter codes can be added to feature 
peakAnno_Sig_df$annotation <- sub(" ", "", peakAnno_Sig_df$annotation)
peakAnno_Sig_df$annotation <- as.factor(peakAnno_Sig_df$annotation)
peakAnno_Sig_df$feature <- as.factor(str_sub(peakAnno_Sig_df$annotation, 1, 1))

#Peaks in INTRONS
intron_peaks <- peakAnno_Sig_df[peakAnno_Sig_df$feature == "I",]
intron_peaks$annotation <- sub(" ", "", intron_peaks$annotation)
intron1_peaks = intron_peaks[grep("intron 1 of", intron_peaks$annotation),] #first intron peaks

#Count genes with peaks in introns and first introns 
length(unique(intron_peaks$geneId)) #temp - 37
length(unique(intron1_peaks$geneId)) #temp - 15 

intron_peaks$geneId = as.factor(intron_peaks$geneId)
intronCov = count(intron_peaks$geneId)
names(intronCov) = c("common_name", "peak_density_Introns")

intron1_peaks$geneId = as.factor(intron1_peaks$geneId)
intron1Cov = count(intron1_peaks$geneId)
names(intron1Cov) = c("common_name", "peak_density_Intron1")

intronsCoverage = merge(intronCov, intron1Cov, by = "common_name", all = TRUE)

#Peaks in EXONS 
exon_peaks <- peakAnno_Sig_df[peakAnno_Sig_df$feature == "E",]
exon_peaks$annotation <- sub(" ", "", exon_peaks$annotation)
exon_peaks$geneId = as.factor(exon_peaks$geneId)
exonCov = count(exon_peaks$geneId)
names(exonCov) = c("common_name", "peak_density_Exons")
length(unique(exon_peaks$geneId)) #temp - 0 

# ChIPseeker graphs of location of significant DA sites 
sigTregionColList = c("#FF0000", #promoter <1kb 
                      "#CC0000", #promoter 1-2kb 
                      "#993300", #promoter 2-3kb 
                       
                      "#3399FF", #3' UTR 
                      
                      
                      "#006600", #1st intron
                      "#003300", #other intron 
                       
                      "#99FFCC" #distal intergenic
)
plotAnnoBar(peakAnno_Sig, y = "18 vs. 14", title = "")
plotAnnoPie(peakAnno_Sig, cex = 1, radius = 0.8, col = sigTregionColList)
#note: have a separate build of the pie chart based on the numbers generated here 
#easier to manipulate in ggplot
#see script atac_sigpeaks_pie.R
plotDistToTSS(peakAnno_Sig, y = "Distance to TSS - temperature", title = "")

##save sig. diff. analysis (with feature codes)
#save(peakAnno_Sig_df, file = "TempATAC_peakAnno_Sig_df.Rdata")


## Marking up ALL peaks 
#note to self: start here, copying code from above but now for .bed of all peaks 
tss_all <- peakAnno_GR[abs(peakAnno_All_df$distanceToTSS) <500]
tss_all_df = data.frame(tss_all)
tss_all_df$feature <- 'P' #add feature line to TSS peaks data frame

length(unique(tss_all$transcriptId)) #temp - 11381

peakCovAll_tss <- count(tss_all_df$geneId) #Sum counts per TSS by geneId 
names(peakCovAll_tss) = c("common_name", "peak_density_TSS")
head(peakCovAll_tss)

### Wrangling rest of peak information
#Clean up annotation column so that single letter codes can be added to feature 
peakAnno_All_df$annotation <- sub(" ", "", peakAnno_All_df$annotation)
peakAnno_All_df$annotation <- as.factor(peakAnno_All_df$annotation)
peakAnno_All_df$feature <- as.factor(str_sub(peakAnno_All_df$annotation, 1, 1))

#Peaks in INTRONS
intron_peaksAll <- peakAnno_All_df[peakAnno_All_df$feature == "I",]
intron_peaksAll$annotation <- sub(" ", "", intron_peaksAll$annotation)
intron1_peaksAll = intron_peaksAll[grep("intron 1 of", intron_peaksAll$annotation),] #first intron peaks

#Count genes with peaks in introns and first introns 
length(unique(intron_peaksAll$geneId)) #temp - 5834
length(unique(intron1_peaksAll$geneId)) #temp - 2200 

intron_peaksAll$geneId = as.factor(intron_peaksAll$geneId)
intronCovAll = count(intron_peaksAll$geneId)
names(intronCovAll) = c("common_name", "peak_density_Introns")

intron1_peaksAll$geneId = as.factor(intron1_peaksAll$geneId)
intron1CovAll = count(intron1_peaksAll$geneId)
names(intron1CovAll) = c("common_name", "peak_density_Intron1")

intronsCoverageAll = merge(intronCovAll, intron1CovAll, by = "common_name", all = TRUE)

#Peaks in EXONS 
exon_peaksAll <- peakAnno_All_df[peakAnno_All_df$feature == "E",]
exon_peaksAll$annotation <- sub(" ", "", exon_peaksAll$annotation)
exon_peaksAll$geneId = as.factor(exon_peaksAll$geneId)
exonCovAll = count(exon_peaksAll$geneId)
names(exonCovAll) = c("common_name", "peak_density_Exons")
length(unique(exon_peaksAll$geneId)) #temp - 5754

# ChIPseeker graphs of location of ALL DA sites (significant or not)
plotAnnoBar(peakAnno_All, y = "18 vs. 14", title = "")
plotDistToTSS(peakAnno_All, y = "Distance to TSS - temperature", title = "")

#save all peak analysis (with feature codes)
#save(peakAnno_All_df, file = "TempATAC_peakAnnoAll_df.Rdata")

#####
## Volcano plots 
DAanno <- read.csv("Temp_DAR_ATACseq.csv")

library('vegan')
library('adegenet')
library('ape')
library('ggplot2')
library('ggrepel')
library('ggprism')
library('ggh4x')
library('data.table')

DAanno$DA <- "No Change"
DAanno$DA[DAanno$Fold > 0 & DAanno$FDR < 0.05] <- "Open"
DAanno$DA[DAanno$Fold < -0 & DAanno$FDR < 0.05] <- "Closed"

DA_colors <- c("blue", "coral2", "azure4")
names(DA_colors) <- c("Closed", "Open", "No change")

volcano <- ggplot(data = DAanno, 
       aes(x = Fold, y = -log10(FDR) , col = DA)) + 
  geom_jitter() + scale_color_manual(values = DA_colors) + 
  geom_vline(xintercept = c(-0.19, 0.19), col = "brown1") + #this is set based on minimum sig fold change
  geom_hline(yintercept = -log10(0.05), col = "brown1") + 
  theme_prism() + 
  xlim (-1, 1) + 
  ylim(0, 7) +
  xlab(expression(bold("log"[2]~Fold~Change))) +
  ylab(expression(bold("-log"[10]~"Adjusted p-value"))) +
  #geom_text_repel(data = subset(DAanno, (DA == "Open" & !(geneId %like% "LOC"))), 
  #                             aes(x = Fold, y = -log10(FDR), label = geneId), 
  #                             size = 3, color = "black", nudge_x = 0.18) +
  #geom_text_repel(data = subset(DAanno, (DA == "Closed" & !(geneId %like% "LOC"))), 
  #          aes(x = Fold, y = -log10(FDR), label = geneId), 
  #          size = 3, color = "black") +
  theme(legend.position = "none") + 
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
print(volcano)

ggsave(plot = volcano, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/atac-temp-volcano.jpg", 
       scale = 1, width = 6, height = 6, units = c("in"), 
       dpi = 300)

#legend
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', 
     xlim=0:1, ylim=0:1) #creates a blank plot
plotSym <- c(16, 16) #shapes list
par(mai = c(0.5,0.2,0.5,0.5)) #bottom, left, top, right
legend("center", 
       legend = c("Increased \naccessibility", 
                  "Decreased \naccessibility", 
                  "No change"), 
       pch = plotSym, pt.cex = 5.1,
       col = c("coral2", "blue", "azure4"), 
       cex = 2.5, bty = 'n', par(mai = c(0.1,0.1,0.1,0.1))
)

#from DiffBind results, min fold change in sig water peaks was 0.19. That's where the x-intercept lines are drawn.

##### 
DAanno <- read.csv("Temp_DAR_ATACseq.csv") #"tempDA_annotated_peaks_ATAC-MDII.csv")
DAanno <- subset(DAanno, select = -c(X, id_start_end, annotation, geneChr, geneStart,
                                     geneEnd, geneLength, geneStrand, transcriptId,
                                     distanceToTSS, width, strand, 
                                     Conc, Conc_18, Conc_14))

DAannoSig <- subset(DAanno, DAanno$FDR <= 0.05)
DAannoSig <- subset(DAannoSig, select = -c(p.value, FDR))
write.csv(DAannoSig, file = "~/Desktop/R/MDII_GO_MWU/temp_Fold_MDATAC.csv", 
          quote = F, row.names = F)

DAanno$pvalue_go[DAanno$p.value <= 0.05] <- 1
DAanno$pvalue_go[DAanno$p.value > 0.05] <- 0
DAanno$pvalue_go = as.factor(DAanno$pvalue_go)
summary(DAanno$pvalue_go)
#0     1 
#64908  5024 
temp_fishers = data.frame(cbind("geneId" = DAanno$geneId, "pvalue" = DAanno$pvalue_go))
write.csv(temp_fishers, file = "~/Desktop/R/MDII_GO_MWU/temp_fishers_MDATAC.csv", 
          quote = F, row.names = F)

#Exporting files with gene conversion names#####
DAanno <- read.csv("Temp_DAR_ATACseq.csv") #"tempDA_annotated_peaks_ATAC-MDII.csv")
conv <- read.csv("~/Desktop/TAMU/Spurp_reference_files/Spurp_genome/AmyGeneNameConversions.csv")
colnames(conv) <- c("names", "description") 

library('vegan')
library('adegenet')
library('ape')
library('ggplot2')
library('ggrepel')

DAanno$DA <- "No Change"
DAanno$DA[DAanno$Fold > 0 & DAanno$FDR < 0.05] <- "Open"
DAanno$DA[DAanno$Fold < -0 & DAanno$FDR < 0.05] <- "Closed"

DAanno <- subset(DAanno, select = -c(X, geneChr,
                                     geneStart, geneEnd, geneLength,
                                     geneStrand, transcriptId,
                                     width, strand, Conc, 
                                     Conc_18, Conc_14
                                     ))

colnames(DAanno)[colnames(DAanno) == "geneId"] <- "names"
DAconv <- merge(DAanno, conv, by = c("names"), all = TRUE) 
DAconv <- DAconv[!is.na(DAconv$Fold),]
write.csv(DAconv, "temp_ATAC-DAsites-all.csv")

convOPEN <- subset(DAconv, DA == "Open")
write.csv(convOPEN, "temp_ATAC-DAsites-OPEN.csv")

convCLOSED <- subset(DAconv, DA == "Closed")
write.csv(convCLOSED, "temp_ATAC-DAsites-CLOSED.csv")

promDA <- subset(DAanno, annotation == "Promoter (<=1kb)")
promDA <- subset(DAanno, select = -c(id_start_end, annotation, distanceToTSS, 
                                     p.value, FDR, DA))
write.csv(promDA, "Temp-allPromoters-Fold.csv", quote = FALSE, row.names = FALSE)