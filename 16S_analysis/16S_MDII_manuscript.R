#Stephanie Hendricks 
#stephaniehendricks@tamu.edu

################### 16S analysis for maroon debut II 
#Amy Tan et al. 2025

#clear working environment
rm(list=ls())

#qiime2R
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")
#faith's div
remotes::install_github("twbattaglia/btools")

#data
library(tidyverse)
library(dplyr)
#comm eco
library(vegan) 
library(phyloseq)
library(qiime2R)
library(btools) #faith's div
library(picante) #faith's div
#plot aesthetics 
library(ggplot2)
library(ggprism)
library(ggh4x) 
#stats 
library(car) 
library(dunn.test) 
#linear regressions
library(gvlma) 

#set working directory 
setwd("~/Desktop/TanV3V4/gzfiles/silva138.1rescript")

#making data into phyloseq object to analyze and plot data  ----
physeq <- qza_to_phyloseq(
  features = "rep_seq_feature_table2.qza",
  tree = "rooted-tree.qza",
  taxonomy = "classified_taxonomy_table.qza",
  metadata = "MDII16S_metadata.tsv") #this is the QIIME2 metadata file
physeq
#otu_table()   OTU Table:         [ 415 taxa and 13 samples ]

#extract metadata from the phyloseq object
metadata <- sample_data(physeq)
#remove larvae sample C12F14_L
metadata <- metadata[rownames(metadata) != "C12F14_L", ]
#change names from sterile and filtered to LMR and HMR
metadata$type <- factor(ifelse(grepl("sterile", metadata$type), "LMR", "HMR"),
  levels = c("LMR", "HMR"))
#create a new treat column in metadata by combining column temperature and type in metadata
metadata$treat <- paste(metadata$temperature, metadata$type, sep = "")
#define colors for each treat group
colorsW <- c("14HMR" = "#336666", 
             "14LMR" = "#99FFFF", 
             "18LMR" = "#FFCC00", 
             "18HMR" = "#CC9933") 
#map colors to treat_c column in metadata
metadata$treat_c <- colorsW[metadata$treat]
#update metadata in phyloseq object
sample_data(physeq) <- sample_data(metadata)
(sample_data(physeq))
nsamples(physeq) #12

#remove chloroplast,  mitochondria, and eukaryota
ps.clean <- subset_taxa(physeq, 
                        !Order %in% "Chloroplast" &
                          !Family %in% "Mitochondria" &
                          !grepl("Eukaryota", Kingdom) &
                          !is.na(Phylum) & Phylum != "")
ps.clean #otu_table()   OTU Table:         [ 390 taxa and 12 samples ]

physeq <- ps.clean

#check - should return FALSE if removal worked
any(tax_table(physeq)[, "Order"] == "Chloroplast", na.rm = TRUE)
any(tax_table(physeq)[, "Family"] == "Mitochondria", na.rm = TRUE)
any(tax_table(physeq)[, "Kingdom"] == "Eukaryota", na.rm = TRUE)

#prepare OTU matrix to plot rarefaction curves
#extract OTU table
otu_data <- otu_table(physeq)
#check dimensions
dim(otu_data) #390  12
#convert to matrix
otu_matrix <- as(otu_data, "matrix")
#transpose matrix
otu_matrix_transposed <- t(otu_matrix)
#remove rows with 0 total counts
otu_matrix_transposed <- otu_matrix_transposed[rowSums(otu_matrix_transposed) > 0, ]
#remove any rows with NA 
otu_matrix_transposed <- na.omit(otu_matrix_transposed)

#plot rarefaction curves
rarecurve(otu_matrix_transposed, step = 50, cex = 0.5)
#confirm value for min sample.size below  
min(sample_sums(physeq)) #[1] 35532
#rarefy samples without replacement to simulate even number of reads per sample
  #rarefaction depth chosen is the minimum sample depth in the dataset
set.seed(1)
ps.rarefied <- rarefy_even_depth(physeq, rngseed = 1, 
                                 sample.size = 35532, 
                                 replace = FALSE)
#1OTUs were removed because they are no longer 
#present in any sample after random subsampling

#number before rarefying
ntaxa(physeq) #390
#number after rarefying
ntaxa(ps.rarefied) #389

#transform to relative abundance
ps.rel <- transform_sample_counts(ps.rarefied, function(x) x / sum(x))

#agglomerate at Class level
ps.class <- tax_glom(ps.rel, taxrank = "Class", NArm = FALSE)
get_taxa_unique(ps.class, "Class") #19 

#calculate max relative abundance per Class
abund_matrix <- otu_table(ps.class)
class_max_ra <- apply(abund_matrix, 1, max)

#extract and modify taxonomy
taxdf <- as.data.frame(tax_table(ps.class))
taxdf$Class <- as.character(taxdf$Class)
taxdf$MaxRA <- class_max_ra

#identify top classes with max RA ≥ 1%
top_classes <- taxdf$Class[taxdf$MaxRA >= 0.01]

#replace other classes with "Other"
taxdf$Class <- ifelse(taxdf$Class %in% top_classes, taxdf$Class, "Other")
tax_table(ps.class) <- tax_table(as.matrix(taxdf))

#re-agglomerate with new Class labels
ps.class.filtered <- tax_glom(ps.class, taxrank = "Class", NArm = FALSE)

#sample data and taxonomy setup
sample_df <- data.frame(sample_data(ps.class.filtered))
tax_df <- as.data.frame(tax_table(ps.class.filtered))

#define class and sample order
class_order <- c(sort(setdiff(unique(tax_df$Class), "Other")), "Other")
tax_df$Class <- factor(tax_df$Class, levels = class_order)
tax_table(ps.class.filtered) <- tax_table(as.matrix(tax_df))

#check to see classes
table(tax_table(ps.class.filtered)[, "Class"])

sample_order <- c("C8St14_W", "C9St14_W", "C10St14_W", 
                  "C8St18.WA", "C9St18.WA", "C10St18.WA",
                  "C11F14_W", "C12F14_W", "C13F14_W", 
                  "C11F18_W", "C12F18_W", "C13F18_W")

sample_df$SampleID <- factor(rownames(sample_df), levels = sample_order)
sample_data(ps.class.filtered) <- sample_data(sample_df)

#define colors
my_colors <- c(
  "#CDC0F8", #light purple
  "#90E0EF", #light blue
  "#F70373", #magenta
  "#FF7F00", #orange
  "#EDC948", #golden yellow
  "#CCCCCC", #grey
  "#018571"  #green
)[seq_along(class_order)]

#plot barplot
class_barplot <- plot_bar(ps.class.filtered, fill = "Class", x = "SampleID") +
  facet_wrap(~ type, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = my_colors) +
  labs(y = "Relative abundance") +
  theme_prism() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 18),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  scale_y_continuous(expand = c(0, 0)) +
  force_panelsizes(rows = unit(4, "in"), cols = unit(3, "in"))
class_barplot
setwd("~/Desktop/TanV3V4/plots")
ggsave(plot = class_barplot, filename = "classPlot.jpg", scale = 1, width = 7.1, height = 5, 
       units = c("in"), dpi = 300)

#tax_glom merges OTUs with the same taxonomy, summing the abundances by family
ps.family <- tax_glom(ps.rel, taxrank = "Family", NArm = FALSE)
ps.family 
tax_table(ps.family)
get_taxa_unique(ps.family, "Family") #100 but NA and Unknown_Family so 98
sample_data(ps.family)$SampleID <- rownames(sample_data(ps.family))
sample_order <- c("C8St14_W", "C9St14_W", "C10St14_W", "C8St18.WA", "C9St18.WA", "C10St18.WA",
                  "C11F14_W", "C12F14_W", "C13F14_W", "C11F18_W", "C12F18_W", "C13F18_W")
sample_data(ps.family)$SampleID <- factor(sample_data(ps.family)$SampleID, levels = sample_order)
levels(sample_data(ps.family)$SampleID)
sample_data(ps.family)
#look at families 
fams <- get_taxa_unique(ps.family, "Family")
fams_clean <- fams[!is.na(fams) & fams != "Unknown_Family"]
length(fams_clean) #should return 98
fams_clean          

#barplot by family for supp fig
family_barplot <- plot_bar(ps.family, fill = "Family", x = "SampleID") + 
  facet_wrap(~ type, scales = "free_x", nrow = 1) +
  labs(x = NULL, y = "Relative abundance") +
  theme_prism() +  
  theme(axis.text.x = element_text(),  
        axis.title.x = element_blank(),  
        axis.text.y = element_text(size = 14, colour = "black"),  
        axis.title.y = element_text(size = 21), 
        axis.ticks.x = element_blank(),  
        legend.position = "none") + #need to stretch fig to right show all families 
  scale_y_continuous(expand = c(0, 0)) +
force_panelsizes(rows = unit(5, "in"), #this is the y-axis height
                 cols = unit(5, "in")) #this is the x-axis width
plot(family_barplot)
setwd("~/Desktop/TanV3V4/plots")
ggsave(plot = family_barplot, filename = "familyPlotSupp.jpg", scale = 1, width = 11.1, height = 6, 
       units = c("in"), dpi = 300)

############## alpha diversity ----
#boxplots of alpha diversity - for supp fig (edited by Amy)
rarefiedPlot <- plot_richness(ps.rarefied, x = "type", measures = c("Observed", "Shannon", "Chao1", "Simpson")) +
  geom_boxplot() +
  theme_prism() +
  labs(x = "Water condition",
       y = "Alpha diversity measure") + 
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(1.25, "in"))
plot(rarefiedPlot)
setwd("~/Desktop/TanV3V4/plots")
ggsave(plot = rarefiedPlot, filename = "rarefiedPlotSupp-type.jpg", scale = 1, width = 7.2, height = 6, 
       units = c("in"), dpi = 300)

rareTemp <- plot_richness(ps.rarefied, x = "temperature", measures = c("Observed", "Shannon", "Chao1", "Simpson")) + 
  geom_boxplot() +
  theme_prism() +
  labs(x = "Temperature (˚C)",
       y = "Alpha diversity measure") + 
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(1.25, "in"))
plot(rareTemp)
setwd("~/Desktop/TanV3V4/plots")
ggsave(plot = rareTemp, filename = "rarefiedPlotSupp-temp.jpg", scale = 1, width = 7.2, height = 6, 
       units = c("in"), dpi = 300)

#reorder treatments
ps.rarefied@sam_data$treat <- factor(ps.rarefied@sam_data$treat, 
                                     levels = c("14LMR", "18LMR", "14HMR", "18HMR"))
rareTreat <- plot_richness(ps.rarefied, x = "treat", 
                           measures = c("Observed", "Shannon", "Chao1", "Simpson")) + 
  geom_boxplot() +
  theme_prism() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.6)) +
  labs(x = "Treatment",
       y = "Alpha diversity measure") + 
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(1.25, "in"))
plot(rareTreat)
setwd("~/Desktop/TanV3V4/plots")
ggsave(plot = rareTreat, filename = "rarefiedPlotSupp-treat.jpg", scale = 1, width = 7.2, height = 6.5, 
       units = c("in"), dpi = 300)

#estimate alpha diversity richness
rich <- estimate_richness(ps.rarefied)

#faiths PD 
fpd <- estimate_pd(ps.rarefied)
rich$Faith_PD <- fpd$PD

#create a columns in metadata from Shannon and Faith_PD column in rich
metadata$Shannon <- rich$Shannon
metadata$Faith_PD <- fpd$PD

#define colors for each treat group
treatColors <- c("14LMR" = "#99FFFF", "14HMR" = "#336666", 
                 "18LMR" = "#FFCC00", "18HMR" = "#CC9933")
treatColors2 <- c("14LMR" = "#33CCCC", "14HMR" = "#003333",
                  "18LMR" = "#CC9900", "18HMR" = "#996633") 
tempShapes <- c(16, 17)

#FPD alpha diversity boxplot for fig 1 ----
dev.new()
FPD_alphadiv_boxplot <- ggplot(metadata, aes(x = type, y = Faith_PD, fill = treat)) +
  geom_boxplot(width = 0.6, position = position_dodge(width = 0.9)) + 
  geom_point(aes(shape = as.factor(temperature), colour = treat), 
             size = 3, alpha = 1,
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9)) +
  theme_prism() +  
  labs(x = NULL, 
       y = "Faith's phylogenetic diversity index") +
  theme(axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5), 
        axis.text.y = element_text(size = 16, colour = "black"),  
        axis.title.y = element_text(size = 20),  
        legend.position = "none", 
        legend.title = element_blank()) +
  scale_fill_manual(values = treatColors) +  
  scale_color_manual(values = treatColors2) +  
  scale_shape_manual(values = tempShapes) + 
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, by = 5)) +  
  force_panelsizes(rows = unit(4, "in"), #this is the y-axis height
                   cols = unit(3, "in")) #this is the x-axis length
FPD_alphadiv_boxplot
setwd("~/Desktop/TanV3V4/plots")
ggsave("FPD_alphadiv_boxplot.jpg", plot = FPD_alphadiv_boxplot, scale = 1, width = 4, height = 5.1,
       units = "in", dpi = 300)

#FPD alpha diversity testing for normality and homogeneity of variances ----
#Shapiro-Wilk test visualization
qqnorm(rich$Faith_PD, main = "Q-Q Plot of Shannon Entropy")
qqline(rich$Faith_PD, col = "red")

#testing type and temp
by(rich$Faith_PD, sample_data(ps.rarefied)$type, shapiro.test)
#W = 0.95038, p-value = 0.7433
#W = 0.90073, p-value = 0.3783
#both NOT sig; data is normally distributed 
by(rich$Faith_PD, sample_data(ps.rarefied)$temperature, shapiro.test)
#W = 0.85713, p-value = 0.1795
#W = 0.80271, p-value = 0.06216
#both NOT sig; data is normally distributed 

#testing each treat 
rich$treat <- sample_data(ps.rarefied)$treat
by(rich$Faith_PD, rich$treat, shapiro.test)
#all NOT sig; data is normally distributed 

#Levene's Test for variances by type
levene_test <- leveneTest(rich$Faith_PD ~ sample_data(ps.rarefied)$type)
print(levene_test)
#group  1  0.8651 0.3742
#no significant differences in variances between types (sterile vs filtered)

#Levene's Test for variances by temp
levene_test <- leveneTest(rich$Faith_PD ~ sample_data(ps.rarefied)$temperature)
print(levene_test)
#group  1  1.7527  0.215
#no significant differences in variances between temp (14 vs 18)

#Levene's Test for variances by treat
levene_test <- leveneTest(rich$Faith_PD ~ sample_data(ps.rarefied)$treat)
print(levene_test)
#group  3  0.4487 0.7251
#no significant differences in variances between treat

#Residuals vs. Fitted Values to check for homoscedasticity
model <- lm(rich$Faith_PD ~ metadata$type)
summary(model)
plot(model$fitted.values, model$residuals,
     main = "Residuals vs. Fitted Values",
     xlab = "Fitted Values",
     ylab = "Residuals",
     pch = 19, col = "blue")  
abline(h = 0, col = "red")

#FPD alpha diversity stats ----
anova_result <- aov(Faith_PD ~ treat, data = rich)
summary(anova_result)
#            Df Sum Sq Mean Sq F value     Pr(>F)    
#treat        3 161.23   53.74   57.89 9.09e-06 ***
#  Residuals    8   7.43    0.93                       
TukeyHSD(anova_result)
#                  diff       lwr       upr     p adj
#18LMR-14LMR -0.9682607 -3.487633  1.551111 0.6263996
#14HMR-14LMR  6.5662719  4.046900  9.085644 0.0001485
#18HMR-14LMR  7.0474885  4.528117  9.566860 0.0000889
#14HMR-18LMR  7.5345326  5.015161 10.053905 0.0000545
#18HMR-18LMR  8.0157492  5.496377 10.535121 0.0000344
#18HMR-14HMR  0.4812165 -2.038155  3.000588 0.9255192

##############################

#Shannon alpha diversity boxplot for fig 1 ----
dev.new()
Shannon_alphadiv_boxplot <- ggplot(metadata, aes(x = type, y = Shannon, fill = treat)) +
  geom_boxplot(width = 0.6, position = position_dodge(width = 0.9)) + 
  geom_point(aes(shape = as.factor(temperature), colour = treat), 
             size = 3, alpha = 1,
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9)) +
  theme_prism() +  
  labs(x = NULL, 
       y = "Shannon diversity index") +
  theme(axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5), 
        axis.text.y = element_text(size = 16, colour = "black"),  
        axis.title.y = element_text(size = 20),  
        legend.position = "none", 
        legend.title = element_blank()) +
  scale_fill_manual(values = treatColors) +  
  scale_color_manual(values = treatColors2) +  
  scale_shape_manual(values = tempShapes) + 
  ylim(min(metadata$Shannon) - 1, max(metadata$Shannon) + 1) +
  force_panelsizes(rows = unit(4, "in"), #this is the y-axis height
                   cols = unit(3, "in")) #this is the x-axis length
Shannon_alphadiv_boxplot
setwd("~/Desktop/TanV3V4/plots")
ggsave("Shannon_alphadiv_boxplot.jpg", plot = Shannon_alphadiv_boxplot, scale = 1, width = 4, height = 5.1,
       units = "in", dpi = 300)

#Shannon alpha diversity testing for normality and homogeneity of variances ----
#Shapiro-Wilk test visualization
qqnorm(rich$Shannon, main = "Q-Q Plot of Shannon Entropy")
qqline(rich$Shannon, col = "red")

#testing type and temp
by(rich$Shannon, sample_data(ps.rarefied)$type, shapiro.test)
#W = 0.97645, p-value = 0.9326
#W = 0.95522, p-value = 0.7823
#both NOT sig; data is normally distributed 
by(rich$Shannon, sample_data(ps.rarefied)$temperature, shapiro.test)
#W = 0.85823, p-value = 0.1832
#W = 0.89398, p-value = 0.3396
#both NOT sig; data is normally distributed 

#testing each treat 
rich$treat <- sample_data(ps.rarefied)$treat
by(rich$Shannon, rich$treat, shapiro.test)
#all NOT sig; data is normally distributed 

#Levene's Test for variances by type
levene_test <- leveneTest(rich$Shannon ~ sample_data(ps.rarefied)$type)
print(levene_test) #group  1  11.852 0.006303 **
#significant differences in variances between types (sterile vs filtered)

#Levene's Test for variances by temp
levene_test <- leveneTest(rich$Shannon ~ sample_data(ps.rarefied)$temperature)
print(levene_test) #group  1   0.872 0.3724
#no significant differences in variances between temp (14 vs 18)

#Levene's Test for variances by treat
levene_test <- leveneTest(rich$Shannon ~ sample_data(ps.rarefied)$treat)
print(levene_test) #group  3  1.9365 0.2023
#no significant differences in variances between treat

#Residuals vs. Fitted Values to check for homoscedasticity
model <- lm(rich$Shannon ~ metadata$type)
summary(model)
plot(model$fitted.values, model$residuals,
     main = "Residuals vs. Fitted Values",
     xlab = "Fitted Values",
     ylab = "Residuals",
     pch = 19, col = "blue")  
abline(h = 0, col = "red")

#Shannon alpha diversity stats ----
anova_shannon <- aov(Shannon ~ treat, data = rich)
summary(anova_shannon)
#Df Sum Sq Mean Sq F value  Pr(>F)   
#treat        3  3.801  1.2670   15.23 0.00114 **
#  Residuals    8  0.665  0.0832                   
TukeyHSD(anova_shannon)
#                   diff        lwr       upr     p adj
#18LMR-14LMR  0.29960279 -0.4544512 1.0536568 0.6031205
#14HMR-14LMR  1.27813650  0.5240825 2.0321905 0.0027762
#18HMR-14LMR  1.23145969  0.4774057 1.9855137 0.0035062
#14HMR-18LMR  0.97853370  0.2244797 1.7325877 0.0135343
#18HMR-18LMR  0.93185690  0.1778029 1.6859109 0.0176330
#18HMR-14HMR -0.04667681 -0.8007308 0.7073772 0.9970101

#define colors for each treat group
treatColors <- c("14LMR" = "#99FFFF", "14HMR" = "#336666", 
                 "18LMR" = "#FFCC00", "18HMR" = "#CC9933")
treatColors2 <- c("14LMR" = "#33CCCC", "14HMR" = "#003333",
                 "18LMR" = "#CC9900", "18HMR" = "#996633") 
tempShapes <- c(16, 17)

############## beta diversity ----
#define colors for each treat group - change 14LMR to be darker blue instead of cyan
colorsW <- c("14HMR" = "#336666", 
             "14LMR" = "#33CCCC", 
             "18LMR" = "#FFCC00", 
             "18HMR" = "#CC9933") 

#weighted unifrac
w_unifrac_dist <- phyloseq::distance(ps.rarefied, method = "unifrac", weighted = TRUE)
ordination <- ordinate(ps.rarefied, method = "PCoA", distance = w_unifrac_dist)

#calculate centroids for each treatment group
ordination_coords <- as.data.frame(ordination$vectors)  
ordination_coords$SampleID <- rownames(ordination_coords)
ordination_coords$treat <- sample_data(ps.rarefied)$treat  
centroids <- aggregate(cbind(Axis.1, Axis.2) ~ treat, data = ordination_coords, FUN = mean)  

w_uni_pcoa_shapes <- plot_ordination(ps.rarefied, ordination, color = "treat") +
  geom_point(aes(shape = temperature), size = 4) +  
  geom_segment(data = ordination_coords, aes(x = Axis.1, y = Axis.2, 
                                             xend = centroids$Axis.1[match(treat, centroids$treat)], 
                                             yend = centroids$Axis.2[match(treat, centroids$treat)], 
                                             colour = treat)) +
  scale_color_manual(values = colorsW) +
  scale_shape_manual(values = c("18" = 17, "14" = 16)) + #triangle = 17, circle = 16
  labs(x = "PCo1 (83%)", 
       y = "PCo2 (7%)") +
  theme_prism() +
  theme(axis.text.x = element_text(size = 20), #diff size than others in fig 1
        axis.text.y = element_text(size = 20), #diff size than others in fig 1
        axis.title.x = element_text(size = 26), #diff size than others in fig 1
        axis.title.y = element_text(size = 26), #diff size than others in fig 1
        legend.position = "none") +
  force_panelsizes(rows = unit(6, "in"), #this is the y-axis height
                   cols = unit(6, "in")) #this is the x-axis width
w_uni_pcoa_shapes
setwd("~/Desktop/TanV3V4/plots")
ggsave("w_uni_pcoaplot_shapes.jpg", plot = w_uni_pcoa_shapes, scale = 1, width = 7.1, height = 7,
       units = "in", dpi = 300)

############## Coracle linear regressions ----
####### length uncorrected ----
setwd("~/Desktop/TanV3V4/Coracle")
length_reg <- read.csv("coracle_length_regression.csv")
#change names from sterile and filtered to SASW and MBE
table(length_reg$type, useNA = "ifany")
length_reg$type <- ifelse(grepl("sterile", length_reg$type, ignore.case = TRUE), "LMR", 
                          ifelse(grepl("filtered", length_reg$type, ignore.case = TRUE), "HMR", NA))
length_reg$type <- factor(length_reg$type, levels = c("LMR", "HMR"))
table(length_reg$type)
str(length_reg$type)

#RhodobacteraceaeSulfitobacter 5 length for fig 2 ----
#check if linear model is appropriate to use
gvmodel <- gvlma(lm(length ~ RhodobacteraceaeSulfitobacter, data = length_reg))
plot(gvmodel)
summary(gvmodel) #all are acceptable !!!
model <- lm(length ~ RhodobacteraceaeSulfitobacter, data = length_reg)
summary(model)
#Residual standard error: 14 on 10 degrees of freedom
#Multiple R-squared:  0.5727,	Adjusted R-squared:   0.53 
#F-statistic:  13.4 on 1 and 10 DF,  p-value: 0.00438
rhodosulfito_length_reg <- ggplot(length_reg, aes(x = RhodobacteraceaeSulfitobacter, y = length)) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", alpha = 0.3) +
  geom_point(aes(shape = type, fill = type), 
             colour = "black", size = 3, shape = 21) +  
  theme_prism() +
  labs(x = "Rhodobacteraceae; Sulfitobacter count", 
       y = "\nMean length (µm)", 
       shape = "Type", 
       fill = "Type") +   # use fill for the legend since shape 21 uses fill
  scale_x_continuous(limits = c(0, 450), breaks = seq(0, 450, by = 100)) +  
  scale_y_continuous(breaks = seq(160, 240, by = 40)) +  
  scale_shape_manual(values = c("LMR" = 21, "HMR" = 21)) +  
  scale_fill_manual(values = c("LMR" = "white", "HMR" = "grey")) + 
  coord_cartesian(ylim = c(160, 280)) +  
  theme(axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),  
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20),
        legend.position = "none",
        legend.title = element_blank()) +
  force_panelsizes(rows = unit(5, "in"), cols = unit(7, "in"))
rhodosulfito_length_reg
setwd("~/Desktop/TanV3V4/plots")
ggsave("rhodosulfito_length_reg.jpg", plot = rhodosulfito_length_reg, height = 6, width = 8.2,
       units = "in", dpi = 300)

####### pigment cell counts corrected ----
#only 11 chambers, not 12
setwd("~/Desktop/TanV3V4/Coracle")
pigment_reg <- read.csv("coracle_pigment_regression.csv")
#change names from sterile and filtered to LMR and HMR
table(pigment_reg$type, useNA = "ifany")
pigment_reg$type <- ifelse(grepl("sterile", pigment_reg$type, ignore.case = TRUE), "LMR", 
                          ifelse(grepl("filtered", pigment_reg$type, ignore.case = TRUE), "HMR", NA))
pigment_reg$type <- factor(pigment_reg$type, levels = c("LMR", "HMR"))
table(pigment_reg$type)
str(pigment_reg$type)

#RhodobacteraceaeSulfitobacter 1 pigment counts for fig 6 ----
#check if linear model is appropriate to use
gvmodel <- gvlma(lm(Pigment_cell_count ~ RhodobacteraceaeSulfitobacter, data = pigment_reg))
plot(gvmodel)
summary(gvmodel) #all are acceptable !!!
model <- lm(Pigment_cell_count ~ RhodobacteraceaeSulfitobacter, data = pigment_reg)
summary(model)
#Residual standard error: 0.07458 on 9 degrees of freedom
#Multiple R-squared:  0.5813,	Adjusted R-squared:  0.5348 
#F-statistic: 12.49 on 1 and 9 DF,  p-value: 0.006367
rhodosulfito_pigment_reg <- ggplot(pigment_reg, aes(x = RhodobacteraceaeSulfitobacter, y = Pigment_cell_count)) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", alpha = 0.3) +
  geom_point(aes(shape = type, fill = type), 
             colour = "black", size = 3, shape = 21) +   
  labs(x = "Rhodobacteraceae; Sulfitobacter count", 
       y = "Mean pigment cell count 
(adjusted for body length)", 
       shape = "Type", 
       colour = "Type") +
  theme_prism() +
  scale_x_continuous(limits = c(0, 450), breaks = seq(0, 450, by = 100)) +  
  scale_y_continuous(breaks = seq(1.5, 1.9, by = 0.2)) +  
  scale_shape_manual(values = c(16, 16)) +  
  scale_fill_manual(values = c("LMR" = "white", "HMR" = "grey")) + 
  coord_cartesian(ylim = c(1.5, 2.1)) +  
  theme(axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),  
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20),
        legend.position = "none",
        legend.title = element_blank()) +
  force_panelsizes(rows = unit(5, "in"), #this is the y-axis height
                   cols = unit(7, "in")) #this is the x-axis width
rhodosulfito_pigment_reg
setwd("~/Desktop/TanV3V4/plots")
ggsave("rhodosulfito_pigment_reg.jpg", plot = rhodosulfito_pigment_reg, width = 8.1, height = 6,
       units = "in", dpi = 300)

####### midgut epithelium thickness corrected ----
#only 11 chambers, not 12
setwd("~/Desktop/TanV3V4/Coracle")
gut_reg <- read.csv("coracle_gutthick_regression.csv")
#change names from sterile and filtered to LMR and HMR
table(gut_reg$type, useNA = "ifany")
gut_reg$type <- ifelse(grepl("sterile", gut_reg$type, ignore.case = TRUE), "LMR", 
                           ifelse(grepl("filtered", gut_reg$type, ignore.case = TRUE), "HMR", NA))
gut_reg$type <- factor(gut_reg$type, levels = c("LMR", "HMR"))
table(gut_reg$type)
str(gut_reg$type)

#RhodobacteraceaeSulfitobacter 1 gut thick for fig 6 ----
#check if linear model is appropriate to use
gvmodel <- gvlma(lm(Mean_gut_thickness ~ RhodobacteraceaeSulfitobacter, data = gut_reg))
plot(gvmodel)
summary(gvmodel) #all are acceptable !!!
model <- lm(Mean_gut_thickness ~ RhodobacteraceaeSulfitobacter, data = gut_reg)
summary(model)
#Residual standard error: 0.06914 on 9 degrees of freedom
#Multiple R-squared:  0.4645,	Adjusted R-squared:  0.405 
#F-statistic: 7.806 on 1 and 9 DF,  p-value: 0.02092
rhodosulfito_gut_reg <- ggplot(gut_reg, aes(x = RhodobacteraceaeSulfitobacter, y = Mean_gut_thickness)) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", alpha = 0.3) +
  geom_point(aes(shape = type, fill = type), 
             colour = "black", size = 3, shape = 21) +   
  labs(x = "Rhodobacteraceae; Sulfitobacter count", 
       y = "Mean midgut epithelium thickness 
(adjusted for body length)", 
       shape = "Type", 
       colour = "Type") +
  theme_prism() +
  scale_x_continuous(limits = c(0, 450), breaks = seq(0, 450, by = 100)) +  
  scale_y_continuous(breaks = seq(0.0, 0.75, by = 0.25)) +  
  scale_shape_manual(values = c(16, 16)) +  
  scale_fill_manual(values = c("LMR" = "white", "HMR" = "grey")) + 
  coord_cartesian(ylim = c(0.0, 1)) +  
  theme(axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),  
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20),
        legend.position = "none",
        legend.title = element_blank()) +
  force_panelsizes(rows = unit(5, "in"), #this is the y-axis height
                   cols = unit(7, "in")) #this is the x-axis width
rhodosulfito_gut_reg
setwd("~/Desktop/TanV3V4/plots")
ggsave("rhodosulfito_gut_reg.jpg", plot = rhodosulfito_gut_reg, width = 8.25, height = 6,
       units = "in", dpi = 300)

############## fin
