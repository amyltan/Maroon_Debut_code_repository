#Stephanie Hendricks 
#stephaniehendricks@tamu.edu

################### 16S analysis for maroon debut II Amy Tan et al.

#clear working environment
rm(list=ls())

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")

#data
library(tidyverse)
library(dplyr)
#comm eco
library(vegan) 
library(phyloseq)
library(qiime2R)
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

#extract metadata from the phyloseq object
metadata <- sample_data(physeq)
#remove larvae sample C12F14_L
metadata <- metadata[-5, ]
#change names from sterile and filtered to SASW and MBE
metadata$type <- factor(ifelse(grepl("sterile", metadata$type), "SASW", "MBE"),
  levels = c("SASW", "MBE"))
#create a new treat column in metadata by combining column temperature and type in metadata
metadata$treat <- paste(metadata$temperature, metadata$type, sep = "")
#define colors for each treat group
colorsW <- c("14MBE" = "#336666", 
             "14SASW" = "#99FFFF", 
             "18SASW" = "#FFCC00", 
             "18MBE" = "#CC9933") 
#map colors to treat_c column in metadata
metadata$treat_c <- colorsW[metadata$treat]
#update metadata in phyloseq object
sample_data(physeq) <- sample_data(metadata)
(sample_data(physeq))
nsamples(physeq) #12

#extract OTU table
otu_data <- otu_table(physeq)
#check dimensions
dim(otu_data) #415  12
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
#rarefy samples without replacement to simulate even number of reads per sample
  #rarefaction depth chosen is the 90% of the minimum sample depth in the dataset
set.seed(1)
ps.rarefied <- rarefy_even_depth(physeq, rngseed = 1, 
                                 sample.size = 0.9*min(sample_sums(physeq)), 
                                 replace = F)
#`set.seed(1)` was used to initialize repeatable random subsampling.
#Please record this for your records so others can reproduce.

#3OTUs were removed because they are no longer 
#present in any sample after random subsampling

#calculate relative abundance (ra) for barplots
ps.rarefied_ra <- transform_sample_counts(ps.rarefied, function(x) x / sum(x))

#tax_glom merges OTUs with the same taxonomy, summing the abundances by class
ps.class <- tax_glom(ps.rarefied_ra, taxrank = "Class", NArm = FALSE)
ps.class
get_taxa_unique(ps.class, "Class") #23 but one is NA so 22

#23 distinct colors for 23 classes
my_colors <- c("#E41A1C", "#377EB8", "#CDC0F8", "#FDD2D4", "#D86FEF", 
               "#90E0EF", "#F70373", "#999999", "#66C2A5", "#FC8D62", 
               "#8DA0CB", "#E78AC3", "#FF7F00", "#EDC948", "#E5C494", 
               "#B3B3B3", "#B15928", "#1F78B4", "#33A02C", "#A6D854", 
               "#018571", "#6A3D9A", "#CCCCCC")

#create a new SampleID column in ps.class based on rownames in ps.class 
sample_data(ps.class)$SampleID <- rownames(sample_data(ps.class))
#order samples 
sample_order <- c("C8St14_W", "C9St14_W", "C10St14_W", "C8St18.WA", "C9St18.WA", "C10St18.WA",
                  "C11F14_W", "C12F14_W", "C13F14_W", "C11F18_W", "C12F18_W", "C13F18_W")
sample_data(ps.class)$SampleID <- factor(sample_data(ps.class)$SampleID, levels = sample_order)
levels(sample_data(ps.class)$SampleID)
sample_data(ps.class)

#barplot by class for fig 1 ----
class_barplot <- plot_bar(ps.class, fill = "Class", x = "SampleID") + 
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
  force_panelsizes(rows = unit(4, "in"), #this is the y-axis height
                   cols = unit(3, "in")) #this is the x-axis width
class_barplot
setwd("~/Desktop/TanV3V4/plots")
ggsave(plot = class_barplot, filename = "classPlot.jpg", scale = 1, width = 7.1, height = 5, 
       units = c("in"), dpi = 300)

#tax_glom merges OTUs with the same taxonomy, summing the abundances by family
ps.family <- tax_glom(ps.rarefied_ra, taxrank = "Family", NArm = FALSE)
ps.family
get_taxa_unique(ps.family, "Family") #102 but one is NA so 101
sample_data(ps.family)$SampleID <- rownames(sample_data(ps.family))
sample_order <- c("C8St14_W", "C9St14_W", "C10St14_W", "C8St18.WA", "C9St18.WA", "C10St18.WA",
                  "C11F14_W", "C12F14_W", "C13F14_W", "C11F18_W", "C12F18_W", "C13F18_W")
sample_data(ps.family)$SampleID <- factor(sample_data(ps.family)$SampleID, levels = sample_order)
levels(sample_data(ps.family)$SampleID)
sample_data(ps.family)

#barplot by family for supp fig
family_barplot <- plot_bar(ps.family, fill = "Family", x = "SampleID") + 
  facet_wrap(~ type, scales = "free_x", nrow = 1) +
  labs(x = NULL, y = "Relative abundance") +
  theme_prism() +  
  theme(axis.text.x = element_text(size = 14, angle = 90),  
        axis.title.x = element_blank(),  
        axis.text.y = element_text(size = 14, colour = "black"),  
        axis.title.y = element_text(size = 21), 
        axis.ticks.x = element_blank(),  
        legend.position = "none") + #need to stretch fig to show all families 
  scale_y_continuous(expand = c(0, 0)) +
force_panelsizes(rows = unit(5, "in"), #this is the y-axis height
                 cols = unit(5, "in")) #this is the x-axis width
plot(family_barplot)
setwd("~/Desktop/TanV3V4/plots")
ggsave(plot = family_barplot, filename = "familyPlotSupp.jpg", scale = 1, width = 11.1, height = 6, 
       units = c("in"), dpi = 300)

############## alpha diversity ----
#boxplots of alpha diversity - for supp fig (edited by Amy)
shannonPlot <- plot_richness(ps.rarefied, x = "type", color = "type", measures = c("Shannon")) +
  theme_prism() + 
  geom_point(size = 3) +
  xlab("") + 
  theme(legend.position = "none") + 
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(5, "in"))
plot(shannonPlot)
setwd("~/Desktop/TanV3V4/plots")
ggsave(plot = shannonPlot, filename = "shannonPlotSupp.jpg", scale = 1, width = 6, height = 6, 
       units = c("in"), dpi = 300)

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
                                     levels = c("14SASW", "18SASW", "14MBE", "18MBE"))
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
#create a new Shannon column in metadata from Shannon column in rich
metadata$Shannon <- rich$Shannon

#alpha diversity testing for normality and homogeneity of variances ----
#Shapiro-Wilk test for normality 
shapiro_test <- shapiro.test(rich$Shannon) 
print(shapiro_test) #W = 0.8576, p-value = 0.04562  
#Shannon diversity data significantly deviates from normality; data is not normally distributed

#Shapiro-Wilk test visualization
qqnorm(rich$Shannon, main = "Q-Q Plot of Shannon Entropy")
qqline(rich$Shannon, col = "red")

#Levene's Test for variances by type
levene_test <- leveneTest(rich$Shannon ~ sample_data(ps.rarefied)$type)
print(levene_test) #group  1  12.892 0.004925 **
#significant differences in variances between types (sterile vs filtered)

#Levene's Test for variances by temp
levene_test <- leveneTest(rich$Shannon ~ sample_data(ps.rarefied)$temperature)
print(levene_test) #group  1   0.619 0.4497
#no significant differences in variances between temp (14 vs 18)

##Levene's Test for variances by treat
levene_test <- leveneTest(rich$Shannon ~ sample_data(ps.rarefied)$treat)
print(levene_test) #group  3    1.96 0.1987
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

#alpha diversity stats ----
#Wilcoxon test - must have exactly 2 levels/groups 
wilcox.test(rich$Shannon ~ metadata$type) #W = 0, p-value = 0.002165
wilcox.test(rich$Shannon ~ metadata$temperature) #W = 19, p-value = 0.9372

#Kruskal-Wallis test - have more than 2 levels/groups 
kruskal.test(rich$Shannon ~ metadata$treat) #Kruskal-Wallis chi-squared = 8.7436, df = 3, p-value = 0.0329
#since KW was significant, use post-hoc Dunn's test with Benjamini-Hochberg correction
dunn.test(rich$Shannon, metadata$treat, method = "bh")
#Kruskal-Wallis chi-squared = 8.7436, df = 3, p-value = 0.03
#Comparison of x by group                            
#(Benjamini-Hochberg)                              

#14MBE-14SASW 0.0382 *
#14MBE-18MBE 0.3428 
#14MBE-18SASW 0.0472 *
#14SASW-18MBE 0.0542 *
#14SASW-18SASW 0.3670 
#18MBE-18SASW 0.0847 

#define colors for each treat group
treatColors <- c("14SASW" = "#99FFFF", "14MBE" = "#336666", 
                 "18SASW" = "#FFCC00", "18MBE" = "#CC9933")
treatColors2 <- c("14SASW" = "#33CCCC", "14MBE" = "#003333",
                 "18SASW" = "#CC9900", "18MBE" = "#996633") 
tempShapes <- c(16, 17)

#alpha diversity boxplot for fig 1 ----
dev.new()
alphadiv_boxplot <- ggplot(metadata, aes(x = type, y = Shannon, fill = treat)) +
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
  force_panelsizes(rows = unit(5, "in"), #this is the y-axis height
                   cols = unit(7, "in")) #this is the x-axis length
alphadiv_boxplot
setwd("~/Desktop/TanV3V4/plots")
ggsave("alphadiv_boxplot.jpg", plot = alphadiv_boxplot, scale = 1, width = 8, height = 6,
        units = "in", dpi = 300)

############## beta diversity ----
bc_dist <- phyloseq::distance(ps.rarefied, method = "bray")
ordination <- ordinate(ps.rarefied, method = "PCoA", distance = bc_dist)

#define colors for each treat group - change 14SASW to be darker blue instead of cyan
colorsW <- c("14MBE" = "#336666", 
             "14SASW" = "#33CCCC", 
             "18SASW" = "#FFCC00", 
             "18MBE" = "#CC9933") 

#calculate centroids for each treatment group
ordination_coords <- as.data.frame(ordination$vectors)  
ordination_coords$SampleID <- rownames(ordination_coords)
ordination_coords$treat <- sample_data(ps.rarefied)$treat  
centroids <- aggregate(cbind(Axis.1, Axis.2) ~ treat, data = ordination_coords, FUN = mean)  

#beta diversity pcoa plot for fig 1 ----
betadiv_pcoa_shapes <- plot_ordination(ps.rarefied, ordination, color = "treat") +
  geom_point(aes(shape = temperature), size = 4) +  
  geom_segment(data = ordination_coords, aes(x = Axis.1, y = Axis.2, 
                                             xend = centroids$Axis.1[match(treat, centroids$treat)], 
                                             yend = centroids$Axis.2[match(treat, centroids$treat)], 
                                             colour = treat)) +
  scale_color_manual(values = colorsW) +
  scale_shape_manual(values = c("18" = 17, "14" = 16)) +  #triangle = 17, circle = 16
  labs(x = "PCoA1 (82%)", 
       y = "PCoA2 (6%)") +
  theme_prism() +
  theme(axis.text.x = element_text(size = 18), #diff size than others in fig 1
        axis.text.y = element_text(size = 18), #diff size than others in fig 1
        axis.title.x = element_text(size = 22), #diff size than others in fig 1
        axis.title.y = element_text(size = 22), #diff size than others in fig 1
        legend.position = "none") +
  force_panelsizes(rows = unit(6, "in"), #this is the y-axis height
                   cols = unit(6, "in")) #this is the x-axis width
betadiv_pcoa_shapes
setwd("~/Desktop/TanV3V4/plots")
ggsave("betadiv_pcoaplot_shapes.jpg", plot = betadiv_pcoa_shapes, scale = 1, width = 7, height = 7,
       units = "in", dpi = 300)

#beta diversity stats ----
set.seed(111)
adonis2(bc_dist ~ sample_data(ps.rarefied)$type, permutations = 999)
#adonis2(formula = bc_dist ~ sample_data(ps.rarefied)$type, permutations = 999)
#Df SumOfSqs      R2      F Pr(>F)   
#Model     1   1.9153 0.81823 45.014  0.003 **
#  Residual 10   0.4255 0.18177                 
#Total    11   2.3408 1.00000                 

set.seed(222)
adonis2(bc_dist ~ sample_data(ps.rarefied)$temperature, permutations = 999)
#adonis2(formula = bc_dist ~ sample_data(ps.rarefied)$temperature, permutations = 999)
#Df SumOfSqs      R2      F Pr(>F)
#Model     1  0.08161 0.03486 0.3612  0.573
#Residual 10  2.25922 0.96514              
#Total    11  2.34083 1.00000 

set.seed(333)
adonis2(bc_dist ~ sample_data(ps.rarefied)$treat, permutations = 999)
#adonis2(formula = bc_dist ~ sample_data(ps.rarefied)$treat, permutations = 999)
#Df SumOfSqs      R2      F Pr(>F)    
#Model     3  2.07967 0.88844 21.236  0.001 ***
#  Residual  8  0.26115 0.11156                  
#Total    11  2.34083 1.00000 

############## Coracle linear regressions ----
####### length uncorrected ----
setwd("~/Desktop/TanV3V4/Coracle")
length_reg <- read.csv("coracle_length_regression.csv")
#change names from sterile and filtered to SASW and MBE
length_reg$type <- factor(
  ifelse(grepl("sterile", length_reg$type), "SASW", "MBE"),
  levels = c("SASW", "MBE"))

#RhodobacteraceaeSulfitobacter 5 length for fig 2 ----
#check if linear model is appropriate to use
gvmodel <- gvlma(lm(length ~ RhodobacteraceaeSulfitobacter, data = length_reg))
plot(gvmodel)
summary(gvmodel) #all are acceptable !!!
model <- lm(length ~ RhodobacteraceaeSulfitobacter, data = length_reg)
summary(model)
#Residual standard error: 12.97 on 10 degrees of freedom
#Multiple R-squared:  0.6333,	Adjusted R-squared:  0.5967 
#F-statistic: 17.27 on 1 and 10 DF,  p-value: 0.001961
rhodosulfito_length_reg <- ggplot(length_reg, aes(x = RhodobacteraceaeSulfitobacter, y = length)) +
  geom_smooth(method = "lm", se = TRUE, colour = "black") +  
  geom_point(aes(shape = as.factor(type), colour = as.factor(type)), size = 3) +  
  theme_prism() +
  labs(x = "Rhodobacteraceae; Sulfitobacter count", 
       y = "\nMean length (µm)", 
       shape = "Type", 
       colour = "Type") +
  scale_x_continuous(limits = c(0, 2000), breaks = seq(0, 2000, by = 500)) +  
  scale_y_continuous(breaks = seq(160, 240, by = 40)) +  
  scale_shape_manual(values = c(16, 16)) +  
  scale_color_manual(values = c("SASW" = "#CDC0F8", "MBE" = "#6A3D9A")) + 
  coord_cartesian(ylim = c(160, 280)) +  
  theme(axis.text.x = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 16, colour = "black"),  
        axis.title.x = element_text(size = 20),  
        axis.title.y = element_text(size = 20),
        legend.position = "none",
        legend.title = element_blank()) +
  force_panelsizes(rows = unit(5, "in"), #this is the y-axis height
                   cols = unit(7, "in")) #this is the x-axis width
rhodosulfito_length_reg
setwd("~/Desktop/TanV3V4/plots")
ggsave("rhodosulfito_length_reg.jpg", plot = rhodosulfito_length_reg, height = 6, width = 8.2,
       units = "in", dpi = 300)

####### pigment cell counts corrected ----
#only 11 chambers, not 12
setwd("~/Desktop/TanV3V4/Coracle")
pigment_reg <- read.csv("coracle_pigment_regression.csv")
#change names from sterile and filtered to SASW and MBE
pigment_reg$type <- factor(
  ifelse(grepl("sterile", pigment_reg$type), "SASW", "MBE"),
  levels = c("SASW", "MBE"))

#RhodobacteraceaeSulfitobacter 1 pigment counts for fig 6 ----
#check if linear model is appropriate to use
gvmodel <- gvlma(lm(Pigment_cell_count ~ RhodobacteraceaeSulfitobacter, data = pigment_reg))
plot(gvmodel)
summary(gvmodel) #all are acceptable !!!
model <- lm(Pigment_cell_count ~ RhodobacteraceaeSulfitobacter, data = pigment_reg)
summary(model)
#Residual standard error: 0.07497 on 9 degrees of freedom
#Multiple R-squared:  0.5769,	Adjusted R-squared:  0.5299 
#F-statistic: 12.27 on 1 and 9 DF,  p-value: 0.006692
rhodosulfito_pigment_reg <- ggplot(pigment_reg, aes(x = RhodobacteraceaeSulfitobacter, y = Pigment_cell_count)) +
  geom_smooth(method = "lm", se = TRUE, colour = "black") +  
  geom_point(aes(shape = as.factor(type), colour = as.factor(type)), size = 3) +  
  labs(x = "Rhodobacteraceae; Sulfitobacter count", 
       y = "Mean pigment cell count 
(adjusted for body length)", 
       shape = "Type", 
       colour = "Type") +
  theme_prism() +
  scale_x_continuous(limits = c(0, 2000), breaks = seq(0, 2000, by = 500)) +  
  scale_y_continuous(breaks = seq(1.5, 1.9, by = 0.2)) +  
  scale_shape_manual(values = c(16, 16)) +  
  scale_color_manual(values = c("SASW" = "#CDC0F8", "MBE" = "#6A3D9A")) +
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
#change names from sterile and filtered to SASW and MBE
gut_reg$type <- factor(
  ifelse(grepl("sterile", gut_reg$type), "SASW", "MBE"),
  levels = c("SASW", "MBE"))

#RhodobacteraceaeSulfitobacter 1 gut thick for fig 6 ----
#check if linear model is appropriate to use
gvmodel <- gvlma(lm(Mean_gut_thickness ~ RhodobacteraceaeSulfitobacter, data = gut_reg))
plot(gvmodel)
summary(gvmodel) #all are acceptable !!!
model <- lm(Mean_gut_thickness ~ RhodobacteraceaeSulfitobacter, data = gut_reg)
summary(model)
#Residual standard error: 0.07398 on 9 degrees of freedom
#Multiple R-squared:  0.3869,	Adjusted R-squared:  0.3188 
#F-statistic:  5.68 on 1 and 9 DF,  p-value: 0.04101
rhodosulfito_gut_reg <- ggplot(gut_reg, aes(x = RhodobacteraceaeSulfitobacter, y = Mean_gut_thickness)) +
  geom_smooth(method = "lm", se = TRUE, colour = "black") +  
  geom_point(aes(shape = as.factor(type), colour = as.factor(type)), size = 3) +  
  labs(x = "Rhodobacteraceae; Sulfitobacter count", 
       y = "Mean midgut epithelium thickness 
(adjusted for body length)", 
       shape = "Type", 
       colour = "Type") +
  theme_prism() +
  scale_x_continuous(limits = c(0, 2000), breaks = seq(0, 2000, by = 500)) +  
  scale_y_continuous(breaks = seq(0.0, 0.75, by = 0.25)) +  
  scale_shape_manual(values = c(16, 16)) +  
  scale_color_manual(values = c("SASW" = "#CDC0F8", "MBE" = "#6A3D9A")) +
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

#plotting ONLY a legend for linear regressions - Amy's code
quartz()
plot(NULL , xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', 
     xlim = 0:1, ylim = 0:1) #creates a blank plot
plotSym <- c(16, 16) #16 is circle 
par(mai = c(0.5,0.5,0.5,0.5))
legend("right", 
       legend = c("SASW", 
                  "MBE"), 
       pch = plotSym, pt.cex = 5.1, #size of points
       col = c("#CDC0F8", "#6A3D9A"), 
       cex = 2.5, 
       bty = 'n', 
       #title = "Treatment",        
       title.adj = 0.1, #moves title name to left or right    
       title.cex = 3)

############## fin
