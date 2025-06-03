setwd("~/Desktop/R/morphology-MDII/")

library("dplyr")
library('stringr')
library('rstatix')
library('ggplot2')
library('GroupStruct')
library('dunn.test')
library('ggprism')
library('ggh4x')

data = read.csv("morphology-MDII-clean-allMetrics.csv")
data <- filter(data, !(Chamber == "C11F18"))
data = subset(data, select = -c(X))

#####
#add OTU and letter column for GroupStruct
data <- data %>%
  mutate(OTU = case_when(
    str_detect(Chamber,"C10St14") ~ "A",
    str_detect(Chamber,"C10St18") ~ "B",
    str_detect(Chamber,"C12F14") ~ "C",
    str_detect(Chamber,"C12F18") ~ "D",
    str_detect(Chamber,"C13F14") ~ "E",
    str_detect(Chamber,"C13F18") ~ "F",
    str_detect(Chamber,"C8St14") ~ "G",
    str_detect(Chamber,"C9St14") ~ "H",
    str_detect(Chamber,"C9St18") ~ "I",
    str_detect(Chamber,"C11F14") ~ "J",
    str_detect(Chamber,"C8St18") ~ "K",
  ), .before = Sample_name)

dataGS <- subset(data, select = -c(Sample_name, Chamber, water, temp, 
                                   treat, Body_area_um2, StDev_CtG))
dataGS <- dataGS %>% 
  relocate(Length.measurement_.micron., .before = Pigment_cell_count)

dataGS <- allom(dataGS, type = "population1")
pc1 <- ez_pca(dataGS)

##adding Chamber names back to corrected dataframe 
dat <- dataGS %>%
  mutate(Chamber = case_when(
    str_detect(OTU,"A") ~ "C10St14",
    str_detect(OTU,"B") ~ "C10St18",
    str_detect(OTU,"C") ~ "C12F14",
    str_detect(OTU,"D") ~ "C12F18",
    str_detect(OTU,"E") ~ "C13F14",
    str_detect(OTU,"F") ~ "C13F18",
    str_detect(OTU,"G") ~ "C8St14",
    str_detect(OTU,"H") ~ "C9St14",
    str_detect(OTU,"I") ~ "C9St18",
    str_detect(OTU,"J") ~ "C11F14",
    str_detect(OTU,"K") ~ "C8St18",
  ), .before = OTU)

dat <- dat %>%
  mutate(water = case_when(
    str_detect(Chamber, "St") ~ "sterile",
    str_detect(Chamber, "F") ~ "filtered" 
  )) %>%
  mutate(temp = case_when(
    str_detect(Chamber, "14") ~ "14",
    str_detect(Chamber, "18") ~ "18"
  ))
dat = mutate(dat, treat = paste(water, temp, sep = "_")) 

## NOTE: how I'm going to proceed with this data
# Cell counts, gut thickness, and CtG have been adjusted based on larval length
# Analyze those three metrics here--corrected
# Analyze length and body area without corrections
# Should be noted that another possible option would be to standardize these 3 metrics by body area instead of length. Not certain which of those would be better.

##Pigment cell count analysis 
shapiro.test(dat$Pigment_cell_count) #YES 0.335 
levene_test(dat, Pigment_cell_count ~ Chamber) #NO 0.00894
levene_test(dat, Pigment_cell_count ~ treat) #YES 0.354 

dunn.test(dat$Pigment_cell_count, g = dat$treat, kw = TRUE, list = TRUE, 
          method = "bh")
#Kruskal-Wallis chi-squared = 57.2327, df = 3, p-value = 0
#filtered_14 - filtered_18 : -1.394590 (0.0816)
#filtered_14 - sterile_14  : -7.179873 (0.0000)*
#filtered_18 - sterile_14  : -4.432384 (0.0000)*
#filtered_14 - sterile_18  : -4.376324 (0.0000)*
#filtered_18 - sterile_18  : -2.383701 (0.0129)*
#sterile_14 - sterile_18   :  2.052717 (0.0241)*

cellsboxplot <-{dat$water <- factor(dat$water, levels = c("sterile", "filtered"))
dat$treat <- factor(dat$treat, levels = c("sterile_14", "filtered_14",
                                            "sterile_18", "filtered_18"))
boxplot(Pigment_cell_count ~ water * temp, data = dat, frame = TRUE,
        ylab = "Pigment Cell Count\n(adjusted for body length)", 
        cex.main = 1.5, cex.lab = 1.9, xlab = "",
        col = c("#99FFFF","#336666","#FFCC00","#CC9933"), xaxt = 'n')
axis(side = 1, tick = FALSE, at = 0:4,
     labels = c("", paste("14ºC\nSASW"),
                paste("14ºC\nMBE"), 
                paste("18ºC\nSASW"), 
                paste("18ºC\nMBE")), 
     cex.axis = 1.4, line = 0.6) 
stripchart(data = dat, Pigment_cell_count ~ treat, method = "jitter", pch = 20,
           vertical = TRUE, add = TRUE, 
           col = c("#33CCCC","#003333","#CC9900","#996633"))
par(mai = c(0.6,1.2,0.1,0.1))
}

#Gut thickness
shapiro.test(dat$Mean_gut_thickness) #YES 0.009398 
levene_test(dat, Mean_gut_thickness ~ Chamber) #NO 0.0135
levene_test(dat, Mean_gut_thickness ~ treat) #NO 0.0174 

kruskal.test(dat, Mean_gut_thickness ~ water) #p-value < 2.2e-16 
kruskal.test(dat, Mean_gut_thickness ~ temp) #p-value < 2.2e-16
dunn.test(dat$Mean_gut_thickness, g = dat$treat, kw = TRUE, list = TRUE, 
          method = "bh")
#Kruskal-Wallis chi-squared = 14.634, df = 3, p-value = 0
#filtered_14 - filtered_18 : -0.102337 (0.4592)
#filtered_14 - sterile_14  : -3.537660 (0.0002)*
#filtered_18 - sterile_14  : -2.765932 (0.0028)*
#filtered_14 - sterile_18  : -1.082112 (0.1396)
#filtered_18 - sterile_18  : -0.813539 (0.2080)
#sterile_14 - sterile_18   :  2.079389 (0.0376)*

gutboxplot <-{dat$water <- factor(dat$water, levels = c("sterile", "filtered"))
dat$treat <- factor(dat$treat, levels = c("sterile_14", "filtered_14",
                                          "sterile_18", "filtered_18"))
boxplot(Mean_gut_thickness ~ water * temp, data = dat, frame = TRUE,
        ylab = "Midgut Epithelium Thickness\n(adjusted for body length)", 
        cex.main = 1.5, cex.lab = 1.9, xlab = "",
        col = c("#99FFFF","#336666","#FFCC00","#CC9933"), xaxt = 'n')
axis(side = 1, tick = FALSE, at = 0:4,
     labels = c("", paste("14ºC\nSASW"),
                paste("14ºC\nMBE"), 
                paste("18ºC\nSASW"), 
                paste("18ºC\nMBE")), 
     cex.axis = 1.4, line = 0.6) 
stripchart(data = dat, Mean_gut_thickness ~ treat, method = "jitter", pch = 20,
           vertical = TRUE, add = TRUE, 
           col = c("#33CCCC","#003333","#CC9900","#996633"))
par(mai = c(0.6,1.4,0.1,0.1))
}

#Cell-to-gut distance
shapiro.test(dat$Mean_CtG) #NO 0.01138 
levene_test(dat, Mean_CtG ~ Chamber) #YES 0.781

ctg_aog <- aov(data = dat, Mean_CtG ~ water * temp)
summary(ctg_aog)
TukeyHSD(ctg_aog) #water and temp both sig 
#interactions sig for F14:S14, S18:F14, F14:F18, F18:S14, S18:F18
#NOT sig for S18:S14

CtGboxplot <-{dat$water <- factor(dat$water, levels = c("sterile", "filtered"))
dat$treat <- factor(dat$treat, levels = c("sterile_14", "filtered_14",
                                          "sterile_18", "filtered_18"))
boxplot(Mean_CtG ~ water * temp, data = dat, frame = TRUE,
        ylab = "Cell-to-Gut distance (adjusted for body length)", xlab = "", 
        main = "11 chambers",
        col = c("#99FFFF","#336666","#FFCC00","#CC9933"), 
        names = c("14ºC sterile","14ºC filtered", "18ºC sterile", "18ºC filtered")) 
stripchart(data = dat, Mean_CtG ~ treat, method = "jitter", pch = 20,
           vertical = TRUE, add = TRUE, 
           col = c("#33CCCC","#003333","#CC9900","#996633"))
}

# checking for correlation between body area and gut size 
ggplot(dat, 
       aes(x=Length.measurement_.micron., y=Mean_gut_thickness)) + 
  geom_point(size = 6, aes(color=treat)) +
  geom_smooth(formula = y ~ x, method = lm, se = F) +
  theme_classic()

#####
##notes to self: 
# try plotting by condition subsets. Might show a more clear relationship there.
# if there's a good relationship, it would indicate there's no obvious inflammation

St14 = subset(dat, treat == "sterile_14")
F14 = subset(dat, treat == "filtered_14")  
St18 = subset(dat, treat == "sterile_18")
F18 = subset(dat, treat == "filtered_18")

ggplot(St14, 
       aes(x=Length.measurement_.micron., y=Mean_gut_thickness)) + 
  geom_point(aes(color = Chamber), size = 6) +
  geom_smooth(formula = y ~ x, method = lm, se = F) +
  theme_classic()

ggplot(F14, 
       aes(x=Length.measurement_.micron., y=Mean_gut_thickness)) + 
  geom_point(aes(color = Chamber), size = 6) +
  geom_smooth(formula = y ~ x, method = lm, se = F) +
  theme_classic()

ggplot(St18, 
       aes(x=Length.measurement_.micron., y=Mean_gut_thickness)) + 
  geom_point(aes(color = Chamber), size = 6) +
  geom_smooth(formula = y ~ x, method = lm, se = F) +
  theme_classic()

ggplot(F18, 
       aes(x=Length.measurement_.micron., y=Mean_gut_thickness)) + 
  geom_point(aes(color = Chamber), size = 6) +
  geom_smooth(formula = y ~ x, method = lm, se = F) +
  theme_classic() 

#####
# trying another method of graphing 
cor.test(St14$Mean_gut_thickness, St14$Length.measurement_.micron., method = "pearson") #NO cor -0.07074384 
cor.test(F14$Mean_gut_thickness, F14$Length.measurement_.micron., method = "pearson") #NO cor = 0.01877806
cor.test(St18$Mean_gut_thickness, St18$Length.measurement_.micron., method = "pearson") #slightly cor = 0.1822534
cor.test(F18$Mean_gut_thickness, F18$Length.measurement_.micron., method = "pearson") #NO cor = 0.0891035

#pigment cell correlations and graphs 
cor.test(dat$Pigment_cell_count, dat$Length.measurement_.micron., method = "pearson") #0.5423416
cor.test(St14$Pigment_cell_count, St14$Length.measurement_.micron., method = "pearson") #NO cor 0.03303216 
cor.test(F14$Pigment_cell_count, F14$Length.measurement_.micron., method = "pearson") #very slight cor = 0.1573511
cor.test(St18$Pigment_cell_count, St18$Length.measurement_.micron., method = "pearson") #yes cor = 0.5492751
cor.test(F18$Pigment_cell_count, F18$Length.measurement_.micron., method = "pearson") #very slightly cor = -0.146104

ggplot(St14, 
       aes(x=Length.measurement_.micron., y=Pigment_cell_count)) + 
  geom_point(aes(color = Chamber), size = 6) +
  geom_smooth(formula = y ~ x, method = lm, se = F) +
  theme_classic()

ggplot(F14, 
       aes(x=Length.measurement_.micron., y=Pigment_cell_count)) + 
  geom_point(aes(color = Chamber), size = 6) +
  geom_smooth(formula = y ~ x, method = lm, se = F) +
  theme_classic()

ggplot(St18, 
       aes(x=Length.measurement_.micron., y=Pigment_cell_count)) + 
  geom_point(aes(color = Chamber), size = 6) +
  geom_smooth(formula = y ~ x, method = lm, se = F) +
  theme_classic()

ggplot(F18, 
       aes(x=Length.measurement_.micron., y=Pigment_cell_count)) + 
  geom_point(aes(color = Chamber), size = 6) +
  geom_smooth(formula = y ~ x, method = lm, se = F) +
  theme_classic() 

ggplot(dat, 
       aes(x=Length.measurement_.micron., y=Pigment_cell_count)) + 
  geom_point(aes(color = treat), size = 6) +
  geom_smooth(formula = y ~ x, method = lm, se = F) +
  theme_classic() 

#pigment cell means 
cellmeans <- aggregate(Pigment_cell_count ~ treat, data = dat, mean)
cellsd <- aggregate(Pigment_cell_count ~ treat, data = dat, sd)
#treat Pigment_cell_count
#1 filtered_14                    1.589047 0.09468654
#2 filtered_18                    1.638984 0.11001893
#3  sterile_14                    1.802265 0.08406487
#4  sterile_18                    1.731553 0.08592712

gutmeans <- aggregate(Mean_gut_thickness ~ treat, data = dat, mean)
gutsd <- aggregate(Mean_gut_thickness ~ treat, data = dat, sd)
#treat Mean_gut_thickness
#1 filtered_14                    0.5388221 0.18104877
#2 filtered_18                    0.5768055 0.09442915
#3  sterile_14                    0.6896927 0.12782026
#4  sterile_18                    0.6084841 0.10931873 

#average of S18, F14, F18 = 0.5747 
#% change - S14 is 16% thicker

aggregate(Mean_gut_thickness ~ water, data = dat, mean)
#water Mean_gut_thickness
#1 filtered          0.5512397
#2  sterile          0.6569243
aggregate(Mean_gut_thickness ~ temp, data = dat, mean)
#temp Mean_gut_thickness
#1   14          0.6131641
#2   18          0.5950207

##### 
# Creating output of all 12 chambers with body length
data = read.csv("morphology-MDII-clean-allMetrics.csv")
data = subset(data, select = -c(X))

dataMeans <- data %>% 
  group_by(Chamber) %>% 
  mutate(Pigment_cell_count = mean(Pigment_cell_count), 
         Length.measurement_.micron. = mean(Length.measurement_.micron.), 
         Body_area_um2 = mean(Body_area_um2), 
         Mean_gut_thickness = mean(Mean_gut_thickness)) 

dataMeans <- dataMeans %>% 
  group_by(Chamber, Pigment_cell_count, Length.measurement_.micron.,
           Body_area_um2, Mean_gut_thickness) %>% 
  summarize('sample_number' = n())

dataMeans <- mutate(dataMeans, Body_area_mm2 = Body_area_um2/1000)
dataMeans <- subset(dataMeans, select = -c(Body_area_um2, sample_number))
write.csv(dataMeans, "~/Desktop/R/MDII_mixOmics/MD-morphologyMeans-columns.csv")

#note: not using length-corrected values right now because GroupStruct needs multiple points 
#also, GroupStruct doesn't like the 0 count samples in C11F18 pigment cells

meansOut <- data.frame(t(dataMeans[-1]))
print(dataMeans$Chamber) #to get order of chamber names
colnames(meansOut) <- c("C10St14", "C10St18", "C11F14", "C11F18", "C12F14", "C12F18",
                        "C13F14", "C13F18", "C8St14", "C8St18", "C9St14", "C9St18")
write.csv(meansOut, "~/Desktop/R/MDII_mixOmics/MD-morphologyMeans.csv")

# generating other versions of plots #####
h2o <- subset(dat, temp == "14")

cells14boxplot <-{h2o$water <- factor(h2o$water, levels = c("sterile", "filtered"))
data$treat <- factor(data$treat, levels = c("sterile_14", "filtered_14"))
boxplot(Pigment_cell_count ~ water, data = h2o, frame = TRUE,
        ylab = "Adjusted Pigment Cell Count", xlab = "", 
        col = c("#99FFFF","#336666"), xaxt = 'n', cex.axis = 1.6, 
        cex.main = 2.4, cex.lab = 2)
axis(side = 1, tick = FALSE, at = 0:2,
     labels = c("", paste("Sterile ASW\n14ºC"),
                paste("Microbiome-exposed\n14ºC")), 
     cex.axis = 2.2, line = 1.9) 
stripchart(data = h2o, Pigment_cell_count ~ water, method = "jitter", pch = 20,
           vertical = TRUE, add = TRUE, 
           col = c("#33CCCC","#003333"))
par(mai = c(0.8,1,0.1,0.1))
}

# X-plots of cell counts and gut thickness
library('Rmisc')
library('ggprism')
#levels = c("sterile_14", "filtered_14", "sterile_18", 
#           "filtered_18")

treatColors = c("#33CCCC", "#336666", "#FFCC00", "#CC9933")
names(treatColors) = c("sterile_14", "filtered_14", 
                       "sterile_18", "filtered_18") 

tempShapes = c(16, 17)
names(waterShapes) = c("14", "18")

cellSumm = summarySE(data = dat, measurevar = "Pigment_cell_count", 
                     groupvars = c("treat", "temp", "water"))
pd <- position_dodge(0.2)

ggplot(cellSumm, aes(x = water, y = Pigment_cell_count, 
                     color = treat)) +
  geom_point(aes(shape= water),size=10,position=pd)+
  #geom_line(aes(group= water, color=temp),position=pd)+
  geom_errorbar(aes(ymin=Pigment_cell_count-sd,
                    ymax=Pigment_cell_count+sd),
                lwd=1.2,width=0.2,position=pd)+
  scale_shape_manual(values= c(17,16))+
  scale_color_manual(values= treatColors)+
  theme_prism()+
  #facet_wrap(~ treat,scales="free_y", ncol=8)+
  theme(legend.text=element_text(size=10)) +
  theme(legend.key = element_blank())+
  theme(legend.position = 'right') + 
  ylab("Adjusted pigment cell count") + 
  xlab(NULL) + 
  xlim("sterile", "filtered")


ggplot(dat, aes(x = temp, y = Pigment_cell_count, 
                     color = treat, shape = water)) +
  geom_boxplot() + 
  geom_jitter() +
  scale_shape_manual(values=c(16,17)) +
  scale_color_manual(values=c("#336666", "#CC9933","#33CCCC", "#FFCC00"))+
  theme_prism()+
  theme(legend.text=element_text(size=10)) +
  theme(legend.key = element_blank())+
  theme(legend.direction = 'horizontal', legend.position = 'top') + 
  ylab("Adjusted pigment cell count") + 
  xlab(NULL)


gutSumm = summarySE(data = dat, measurevar = "Mean_gut_thickness", 
                     groupvars = c("treat", "temp", "water"))
pd <- position_dodge(0.1)

ggplot(gutSumm, aes(x = temp, y = Mean_gut_thickness, 
                     color = water)) +
  geom_point(aes(shape= water),size=5,position=pd)+
  #geom_line(aes(group= water, color=temp),position=pd)+
  geom_errorbar(aes(ymin=Mean_gut_thickness-sd,
                    ymax=Mean_gut_thickness+sd),
                lwd=0.4,width=0.3,position=pd)+
  scale_shape_manual(values=c(16,17))+
  scale_color_manual(values=c("#FFCC00", "#CC9933", "#336666", "#33CCCC"))+
  theme_prism()+
  #facet_wrap(~ treat,scales="free_y", ncol=8)+
  theme(legend.text=element_text(size=10)) +
  theme(legend.key = element_blank())+
  theme(legend.direction = 'horizontal', legend.position = 'top') + 
  ylab("Adjusted mid-gut epithelial thickness") + 
  xlab(NULL)


ggplot(dat, aes(x = temp, y = Mean_gut_thickness, 
                color = treat, shape = water)) +
  geom_boxplot() + 
  geom_jitter() +
  scale_shape_manual(values=c(16,17)) +
  scale_color_manual(values=c("#336666", "#CC9933","#33CCCC", "#FFCC00"))+
  theme_prism()+
  theme(legend.text=element_text(size=10)) +
  theme(legend.key = element_blank())+
  theme(legend.direction = 'horizontal', legend.position = 'top') + 
  ylab("Adjusted gut thickness") + 
  xlab(NULL)

# new boxplots #####
library('patchwork')

treatColors = c("#99FFFF", "#336666", "#FFCC00", "#CC9933")
names(treatColors) = c("sterile_14", "filtered_14", 
                       "sterile_18", "filtered_18") 

treatColors2 = c("#33CCCC","#003333","#CC9900","#996633") 
names(treatColors2) = c("sterile_14", "filtered_14", 
                        "sterile_18", "filtered_18") 

tempShapes = c(16, 17)
names(tempShapes) = c("14", "18")

orderList = c("sterile", "filtered") 
names(orderList) = c("SASW", "MBE")


cellPlot <- ggplot() + 
  theme_prism() +
  geom_boxplot(data = dat, aes(
    x = water, y = Pigment_cell_count, fill = treat), 
    show.legend = FALSE, width = 0.5, 
    position = position_dodge(width = 0.7))  + 
  geom_point(data = dat, aes(
    x = water, y = Pigment_cell_count, 
    shape = temp, color = treat),
    show.legend = FALSE, size = 3,
    position = position_jitterdodge(dodge.width = 0.7, 
                                    jitter.width = 0.3)) + 
  scale_fill_manual(values = treatColors) + 
  scale_color_manual(values = treatColors2) + 
  scale_shape_manual(values = tempShapes) + 
  scale_x_discrete(limits = orderList, labels = names(orderList)) +
  xlab(NULL) + 
  ylab("Pigment cell count\n(adjusted for body length)") + #"Pigment cell count\n(adjusted for body length)"
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) + 
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(7, "in"))
print(cellPlot)
ggsave(plot = cellPlot, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/cellPlot.jpg", 
       scale = 1, width = 8.1, height = 6, units = c("in"), 
       dpi = 300)

gutPlot <- ggplot() + 
  theme_prism() +
  geom_boxplot(data = dat, aes(
    x = water, y = Mean_gut_thickness, fill = treat), 
    show.legend = FALSE, width = 0.5, 
    position = position_dodge(width = 0.7))  + 
  geom_point(data = dat, aes(
    x = water, y = Mean_gut_thickness, 
    shape = temp, color = treat),
    show.legend = FALSE, size = 3,
    position = position_jitterdodge(dodge.width = 0.7, 
                                    jitter.width = 0.3)) + 
  scale_fill_manual(values = treatColors) + 
  scale_color_manual(values = treatColors2) + 
  scale_shape_manual(values = tempShapes) + 
  scale_x_discrete(limits = orderList, labels = names(orderList)) +
  xlab(NULL) + 
  ylab("Midgut epithelium thickness\n(adjusted for body length)") + #"Midgut epithelium thickness\n(adjusted for body length)"
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) + 
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(7, "in"))
print(gutPlot)

ggsave(plot = gutPlot, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/gutPlot.jpg", 
       scale = 1, width = 8.25, height = 6, units = c("in"), 
       dpi = 300)
