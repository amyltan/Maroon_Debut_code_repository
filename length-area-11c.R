setwd("~/Desktop/R/morphology-MDII/")

library("dplyr")
library('stringr')
library('rstatix')
library('ggplot2')
library('ggprism')
library('ggh4x')

data = read.csv("morphology-MDII-clean-allMetrics.csv")
data = subset(data, select = -c(X)) 
data <- filter(data, !(Chamber == "C11F18"))
data$water <- as.factor(data$water)
data$temp <- as.factor(data$temp)

## Length analysis
shapiro.test(data$Length.measurement_.micron.) #YES 0.96 
levene_test(data, Length.measurement_.micron. ~ Chamber) #YES 0.353

length_aov <- aov(data = data, Length.measurement_.micron. ~ water * temp)
summary(length_aov)
TukeyHSD(length_aov) #water is significant, temp is also (why did this say not earlier?) 
#interactions sig: S14:F14, F14:F18, F14:S18, S14:F18, S18:F18 
#not sig: S14:S18 (but not S18:F18)

lengthboxplot <-{data$water <- factor(data$water, levels = c("sterile", "filtered"))
data$treat <- factor(data$treat, levels = c("sterile_14", "filtered_14",
                                            "sterile_18", "filtered_18"))
boxplot(Length.measurement_.micron. ~ water * temp, data = data, frame = TRUE,
        ylab = "Length (µm)", xlab = "", cex.main = 1.5, cex.lab = 2,
        col = c("#99FFFF","#336666","#FFCC00","#CC9933"), xaxt = 'n')
        axis(side = 1, tick = FALSE, at = 0:4,
             labels = c("", paste("14ºC\nSASW"),
                                               paste("14ºC\nMBE"), 
                                               paste("18ºC\nSASW"), 
                                               paste("18ºC\nMBE")),
             cex.axis = 1.4, line = 0.6) 
stripchart(data = data, Length.measurement_.micron. ~ treat, method = "jitter", pch = 20,
           vertical = TRUE, add = TRUE, 
           col = c("#33CCCC","#003333","#CC9900","#996633"))
par(mai = c(0.6,1,0.1,0.1))
}


## Area analysis
#convert um2 to mm2 
data <- mutate(data, Body_area_mm2 = Body_area_um2/1000)

shapiro.test(data$Body_area_um2) #NO 0.002138 
levene_test(data, Body_area_um2 ~ Chamber) #YES 0.0578

area_aov <- aov(data = data, Body_area_um2 ~ water * temp)
summary(area_aov)
TukeyHSD(area_aov) #water and temp are both sig
#interactions sig: S14:F14, S14:F18, F14:S18, F14:F18, S18:F18
#not sig: S14:S18 #note: same as 9 chambers

areaboxplot <-{data$water <- factor(data$water, levels = c("sterile", "filtered"))
data$treat <- factor(data$treat, levels = c("sterile_14", "filtered_14",
                                            "sterile_18", "filtered_18"))
boxplot(Body_area_mm2 ~ water * temp, data = data, frame = TRUE,
        ylab = expression(Body~Area~"(mm"^"2"~")"), cex.main = 1.5, 
        cex.lab = 2, xlab = "", 
        col = c("#99FFFF","#336666","#FFCC00","#CC9933"), xaxt = 'n')
axis(side = 1, tick = FALSE, at = 0:4,
     labels = c("", paste("14ºC\nSASW"),
                paste("14ºC\nMBE"), 
                paste("18ºC\nSASW"), 
                paste("18ºC\nMBE")), 
     cex.axis = 1.4, line = 0.6) 
stripchart(data = data, Body_area_mm2 ~ treat, method = "jitter", pch = 20,
           vertical = TRUE, add = TRUE, 
           col = c("#33CCCC","#003333","#CC9900","#996633"))
par(mai = c(0.6,1.4,0.1,0.1))
}
#mar/mai is bottom, L, top, R
#mai is margin inches, mar is margin lines

#means and sd
St14 = subset(data, water == "sterile" & temp == "14") 
St18 = subset(data, water == "sterile" & temp == "18") 
F14 = subset(data, water == "filtered" & temp == "14") 
F18 = subset(data, water == "filtered" & temp == "18") 

sizemeans <- aggregate(Length.measurement_.micron. ~ treat, data = data, mean)
sizesd <- aggregate(Length.measurement_.micron. ~ treat, data = data, sd)
#treat Length.measurement_.micron.
#1 filtered_14                    185.2403 20.92
#2 filtered_18                    201.9052 11.55
#3  sterile_14                    225.6127 21.31
#4  sterile_18                    226.1085 19.22

sizemeansW <- aggregate(Length.measurement_.micron. ~ water, data = data, mean)
sizesdW <- aggregate(Length.measurement_.micron. ~ water, data = data, sd)
#water Length.measurement_.micron.
#1 filtered                    190.6885 19.46
#2  sterile                    225.8128 20.32

sizemeansT <- aggregate(Length.measurement_.micron. ~ temp, data = data, mean)
sizesdT <- aggregate(Length.measurement_.micron. ~ temp, data = data, sd)
#temp Length.measurement_.micron.
#1   14                    205.1340 28.98
#2   18                    215.8221 20.24

areameans <- aggregate(Body_area_mm2 ~ treat, data = data, mean)
areasd <- aggregate(Body_area_mm2 ~ treat, data = data, sd)
#treat Body_area_mm2
#1 filtered_14      16.47343 1.87
#2 filtered_18      19.18706 1.25
#3  sterile_14      25.75324 3.37
#4  sterile_18      27.15826 3.50

areameansW <- aggregate(Body_area_mm2 ~ water, data = data, mean)
areasdW <- aggregate(Body_area_mm2 ~ water, data = data, sd)
#water Body_area_mm2
#1 filtered      17.36058 2.11
#2  sterile      26.32018 3.56

areameansT <- aggregate(Body_area_mm2 ~ temp, data = data, mean)
areasdT <- aggregate(Body_area_mm2 ~ temp, data = data, sd)
#temp Body_area_mm2
#1   14      21.04609 5.39
#2   18      23.77050 4.84

# variations of boxplots #####
head(data)
h2o <- subset(data, temp == "14")

length14boxplot <-{h2o$water <- factor(h2o$water, levels = c("sterile", "filtered"))
data$treat <- factor(data$treat, levels = c("sterile_14", "filtered_14"))
boxplot(Length.measurement_.micron. ~ water, data = h2o, frame = TRUE,
        ylab = "Length (µm)", xlab = "", 
        col = c("#99FFFF","#336666"), xaxt = 'n', cex.axis = 1.6, 
        cex.main = 2.4, cex.lab = 2)
axis(side = 1, tick = FALSE, at = 0:2, 
     labels = c("", paste("Sterile ASW\n14ºC"),
                paste("Microbiome-exposed\n14ºC")), 
     cex.axis = 2.2, line = 1.9) 
stripchart(data = h2o, Length.measurement_.micron. ~ water, method = "jitter", pch = 20,
           vertical = TRUE, add = TRUE, 
           col = c("#33CCCC","#003333"))
par(mai = c(0.8,1,0.1,0.1))
}

# x plots for size and area #### 
library('Rmisc')
treatColors = c("#33CCCC", "#336666", "#FFCC00", "#CC9933")
names(treatColors) = c("sterile_14", "filtered_14", 
                       "sterile_18", "filtered_18") 

waterShapes = c(16, 15)
names(waterShapes) = c("sterile", "filtered")

summ = summarySE(data = data, measurevar = "Length.measurement_.micron.", 
                 groupvars = c("treat", "temp", "water"))
pd <- position_dodge(0.2) 

ggplot(summ, aes(x = temp, y = Length.measurement_.micron., 
                 color = treat)) +
  geom_point(aes(shape= water),size=10,position=pd)+
  #geom_line(aes(group= water, color=temp),position=pd)+
  geom_errorbar(aes(ymin=Length.measurement_.micron.- sd,
                    ymax=Length.measurement_.micron.+ sd),
                lwd=1.2,width=0.2,position=pd)+
  scale_shape_manual(values = waterShapes)+
  scale_color_manual(values = treatColors)+
  theme_prism()+
  #facet_wrap(~ treat,scales="free_y", ncol=8)+
  theme(legend.text=element_text(size=10)) +
  theme(legend.key = element_blank())+
  theme(legend.position = 'right') + 
  ylab("Length (µm)") + 
  xlab(NULL) 

# new boxplot and key code #####
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

lengthPlot <- ggplot() + 
  theme_prism() +
  geom_boxplot(data = data, aes(
    x = water, y = Length.measurement_.micron., fill = treat), 
    show.legend = FALSE, width = 0.5, 
    position = position_dodge(width = 0.7))  + 
  geom_point(data = data, aes(
    x = water, y = Length.measurement_.micron., 
    shape = temp, color = treat),
    show.legend = FALSE, size = 3,
    position = position_jitterdodge(dodge.width = 0.7, 
                                    jitter.width = 0.5)) + 
  scale_fill_manual(values = treatColors) + 
  scale_color_manual(values = treatColors2) + 
  scale_shape_manual(values = tempShapes) + 
  scale_x_discrete(limits = orderList, labels = names(orderList)) +
  xlab(NULL) + 
  ylab("\nLength (µm)") + 
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) + 
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(7, "in"))
print(lengthPlot)

ggsave(plot = lengthPlot, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/lengthPlot.jpg", 
       scale = 1, width = 8, height = 6, units = c("in"), 
       dpi = 300)

areaPlot <- ggplot() + 
  theme_prism() +
  geom_boxplot(data = data, aes(
    x = water, y = Body_area_mm2, fill = treat), 
    show.legend = FALSE, width = 0.5, 
    position = position_dodge(width = 0.7))  + 
  geom_point(data = data, aes(
    x = water, y = Body_area_mm2, 
    shape = temp, color = treat),
    show.legend = FALSE, size = 3,
    position = position_jitterdodge(dodge.width = 0.7, 
                                    jitter.width = 0.5)) + 
  scale_fill_manual(values = treatColors) + 
  scale_color_manual(values = treatColors2) + 
  scale_shape_manual(values = tempShapes) + 
  scale_x_discrete(limits = orderList, labels = names(orderList)) +
  xlab(NULL) + 
  ylab("\nArea (mm\u00B2)") + #unicode for superscript 2 
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, colour = "black"),  # y-axis numbers size
        axis.title.y = element_text(size = 20),  # y-axis title size
        legend.position = "none",
        legend.title = element_blank()) + 
  force_panelsizes(rows = unit(5, "in"), 
                   cols = unit(7, "in"))
print(areaPlot)

ggsave(plot = areaPlot, 
       filename = "~/Desktop/Maroon_Debut/paper_fig_final/areaPlot.jpg", 
       scale = 1, width = 8, height = 6, units = c("in"), 
       dpi = 300)
  
#plotting only legend
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
       col = c("#99FFFF", "#FFCC00", "#336666", "#CC9933"), 
       cex = 2.5, bty = 'n', par(mai = c(0.1,0.1,0.1,0.1))
)   
