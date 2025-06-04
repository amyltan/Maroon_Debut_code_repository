setwd("~/Desktop/R/morphology-MDII/")

library("dplyr")
library('stringr')
library('rstatix')
library('ggplot2')

area = read.csv("Rsheet_Morphology_MDII_RawData_6May24/2D_body_area-Table_1.csv")
CtG = read.csv("Rsheet_Morphology_MDII_RawData_6May24/Cell-to-gut_distance-Table_1.csv")
gut = read.csv("Rsheet_Morphology_MDII_RawData_6May24/Gut_wall_thickness-Table_1.csv")
lenCell = read.csv("Rsheet_Morphology_MDII_RawData_6May24/Length_pigment-cells-Table_1.csv")

###make dataframes match and combining into a single dataframe
area <- area %>% 
  rename(
    Sample_name = Photo_ID
  )
area$Sample_name <- sub("-_", "-", area$Sample_name)

CtG <- CtG %>% 
  rename(
    Sample_name = Photo_ID,
    Mean_CtG = Mean_distance,
    StDev_CtG = StDev
  )
CtG$Sample_name <- sub("-_", "-", CtG$Sample_name)

gut <- gut %>% 
  rename(
    Sample_name = Photo_ID
  )
gut$Sample_name <- sub("-_", "-", gut$Sample_name)

lenCell$Sample_name <- sub("20x", "25x", lenCell$Sample_name)
#necessary because of naming difference between the sheets. Otherwise you lose C9's.

count(area) #165
count(lenCell) #337
count(CtG) #165
count(gut) #165
lenCellArea <- merge(lenCell, area, by = "Sample_name", all = FALSE)
count(lenCellArea) #157

LCAgut <- merge (lenCellArea, gut, by = "Sample_name", all = FALSE)
count(LCAgut) #157

data <- merge(LCAgut, CtG, by = "Sample_name", all = FALSE)
count(data) #128

#Note: 128 overlap between the cell count and lengths files and the other 3 measurements
#Might raise question of whether we should merge and then do outliers or other way round 

data <- subset(data, select = -c(Chamber.x, Slice.x, Chamber.y, Slice.y, Slice, 
                                 CtG1, CtG2, CtG3, CtG4, CtG5, CtG6, 
                                 CtG7, CtG8, CtG9, CtG10, CtG11, CtG12, CtG13, 
                                 CtG14, CtG15, CtG16, CtG17, CtG18, CtG19, 
                                 CtG20, CtG21, CtG22, CtG23, CtG24, 
                                 Gut_thickness_.1, Gut_thickness_.2,
                                 Gut_thickness_.3))

data <- data %>%
  mutate(water = case_when(
    str_detect(Chamber, "St") ~ "sterile",
    str_detect(Chamber, "F") ~ "filtered" 
  )) %>%
  mutate(temp = case_when(
    str_detect(Chamber, "14") ~ "14",
    str_detect(Chamber, "18") ~ "18"
  ))
data = mutate(data, treat = paste(water, temp, sep = "_"))

unique(data$Chamber) #confirm that all 12 chambers are present

cell_summary <- data %>% 
  group_by(Chamber) %>% 
  get_summary_stats(Pigment_cell_count, type = "mean_sd")
#Chamber variable               n  mean    sd
#<chr>   <fct>              <dbl> <dbl> <dbl>
#  1 C10St14 Pigment_cell_count    10  50.9  5.51
#2 C10St18 Pigment_cell_count     6  46.5  5.36
#3 C11F14  Pigment_cell_count    15  47.8  8.78
#4 C11F18  Pigment_cell_count    14  12.9 17.2 
#5 C12F14  Pigment_cell_count    14  38.4  7.93
#6 C12F18  Pigment_cell_count    14  40.3  7.00
#7 C13F14  Pigment_cell_count    15  39.5  4.64
#8 C13F18  Pigment_cell_count    15  56.7 10.6 
#9 C8St14  Pigment_cell_count    13  60.7  9.23
#10 C8St18  Pigment_cell_count    14  54.7  9.04
#11 C9St14  Pigment_cell_count    15  65.6 19.2 
#12 C9St18  Pigment_cell_count    12  54.2 11.2

length_summary <- data %>% 
  group_by(Chamber) %>% 
  get_summary_stats(Length.measurement_.micron., type = "mean_sd")

area_summary <- data %>% 
  group_by(Chamber) %>% 
  get_summary_stats(Body_area_um2, type = "mean_sd")

gut_summary <- data %>% 
  group_by(Chamber) %>% 
  get_summary_stats(Mean_gut_thickness, type = "mean_sd")

CtG_summary <- data %>% 
  group_by(Chamber) %>% 
  get_summary_stats(Mean_CtG, type = "mean_sd")

### Determining outliers from each stat
gut_out <- data %>% 
  group_by(Chamber) %>% identify_outliers(Mean_gut_thickness)
count(gut_out) #5 
#[1] "C11F18-5dpf-MDII-25x-4u-9"  "C12F18-5dpf-MDII-25x-4u-15"
#[3] "C12F18-5dpf-MDII-25x-4u-8"  "C13F14-5dpf-MDII-25x-4u-5" 
#[5] "C13F18-5dpf-MDII-25x-4u-12"

pigment_out <- data %>% 
  group_by(Chamber) %>% identify_outliers(Pigment_cell_count)
count(pigment_out) #6 
#[1] "C13F14-5dpf-MDII-25x-4u-9"  "C13F18-5dpf-MDII-25x-4u-11"
#[3] "C13F18-5dpf-MDII-25x-4u-14" "C8St14-5dpf-MDII-25x-4u-8" 
#[5] "C8St18-5dpf-MDII-25x-4u-14" "C8St18-5dpf-MDII-25x-4u-4"

length_out <- data %>% 
  group_by(Chamber) %>% identify_outliers(Length.measurement_.micron.)
count(length_out) #6 
#[1] "C10St14-5dpf-MDII-25x-4u-8" "C11F14-5dpf-MDII-25x-4u-12"
#[3] "C11F14-5dpf-MDII-25x-4u-15" "C12F18-5dpf-MDII-25x-4u-7" 
#[5] "C8St18-5dpf-MDII-25x-4u-10" "C9St18-5dpf-MDII-25x-4u-15"

CtG_out <- data %>% 
  group_by(Chamber) %>% identify_outliers(Mean_CtG)
count(CtG_out) #9 
#[1] "C12F14-5dpf-MDII-25x-4u-12" "C12F14-5dpf-MDII-25x-4u-2" 
#[3] "C13F18-5dpf-MDII-25x-4u-2"  "C13F18-5dpf-MDII-25x-4u-4" 
#[5] "C13F18-5dpf-MDII-25x-4u-5"  "C8St18-5dpf-MDII-25x-4u-2" 
#[7] "C8St18-5dpf-MDII-25x-4u-9"  "C9St14-5dpf-MDII-25x-4u-12"
#[9] "C9St14-5dpf-MDII-25x-4u-14"

area_out <- data %>% 
  group_by(Chamber) %>% identify_outliers(Body_area_um2)
count(area_out) #12
#[1] "C11F14-5dpf-MDII-25x-4u-12" "C11F14-5dpf-MDII-25x-4u-15"
#[3] "C11F14-5dpf-MDII-25x-4u-2"  "C11F14-5dpf-MDII-25x-4u-8" 
#[5] "C12F14-5dpf-MDII-25x-4u-1"  "C13F18-5dpf-MDII-25x-4u-14"
#[7] "C13F18-5dpf-MDII-25x-4u-15" "C13F18-5dpf-MDII-25x-4u-16"
#[9] "C13F18-5dpf-MDII-25x-4u-6"  "C8St18-5dpf-MDII-25x-4u-12"
#[11] "C9St18-5dpf-MDII-25x-4u-11" "C9St18-5dpf-MDII-25x-4u-3"

out1 <- merge(gut_out, pigment_out, by = "Sample_name", all = TRUE)
out1 <- subset(out1, select = c(Sample_name))
out2 <- merge(out1, length_out, by = "Sample_name", all = TRUE)
out2 <- subset(out2, select = c(Sample_name))
out3 <- merge(out2, CtG_out, by = "Sample_name", all = TRUE)
out3 <- subset(out3, select = c(Sample_name))
outliers <- merge(out3, area_out, by = "Sample_name", all = TRUE)
outliers <- subset(outliers, select = c(Sample_name))

#lines to confirm that all outlier names got moved to unified list
#returns 1 if present, 0 if absent
as.integer(gut_out$Sample_name %in% outliers$Sample_name) #5
as.integer(pigment_out$Sample_name %in% outliers$Sample_name) #6
as.integer(length_out$Sample_name %in% outliers$Sample_name) #6
as.integer(CtG_out$Sample_name %in% outliers$Sample_name) #9
as.integer(area_out$Sample_name %in% outliers$Sample_name) #12

#remove outliers from datasheet
cleandata <- data[!(data$Sample_name %in% outliers$Sample_name),] #122 remain 

#export clean data as .csv 
write.csv(cleandata, "morphology-MDII-clean-allMetrics.csv")

## ALL data (including outliers) for comparison to clean data 
### Visual checks
#pigment cell counts
data %>% 
  ggplot(aes(x = treat, y = Pigment_cell_count, colour = Chamber)) + 
  geom_jitter(size = 2) + 
  scale_color_brewer(palette = "Paired") +
  ggtitle("ALL data")

data %>% 
  ggplot(aes(x = treat, y = Length.measurement_.micron., colour = Chamber)) + 
  geom_jitter(size = 2) + 
  scale_color_brewer(palette = "Paired") + 
  ggtitle("ALL data")

data %>% 
  ggplot(aes(x = treat, y = Body_area_um2, colour = Chamber)) + 
  geom_jitter(size = 2) + 
  scale_color_brewer(palette = "Paired") +
  ggtitle("ALL data")

data %>% 
  ggplot(aes(x = treat, y = Mean_gut_thickness, colour = Chamber)) + 
  geom_jitter(size = 2) + 
  scale_color_brewer(palette = "Paired") + 
  ggtitle("ALL data")

data %>% 
  ggplot(aes(x = treat, y = Mean_CtG, colour = Chamber)) + 
  geom_jitter(size = 2) + 
  scale_color_brewer(palette = "Paired") +
  ggtitle("ALL data")

## Clean data - checking for large differences between groups 
### Visual checks
#pigment cell counts
cleandata %>% 
  ggplot(aes(x = treat, y = Pigment_cell_count, colour = Chamber)) + 
  geom_jitter(size = 2) + 
  scale_color_brewer(palette = "Paired")
#C11F18 has two distinct groupings 
#C9St14 has large spread

cleandata %>% 
  ggplot(aes(x = treat, y = Length.measurement_.micron., colour = Chamber)) + 
  geom_jitter(size = 2) + 
  scale_color_brewer(palette = "Paired")
#everything looks pretty equally spread 

cleandata %>% 
  ggplot(aes(x = treat, y = Body_area_um2, colour = Chamber)) + 
  geom_jitter(size = 2) + 
  scale_color_brewer(palette = "Paired")
#even spread

cleandata %>% 
  ggplot(aes(x = treat, y = Mean_gut_thickness, colour = Chamber)) + 
  geom_jitter(size = 2) + 
  scale_color_brewer(palette = "Paired")
#C12F14 might have a larger variance than other F14 

cleandata %>% 
  ggplot(aes(x = treat, y = Mean_CtG, colour = Chamber)) + 
  geom_jitter(size = 2) + 
  scale_color_brewer(palette = "Paired")
#C11F18 is probably an issue again

# Sterile 14 -- keep all
St14 = subset(cleandata, treat == "sterile_14")
#Shapiro-Wilk normality: p>0.5 is NOT sig diff from normal
shapiro.test(St14$Pigment_cell_count) #0.01618 NO 
shapiro.test(St14$Length.measurement_.micron.) #0.6165 YES 
shapiro.test(St14$Body_area_um2) #0.4199 YES
shapiro.test(St14$Mean_gut_thickness) #0.01662 NO 
shapiro.test(St14$Mean_CtG) #0.8959 YES

#Levene's variance:p>0.05 is NOT sig diff (ARE homogeneous)
levene_test(St14, Pigment_cell_count ~ Chamber) #0.00223 NO
levene_test(St14, Length.measurement_.micron. ~ Chamber) #0.439 YES 
levene_test(St14, Body_area_um2 ~ Chamber) #0.693 YES
levene_test(St14, Mean_gut_thickness ~ Chamber) #0.885 YES 
levene_test(St14, Mean_CtG ~ Chamber) #0.361 YES

St14pig <- aov(data = St14, Pigment_cell_count ~ Chamber)
summary(St14pig)
TukeyHSD(St14pig)
#C9 and C10 are sig diff from each other, but neither from C8
#probably keep all

St14len <- aov(data = St14, Length.measurement_.micron. ~ Chamber)
summary(St14len)
TukeyHSD(St14len)
#all same

St14area <- aov(data = St14, Body_area_um2 ~ Chamber)
summary(St14area)
TukeyHSD(St14area) 
#C8 and C10 different but neither from C9 
#probably keep all

St14gut <- aov(data = St14, Mean_gut_thickness ~ Chamber)
summary(St14gut)
TukeyHSD(St14gut) 
#all same

St14ctg <- aov(data = St14, Mean_CtG ~ Chamber)
summary(St14ctg)
TukeyHSD(St14ctg)
#all same

# Sterile 18 --probably remove C8St18
St18 = subset(cleandata, treat == "sterile_18")
#Shapiro-Wilk normality: p>0.5 is NOT sig diff from normal
shapiro.test(St18$Pigment_cell_count) #0.05105 YES 
shapiro.test(St18$Length.measurement_.micron.) #0.3555 YES 
shapiro.test(St18$Body_area_um2) #0.5788 YES
shapiro.test(St18$Mean_gut_thickness) #0.756 YES 
shapiro.test(St18$Mean_CtG) #0.7482 YES

#Levene's variance:p>0.05 is NOT sig diff (ARE homogeneous)
levene_test(St18, Pigment_cell_count ~ Chamber) #0.178 YES
levene_test(St18, Length.measurement_.micron. ~ Chamber) #0.588 YES 
levene_test(St18, Body_area_um2 ~ Chamber) #0.973 YES
levene_test(St18, Mean_gut_thickness ~ Chamber) #0.0914 YES 
levene_test(St18, Mean_CtG ~ Chamber) #0.690 YES

St18pig <- aov(data = St18, Pigment_cell_count ~ Chamber)
summary(St18pig)
TukeyHSD(St18pig)
#all same

St18len <- aov(data = St18, Length.measurement_.micron. ~ Chamber)
summary(St18len)
TukeyHSD(St18len)
#C8 sig diff from both

St18area <- aov(data = St18, Body_area_um2 ~ Chamber)
summary(St18area)
TukeyHSD(St18area) 
#C8 sig diff from both

St18gut <- aov(data = St18, Mean_gut_thickness ~ Chamber)
summary(St18gut)
TukeyHSD(St18gut) 
#all same

St18ctg <- aov(data = St18, Mean_CtG ~ Chamber)
summary(St18ctg)
TukeyHSD(St18ctg)
#C8 sig diff from both

#F14 --maybe remove C11F14?
F14 = subset(cleandata, treat == "filtered_14")
#Shapiro-Wilk normality: p>0.5 is NOT sig diff from normal
shapiro.test(F14$Pigment_cell_count) #0.3633 YES 
shapiro.test(F14$Length.measurement_.micron.) #0.7829 YES 
shapiro.test(F14$Body_area_um2) #0.1437 YES
shapiro.test(F14$Mean_gut_thickness) #0.8876 YES 
shapiro.test(F14$Mean_CtG) #0.04932 NO

#Levene's variance:p>0.05 is NOT sig diff (ARE homogeneous)
levene_test(F14, Pigment_cell_count ~ Chamber) #0.298 YES
levene_test(F14, Length.measurement_.micron. ~ Chamber) #0.675 YES 
levene_test(F14, Body_area_um2 ~ Chamber) #0.0957 YES
levene_test(F14, Mean_gut_thickness ~ Chamber) #0.0667 YES 
levene_test(F14, Mean_CtG ~ Chamber) #0.321 YES

F14pig <- aov(data = F14, Pigment_cell_count ~ Chamber)
summary(F14pig)
TukeyHSD(F14pig)
#C11 sig diff from both

F14len <- aov(data = F14, Length.measurement_.micron. ~ Chamber)
summary(F14len)
TukeyHSD(F14len)
#all same

F14area <- aov(data = F14, Body_area_um2 ~ Chamber)
summary(F14area)
TukeyHSD(F14area) 
#C12 and C11 diff but neither from C13 
#probably keep all

F14gut <- aov(data = F14, Mean_gut_thickness ~ Chamber)
summary(F14gut)
TukeyHSD(F14gut) 
#C13 and C11 diff but neither from C12 
#probably keep all

F14ctg <- aov(data = F14, Mean_CtG ~ Chamber)
summary(F14ctg)
TukeyHSD(F14ctg)
#all same

#F18 --should possibly remove C11F18
F18 = subset(cleandata, treat == "filtered_18")
#Shapiro-Wilk normality: p>0.5 is NOT sig diff from normal
shapiro.test(F18$Pigment_cell_count) #0.005821 NO 
shapiro.test(F18$Length.measurement_.micron.) #0.2665 YES 
shapiro.test(F18$Body_area_um2) #0.4367 YES
shapiro.test(F18$Mean_gut_thickness) #0.7267 YES 
shapiro.test(F18$Mean_CtG) #0.0006887 NO

#Levene's variance:p>0.05 is NOT sig diff (ARE homogeneous)
levene_test(F18, Pigment_cell_count ~ Chamber) #0.225 YES
levene_test(F18, Length.measurement_.micron. ~ Chamber) #0.415 YES 
levene_test(F18, Body_area_um2 ~ Chamber) #0.161 YES
levene_test(F18, Mean_gut_thickness ~ Chamber) #0.879 YES 
levene_test(F18, Mean_CtG ~ Chamber) #0.00609 NO

F18pig <- aov(data = F18, Pigment_cell_count ~ Chamber)
summary(F18pig)
TukeyHSD(F18pig)
#C11 diff from both (C12 and 13 also diff, actually, but way closer)

F18len <- aov(data = F18, Length.measurement_.micron. ~ Chamber)
summary(F18len)
TukeyHSD(F18len)
#all same

F18area <- aov(data = F18, Body_area_um2 ~ Chamber)
summary(F18area)
TukeyHSD(F18area) 
#all same

F18gut <- aov(data = F18, Mean_gut_thickness ~ Chamber)
summary(F18gut)
TukeyHSD(F18gut) 
#all same

F18ctg <- aov(data = F18, Mean_CtG ~ Chamber)
summary(F18ctg)
TukeyHSD(F18ctg)
#C11 sig diff from both

##Remove chambers
badC = c("C8St18", "C11F14", "C11F18")
altdata <- cleandata[!(cleandata$Chamber %in% badC),]

altdata %>% 
  ggplot(aes(x = treat, y = Pigment_cell_count, colour = Chamber)) + 
  geom_jitter(size = 2) + 
  scale_color_brewer(palette = "Paired") +
  ggtitle("9 chambers")

altdata %>% 
  ggplot(aes(x = treat, y = Length.measurement_.micron., colour = Chamber)) + 
  geom_jitter(size = 2) + 
  scale_color_brewer(palette = "Paired") +
  ggtitle("9 chambers")

altdata %>% 
  ggplot(aes(x = treat, y = Body_area_um2, colour = Chamber)) + 
  geom_jitter(size = 2) + 
  scale_color_brewer(palette = "Paired") +
  ggtitle("9 chambers")

altdata %>% 
  ggplot(aes(x = treat, y = Mean_gut_thickness, colour = Chamber)) + 
  geom_jitter(size = 2) + 
  scale_color_brewer(palette = "Paired") +
  ggtitle("9 chambers")

altdata %>% 
  ggplot(aes(x = treat, y = Mean_CtG, colour = Chamber)) + 
  geom_jitter(size = 2) + 
  scale_color_brewer(palette = "Paired") +
  ggtitle("9 chambers")

length(St14$Sample_name) #n=34
length(F14$Sample_name) #n=35
length(St18$Sample_name) #n=23 
length(F18$Sample_name) #n=30 

St14_alt = subset(altdata, treat == "sterile_14")
length(St14_alt$Sample_name) #n=34
F14_alt = subset(altdata, treat == "filtered_14")
length(F14_alt$Sample_name) #24
St18_alt = subset(altdata, treat == "sterile_18")
length(St18_alt$Sample_name) #15
F18_alt = subset(altdata, treat == "filtered_18")
length(F18_alt$Sample_name) #n=17

###NOTE: my inclination is to use the 9 chambers (or maybe 11?) 
#but then lump the other chambers from each condition and NOT use chamber as a variable

##save version of file with 9 chambers 
write.csv(altdata, "morphology-9chambers.csv")
