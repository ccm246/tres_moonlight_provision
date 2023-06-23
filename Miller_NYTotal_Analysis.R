##----------------------------------------------##
#   TOTAL FORAGING TIME & ILLUMINATION NY TRES   #
#     Colleen Miller, Cornell University         #
#                   2023                         #
##----------------------------------------------##

#### Data Set Up ####
# library packages that may be useful later on
library(lme4)
library(ggplot2)
library(sjPlot)
library(data.table)
library(MuMIn)
library(wesanderson)
library(nlme)
library(lubridate)
library(dplyr)
library(gapminder)
library(stringr)
library(data.table)
library(car)
library(emmeans)
library(pacman)

# N.B. Illumination data is in percentage and comes from timeanddate.com. 
# N.B. Cloud cover data is confined from 0 to 1 and comes from the Northeast Climate Center. 

setwd("/Users/colleenmiller/Dropbox/TRES LUNAR/Cloud Data")

LFM.ag.unique <- read.csv('ForagingTime_Male_Aggregated.csv')
LFF.ag.unique <- read.csv('ForagingTime_Female_Aggregated.csv')

#CLEAN UP COLUMNS
LFM.ag.unique$Sky.Cov.x <- NULL
LFF.ag.unique$Sky.Cov.x <- NULL
LFM.ag.unique$X <- NULL
LFF.ag.unique$X <- NULL

#UNIQUE DATA
LFF.ag.unique <- LFF.ag.unique[!duplicated(LFF.ag.unique),]
LFM.ag.unique <- LFM.ag.unique[!duplicated(LFM.ag.unique),]

#CHARACTER VARIABLE FOR YEAR AND SCALE REST OF VARIABLES
LFF.ag.unique$Ch.Year <- as.character(LFF.ag.unique$Sc.Year)
LFM.ag.unique$Ch.Year <- as.character(LFM.ag.unique$Sc.Year)
LFF.ag.unique$sc.cloud <- scale(LFF.ag.unique$cloud)
LFM.ag.unique$sc.cloud <- scale(LFM.ag.unique$cloud)
LFF.ag.unique$sc.Illumination <- scale(LFF.ag.unique$Illumination)
LFM.ag.unique$sc.Illumination <- scale(LFM.ag.unique$Illumination)

#### RUN MODELS FEMALES ####

#model suite
red1 <- lme(rank(totalTime) ~ 1, random=list(~1|RFID), data=LFF.ag.unique)
red2 <- lme(rank(totalTime) ~ Ch.Year + Sc.DOY  + sc.cloud, random=list(~1|RFID), data=LFF.ag.unique)
red3 <- lme(rank(totalTime) ~ Ch.Year + Sc.DOY  + sc.cloud + sc.Illumination, random=list(~1|RFID), data=LFF.ag.unique)
red4 <- lme(rank(totalTime) ~ Ch.Year + Sc.DOY  + sc.cloud + sc.Illumination + sc.cloud * sc.Illumination, random=list(~1|RFID), data=LFF.ag.unique)

#review
model.sel(red1, red2, red3, red4)  
summary(red4)
intervals(red4, level = 0.95)

#### PLOTTING TOTAL TIME FEMALES ####

#creating color scaled for illumination plotting
illa <- mean(LFF.ag.unique$sc.Illumination) + sd(LFF.ag.unique$sc.Illumination)
ill <- mean(LFF.ag.unique$sc.Illumination)
illb <- mean(LFF.ag.unique$sc.Illumination) - sd(LFF.ag.unique$sc.Illumination)

illa <- mean(LFF.ag.unique$Illumination) + sd(LFF.ag.unique$Illumination)
ill <- mean(LFF.ag.unique$Illumination)
illb <- mean(LFF.ag.unique$Illumination) - sd(LFF.ag.unique$Illumination)

#adjusting cloud cover
mylist <- list(sc.cloud = seq(-1.34, 2.17, by = 0.4), sc.Illumination = c(illa, ill, illb))
mylist <- list(cloud = seq(0, 1, by = 0.15), Illumination = c(illa, ill, illb))

#create colors
cols = c('navyblue', 'cornflowerblue', 'grey')

#base model
red4.1 <- lme(rank(totalTime) ~ Ch.Year + JDate  + cloud + Illumination + cloud * Illumination, random=list(~1|RFID), data=LFF.ag.unique)

#plot
plotmoon <- plot_model(red4.1, type = 'eff', terms = c('cloud', 'Illumination'), color=c('navyblue', 'cornflowerblue', 'grey')) + 
  theme_classic() + theme(text=element_text(size = 20)) + xlab('Proportion Cloud Cover') + ylab('Ranked Total Activity Time') 
p2 <- plotmoon + ggtitle('')
p2 

#### RUN MALE MODELS ####

#LOG TRANFORM AND FLIPPING DATA
LFM.ag.unique$flip_totaltime <- max(LFM.ag.unique$totalTime)+1-LFM.ag.unique$totalTime
LFM.ag.unique$logflip_totaltime <- log(LFM.ag.unique$flip_totaltime)
LFM.ag.unique$scale_logflip_totaltime <- scale(LFM.ag.unique$logflip_totaltime)

#model suite
red1 <- lme(rank(totalTime) ~ 1, random=list(~1|RFID), data=LFM.ag.unique)
red2 <- lme(rank(totalTime) ~ Ch.Year + Sc.DOY  + sc.cloud, random=list(~1|RFID), data=LFM.ag.unique)
red3 <- lme(rank(totalTime) ~ Ch.Year + Sc.DOY  + sc.cloud + sc.Illumination, random=list(~1|RFID), data=LFM.ag.unique)
red4 <- lme(rank(totalTime) ~ Ch.Year + Sc.DOY  + sc.cloud + sc.Illumination + sc.cloud * sc.Illumination, random=list(~1|RFID), data=LFM.ag.unique)

#review
model.sel(red1, red2, red3, red4)  
summary(red4)
intervals(red4, level = 0.95)

