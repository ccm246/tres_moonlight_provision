##----------------------------------------------##
# Onset of Provisioning & Moonlight NY TRES      #
#     Colleen Miller, Cornell University         #
#                   2023                         #
##----------------------------------------------##

#### SET UP DATA ####
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
library(tidyverse)
library(cowplot)
library(sjPlot)
library(effects)
library(sjstats)
library(lmerTest)
library(lubridate)
library(chron)
library(emmeans)

setwd("/Users/colleenmiller/Dropbox/TRES LUNAR/Cloud Data")

#read in data
onsetfem.ag <- read.csv('Onset_Female_Aggregated.csv')

#CREATE LOG TRANSFORMED RESPONSE VARIABLE AND FINISH SCALING
onsetfem.ag$log_sec<- log(onsetfem.ag$Seconds)
onsetfem.ag$sc.cloud <- scale(onsetfem.ag$cloud)
onsetfem.ag$sc.illum <- scale(onsetfem.ag$Illumination)
onsetfem.ag$Ch.Year <- as.character(onsetfem.ag$Sc.Year)

# N.B. Illumination data is in percentage and comes from timeanddate.com. 
# N.B. Cloud cover data is confined from 0 to 1 and comes from the Northeast Climate Center. 

#### MODELING FEMALES ####

#models for females onset
red1 <- lme(log_sec ~ 1, random=list(~1|RFID), data=onsetfem.ag, method = 'ML')
red2 <- lme(log_sec ~ Ch.Year + Sc.DOY  + sc.cloud, random=list(~1|RFID), data=onsetfem.ag, method = 'ML')
red3 <- lme(log_sec ~ Ch.Year + Sc.DOY  + sc.cloud + sc.illum, random=list(~1|RFID), data=onsetfem.ag, method = 'ML')
red4 <- lme(log_sec ~ Ch.Year + Sc.DOY  + sc.cloud * sc.illum, random=list(~1|RFID), data=onsetfem.ag, method = 'ML')

model.sel(red1, red2, red3, red4)  
summary(red4)
intervals(red4, level = 0.95)

#models for females onset, DOY^2
red1 <- lme(log_sec ~ 1, random=list(~1|RFID), data=onsetfem.ag, method = 'ML')
red2 <- lme(log_sec ~ Ch.Year + Sc.DOY^2  + sc.cloud, random=list(~1|RFID), data=onsetfem.ag, method = 'ML')
red3 <- lme(log_sec ~ Ch.Year + Sc.DOY^2  + sc.cloud + sc.illum, random=list(~1|RFID), data=onsetfem.ag, method = 'ML')
red4 <- lme(log_sec ~ Ch.Year + Sc.DOY^2  + sc.cloud * sc.illum, random=list(~1|RFID), data=onsetfem.ag, method = 'ML')

model.sel(red1, red2, red3, red4)  
summary(red4)
intervals(red4, level = 0.95)

summary(red3)
intervals(red3, level = 0.95)


#### SET UP PLOTS ####

#to create sd points
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#create broader illumination value 
onsetfem.ag$ill5[onsetfem.ag$Illumination < 6] <- 5
onsetfem.ag$ill5[onsetfem.ag$Illumination > 5 & onsetfem.ag$Illumination < 10] <- 10
onsetfem.ag$ill5[onsetfem.ag$Illumination > 10 & onsetfem.ag$Illumination < 15] <- 15
onsetfem.ag$ill5[onsetfem.ag$Illumination > 15 & onsetfem.ag$Illumination < 20] <- 20
onsetfem.ag$ill5[onsetfem.ag$Illumination > 20 & onsetfem.ag$Illumination < 25] <- 25
onsetfem.ag$ill5[onsetfem.ag$Illumination > 25 & onsetfem.ag$Illumination < 30] <- 30
onsetfem.ag$ill5[onsetfem.ag$Illumination > 30 & onsetfem.ag$Illumination < 35] <- 35
onsetfem.ag$ill5[onsetfem.ag$Illumination > 35 & onsetfem.ag$Illumination < 40] <- 40
onsetfem.ag$ill5[onsetfem.ag$Illumination > 40 & onsetfem.ag$Illumination < 45] <- 45
onsetfem.ag$ill5[onsetfem.ag$Illumination > 45 & onsetfem.ag$Illumination < 50] <- 50
onsetfem.ag$ill5[onsetfem.ag$Illumination > 50 & onsetfem.ag$Illumination < 55] <- 55
onsetfem.ag$ill5[onsetfem.ag$Illumination > 55 & onsetfem.ag$Illumination < 60] <- 60
onsetfem.ag$ill5[onsetfem.ag$Illumination > 60 & onsetfem.ag$Illumination < 65] <- 65
onsetfem.ag$ill5[onsetfem.ag$Illumination > 65 & onsetfem.ag$Illumination < 70] <- 70
onsetfem.ag$ill5[onsetfem.ag$Illumination > 70 & onsetfem.ag$Illumination < 75] <- 75
onsetfem.ag$ill5[onsetfem.ag$Illumination > 75 & onsetfem.ag$Illumination < 80] <- 80
onsetfem.ag$ill5[onsetfem.ag$Illumination > 80 & onsetfem.ag$Illumination < 85] <- 85
onsetfem.ag$ill5[onsetfem.ag$Illumination > 85 & onsetfem.ag$Illumination < 90] <- 90
onsetfem.ag$ill5[onsetfem.ag$Illumination > 90 & onsetfem.ag$Illumination < 95] <- 95
onsetfem.ag$ill5[onsetfem.ag$Illumination > 95] <- 100
head(onsetfem.ag)
unique(onsetfem.ag$ill5)


#### PLOT ONSET FOR FEMALES ####
#output plot seconds with points, fit, se
summonsetfem <- summarySE(onsetfem.ag, measurevar='Seconds', groupvars = c('ill5'))
head(summonsetfem)
summonsetfem <- summonsetfem[summonsetfem$Seconds < 24000,]

pd <- position_dodge(0.1) # move them .05 to the left and right

#base model seconds
red3.1 <- lme(Seconds ~ Year + JDate^2  + cloud + Illumination, random=list(~1|RFID), data=onsetfem.ag)

#output plot seconds
ggplot(summonsetfem, aes(x=ill5, y=Seconds, colour='black')) + theme_classic() +
  geom_smooth(data = cbind(onsetfem.ag, pred=predict(red3.1)), method=lm, color='black', fill='cornflowerblue') +
  geom_errorbar(aes(ymin=Seconds-se, ymax=Seconds+se), width=.1, position=pd, color=alpha('black', .3)) +
  geom_point(position=pd, color=alpha('black', .3), cex=2) +
  labs(x="Illumination (%)", y="Onset of Activity (Seconds after Midnight)") + 
  theme(text = element_text(size = 20))

#output plot hours
summonsetfem <- summarySE(onsetfem.ag, measurevar='Hours', groupvars = c('ill5'))
head(summonsetfem)
summonsetfem <- summonsetfem[summonsetfem$Hours < 6.5,]

pd <- position_dodge(0.1) # move them .05 to the left and right

#base model hours
red3.1 <- lme(Hours ~ Year + JDate  + cloud + Illumination, random=list(~1|RFID), data=onsetfem.ag)

#output hours with points, fit, se
ggplot(summonsetfem, aes(x=ill5, y=Hours, colour='black')) + theme_classic() +
  geom_smooth(data = cbind(onsetfem.ag, pred=predict(red3.1)), method=lm, color='black', fill='cornflowerblue') +
  geom_errorbar(aes(ymin=Hours-se, ymax=Hours+se), width=.1, position=pd, color=alpha('black', .3)) +
  geom_point(position=pd, color=alpha('black', .3), cex=2) +
  labs(x="Illumination (%)", y="Onset of Activity (Hours after Midnight)") + 
  theme(text = element_text(size = 20))

#### MODELING MALES ####
rm(list=ls())

#set up for males

#setwd("~/Dropbox/TRES LUNAR/Cloud Data")

onsetmal.ag <- read.csv('Onset_Male_Aggregated.csv')
onsetmal.ag$Ch.Year <- as.character(onsetmal.ag$Sc.Year)

#CREATE LOG TRANSFORMED RESPONSE AND FINISH SCALING DATA
onsetmal.ag$log_sec <- log(onsetmal.ag$Seconds)
onsetmal.ag$sc.cloud <- scale(onsetmal.ag$cloud)
onsetmal.ag$sc.illum <- scale(onsetmal.ag$Illumination)

#EVALUATE ONLY SIMPLIFIED MODEL SUITE AND RUN MODELS
red1 <- lme(log_sec ~ 1, random=list(~1|RFID), data=onsetmal.ag, method = 'ML')
red2 <- lme(log_sec ~ Ch.Year + Sc.DOY  + sc.cloud, random=list(~1|RFID), data=onsetmal.ag, method = 'ML')
red3 <- lme(log_sec ~ Ch.Year + Sc.DOY  + sc.cloud + sc.illum, random=list(~1|RFID), data=onsetmal.ag, method = 'ML')
red4 <- lme(log_sec ~ Ch.Year + Sc.DOY  + sc.cloud * sc.illum, random=list(~1|RFID), data=onsetmal.ag, method = 'ML')

model.sel(red1, red2, red3, red4)  
summary(red3)
intervals(red3, level = 0.95)

#EVALUATE ONLY SIMPLIFIED MODEL SUITE AND RUN MODELS, DOY^2
red1 <- lme(log_sec ~ 1, random=list(~1|RFID), data=onsetmal.ag, method = 'ML')
red2 <- lme(log_sec ~ Ch.Year + Sc.DOY^2  + sc.cloud, random=list(~1|RFID), data=onsetmal.ag, method = 'ML')
red3 <- lme(log_sec ~ Ch.Year + Sc.DOY^2  + sc.cloud + sc.illum, random=list(~1|RFID), data=onsetmal.ag, method = 'ML')
red4 <- lme(log_sec ~ Ch.Year + Sc.DOY^2  + sc.cloud * sc.illum, random=list(~1|RFID), data=onsetmal.ag, method = 'ML')

model.sel(red1, red2, red3, red4)  
summary(red3)
intervals(red3, level = 0.95)

#### CLUTCH INITIATION ####

#using female data
head(onsetfem.ag)

#quick correlation
model <- lm(HDJ ~ Illumination, onsetfem.ag)
summary(model)

#quick plot
plot(HDJ ~ Illumination, onsetfem.ag, pch=19, col = alpha('black', .7), xlab = 'Illumination (%)', ylab = 'Nest Hatch Date (DOY)')

####### MOONLIGHT YEAR VARIATION PLOT ######

head(onsetfem.ag)

ggplot(data = onsetfem.ag, aes(log_sec*Illumination*Year)) +
  geom_point(color='steelblue')

tiff("moonyearsuppfigure.tif", units="in", width=5, height=4, res=600)
q <- ggplot(onsetfem.ag, aes(x=Illumination, y=log_sec)) + geom_point(color = alpha('steelblue', .5)) + theme_classic() + theme(text = element_text(size=8))+ xlab('Illumination (%)') + ylab('Provisioning Onset (Log Seconds from Midnight)')
q + facet_wrap(~Year) + theme(text = element_text(size=8), axis.text.x = element_text(angle = 50, hjust = 1))
dev.off() 
