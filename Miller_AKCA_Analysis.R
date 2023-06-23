##----------------------------------------------##
#     ANALYSIS OF CA AND AK DATA SET             #
#     Colleen Miller, Cornell University         #
#                 2023                           #
##----------------------------------------------##

#### READ IN DATA ####

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
library(raster)
library(prism)
library(rnoaa)
library(lmerTest)
library(lubridate)
library(chron)
library(emmeans)
library(dabestr)

#setwd("/Users/colleenmiller/Dropbox/TRES LUNAR/Cloud Data/")

akcatot <- read.csv('total_alaskacali.csv', header=T)

#create scaled and tranformed variables
akcatot$Illumination.sc <- scale(akcatot$Illumination)
akcatot$dayofyear.sc <- scale(akcatot$dayofyear)
akcatot$starttime.sqrt <- sqrt(akcatot$starttime)
akcatot$starttime.sqrtsc <- scale(akcatot$starttime.sqrt)
akcatot <- na.omit(akcatot)
akcatot$Sunrisetime.sc <- scale(akcatot$Sunrise)
akcatot$Sunsettime.sc <-scale(akcatot$Sunset)

head(akcatot)

#### MODEL ONSET DATA ####

#run series of models
m1 <- lmer(starttime.sqrtsc ~ 1 + (1|NestNumber), data=akcatot)
m2 <- lmer(starttime.sqrtsc ~ site + year + Sunrisetime.sc  + (1|NestNumber), data=akcatot)
m3 <- lmer(starttime.sqrtsc ~ site + year + Sunrisetime.sc   + Illumination.sc + (1|NestNumber), data=akcatot)
m4 <- lmer(starttime.sqrtsc ~ site + year + Sunrisetime.sc   + Illumination.sc + site*Illumination.sc + (1|NestNumber), data=akcatot)

#look into results
model.sel(m1, m2, m3, m4) 
summary(m4)
anova(m4)
confint(m4, level=0.95)
summary(m2)
anova(m2)
confint(m2, level=0.95)

#### SET UP ONSET PLOTS ####

#to create plotted sd points
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
akcatot$ill5[akcatot$Illumination < 6] <- 5
akcatot$ill5[akcatot$Illumination > 5 & akcatot$Illumination < 10] <- 10
akcatot$ill5[akcatot$Illumination > 10 & akcatot$Illumination < 15] <- 15
akcatot$ill5[akcatot$Illumination > 15 & akcatot$Illumination < 20] <- 20
akcatot$ill5[akcatot$Illumination > 20 & akcatot$Illumination < 25] <- 25
akcatot$ill5[akcatot$Illumination > 25 & akcatot$Illumination < 30] <- 30
akcatot$ill5[akcatot$Illumination > 30 & akcatot$Illumination < 35] <- 35
akcatot$ill5[akcatot$Illumination > 35 & akcatot$Illumination < 40] <- 40
akcatot$ill5[akcatot$Illumination > 40 & akcatot$Illumination < 45] <- 45
akcatot$ill5[akcatot$Illumination > 45 & akcatot$Illumination < 50] <- 50
akcatot$ill5[akcatot$Illumination > 50 & akcatot$Illumination < 55] <- 55
akcatot$ill5[akcatot$Illumination > 55 & akcatot$Illumination < 60] <- 60
akcatot$ill5[akcatot$Illumination > 60 & akcatot$Illumination < 65] <- 65
akcatot$ill5[akcatot$Illumination > 65 & akcatot$Illumination < 70] <- 70
akcatot$ill5[akcatot$Illumination > 70 & akcatot$Illumination < 75] <- 75
akcatot$ill5[akcatot$Illumination > 75 & akcatot$Illumination < 80] <- 80
akcatot$ill5[akcatot$Illumination > 80 & akcatot$Illumination < 85] <- 85
akcatot$ill5[akcatot$Illumination > 85 & akcatot$Illumination < 90] <- 90
akcatot$ill5[akcatot$Illumination > 90 & akcatot$Illumination < 95] <- 95
akcatot$ill5[akcatot$Illumination > 95] <- 100

head(akcatot)


#### ACTUAL ONSET PLOTS ####

#output point  plot AK
summonsetak <- summarySE(akcatot[akcatot$site =='AK',], measurevar='starttime', groupvars = c('ill5'))
head(summonsetak)

#adding in a dodge
pd <- position_dodge(0.1) # move them .05 to the left and right

#create point/se plot
akonspts <- ggplot(summonsetak, aes(x=ill5, y=starttime, colour='black')) + theme_classic() +
  coord_cartesian(xlim =c(0, 100), ylim = c(17000, 26000)) + 
  geom_errorbar(aes(ymin=starttime-se, ymax=starttime+se), width=.1, position=pd, color=alpha('navy', .4)) +
  geom_point(position=pd, color=alpha('navy', .4), cex=2) +
  labs(x="Illumination (%)", y="Onset of Activity (Seconds after Midnight)") + 
  theme(text = element_text(size = 20))+
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  ) 

#output
akonspts
ggsave('onsetakpts.png', akonspts, bg='transparent', dpi=600)

#output point plot CA
summonsetca <- summarySE(akcatot[akcatot$site =='CA',], measurevar='starttime', groupvars = c('ill5'))
head(summonsetca)
summonsetca <- na.omit(summonsetca)

#adding in a dodge
pd <- position_dodge(0.1) # move them .05 to the left and right

#create point/se plot
caonspts <- ggplot(summonsetca, aes(x=ill5, y=starttime, colour='black')) + theme_classic() +
  coord_cartesian(xlim =c(0, 100), ylim = c(17000, 26000)) + 
  geom_errorbar(aes(ymin=starttime-se, ymax=starttime+se), width=.1, position=pd, color=alpha('chocolate2', .4)) +
  geom_point(position=pd, color=alpha('chocolate2', .4), cex=2) +
  labs(x="Illumination (%)", y="Onset of Activity (Seconds after Midnight)") + 
  theme(text = element_text(size = 20))+
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  ) 

#output
caonspts
ggsave('onsetcapts.png', caonspts, bg='transparent', dpi=600)

#interaction plot overlay

#base model
m5.1 <- lmer(starttime ~ site + year + Sunrise + Illumination + site*Illumination + (1|NestNumber), data=akcatot)

#initial plot
plotmoon <- plot_model(m5.1, type = 'pred', terms = c('Illumination', 'site'), color=c('navyblue','chocolate1')) + 
  theme_classic() + theme(text=element_text(size = 20)) + xlab('Illumination (%)') + ylab('Onset of Activity (Seconds From Midnight)') 

#customization
p2 <- plotmoon +  scale_color_manual(values = c('navyblue', 'chocolate1'), name='Location') +  scale_fill_manual(values = c('navyblue', 'chocolate1'), name='Location') + ggtitle('') + theme(legend.position = "none") +
  coord_cartesian(xlim =c(0, 100), ylim = c(17000, 26000)) + 
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )

#output
p2
ggsave('onsetakcaint.png', p2, bg='transparent', dpi=600)

#### MODEL TOTAL DATA ####

#recreate variable
akcatot$totaltime <- akcatot$endtime - akcatot$starttime

#create transformed and scaled variables
akcatot$totaltime.sq <- (akcatot$totaltime)^2
akcatot$totaltime.sqsc  <- scale(akcatot$totaltime.sq)
akcatot$totaltime.log <- log(akcatot$totaltime)
akcatot$totaltime.sc <- scale(akcatot$totaltime)

#mixed model suite
m1 <- lmer(totaltime.sc ~ 1 + (1|NestNumber), data=akcatot)
m2 <- lmer(totaltime.sc ~ site + year + Sunrisetime.sc  + (1|NestNumber), data=akcatot)
m3 <- lmer(totaltime.sc ~ site + year + Sunrisetime.sc   + Illumination.sc + (1|NestNumber), data=akcatot)
m4 <- lmer(totaltime.sc ~ site + year + Sunrisetime.sc   + Illumination.sc + site*Illumination.sc + (1|NestNumber), data=akcatot)

#check results
model.sel(m1, m2, m3, m4)
summary(m4)#top 
confint(m4)
anova(m4)
summary(m2)#top 
confint(m2)

#### ACTUAL TOTAL PLOTS ####

#point plot AK
summonsetak <- summarySE(akcatot[akcatot$site=='AK',], measurevar='totaltime', groupvars = c('ill5'))
head(summonsetak)

#add in a dodge
pd <- position_dodge(0.1) # move them .05 to the left and right

#plotting AK points/SE
akpts <- ggplot(summonsetak, aes(x=ill5, y=totaltime, colour='black')) + theme_classic() +
  coord_cartesian(xlim =c(0, 100), ylim = c(50000, 65000)) + 
  geom_errorbar(aes(ymin=totaltime-se, ymax=totaltime+se), width=.1, position=pd, color=alpha('navy', .4)) +
  geom_point(position=pd, color=alpha('navy', .4), cex=2) +
  labs(x="Illumination (%)", y="Onset of Activity (Seconds after Midnight)") + 
  theme(text = element_text(size = 20))+
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )

#output
akpts
ggsave('totalakpts.png', akpts, bg='transparent', dpi=600)

#output point plot CA
summonsetca <- summarySE(akcatot[akcatot$site=='CA',], measurevar='totaltime', groupvars = c('ill5'))
head(summonsetca)

#add in a dodge
pd <- position_dodge(0.1) # move them .05 to the left and right

#point/SE plot CA
capts <- ggplot(summonsetca, aes(x=ill5, y=totaltime, colour='black')) + theme_classic() +
  coord_cartesian(xlim =c(0, 100), ylim = c(50000, 65000)) + 
  geom_errorbar(aes(ymin=totaltime-se, ymax=totaltime+se), width=.1, position=pd, color=alpha('chocolate2', .4)) +
  geom_point(position=pd, color=alpha('chocolate2', .4), cex=2) +
  labs(x="Illumination (%)", y="Onset of Activity (Seconds after Midnight)") + 
  theme(text = element_text(size = 20)) +
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  ) 

#output
capts
ggsave('totalcapts.png', capts, bg='transparent', dpi=600)

#interaction plot overlay

#base model
m5.2 <- lmer(totaltime ~ site + year + Sunrise + Illumination + site*Illumination + (1|NestNumber), data=akcatot)

#initial plot
plotmoon <- plot_model(m5.2, type = 'pred', terms = c('Illumination', 'site'), color=c('navyblue','chocolate1')) + 
  theme_classic() + theme(text=element_text(size = 20)) + xlab('Illumination (%)') + ylab('Total Activity Time (Seconds)') 

#customization
p2 <- plotmoon +  scale_color_manual(values = c('navyblue', 'chocolate1'), name='Location') +  scale_fill_manual(values = c('navyblue', 'chocolate1'), name='Location') + ggtitle('') + theme(legend.position = "none") +
  coord_cartesian(xlim =c(0, 100), ylim = c(50000, 65000)) + 
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )

#output
p2
ggsave('totalakcaint.png', p2, bg='transparent', dpi=600)




