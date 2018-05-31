
PARdata <- read.csv("PAR_rawdata.csv")
head(PARdata)

library(tidyverse)

# Average PAR data by date
PARdata.averaged <- PARdata %>%
  group_by(Date) %>%
  summarize(PAR = mean(PAR.mmol.m2.min))

PARdata.averaged$Date <- as.Date(PARdata.averaged$Date, "%m/%d/%y")
PARdata.averaged <- arrange(PARdata.averaged, Date)

# Convert PAR from mmol m-2 min-1 to umol m-2 s-1
PARdata.averaged$PAR <- PARdata.averaged$PAR * 1000/60

# Compute total radiation from PAR using factor of 0.45
# Source: Meek et al. 1984
PARdata.averaged <- mutate (PARdata.averaged, totalrad = PAR/0.45)

# Convert total radiation from umol m-2 s-1 to MJ m-2 d-1
# 1 J m-2 s-1 = 4.6 umol m-2 s-1
PARdata.averaged$totalrad <- PARdata.averaged$totalrad * (60*60*24) / (4.6*1000000)

# Find mean for each Julian day 
PARdata.averaged <- mutate(PARdata.averaged, JulianDay = format(PARdata.averaged$Date, "%j"))
MeanTotalRad <- aggregate(PARdata.averaged, list(JulianDay = PARdata.averaged$JulianDay), mean)

# Create list of all dates in modeled period
alldates <- data.frame(Date = seq(as.Date("1969-06-27"), as.Date("2016-12-31"), by = "days"))

TotalRad.df <- merge(x = alldates, y = PARdata.averaged, by = "Date", all.x = TRUE)
TotalRad.df <- TotalRad.df[, -c(2, 4)]
TotalRad.df <- mutate(TotalRad.df, JulianDay = format(TotalRad.df$Date, "%j"))


daily$daymonth <- strftime(daily$Date, format="%m-%d") 
TotalRad.df <- transform(TotalRad.df, meanrad =ave(totalrad, JulianDay,
                                  FUN=function(x) mean(x,na.rm=TRUE) ))
TotalRad.df <- transform(TotalRad.df, totalrad =ifelse(is.na(totalrad),meanrad,totalrad))

 
# Create new csv file with total radiation by date
write.csv(TotalRad.df[, 1:2], "TotalRadiation.csv", row.names = F)

