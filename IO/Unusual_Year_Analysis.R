#### Lake 227 observation analysis ----
# This script examines patterns in observations for ELA Lake 227 data and associated weather and inflow data.
# The purpose is to identify years with highest and lowest cumulative PP observations 
# and to examine potential connections with unusual weather and inflow patterns.

#### Setup ----

library(ggplot2)
library(tidyverse)
library(zoo)
library(dplyr)
library(plyr)
theme_std <- function (base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.ticks = element_line(colour = "black", size = 1), 
          legend.key = element_rect(colour = "white"), 
          panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA, colour = NA),
          axis.line = element_line(size = 0.5, colour = "black"),
          panel.grid.major = element_line(NA), 
          panel.grid.minor = element_line(NA), 
          strip.background = element_rect(fill = "grey80", colour = "grey50", size = 0.2),
          axis.text  = element_text(size=rel(0.9)),
          axis.title.x = element_text(margin = unit(c(3, 0, 0, 0), "mm"),size=rel(1)),
          axis.title.y = element_text(margin = unit(c(0, 3, 0, 0), "mm"),size=rel(1), angle = 90),
          strip.text = element_text(size = rel(1.15), colour = "black", face = "bold"),
          plot.margin=unit(c(10,10,10,10),"pt"))
}
theme_set(theme_std())

# inflow and weather spreadsheet
inflowdata <- read.csv("Inflow_for_interpolation.csv")
head(inflowdata)
inflowdata <- inflowdata %>% unite(Date, Year, Month, Day, sep = '-')
inflowdata$Date <- as.Date(inflowdata$Date, format = "%Y-%m-%d")
inflowdata <- mutate(inflowdata, Year = as.numeric(format(inflowdata$Date, "%Y")))
mutate(Day = as.numeric(format(mod.match.daily$Date, "%d")))  %>% 
  

# lake spreadsheet
obsinit <- read.csv("../Postproc_code/L227/Observed_IntegratedEpi_icefree_lateMay.csv")
colnames(obsinit) <- c("org.date","Year","Month","Day", "obs.TP","obs.chla","obs.TDP", "obs.PP")
obs <- obsinit %>% 
  unite(Date, Year, Month, Day, sep = '-') #%>%
obs <- data.frame(obs, obsinit$Month, obsinit$Year)
obs <- na.omit(obs)
obs <- obs[,-1]
colnames(obs) <- c("Date", "obs.TP","obs.chla","obs.TDP", "obs.PP", "Month", "Year")
obs$Date <- as.Date(obs$Date, format = "%Y-%m-%d") #Converts data to date structure

mod.match.daily <- right_join(obs, inflowdata, by = "Date")
mod.match.daily$Month <- as.numeric(format(mod.match.daily$Date, "%m"))
mod.match.daily$Year <- format(mod.match.daily$Date, "%Y")
mod.match.daily <- mod.match.daily %>% 
  mutate(Day = as.numeric(format(mod.match.daily$Date, "%d")))  %>% 
  filter(Year != 1996) %>% #1996 was the pike year
  filter(Month > 4 & Month < 11)  %>% 
  filter(Month == 5 & Day > 15 | Month > 5 )

mod.match.daily$obs.PP[mod.match.daily$Month == 5 & mod.match.daily$Day == 16] <- 0
obs.PP.interpolated <- na.approx(mod.match.daily$obs.PP)
# remove dates before and after the first and last measured values, respectively
mod.match.daily <- mod.match.daily[-c(1:6, 6694:6718),]
mod.match.daily$obs.PP <- obs.PP.interpolated
mod.match.cumulative <- ddply(mod.match.daily, .(Year), transform,  Cum.Sum.obs.PP = cumsum(obs.PP))
mod.match.cumulative <- filter(mod.match.cumulative, Year > 1969)
# calculate end-of season comparison between cumulative modeled vs. observed PP
end.of.season.cumulative <- as.data.frame(mod.match.cumulative) %>%
  filter(Day == 31 & Month == 10) %>%
  select(Date, Month, Year, Day, Cum.Sum.obs.PP)

end.of.season.cumulative <- arrange(end.of.season.cumulative, Cum.Sum.obs.PP)
summary(end.of.season.cumulative$Cum.Sum.obs.PP)
boxplot(end.of.season.cumulative$Cum.Sum.obs.PP)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3216    3877    4460    5202    5698   10638 

# years above 3rd quartile: 
  # 1987, 1980, 1984, 1979, 1973, 1998, 1972, 1981, 1990, 1985

# years below 1st quartile:
  # 1971, 2008, 1976, 1975, 1977, 1992, 1983, 1970, 2003, 2001


#### Plots of variables over time ----

#Air temperature
ggplot(inflowdata, aes(x = Date, y = AirTemp, color = Year)) +
  geom_point(size = 0.5) + 
  scale_color_manual(color = c("black", "blue", "red"), breaks = c("1988", "1989", "1990"))

#Inflow temperature
ggplot(inflowdata, aes(x = Date, y = InflowTemp)) +
  geom_point(size = 0.5)

#Wind speed
ggplot(inflowdata, aes(x = Date, y = WindSpeed)) +
  geom_point(size = 0.5)

#TP
ggplot(inflowdata, aes(x = Date, y = TP)) +
  geom_point(size = 0.5)

#DOC
ggplot(inflowdata, aes(x = Date, y = DOC)) +
  geom_point(size = 0.5)

#NO3
ggplot(inflowdata, aes(x = Date, y = NO3)) +
  geom_point(size = 0.5)

#NH4
ggplot(inflowdata, aes(x = Date, y = NH4)) +
  geom_point(size = 0.5)