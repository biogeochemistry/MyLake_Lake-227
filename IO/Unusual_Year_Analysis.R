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
  
#### inflow and weather spreadsheet ----
inflowdata <- read.csv("Inflow_for_interpolation.csv")
head(inflowdata)
inflowdata <- inflowdata %>% unite(Date, Year, Month, Day, sep = '-')
inflowdata$Date <- as.Date(inflowdata$Date, format = "%Y-%m-%d")
inflowdata <- mutate(inflowdata, Year = as.numeric(format(inflowdata$Date, "%Y")))  %>% 
  mutate(inflowdata, Class = 0)

inflowdata$Class <- ifelse(inflowdata$Year == 1987 | inflowdata$Year == 1980 | 
                             inflowdata$Year == 1984 | inflowdata$Year == 1979 | 
                             inflowdata$Year == 1973 | inflowdata$Year == 1998 | 
                             inflowdata$Year == 1972 | inflowdata$Year == 1981 | 
                             inflowdata$Year == 1990 | inflowdata$Year == 1985, 
                           "High", 
                           ifelse(inflowdata$Year == 1971 | inflowdata$Year == 2008 | 
                                    inflowdata$Year == 1976| inflowdata$Year == 1975 | 
                                    inflowdata$Year == 1977 | inflowdata$Year == 1992 | 
                                    inflowdata$Year == 1983 | inflowdata$Year == 1970 | 
                                    inflowdata$Year == 2003 | inflowdata$Year == 2001, 
                                  "Low", "Typical"))  
# lake spreadsheet
obsinit <- read.csv("../Postproc_code/L227/Observed_IntegratedEpi_icefree_lateMay.csv")
colnames(obsinit) <- c("org.date","Year","Month","Day", "obs.TP","obs.chla","obs.TDP", "obs.PP", "obs.DOC")
obs <- obsinit %>% 
  unite(Date, Year, Month, Day, sep = '-') #%>%
obs <- data.frame(obs, obsinit$Month, obsinit$Year)
obs <- na.omit(obs)
obs <- obs[,-1]
colnames(obs) <- c("Date", "obs.TP","obs.chla","obs.TDP", "obs.PP", "obs.DOC", "Month", "Year")
obs$Date <- as.Date(obs$Date, format = "%Y-%m-%d") #Converts data to date structure

mod.match.daily <- right_join(obs, inflowdata, by = "Date")
mod.match.daily$Month <- as.numeric(format(mod.match.daily$Date, "%m"))
mod.match.daily$Year <- format(mod.match.daily$Date, "%Y")
mod.match.daily <- mod.match.daily %>% 
  mutate(Day = as.numeric(format(mod.match.daily$Date, "%d")))  %>% 
  filter(Year != 1996) %>% #1996 was the pike year
  filter(Month > 4 & Month < 11)  %>% 
  filter(Month == 5 & Day > 15 | Month > 5 )

# enable October of 2009 to be part of interpolation
mod.match.daily[6718, 5] <- 0

mod.match.daily$obs.PP[mod.match.daily$Month == 5 & mod.match.daily$Day == 16] <- 0
obs.PP.interpolated <- na.approx(mod.match.daily$obs.PP)
# remove dates before and after the first and last measured values, respectively
mod.match.daily <- mod.match.daily[-c(1:6),]
mod.match.daily$obs.PP <- obs.PP.interpolated
mod.match.cumulative <- ddply(mod.match.daily, .(Year), transform,  Cum.Sum.obs.PP = cumsum(obs.PP))
mod.match.cumulative <- filter(mod.match.cumulative, Year > 1969)
# calculate end-of season comparison between cumulative modeled vs. observed PP
end.of.season.cumulative <- as.data.frame(mod.match.cumulative) %>%
  filter(Day == 31 & Month == 10) %>%
  select(Date, Month, Year, Day, Cum.Sum.obs.PP)

end.of.season.cumulative <- arrange(end.of.season.cumulative, Cum.Sum.obs.PP)
summary(end.of.season.cumulative$Cum.Sum.obs.PP)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3182    3833    4434    5150    5545   10638  

# years above 3rd quartile: 
  # 1987, 1980, 1984, 1979, 1973, 1998, 1972, 1981, 1990, 1985

# years below 1st quartile:
  # 1971, 2008, 1976, 1975, 1977, 1992, 1983, 1970, 2003, 2001

end.of.season.cumulative <- mutate(end.of.season.cumulative, Class = 0)
end.of.season.cumulative$Class[1:10] <- "Low"
end.of.season.cumulative$Class[11:29] <- "Typical"
end.of.season.cumulative$Class[30:39] <- "High"

mod.match.cumulative <- mutate(mod.match.cumulative, Class = 0)
mod.match.cumulative$Class <- ifelse(mod.match.cumulative$Year == 1987 | mod.match.cumulative$Year == 1980 | 
                                       mod.match.cumulative$Year == 1984 | mod.match.cumulative$Year == 1979 | 
                                       mod.match.cumulative$Year == 1973 | mod.match.cumulative$Year == 1998 | 
                                       mod.match.cumulative$Year == 1972 | mod.match.cumulative$Year == 1981 | 
                                       mod.match.cumulative$Year == 1990 | mod.match.cumulative$Year == 1985, 
                           "High", 
                           ifelse(mod.match.cumulative$Year == 1971 | mod.match.cumulative$Year == 2008 | 
                                    mod.match.cumulative$Year == 1976| mod.match.cumulative$Year == 1975 | 
                                    mod.match.cumulative$Year == 1977 | mod.match.cumulative$Year == 1992 | 
                                    mod.match.cumulative$Year == 1983 | mod.match.cumulative$Year == 1970 | 
                                    mod.match.cumulative$Year == 2003 | mod.match.cumulative$Year == 2001, 
                                  "Low", "Typical"))  


  

#### Plots of variables over time ----

InflowVolumeplot <- 
ggplot(inflowdata, aes(x = Date, y = InflowVolume, color = Class)) +
  geom_point(size = 0.5) + 
  scale_color_manual(values = c("red", "blue", "black")) + 
  theme(legend.position = "none")

AirTempplot <-
ggplot(inflowdata, aes(x = Date, y = AirTemp, color = Class)) +
  geom_point(size = 0.5) + 
  scale_color_manual(values = c("red", "blue", "black")) + 
  theme(legend.position = "none")

InflowTempplot <- 
ggplot(inflowdata, aes(x = Date, y = InflowTemp, color = Class)) +
  geom_point(size = 0.5) + 
  scale_color_manual(values = c("red", "blue", "black")) + 
  theme(legend.position = "none")

WindSpeedplot <-
ggplot(inflowdata, aes(x = Date, y = WindSpeed, color = Class)) +
  geom_point(size = 0.5) + 
  scale_color_manual(values = c("red", "blue", "black")) + 
  theme(legend.position = "none")

TPplot <-
ggplot(inflowdata, aes(x = Date, y = TP, color = Class)) +
  geom_point(size = 0.5) + 
  scale_color_manual(values = c("red", "blue", "black")) + 
  theme(legend.position = "none")

DOCplot <-
ggplot(inflowdata, aes(x = Date, y = DOC, color = Class)) +
  geom_point(size = 0.5) + 
  scale_color_manual(values = c("red", "blue", "black")) + 
  theme(legend.position = "none")

#NO3
ggplot(inflowdata, aes(x = Date, y = NO3, color = Class)) +
  geom_point(size = 0.5) + 
  scale_color_manual(values = c("red", "blue", "black")) + 
  theme(legend.position = "none")

#NH4
ggplot(inflowdata, aes(x = Date, y = NH4, color = Class)) +
  geom_point(size = 0.5) + 
  scale_color_manual(values = c("red", "blue", "black")) + 
  theme(legend.position = "none")

LakePPplot <-
ggplot(mod.match.cumulative, aes(x = Date, y = Cum.Sum.obs.PP, color = Class)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("red", "blue", "black")) + 
  ylab(expression(Cumulative ~ PP ~ (mu*g / L))) + 
  theme(legend.position = "none")

grid.arrange(LakePPplot, AirTempplot, InflowTempplot, InflowVolumeplot, WindSpeedplot, DOCplot, ncol = 3)
