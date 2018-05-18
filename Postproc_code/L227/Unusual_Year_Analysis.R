#### Lake 227 observation analysis ----
# This script examines patterns in observations for ELA Lake 227 data and associated weather and inflow data.
# The purpose is to identify years with highest and lowest cumulative PP observations 
# and to examine potential connections with unusual weather and inflow patterns.

#### Setup ----

library(tidyverse)
library(zoo)
library(dplyr)
library(plyr)
library(lubridate)
library(chron)
library(trend)
library(gridExtra)
library(pander)

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

#### lake spreadsheet ----
obsinit <- read.csv("Observed_IntegratedEpi_icefree_lateMay.csv")
colnames(obsinit) <- c("org.date","Year","Month","Day", "obs.TP","obs.chla","obs.TDP", "obs.PP", "obs.DOC")
obs <- obsinit %>% 
  unite(Date, Year, Month, Day, sep = '-') #%>%
obs <- data.frame(obs, obsinit$Month, obsinit$Year)
#obs <- na.omit(obs)
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

# enable October of 2016 to be part of interpolation
mod.match.daily[7922, 5] <- 0

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
# 3190    3877    4404    5002    5095   10638 

# years above 3rd quartile: 
  # 1987, 1980, 1984, 1979, 1973, 1998, 1972, 1981, 1990, 1985, 1982, 1997

# years below 1st quartile:
  # 2009, 2008, 1971, 1976, 1975, 2016, 1977, 1992, 1983, 2003, 2015, 2001

end.of.season.cumulative <- mutate(end.of.season.cumulative, Class = 0)
# end.of.season.cumulative$Class[1:10] <- "Low"
# end.of.season.cumulative$Class[11:29] <- "Typical"
# end.of.season.cumulative$Class[30:39] <- "High"

mod.match.cumulative <- mutate(mod.match.cumulative, Class = 0)
mod.match.cumulative$Class <- ifelse(mod.match.cumulative$Year == 1987 | mod.match.cumulative$Year == 1980 | 
                                       mod.match.cumulative$Year == 1984 | mod.match.cumulative$Year == 1979 | 
                                       mod.match.cumulative$Year == 1973 | mod.match.cumulative$Year == 1998 | 
                                       mod.match.cumulative$Year == 1972 | mod.match.cumulative$Year == 1981 | 
                                       mod.match.cumulative$Year == 1990 | mod.match.cumulative$Year == 1985 |
                                       mod.match.cumulative$Year == 1982 | mod.match.cumulative$Year == 1997, 
                           "High", 
                           ifelse(mod.match.cumulative$Year == 2009 | mod.match.cumulative$Year == 2008 |
                                    mod.match.cumulative$Year == 1971 |  mod.match.cumulative$Year == 1976| 
                                    mod.match.cumulative$Year == 1975 | mod.match.cumulative$Year == 2016 |
                                    mod.match.cumulative$Year == 1977 | mod.match.cumulative$Year == 1992 | 
                                    mod.match.cumulative$Year == 1983 | mod.match.cumulative$Year == 2003 |
                                    mod.match.cumulative$Year == 2015 | mod.match.cumulative$Year == 2001, 
                                  "Low", "Typical"))  

#### inflow and weather spreadsheet ----
inflowdata <- read.csv("../../IO/Inflow_for_interpolation.csv")
head(inflowdata)
inflowdata <- inflowdata %>% unite(Date, Year, Month, Day, sep = '-')
inflowdata$Date <- as.Date(inflowdata$Date, format = "%Y-%m-%d")
inflowdata <- mutate(inflowdata, Year = as.numeric(format(inflowdata$Date, "%Y")))  %>% 
  mutate(inflowdata, Class = 0)

inflowdata$Class <- ifelse(inflowdata$Year == 1987 | inflowdata$Year == 1980 | 
                             inflowdata$Year == 1984 | inflowdata$Year == 1979 | 
                             inflowdata$Year == 1973 | inflowdata$Year == 1998 | 
                             inflowdata$Year == 1972 | inflowdata$Year == 1981 | 
                             inflowdata$Year == 1990 | inflowdata$Year == 1985|
                             inflowdata$Year == 1982 | inflowdata$Year == 1997, 
                           "High", 
                           ifelse(inflowdata$Year == 2009 | inflowdata$Year == 2008 |
                                    inflowdata$Year == 1971 |  inflowdata$Year == 1976| 
                                    inflowdata$Year == 1975 | inflowdata$Year == 2016 |
                                    inflowdata$Year == 1977 | inflowdata$Year == 1992 | 
                                    inflowdata$Year == 1983 | inflowdata$Year == 2003 |
                                    inflowdata$Year == 2015 | inflowdata$Year == 2001, 
                                  "Low", "Typical"))   

#### Phytoplankton diversity data ----  
ice.dates <- read.csv("./PhytoDiversity/IISDELA_ice_1969-2017.csv")

ice.dates <- ice.dates %>%
  mutate(DOY.IceOff = yday(as.Date(ice.dates$Ice.Off.Date, format = "%Y-%m-%d")),
         DOY.IceOn = yday(as.Date(ice.dates$Ice.On.Date, format = "%Y-%m-%d")),
         year = as.numeric(format(as.Date(ice.dates$Ice.Off.Date, format = "%Y-%m-%d"),"%Y"))) %>%
  dplyr::select(Ice.Off.Date, Ice.On.Date, year, DOY.IceOff, DOY.IceOn) %>%
  mutate(NoIceFree = DOY.IceOn - DOY.IceOff) 

fert.dates <- read.csv("./PhytoDiversity/L227_fertilizerAdditions.csv")

fert.dates <- fert.dates %>%
  mutate(date = as.Date(fert.dates$start_date, format = "%d-%b-%y"),
         year = as.numeric(format(date, "%Y")),
         DOY.fert = yday(date)) %>%
  group_by(year) %>%
  filter(DOY.fert == min(DOY.fert)) %>%
  dplyr::select(date, year, DOY.fert) %>%
  arrange(year) %>%
  distinct(date, year, DOY.fert)
fert.dates <- as.data.frame(fert.dates)

## Updated fertilzation data from 2013 - 2017
fert.dates.add <- as.matrix(read.csv("./PhytoDiversity/L227_fertilizerAdditions_2013-2017.csv", header = T))
colnames(fert.dates.add) <- fert.dates.add[1,]
fert.dates.add <- fert.dates.add[-1,]
fert.dates.add <- as.data.frame(fert.dates.add)

fert.dates.add <- fert.dates.add %>%
  mutate(date = as.Date(fert.dates.add$start_date, format = "%Y-%m-%d"),
         year = as.numeric(format(date, "%Y")),
         DOY.fert = yday(date)) %>%
  group_by(year) %>%
  filter(DOY.fert == min(DOY.fert)) %>%
  dplyr::select(date, year, DOY.fert) %>%
  arrange(year) %>%
  distinct(date, year, DOY.fert)
fert.dates.add <- as.data.frame(fert.dates.add)

# merge the data 
fert.dates <- rbind(fert.dates, fert.dates.add)

# Read in data from 1968-2012
## Import data
dat.tax1 <- read.csv("./PhytoDiversity/L227-phytocounts_19682012.csv", header = T) #1968 - 2012
#head(dat.tax1)

## change column names
colnames(dat.tax1) <- c("year","data_set","sample.id","sequence.number","location","sublocation","station","station.type","start.date","start.time",
                        "stratum","stratum.volume","sample.type","start.depth","end.depth","taxon.code","old.cell.volume","cell.volume","correction.factor","cell.count",
                        "cell.count2","cell.length","cell.width","field.ref","cell.density","biomass","collection.authority","collection.method","analyst")

## set date structure
dat.tax1$start.date <- as.Date(dat.tax1$start.date, format = "%Y-%m-%d")

## Remove all 1996 data - these are incomplete and will be amended in the following section
dat.tax1.minus96 <- dat.tax1[dat.tax1$year != 1996,]

# Read in and clean data from 2013-2015
## Import data
dat.tax3 <- read.csv("./PhytoDiversity/L227-phytocounts_20132015.csv", header = F)
#head(dat.tax3)

## Subset data
dat.tax3 <- dat.tax3[3:1745,1:(ncol(dat.tax3)-1)]

## change column names
colnames(dat.tax3) <- c("location","start.date","stratum","start.depth","end.depth","taxon.code","species.name","subsample.volume","cell.count",
                        "correction.factor","cell.length","cell.width","cell.volume","cell.density","biomass")

dat.tax3$stratum.volume <- c(NA)
dat.tax3$start.date <- as.Date(dat.tax3$start.date, format = "%Y-%m-%d") #Converts data to date structure
dat.tax3$taxon.code <- paste0("P",dat.tax3$taxon.code) # Add P tp species code

subset.cols <- c("start.date","stratum","stratum.volume","start.depth","end.depth","taxon.code","correction.factor","cell.count","cell.density","cell.volume","biomass")

## select columns from each data frame
dat.tax1.subset <- dplyr::select(dat.tax1.minus96, subset.cols)
dat.tax3.subset <- dplyr::select(dat.tax3, subset.cols)

phyto.species <- rbind(dat.tax1.subset,dat.tax3.subset)

#head(phyto.species)
#str(phyto.species)

# Convert specific columns to numeric class
phyto.species[,c(3:5,7:11)]<- lapply(phyto.species[,c(3:5,7:11)], as.numeric)

# fix the stratum layers
phyto.species$stratum <- as.factor(gsub("META","MET",phyto.species$stratum)) 

# add columns for day and year; recalculate biomass and density
phyto.species <- phyto.species %>%
  mutate(year = as.numeric(format(start.date,'%Y')),
         DOY = yday(start.date)) %>%
  left_join(ice.dates, by = "year") %>%
  left_join(fert.dates, by = "year") %>%
  mutate(DOY.IceAdj = DOY-DOY.IceOff,
         DOY.FertAdj = DOY-DOY.fert,
         cell.density2 = correction.factor * cell.count,
         biomass2 = cell.volume*correction.factor*cell.count*0.000001)

phyto.species <- droplevels(phyto.species)

sampling.dates <- phyto.species %>% 
  filter(stratum == "EPI") %>%
  arrange(start.date) %>%
  distinct(start.date) %>%
  mutate(yr = format(start.date,"%Y"))

# create sampling sequence for each year
yr.list <- c()
samp.no.list <- c()

for(i in unique(sampling.dates$yr)){
  yr <- rep(i, length(sampling.dates$start.date[sampling.dates$yr == i]))
  samp.no <- seq(1,length(sampling.dates$start.date[sampling.dates$yr == i]),1)
  
  yr.list <- append(yr, yr.list)
  samp.no.list <- append(samp.no, samp.no.list)
}

# Create a new list of sampling numbers
samp.list <- data.frame(yr.list, samp.no.list) %>%
  mutate(samp.no = paste0(yr.list,"_",samp.no.list)) %>%
  arrange(yr.list,samp.no.list) 

beta.list <- samp.list
samp.list <- data.frame(samp.list, sampling.dates$start.date)[,3:4]

# use list to create species matrix
epi.species.tab <- phyto.species %>% 
  filter(stratum == "EPI") %>%    # Select EPI samples
  arrange(start.date) %>%      # order dates 
  left_join(samp.list, by = c("start.date"="sampling.dates.start.date")) %>% 
  dplyr::select(samp.no, taxon.code, cell.density2) # select appropriate columns; diversity uses count rather than density or biomass 

epi.species.tab$samp.no <- as.factor(epi.species.tab$samp.no)

# Convert long form data into matrix by sampling number
library(reshape2)
epi.species.mat <- dcast(data = epi.species.tab,  samp.no ~ taxon.code, value.var = "cell.density2", fun.aggregate = sum)

row.names(epi.species.mat) <- epi.species.mat[,1]
epi.species.mat <- epi.species.mat[,-1]

library(vegan)

epi.species.mat[is.na(epi.species.mat)] <- 0

shan <- as.data.frame(diversity(epi.species.mat, index = "shannon"))
colnames(shan) <- c("shannon")
shan$site <- rownames(shan)
shan$year <- substr(shan$site, 1, 4)
shan$year <- as.factor(shan$year)

shan$sampling <- substring(shan$site, 6)

MeanShannon <- aggregate(shan, list(Year = shan$year), mean)
MeanShannon <- select(MeanShannon, Year:shannon)
MeanShannon <- mutate(MeanShannon, Date = as.Date(MeanShannon$Year, "%Y"), Class = 0)

MeanShannon$Year <- as.numeric(format(MeanShannon$Year, "%Y"))
inflowdata <- mutate(inflowdata, Year = as.numeric(format(inflowdata$Date, "%Y")))  %>% 
  mutate(inflowdata, Class = 0)
MeanShannon$Year <- format(MeanShannon$Year , "%Y"))

MeanShannon$Class <- ifelse(MeanShannon$Year == 1987 | MeanShannon$Year == 1980 | 
                              MeanShannon$Year == 1984 | MeanShannon$Year == 1979 | 
                              MeanShannon$Year == 1973 | MeanShannon$Year == 1998 | 
                              MeanShannon$Year == 1972 | MeanShannon$Year == 1981 | 
                              MeanShannon$Year == 1990 | MeanShannon$Year == 1985 |
                              MeanShannon$Year == 1982 | MeanShannon$Year == 1997, 
                                     "High", 
                            ifelse(MeanShannon$Year == 2009 | MeanShannon$Year == 2008 |
                                     MeanShannon$Year == 1971 |  MeanShannon$Year == 1976| 
                                     MeanShannon$Year == 1975 | 
                                     MeanShannon$Year == 1977 | MeanShannon$Year == 1992 | 
                                     MeanShannon$Year == 1983 | MeanShannon$Year == 2003 |
                                     MeanShannon$Year == 2015 | MeanShannon$Year == 2001, 
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

PhytoDiversityplot <-
  ggplot(MeanShannon, aes(x = Date, y = shannon, color = Class)) +
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

grid.arrange(LakePPplot, AirTempplot, WindSpeedplot, InflowVolumeplot, DOCplot, PhytoDiversityplot, ncol = 3)
