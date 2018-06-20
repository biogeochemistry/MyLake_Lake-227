#### Information ----

#### Setup ----
library(tidyverse)
library(knitr)
knitr::opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE)
library(magrittr)
library(plyr)
library(hydroGOF)
library(gridExtra)
library(zoo)
library(colormap)
library(trend)
library(png)
library(reshape2)
library(vegan)

theme_std <- function (base_size = 18, base_family = "") {
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
          plot.margin=unit(c(10,10,10,10),"pt"))}
theme_set(theme_std())

scales::show_col(colormap(colormap = colormaps$plasma, nshades=15))
scales::show_col(colormap(colormap = colormaps$inferno, nshades=15))
scales::show_col(colormap(colormap = colormaps$viridis, nshades=15))

#### Historical Data Analysis ----
# Figures 1 and S1 for manuscript
# Code taken from 
#### Spreadsheet Import ====
NtoPLake <- read.csv("./Stoichiometry/NP_Stoichiometry_L227.csv")
attach(NtoPLake)
NtoPInflow <- read.csv("./Stoichiometry/NP_Stoichiometry_Inflow.csv")
attach(NtoPInflow)
CyanoPercent <- read.csv("./Phytoplankton/L227_cyanopercent.csv")
InflowData <- read.csv("../../IO/Inflow_for_interpolation.csv")
InflowData <- InflowData %>% unite(Date, Year, Month, Day, sep = '-')
dat.tax1 <- read.csv("./PhytoDiversity/L227-phytocounts_19682012.csv", header = T) #1968 - 2012
dat.tax3 <- read.csv("./PhytoDiversity/L227-phytocounts_20132015.csv", header = F)
Radiation <- read.csv("../../IO/TotalRadiation.csv")

# Convert dates from factor to date format, add month or year
Datelake <- as.Date(NtoPLake$Date, "%m/%d/%y")
Dateinflow <- as.Date(NtoPInflow$Date, "%m/%d/%y")
Monthlake <- as.numeric(format(Datelake, "%m"))
Monthinflow <- as.numeric(format(Dateinflow, "%m"))
CyanoPercent$Date <- as.Date(CyanoPercent$Date, "%m/%d/%y")
InflowData$Date <- as.Date(InflowData$Date, format = "%Y-%m-%d")
InflowData <- mutate(InflowData, Year = as.numeric(format(InflowData$Date, "%Y")))  %>% 
  mutate(InflowData, Class = 0)
Radiation$Date <- as.Date(Radiation$Date, "%Y-%m-%d")

#### Stoichiometry ==== 
#### Data wrangling ####
# Change concentrations from mass (ug/L or mg/d) to molar (umol/L or mmol/d)
Fert_TP_molar <- Fert_TP_mg.d/30.97
Inflow_TP_molar <- Inflow_TP_mg.d/30.97
Fert_TN_molar <- Fert_TN_mg.d/14.01
Inflow_TN_molar <- Inflow_TN_mg.d/14.01
NO3_molar <- NO3/14.01
NH4_molar <- NH4/14.01
DIN_molar <- NO3_molar + NH4_molar
PN_molar <- PN/14.01
TDN_molar <- TDN/14.01
PP_molar <- PP/30.97
TDP_molar <- TDP/30.97
TN_molar <- TDN_molar + PN_molar
TP_molar <- TDP_molar + PP_molar
chl <- NtoPLake$Chl

# Combine inflows + fertilization
Input_TN_molar <- Fert_TN_molar + Inflow_TN_molar
Input_TP_molar <- Fert_TP_molar + Inflow_TP_molar

# Create N:P stoichiometric ratios in lake and inputs
Fert_NtoP <- Fert_TN_molar/Fert_TP_molar
Inflow_NtoP <- Inflow_TN_molar/Inflow_TP_molar
Input_NtoP <- Input_TN_molar/Input_TP_molar
DINtoTDP <- DIN_molar/TDP_molar
TDNtoTDP <- TDN_molar/TDP_molar
PNtoPP <- PN_molar/PP_molar
TNtoTP <- TN_molar/TP_molar

# Create data frames for N:P stoichiometry in lake and inflows for May-October
NtoPinsitu <- data.frame(Datelake, Monthlake, TNtoTP, PNtoPP, DINtoTDP, TDNtoTDP, DIN_molar, PP, chl)
NtoPinsitu <- filter(NtoPinsitu, Monthlake > 4 & Monthlake < 11)
NtoPinput <- data.frame(Dateinflow, Monthinflow, Fert_NtoP, Inflow_NtoP, Input_NtoP, Input_TN_molar)
NtoPinput <- filter(NtoPinput, Monthinflow > 4 & Monthinflow < 11)

# Add Year to datasets
NtoPinsitu <- mutate(NtoPinsitu, Year = format(NtoPinsitu$Datelake, "%Y"))
NtoPinput <- mutate(NtoPinput, Year = format(NtoPinput$Dateinflow, "%Y"))

# Test for normality
shapiro.test(NtoPinsitu$TNtoTP)
shapiro.test(NtoPinsitu$TDNtoTDP)
shapiro.test(NtoPinsitu$PP)
shapiro.test(NtoPinsitu$chl)

# Separate datasets and eliminate NAs
TNtoTPdataset <- NtoPinsitu %>%
  select(Datelake, TNtoTP) %>%
  drop_na()
rownames(TNtoTPdataset) <- NULL

TDNtoTDPdataset <- NtoPinsitu %>%
  select(Datelake, TDNtoTDP) %>%
  drop_na()
rownames(TDNtoTDPdataset) <- NULL

PPdataset <- NtoPinsitu %>%
  select(Datelake, PP) %>%
  drop_na()
rownames(PPdataset) <- NULL

chldataset <- NtoPinsitu %>%
  select(Datelake, chl) %>%
  drop_na()
rownames(chldataset) <- NULL

# Test differences among fertilization periods
var.test(NtoPinsitu$PP[NtoPinsitu$Year < 1975], NtoPinsitu$PP[NtoPinsitu$Year > 1974 & NtoPinsitu$Year < 1990])
var.test(NtoPinsitu$PP[NtoPinsitu$Year < 1975], NtoPinsitu$PP[NtoPinsitu$Year > 1989])
var.test(NtoPinsitu$PP[NtoPinsitu$Year > 1974 & NtoPinsitu$Year < 1990], NtoPinsitu$PP[NtoPinsitu$Year > 1989])

# Variance of PP   Period 3 < Period 1 < Period 2

wilcox.test(NtoPinsitu$PP[NtoPinsitu$Year < 1975], NtoPinsitu$PP[NtoPinsitu$Year > 1974 & NtoPinsitu$Year < 1990])
wilcox.test(NtoPinsitu$PP[NtoPinsitu$Year < 1975], NtoPinsitu$PP[NtoPinsitu$Year > 1989])
wilcox.test(NtoPinsitu$PP[NtoPinsitu$Year > 1974 & NtoPinsitu$Year < 1990], NtoPinsitu$PP[NtoPinsitu$Year > 1989])

# PP   Period 1 < Period 3 < Period 2

var.test(NtoPinsitu$chl[NtoPinsitu$Year < 1975], NtoPinsitu$chl[NtoPinsitu$Year > 1974 & NtoPinsitu$Year < 1990])
var.test(NtoPinsitu$chl[NtoPinsitu$Year < 1975], NtoPinsitu$chl[NtoPinsitu$Year > 1989])
var.test(NtoPinsitu$chl[NtoPinsitu$Year > 1974 & NtoPinsitu$Year < 1990], NtoPinsitu$chl[NtoPinsitu$Year > 1989])

# Variance of chl   Period 1 > Period 2 and 3

wilcox.test(NtoPinsitu$chl[NtoPinsitu$Year < 1975], NtoPinsitu$chl[NtoPinsitu$Year > 1974 & NtoPinsitu$Year < 1990])
wilcox.test(NtoPinsitu$chl[NtoPinsitu$Year < 1975], NtoPinsitu$chl[NtoPinsitu$Year > 1989])
wilcox.test(NtoPinsitu$chl[NtoPinsitu$Year > 1974 & NtoPinsitu$Year < 1990], NtoPinsitu$chl[NtoPinsitu$Year > 1989])

# chl no significant differences

#### Pettitt's test ####
# Pettitt's test: detects a single changepoint in hydrological/climate series with continuous data
TNtoTPvec <- as.vector(TNtoTPdataset$TNtoTP)
pettitt.test(TNtoTPvec)
    # data:  TNtoTPvec
    # U* = 30578, p-value = 3.796e-06
    # alternative hypothesis: two.sided
    # sample estimates:
    #   probable change point at time K 
    # 470 

TDNtoTDPvec <- as.vector(TDNtoTDPdataset$TDNtoTDP)
pettitt.test(TDNtoTDPvec)
    # data:  TDNtoTDPvec
    # U* = 27951, p-value = 0.0001079
    # alternative hypothesis: two.sided
    # sample estimates:
    #   probable change point at time K 
    # 120 

PPvec <- as.vector(PPdataset$PP)
pettitt.test(PPvec)
    # data:  PPvec
    # U* = 31464, p-value = 1.245e-05
    # alternative hypothesis: two.sided
    # sample estimates:
    #   probable change point at time K 
    # 179 

chlvec <- as.vector(chldataset$chl)
pettitt.test(chlvec)
    # data:  chlvec
    # U* = 33080, p-value = 0.0003617
    # alternative hypothesis: two.sided
    # sample estimates:
    #   probable change point at time K 
    # 435 

#### Mann-Kendall test ####
# Mann-Kendall test: detect monotonic trends in series of environmental/climate/hydrological data
mk.test(TNtoTPvec[1:469])
# data:  TNtoTPvec[1:469]
# z = 4.0215, n = 469, p-value = 5.782e-05
# alternative hypothesis: true S is not equal to 0
# sample estimates:
#   S         varS          tau 
# 1.363800e+04 1.149894e+07 1.242727e-01 

mk.test(TNtoTPvec[470:752])
# data:  TNtoTPvec[470:752]
# z = 0.77053, n = 283, p-value = 0.441
# alternative hypothesis: true S is not equal to 0
# sample estimates:
#   S         varS          tau 
# 1.227000e+03 2.531613e+06 3.075342e-02 

mk.test(TDNtoTDPvec[1:119])
# data:  TDNtoTDPvec[1:119]
# z = 1.4424, n = 119, p-value = 0.1492
# alternative hypothesis: true S is not equal to 0
# sample estimates:
#   S         varS          tau 
# 6.290000e+02 1.895640e+05 8.960752e-02 

mk.test(TDNtoTDPvec[120:781])
# data:  TDNtoTDPvec[120:781]
# z = -3.5267, n = 662, p-value = 0.0004207
# alternative hypothesis: true S is not equal to 0
# sample estimates:
#   S          varS           tau 
# -2.004700e+04  3.230782e+07 -9.167299e-02 

mk.test(PPvec[1:178])
# data:  PPvec[1:178]
# z = 1.3246, n = 178, p-value = 0.1853
# alternative hypothesis: true S is not equal to 0
# sample estimates:
#   S         varS          tau 
# 1.053000e+03 6.307977e+05 6.781124e-02 

mk.test(PPvec[179:791])
# Mann-Kendall trend test
# 
# data:  PPvec[179:791]
# z = -6.7774, n = 613, p-value = 1.223e-11
# alternative hypothesis: true S is not equal to 0
# sample estimates:
#   S          varS           tau 
# -3.431400e+04  2.563227e+07 -1.852262e-01 

mk.test(chlvec[1:434])
    # Mann-Kendall trend test
    # 
    # data:  chlvec[1:434]
    # z = -5.0183, n = 434, p-value = 5.212e-07
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S          varS           tau 
    # -1.515100e+04  9.113908e+06 -1.614472e-01 

mk.test(chlvec[435:913])
    # Mann-Kendall trend test
    # 
    # data:  chlvec[435:913]
    # z = -3.7101, n = 479, p-value = 0.0002072
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S          varS           tau 
    # -1.298600e+04  1.224942e+07 -1.134540e-01 

#### Plots ####
TNtoTPearly <- TNtoTPdataset[1:469,]
TNtoTPlate <- TNtoTPdataset[470:752,]
TNtoTPplot <- 
  ggplot(TNtoTPearly) +
  geom_rect(xmin = -Inf, xmax = as.numeric(as.Date("1975-01-01")), ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_rect(xmin = as.numeric(as.Date("1990-01-01")), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_point(data = TNtoTPearly, aes(x = Datelake, y = TNtoTP), size = 0.5, color = "#f99d15ff") + #significant positive slope
  geom_point(data = TNtoTPlate, aes(x = Datelake, y = TNtoTP), size = 0.5, color = "gray40") + #non-significant slope
  ylim(0,200) +
  ylab(expression(TN:TP)) +
  xlab(" ") +
  #theme(axis.text.x=element_blank()) + 
  geom_vline(xintercept = as.numeric(TNtoTPdataset$Datelake[470]), lty = 5)
print(TNtoTPplot)

TDNtoTDPearly <- TDNtoTDPdataset[1:119,]
TDNtoTDPlate <- TDNtoTDPdataset[120:781,]
TDNtoTDPplot <- 
  ggplot(TDNtoTDPearly) +
  geom_rect(xmin = -Inf, xmax = as.numeric(as.Date("1975-01-01")), ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_rect(xmin = as.numeric(as.Date("1990-01-01")), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_point(data = TDNtoTDPearly, aes(x = Datelake, y = TDNtoTDP), size = 0.5, color = "gray40") + #non-significant slope
  geom_point(data = TDNtoTDPlate, aes(x = Datelake, y = TDNtoTDP), size = 0.5, color = "#240c4cff") + #significant negative slope
  ylim(0,300) +
  ylab(expression(TDN:TDP)) +
  xlab(" ") +
  #theme(axis.text.x=element_blank()) + 
  geom_vline(xintercept = as.numeric(TDNtoTDPdataset$Datelake[120]), lty = 5)
print(TDNtoTDPplot)

PPearly <- PPdataset[1:178,]
PPlate <- PPdataset[179:791,]
PPplot <- 
  ggplot(PPearly) +
  geom_rect(xmin = -Inf, xmax = as.numeric(as.Date("1975-01-01")), ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_rect(xmin = as.numeric(as.Date("1990-01-01")), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_point(data = PPearly, aes(x = Datelake, y = PP), size = 0.5, color = "gray40") + #non-significant slope
  geom_point(data = PPlate, aes(x = Datelake, y = PP), size = 0.5, color = "#240c4cff") + #significant negative slope
  ylab(expression(PP)) +
  xlab(" ") +
  theme(axis.text.x=element_blank()) + 
  geom_vline(xintercept = as.numeric(PPdataset$Datelake[179]), lty = 5)
print(PPplot)

chlearly <- chldataset[1:434,]
chllate <- chldataset[435:913,]
chlplot <- 
  ggplot(chlearly) +
  geom_rect(xmin = -Inf, xmax = as.numeric(as.Date("1975-01-01")), ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_rect(xmin = as.numeric(as.Date("1990-01-01")), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_point(data = chlearly, aes(x = Datelake, y = chl), size = 0.5, color = "#240c4cff") + #significant negative slope
  geom_point(data = chllate, aes(x = Datelake, y = chl), size = 0.5, color = "#240c4cff") + #significant negative slope
  ylab(expression(chl)) +
  xlab(" ") +
  theme() + 
  geom_vline(xintercept = as.numeric(chldataset$Datelake[435]), lty = 5)
print(chlplot)

FertNtoPplot <-
  ggplot(NtoPinput, aes(x = Dateinflow, y = Fert_NtoP)) +
  geom_rect(xmin = -Inf, xmax = as.numeric(as.Date("1975-01-01")), ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_rect(xmin = as.numeric(as.Date("1990-01-01")), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_point(size = 0.5) +
  #ylim(0,40) +
  ylab(expression(FertN:FertP)) +
  xlab(" ") +
  annotate("text", x = as.Date("1968-01-01"), y = 38, hjust = 0, label = "High N:P", fontface = "bold") +
  annotate("text", x = as.Date("1980-01-01"), y = 38, hjust = 0, label = "Low N:P", fontface = "bold") +
  annotate("text", x = as.Date("2003-01-01"), y = 38, hjust = 0, label = "P Only", fontface = "bold") +
  annotate("point", x = as.Date("1980-01-01"), y = 30, pch = 14, size = 2) +
  annotate("point", x = as.Date("1996-01-01"), y = 30, pch = 13, size = 2) +
  theme(axis.text.x=element_blank()) 
print(FertNtoPplot)

#### test for equal sample number by year ####
tapply(NtoPinsitu$TNtoTP, NtoPinsitu$Year, length)
# range from 6 to 44
tapply(NtoPinsitu$TDNtoTDP, NtoPinsitu$Year, length)
# range from 6 to 44
tapply(NtoPinsitu$PP, NtoPinsitu$Year, length)
# range from 6 to 44
tapply(NtoPinsitu$chl, NtoPinsitu$Year, length)
# range from 6 to 44

#### Cyano % of biomass ====
# test for normality
shapiro.test(CyanoPercent$Nfixcyano.percent)

pettitt.test(CyanoPercent$Nfixcyano.percent)
    # data:  CyanoPercent$Nfixcyano.percent
    # U* = 40536, p-value < 2.2e-16
    # alternative hypothesis: two.sided
    # sample estimates:
    #   probable change point at time K 
    # 87 

mk.test(CyanoPercent$Nfixcyano.percent[1:86])
    # data:  CyanoPercent$Nfixcyano.percent[1:86]
    # z = 4.1489, n = 86, p-value = 3.34e-05
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S         varS          tau 
    # 4.910000e+02 1.394833e+04 3.650356e-01 

mk.test(CyanoPercent$Nfixcyano.percent[87:611])
    # data:  CyanoPercent$Nfixcyano.percent[87:611]
    # z = -1.4936, n = 525, p-value = 0.1353
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S          varS           tau 
    # -5.998000e+03  1.612026e+07 -4.370521e-02 

cyanoplot <- 
  ggplot(CyanoPercent) +
  geom_rect(xmin = -Inf, xmax = as.numeric(as.Date("1975-01-01")), ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_rect(xmin = as.numeric(as.Date("1990-01-01")), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_point(data = CyanoPercent, aes(x = Date, y = Nfixcyano.percent), size = 1) +
  ylab(expression("Prop. Diazotrophs")) + 
  xlab("") #+ 
 # theme(axis.text.x=element_blank())
print(cyanoplot)

#ggsave("PropCyano.jpg", cyanoplot, dpi = 300)

#### Phytoplankton ====
#### Data wrangling ####
# Read in data from 1968-2012
colnames(dat.tax1) <- c("year","data_set","sample.id","sequence.number","location","sublocation","station","station.type","start.date","start.time",
                        "stratum","stratum.volume","sample.type","start.depth","end.depth","taxon.code","old.cell.volume","cell.volume","correction.factor","cell.count",
                        "cell.count2","cell.length","cell.width","field.ref","cell.density","biomass","collection.authority","collection.method","analyst")
dat.tax1$start.date <- as.Date(dat.tax1$start.date, format = "%Y-%m-%d")
dat.tax1.minus96 <- dat.tax1[dat.tax1$year != 1996,] #1996 was the pike year

# Read in data from 2013-2015
dat.tax3 <- dat.tax3[3:1745,1:(ncol(dat.tax3)-1)]
colnames(dat.tax3) <- c("location","start.date","stratum","start.depth","end.depth","taxon.code","species.name","subsample.volume","cell.count",
                        "correction.factor","cell.length","cell.width","cell.volume","cell.density","biomass")
dat.tax3$stratum.volume <- c(NA)
dat.tax3$start.date <- as.Date(dat.tax3$start.date, format = "%Y-%m-%d") #Converts data to date structure
dat.tax3$taxon.code <- paste0("P",dat.tax3$taxon.code) # Add P tp species code
subset.cols <- c("start.date","stratum","stratum.volume","start.depth","end.depth","taxon.code","correction.factor","cell.count","cell.density","cell.volume","biomass")

# select columns from each data frame
dat.tax1.subset <- dplyr::select(dat.tax1.minus96, subset.cols)
dat.tax3.subset <- dplyr::select(dat.tax3, subset.cols)
phyto.species <- rbind(dat.tax1.subset,dat.tax3.subset)

# data wrangling
phyto.species[,c(3:5,7:11)]<- lapply(phyto.species[,c(3:5,7:11)], as.numeric)
phyto.species$stratum <- as.factor(gsub("META","MET",phyto.species$stratum)) 
phyto.species <- phyto.species %>%
  mutate(year = as.numeric(format(start.date,'%Y'))) 
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
  dplyr::select(samp.no, taxon.code, cell.count) # select appropriate columns; diversity uses count rather than density or biomass 

epi.species.tab$samp.no <- as.factor(epi.species.tab$samp.no)

# Convert long form data into matrix by sampling number
epi.species.mat <- dcast(data = epi.species.tab,  samp.no ~ taxon.code, value.var = "cell.count", fun.aggregate = sum)
row.names(epi.species.mat) <- epi.species.mat[,1]
epi.species.mat <- epi.species.mat[,-1]
epi.species.mat[is.na(epi.species.mat)] <- 0

shan <- as.data.frame(diversity(epi.species.mat, index = "shannon"))
colnames(shan) <- c("shannon")
shan$site <- rownames(shan)
shan$year <- substr(shan$site, 1, 4)
shan$year <- as.factor(shan$year)
shan$sampling <- substring(shan$site, 6)

MeanShannon <- aggregate(shan, list(Year = shan$year), mean)
MeanShannon <- select(MeanShannon, Year:shannon)
SDShannon <- aggregate(shan, list(Year = shan$year), sd)
SDShannon <- select(SDShannon, Year:shannon)
MeanShannon <- mutate(MeanShannon, sdshannon = SDShannon$shannon)
MeanShannon$Year <- as.Date(MeanShannon$Year, "%Y")

shan$sampling <- as.numeric(shan$sampling)
shan <- shan %>%
  arrange(year, sampling)
Shannondataset <- as.data.frame(sampling.dates[, "start.date"])
Shannondataset <- mutate(Shannondataset, ShannonIndex = shan$shannon)
colnames(Shannondataset) <- c("Date", "ShannonIndex")

#### Pettitt's test ####
Shannonvec <- as.vector(Shannondataset$ShannonIndex)
pettitt.test(Shannonvec)
    # data:  Shannonvec
    # U* = 35268, p-value < 2.2e-16
    # alternative hypothesis: two.sided
    # sample estimates:
    #   probable change point at time K 
    # 232 

Shannonmeanvec <- as.vector(MeanShannon$shannon)
pettitt.test(Shannonmeanvec)
    # data:  Shannonmeanvec
    # U* = 380, p-value = 0.0001826
    # alternative hypothesis: two.sided
    # sample estimates:
    #   probable change point at time K 
    # 17 

#### Mann-Kendall test ####
mk.test(Shannonvec[1:231])
    # data:  Shannonvec[1:231]
    # z = 3.5535, n = 231, p-value = 0.0003802
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S         varS          tau 
    # 4.173000e+03 1.378428e+06 1.570864e-01 

mk.test(Shannonvec[232:581])
    # data:  Shannonvec[232:581]
    # z = 1.4566, n = 350, p-value = 0.1452
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S         varS          tau 
    # 3.187000e+03 4.784208e+06 5.218174e-02 

mk.test(Shannonmeanvec[1:16])
    # data:  Shannonmeanvec[1:16]
    # z = 1.3957, n = 16, p-value = 0.1628
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S        varS         tau 
    # 32.0000000 493.3333333   0.2666667 

mk.test(Shannonmeanvec[17:45])
    # data:  Shannonmeanvec[17:45]
    # z = 0.50647, n = 29, p-value = 0.6125
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S         varS          tau 
    # 2.800000e+01 2.842000e+03 6.896552e-02 
#### Statistics ####
shapiro.test(Shannonmeanvec[1:16])
shapiro.test(Shannonmeanvec[17:45])
var.test(Shannonmeanvec[1:16], Shannonmeanvec[17:45])
# mean annual Shannon indices are approximated by a normal distribution and have equal variance

t.test(Shannonmeanvec[1:16], Shannonmeanvec[17:45])
    # data:  Shannonmeanvec[1:16] and Shannonmeanvec[17:45]
    # t = 5.4592, df = 24.82, p-value = 1.169e-05
    # alternative hypothesis: true difference in means is not equal to 0
    # 95 percent confidence interval:
    #   0.3183995 0.7044030
    # sample estimates:
    #   mean of x mean of y 
    # 1.854406  1.343005 

Shannondataset <- mutate(Shannondataset, Year = format(Shannondataset$Date, "%Y"))
tapply(Shannondataset$ShannonIndex, Shannondataset$Year, length)
# range from 5 to 20

#### Plot ####
Shannonearly <- Shannondataset[1:231,]
Shannonlate <- Shannondataset[232:581,]
Shannonearly$ShannonIndex <- round(Shannonearly$ShannonIndex, digits = 2)
Shannonlate$ShannonIndex <- round(Shannonlate$ShannonIndex, digits = 2)

Shannonplot <- 
  ggplot(Shannonearly) +
  geom_rect(xmin = -Inf, xmax = as.numeric(as.Date("1975-01-01")), ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_rect(xmin = as.numeric(as.Date("1990-01-01")), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_point(data = Shannonearly, aes(x = Date, y = ShannonIndex), size = 0.5, color = "#f99d15ff") + #significant positive slope
  geom_point(data = Shannonlate, aes(x = Date, y = ShannonIndex), size = 0.5, color = "gray40") + #non-significant slope
  geom_point(data = MeanShannon, aes(x = Year, y = shannon), pch = 1) +
  ylab(expression("Shannon Index")) +
  xlab(" ") +
  #theme(axis.text.x=element_blank()) + 
  geom_vline(xintercept = as.numeric(Shannondataset$Date[232]), lty = 5)
print(Shannonplot)


#### Long-term drivers ====
#### Data Wrangling ####
InflowData <- mutate(InflowData, DIN = NO3 + NH4)
Inflowlimited <- sample_n(InflowData, 5000)

Radiation <- Radiation %>%
  filter(Date > "1973-04-14")  %>%
  filter(Date < "2005-12-07" | Date > "2008-10-15")
Radiationlimited <- sample_n(Radiation, 5000)

# test for normality
shapiro.test(Inflowlimited$AirTemp)
shapiro.test(Inflowlimited$WindSpeed)
shapiro.test(Inflowlimited$Precipitation)
shapiro.test(InflowData$TP)
shapiro.test(InflowData$DIN)
shapiro.test(Radiationlimited$totalrad)

# Separate datasets and eliminate NAs
AirTempdataset <- InflowData %>%
  select(Date, AirTemp) %>%
  drop_na()
rownames(AirTempdataset) <- NULL

WindSpeeddataset <- InflowData %>%
  select(Date, WindSpeed) %>%
  drop_na()
rownames(WindSpeeddataset) <- NULL

Precipdataset <- InflowData %>%
  select(Date, Precipitation) %>%
  drop_na()
rownames(Precipdataset) <- NULL

TPinflowdataset <- InflowData %>%
  select(Date, TP) %>%
  drop_na()
rownames(TPinflowdataset) <- NULL

DINinflowdataset <- InflowData %>%
  select(Date, DIN) %>%
  drop_na()
rownames(DINinflowdataset) <- NULL

#### Pettitt's test ####
# Pettitt's test: detects a single changepoint in hydrological/climate series with continuous data
AirTempvec <- as.vector(AirTempdataset$AirTemp)
pettitt.test(AirTempvec)
    # data:  AirTempvec
    # U* = 3696400, p-value = 2.699e-07
    # alternative hypothesis: two.sided
    # sample estimates:
    #   probable change point at time K 
    # 10121 

WindSpeedvec <- as.vector(WindSpeeddataset$WindSpeed)
pettitt.test(WindSpeedvec)
    # data:  WindSpeedvec
    # U* = 7456000, p-value < 2.2e-16
    # alternative hypothesis: two.sided
    # sample estimates:
    #   probable change point at time K 
    # 4430 

Precipvec <- as.vector(Precipdataset$Precipitation)
pettitt.test(Precipvec)
    # data:  Precipvec
    # U* = 3790800, p-value = 1.374e-07
    # alternative hypothesis: two.sided
    # sample estimates:
    #   probable change point at time K 
    # 11270 

TPinflowvec <- as.vector(TPinflowdataset$TP)
pettitt.test(TPinflowvec)
    # data:  TPinflowvec
    # U* = 153180, p-value < 2.2e-16
    # alternative hypothesis: two.sided
    # sample estimates:
    #   probable change point at time K 
    # 677 

DINinflowvec <- as.vector(DINinflowdataset$DIN)
pettitt.test(DINinflowvec)
    # data:  DINinflowvec
    # U* = 84911, p-value = 1.378e-10
    # alternative hypothesis: two.sided
    # sample estimates:
    #   probable change point at time K 
    # 308 

Radiationvec <- as.vector(Radiation$totalrad)
pettitt.test(Radiationvec)
    # data:  Radiationvec
    # U* = 7430100, p-value < 2.2e-16
    # alternative hypothesis: two.sided
    # sample estimates:
    #   probable change point at time K 
    # 11112 

#### Mann-Kendall test ####
# Mann-Kendall test: detect monotonic trends in series of environmental/climate/hydrological data
mk.test(AirTempvec[1:10120])
    # data:  AirTempvec[1:10120]
    # z = -0.79636, n = 10120, p-value = 0.4258
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S          varS           tau 
    # -2.702660e+05  1.151746e+11 -5.286991e-03 

mk.test(AirTempvec[10121:17305])
    # data:  AirTempvec[10121:17305]
    # z = -1.0596, n = 7185, p-value = 0.2893
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S          varS           tau 
    # -2.151250e+05  4.122037e+10 -8.358701e-03 

mk.test(WindSpeedvec[1:4429])
    # data:  WindSpeedvec[1:4429]
    # z = 4.7796, n = 4429, p-value = 1.757e-06
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S         varS          tau 
    # 4.686810e+05 9.615581e+09 4.928511e-02 

mk.test(WindSpeedvec[4430:15453])
    # data:  WindSpeedvec[4430:15453]
    # z = -20.828, n = 11024, p-value < 2.2e-16
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S          varS           tau 
    # -8.022001e+06  1.483467e+11 -1.358043e-01 

mk.test(Precipvec[1:11269])
    # data:  Precipvec[1:11269]
    # z = 0.29164, n = 11269, p-value = 0.7706
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S         varS          tau 
    # 1.027980e+05 1.242402e+11 2.034715e-03 

mk.test(Precipvec[11270:17355])
    # data:  Precipvec[11270:17355]
    # z = 1.496, n = 6086, p-value = 0.1347
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S         varS          tau 
    # 2.173940e+05 2.111648e+10 1.399441e-02 

mk.test(TPinflowvec[1:676])
    # data:  TPinflowvec[1:676]
    # z = 1.9682, n = 676, p-value = 0.04905
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S         varS          tau 
    # 1.154200e+04 3.438425e+07 5.099094e-02 

mk.test(TPinflowvec[677:1276])
    # data:  TPinflowvec[677:1276]
    # z = -0.12095, n = 600, p-value = 0.9037
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S          varS           tau 
    # -5.940000e+02  2.403618e+07 -3.343057e-03 

mk.test(DINinflowvec[1:307])
    # data:  DINinflowvec[1:307]
    # z = 1.3065, n = 307, p-value = 0.1914
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S         varS          tau 
    # 2.349000e+03 3.229839e+06 5.022712e-02 

mk.test(DINinflowvec[308:1227])
    # data:  DINinflowvec[308:1227]
    # z = -8.7875, n = 920, p-value < 2.2e-16
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S          varS           tau 
    # -8.180400e+04  8.665817e+07 -1.938433e-01 

mk.test(Radiationvec[1:1111])
    # data:  Radiationvec[1:1111]
    # z = 0.031331, n = 1111, p-value = 0.975
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S         varS          tau 
    # 3.880000e+02 1.525748e+08 6.295335e-04 

mk.test(Radiationvec[11112:14923])
    # data:  Radiationvec[11112:14923]
    # z = -0.60987, n = 3812, p-value = 0.542
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S          varS           tau 
    # -4.785600e+04  6.157252e+09 -6.588474e-03 

#### test for equal sample number by year ####
tapply(InflowData$AirTemp, InflowData$Year, length)
# range from 365 to 366
tapply(InflowData$WindSpeed, InflowData$Year, length)
# range from 365 to 366
tapply(InflowData$Precipitation, InflowData$Year, length)
# range from 365 to 366
tapply(InflowData$TP, InflowData$Year, length)
# range from 365 to 366
tapply(InflowData$DIN, InflowData$Year, length)
# range from 365 to 366

#### Plots ####
AirTempearly <- AirTempdataset[1:10120,]
AirTemplate <- AirTempdataset[10121:17305,]
AirTempplot <-
  ggplot(AirTempearly) +
  geom_rect(xmin = -Inf, xmax = as.numeric(as.Date("1975-01-01")), ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_rect(xmin = as.numeric(as.Date("1990-01-01")), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_point(data = AirTempearly, aes(x = Date, y = AirTemp), size = 0.5, color = "gray40") + #non-significant slope
  geom_point(data = AirTemplate, aes(x = Date, y = AirTemp), size = 0.5, color = "gray40") + #non-significant slope
  ylab(expression(Air ~ Temp ~(degree*C))) +
  xlab(" ") #+
  #geom_vline(xintercept = as.numeric(AirTempdataset$Date[10121]), lty = 5) # took out because both periods have no trend
print(AirTempplot)
#ggsave("AirTemp.jpg", AirTempplot, dpi = 300)

WindSpeedearly <- WindSpeeddataset[1:4429,]
WindSpeedlate <- WindSpeeddataset[4430:15453,]
WindSpeedplot <-
  ggplot(WindSpeedearly) +
  geom_rect(xmin = -Inf, xmax = as.numeric(as.Date("1975-01-01")), ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_rect(xmin = as.numeric(as.Date("1990-01-01")), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_point(data = WindSpeedearly, aes(x = Date, y = WindSpeed), size = 0.5, color = "#f99d15ff") + #significant positive slope
  geom_point(data = WindSpeedlate, aes(x = Date, y = WindSpeed), size = 0.5, color = "#240c4cff") + #significant negative slope
  ylab(expression(Wind ~ Speed ~ (m/s))) +
  xlab(" ") +
  geom_vline(xintercept = as.numeric(WindSpeeddataset$Date[4430]), lty = 5)
print(WindSpeedplot)
#ggsave("WindSpeed.jpg", WindSpeedplot, dpi = 300)

Precipearly <- Precipdataset[1:11269,]
Preciplate <- Precipdataset[11270:17355,]
Precipplot <-
  ggplot(WindSpeedearly) +
  geom_rect(xmin = -Inf, xmax = as.numeric(as.Date("1975-01-01")), ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_rect(xmin = as.numeric(as.Date("1990-01-01")), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_point(data = Precipearly, aes(x = Date, y = Precipitation), size = 0.5, color = "gray40") + #non-significant slope
  geom_point(data = Preciplate, aes(x = Date, y = Precipitation), size = 0.5, color = "gray40") + #non-significant slope
  ylab(expression(Precipitation ~ (mm/d))) +
  xlab(" ") #+
  #geom_vline(xintercept = as.numeric(Precipdataset$Date[11270]), lty = 5) # took out because both periods have no trend
print(Precipplot)
#ggsave("Precip.jpg", Precipplot, dpi = 300)

TPinflowearly <- TPinflowdataset[1:676,]
TPinflowlate <- TPinflowdataset[677:1276,]
TPinflowplot <-
  ggplot(TPinflowearly) +
  geom_rect(xmin = -Inf, xmax = as.numeric(as.Date("1975-01-01")), ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_rect(xmin = as.numeric(as.Date("1990-01-01")), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_point(data = TPinflowearly, aes(x = Date, y = TP), size = 0.5, color = "#f99d15ff") + #significant positive slope
  geom_point(data = TPinflowlate, aes(x = Date, y = TP), size = 0.5, color = "gray40") + #non-significant slope
  ylab(expression(TP)) +
  xlab(" ") +
  geom_vline(xintercept = as.numeric(TPinflowdataset$Date[677]), lty = 5) # took out because both periods have no trend
print(TPinflowplot)

DINinflowearly <- DINinflowdataset[1:307,]
DINinflowlate <- DINinflowdataset[308:1227,]
DINinflowplot <-
  ggplot(DINinflowearly) +
  geom_rect(xmin = -Inf, xmax = as.numeric(as.Date("1975-01-01")), ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_rect(xmin = as.numeric(as.Date("1990-01-01")), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_point(data = DINinflowearly, aes(x = Date, y = DIN), size = 0.5, color = "gray40") + #non-significant slope
  geom_point(data = DINinflowlate, aes(x = Date, y = DIN), size = 0.5, color = "#240c4cff") + #significant negative slope
  ylab(expression(DIN)) +
  xlab(" ") +
  geom_vline(xintercept = as.numeric(DINinflowdataset$Date[308]), lty = 5) 
print(DINinflowplot)

Radiationearly <- Radiation[1:11111,]
Radiationlate <- Radiation[11112:14923,]
Radiationplot <-
  ggplot(Radiationearly) +
  geom_rect(xmin = -Inf, xmax = as.numeric(as.Date("1975-01-01")), ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_rect(xmin = as.numeric(as.Date("1990-01-01")), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_point(data = Radiationearly, aes(x = Date, y = totalrad), size = 0.5, color = "gray40") + #non-significant slope
  geom_point(data = Radiationlate, aes(x = Date, y = totalrad), size = 0.5, color = "gray40") + #non-significant slope
  ylab(expression(Radiation ~ (MJ ~ m^-2 ~ d^-1))) +
  xlab(" ") #+
#geom_vline(xintercept = as.numeric(AirTempdataset$Date[10121]), lty = 5) # took out because both periods have no trend
print(Radiationplot)
#ggsave("Radiation.jpg", Radiationplot, dpi = 300)

grid.arrange(AirTempplot, WindSpeedplot, Precipplot, TPinflowplot, DINinflowplot, ncol = 1)


#### Combined plots ====
grid.arrange(FertNtoPplot, cyanoplot, Shannonplot, TNtoTPplot, TDNtoTDPplot, PPplot, chlplot,  ncol = 1)
grid.arrange(PPplot, chlplot, ncol = 1)

#### Model Performance Analysis ----
# Code taken from ModelIterationReport
#### Spreadsheet Import ====
# Read in data
obsinit <- read.csv("Observed_IntegratedEpi_icefree_lateMay.csv", header = T)
colnames(obsinit) <- c("org.date","Year","Month","Day", "obs.TP","obs.chla","obs.TDP", "obs.PP", "obs.DOC")
obsinit$obs.DOC <- obsinit$obs.DOC * 12 # concentrations are in umol/L and need to be in ug/L)
obs_temp <- read.csv("Observed_Temperature.csv")
obs_O2 <- read.csv("Observed_Oxygen.csv")
mod <- read.csv("Output_IntegratedEpi.csv", header = F)
colnames(mod) <- c("Year", "Month", "Day", "mod.TDP", "mod.PP", "mod.DOC")
mod2 <- read.csv("Output_Depths.csv", header = F)
colnames(mod2) <- c("Year", "Month", "Day", "mod.Temp1m", "mod.Temp4m", "mod.Temp9m", "mod.Oxy2m", "mod.Oxy3m", "mod.Oxy4m", "mod.Oxy6m", "mod.Oxy8m", "mod.Oxy10m", "mod.Fe4m", "mod.Fe6m", "mod.Fe8m", "mod.Fe10m")
obs.ice <- read.csv("Lake_239_ice-on_ice-off.csv", header = T)
out.ice <- read.csv("Output_Ice.csv", header = F)

# Tidy
obs <- obsinit %>% 
  unite(date, Year, Month, Day, sep = '-') #%>%
obs <- data.frame(obs, obsinit$Month, obsinit$Year)
#obs <- na.omit(obs)
obs <- obs[,-1]
colnames(obs) <- c("date", "obs.TP","obs.chla","obs.TDP", "obs.PP", "obs.DOC", "Month", "Year")
mod <- mod %>% unite(date, Year, Month, Day, sep = '-')
mod2 <- mod2 %>% unite(date, Year, Month, Day, sep = '-')
colnames(out.ice) <- c("Year.Break", "Month.Break", "Day.Break", "Year.Freeze", "Month.Freeze", "Day.Freeze")
out.ice.edit <- out.ice %>% unite(Ice.Off.Date, Year.Break, Month.Break, Day.Break, sep = '-')
out.ice.edit <- out.ice.edit %>% unite(Ice.On.Date, Year.Freeze, Month.Freeze, Day.Freeze, sep = '-')

# Convert date structure
obs$date <- as.Date(obs$date, format = "%Y-%m-%d") 
obs_temp$date <- as.Date(obs_temp$date, format = "%m/%d/%y") 
obs_O2$date <- as.Date(obs_O2$date, format = "%d/%m/%y")
mod$date <- as.Date(mod$date, format = "%Y-%m-%d") 
mod2$date <- as.Date(mod2$date, format = "%Y-%m-%d") 
obs.ice$Ice.Off.Date <- as.Date(obs.ice$Ice.Off.Date, format = "%m/%d/%y")
obs.ice$Ice.On.Date <- as.Date(obs.ice$Ice.On.Date, format = "%m/%d/%y")
obs.daybreak <- strftime(obs.ice$Ice.Off.Date, format = "%j")
obs.dayfreeze <- strftime(obs.ice$Ice.On.Date, format = "%j")
obs.year <- format(obs.ice$Ice.On.Date, "%Y")
out.ice.edit$Ice.Off.Date <- as.Date(out.ice.edit$Ice.Off.Date, format = "%Y-%m-%d")
out.ice.edit$Ice.On.Date <- as.Date(out.ice.edit$Ice.On.Date, format = "%Y-%m-%d")
out.daybreak <- strftime(out.ice.edit$Ice.Off.Date, format = "%j")
out.dayfreeze <- strftime(out.ice.edit$Ice.On.Date, format = "%j")
out.year <- as.character(out.ice$Year.Freeze)

# Tidy
mod.match <- inner_join(obs,mod, by = "date") 
mod.match <- filter(mod.match, obs.TDP < 50)
mod2.match <- inner_join(obs_temp, mod2, by = "date")
mod3.match <- inner_join(obs_O2, mod2, by = "date")
obs.ice.edit <- cbind(obs.year, obs.ice[4], obs.ice[5], obs.daybreak, obs.dayfreeze)
colnames(obs.ice.edit) <- c("Year", "Ice.Off.Date", "Ice.On.Date", "obs.daybreak", "obs.dayfreeze")
out.ice.edit2 <- cbind(out.year, out.ice.edit, out.daybreak, out.dayfreeze)
colnames(out.ice.edit2) <- c("Year", "Ice.Off.Date", "Ice.On.Date", "out.daybreak", "out.dayfreeze")
match.ice <- inner_join(obs.ice.edit, out.ice.edit2, by = "Year") 
match.ice$obs.daybreak <- as.numeric(as.character(match.ice$obs.daybreak))
match.ice$out.daybreak <- as.numeric(as.character(match.ice$out.daybreak))
match.ice$obs.dayfreeze <- as.numeric(as.character(match.ice$obs.dayfreeze))
match.ice$out.dayfreeze <- as.numeric(as.character(match.ice$out.dayfreeze))
match.ice$Year <- as.numeric(match.ice$Year)
match.ice$obs.daybreak[1] <- NA
mod$mod.DOC <- mod$mod.DOC/1000 #convert DOC to mg/L instead of ug/L
mod.match$obs.DOC <- mod.match$obs.DOC/1000 #convert DOC to mg/L instead of ug/L
mod.match$mod.DOC <- mod.match$mod.DOC/1000 #convert DOC to mg/L instead of ug/L
mod.match <- mutate(mod.match, Year = format(mod.match$date, "%Y"))
mod2.match <- mutate(mod2.match, Year = format(mod2.match$date, "%Y"))
mod3.match <- mutate(mod3.match, Year = format(mod3.match$date, "%Y"))
mod.match <- filter(mod.match, Year != 1969)
mod2.match <- filter(mod2.match, Year != 1969)
mod3.match <- filter(mod3.match, Year != 1969)


#### Ice Fit metrics ====
icebreakregression.period1 <- lm (match.ice$obs.daybreak[match.ice$Year < 1975] ~ match.ice$out.daybreak[match.ice$Year < 1975])
summary(icebreakregression.period1)$adj.r.squared
msebreak.period1 <- mean(residuals(icebreakregression.period1)^2); rmsebreak.period1 <- sqrt(msebreak.period1); rmsebreak.period1
rmsebreak.period1/sd(na.omit(match.ice$obs.daybreak[match.ice$Year < 1975]))
NashSutcliffe.icebreak.period1 <- NSE(match.ice$out.daybreak[match.ice$Year < 1975], match.ice$obs.daybreak[match.ice$Year < 1975]); NashSutcliffe.icebreak.period1

icefreezeregression.period1 <- lm (match.ice$obs.dayfreeze[match.ice$Year < 1975] ~ match.ice$out.dayfreeze[match.ice$Year < 1975])
summary(icefreezeregression.period1)$adj.r.squared
msefreeze.period1 <- mean(residuals(icefreezeregression.period1)^2); rmsefreeze.period1 <- sqrt(msefreeze.period1); rmsefreeze.period1
rmsefreeze.period1/sd(na.omit(match.ice$obs.dayfreeze[match.ice$Year < 1975]))
NashSutcliffe.icefreeze.period1 <- NSE(match.ice$out.dayfreeze[match.ice$Year < 1975], match.ice$obs.dayfreeze[match.ice$Year < 1975]); NashSutcliffe.icefreeze.period1

icebreakregression.period2 <- lm (match.ice$obs.daybreak[match.ice$Year > 1974 & match.ice$Year < 1990] ~ match.ice$out.daybreak[match.ice$Year > 1974 & match.ice$Year < 1990])
summary(icebreakregression.period2)$adj.r.squared
msebreak.period2 <- mean(residuals(icebreakregression.period2)^2); rmsebreak.period2 <- sqrt(msebreak.period2); rmsebreak.period2
rmsebreak.period2/sd(na.omit(match.ice$obs.daybreak[match.ice$Year > 1974 & match.ice$Year < 1990]))
NashSutcliffe.icebreak.period2 <- NSE(match.ice$out.daybreak[match.ice$Year > 1974 & match.ice$Year < 1990], match.ice$obs.daybreak[match.ice$Year > 1974 & match.ice$Year < 1990]); NashSutcliffe.icebreak.period2

icefreezeregression.period2 <- lm (match.ice$obs.dayfreeze[match.ice$Year > 1974 & match.ice$Year < 1990] ~ match.ice$out.dayfreeze[match.ice$Year > 1974 & match.ice$Year < 1990])
summary(icefreezeregression.period2)$adj.r.squared
msefreeze.period2 <- mean(residuals(icefreezeregression.period2)^2); rmsefreeze.period2 <- sqrt(msefreeze.period2); rmsefreeze.period2
rmsefreeze.period2/sd(na.omit(match.ice$obs.dayfreeze[match.ice$Year > 1974 & match.ice$Year < 1990]))
NashSutcliffe.icefreeze.period2 <- NSE(match.ice$out.dayfreeze[match.ice$Year > 1974 & match.ice$Year < 1990], match.ice$obs.dayfreeze[match.ice$Year > 1974 & match.ice$Year < 1990]); NashSutcliffe.icefreeze.period2

icebreakregression.period3 <- lm (match.ice$obs.daybreak[match.ice$Year > 1989] ~ match.ice$out.daybreak[match.ice$Year > 1989])
summary(icebreakregression.period3)$adj.r.squared
msebreak.period3 <- mean(residuals(icebreakregression.period3)^2); rmsebreak.period3 <- sqrt(msebreak.period3); rmsebreak.period3
rmsebreak.period3/sd(na.omit(match.ice$obs.daybreak[match.ice$Year > 1989]))
NashSutcliffe.icebreak.period3 <- NSE(match.ice$out.daybreak[match.ice$Year > 1989], match.ice$obs.daybreak[match.ice$Year > 1989]); NashSutcliffe.icebreak.period3

icefreezeregression.period3 <- lm (match.ice$obs.dayfreeze[match.ice$Year > 1989] ~ match.ice$out.dayfreeze[match.ice$Year > 1989])
summary(icefreezeregression.period3)$adj.r.squared
msefreeze.period3 <- mean(residuals(icefreezeregression.period3)^2); rmsefreeze.period3 <- sqrt(msefreeze.period3); rmsefreeze.period3
rmsefreeze.period3/sd(na.omit(match.ice$obs.dayfreeze[match.ice$Year > 1989]))
NashSutcliffe.icefreeze.period3 <- NSE(match.ice$out.dayfreeze[match.ice$Year > 1989], match.ice$obs.dayfreeze[match.ice$Year > 1989]); NashSutcliffe.icefreeze.period3

#### Ice Plot ====
icedateplot <- 
  ggplot(data = match.ice, aes(x = Year)) +
  geom_rect(xmin = -Inf, xmax = 1975, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_rect(xmin = 1990, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_line(aes(y = out.daybreak, col = "Ice Break"), size = 0.5) +
  geom_point(aes(y = obs.daybreak, col = "Ice Break"), size = 0.5) +
  geom_line(aes(y = out.dayfreeze - 160, col = "Ice Freeze"), size = 0.5) +
  geom_point(aes(y = obs.dayfreeze - 160, col = "Ice Freeze"), size = 0.5) + 
  geom_linerange(aes(ymax = obs.daybreak-5, ymin = obs.daybreak, col = "Ice Break")) +
  geom_linerange(aes(ymax = obs.dayfreeze-175, ymin = obs.dayfreeze - 160, col = "Ice Freeze")) +
  scale_y_continuous(sec.axis = sec_axis(~.+160)) + 
  ylab("Julian Day") +
  xlab(" ") +
  scale_colour_manual("", breaks = c("Ice Break", "Ice Freeze"), values = c("#f99d15ff", "#240c4cff")) +
  geom_vline(xintercept = 1975, color = "white", alpha = 0) +
  geom_vline(xintercept = 1990, color = "white", alpha = 0) +
  theme(legend.position = "right")#, axis.text.x = element_blank())
print(icedateplot)

#ggsave("Ice.jpg", icedateplot, dpi = 300)

#### Temperature Fit metrics ====
temp1m.regression.period1 <- lm(mod2.match$obs.Temp1m[mod2.match$Year < 1975] ~ mod2.match$mod.Temp1m[mod2.match$Year < 1975])
summary(temp1m.regression.period1)$adj.r.squared
mse.temp1m.period1 <- mean(residuals(temp1m.regression.period1)^2); rmse.temp1m.period1 <- sqrt(mse.temp1m.period1); rmse.temp1m.period1
rmse.temp1m.period1/sd(na.omit(mod2.match$obs.Temp1m[mod2.match$Year < 1975]))
NashSutcliffe.temp1m.period1 <- NSE(mod2.match$mod.Temp1m[mod2.match$Year < 1975], mod2.match$obs.Temp1m[mod2.match$Year < 1975]); NashSutcliffe.temp1m.period1

temp4m.regression.period1 <- lm(mod2.match$obs.Temp4m[mod2.match$Year < 1975] ~ mod2.match$mod.Temp4m[mod2.match$Year < 1975])
summary(temp4m.regression.period1)$adj.r.squared
mse.temp4m.period1 <- mean(residuals(temp4m.regression.period1)^2); rmse.temp4m.period1 <- sqrt(mse.temp4m.period1); rmse.temp4m.period1
rmse.temp4m.period1/sd(na.omit(mod2.match$obs.Temp4m[mod2.match$Year < 1975]))
NashSutcliffe.temp4m.period1 <- NSE(mod2.match$mod.Temp4m[mod2.match$Year < 1975], mod2.match$obs.Temp4m[mod2.match$Year < 1975]); NashSutcliffe.temp4m.period1

temp9m.regression.period1 <- lm(mod2.match$obs.Temp9m[mod2.match$Year < 1975] ~ mod2.match$mod.Temp9m[mod2.match$Year < 1975])
summary(temp9m.regression.period1)$adj.r.squared
mse.temp9m.period1 <- mean(residuals(temp9m.regression.period1)^2); rmse.temp9m.period1 <- sqrt(mse.temp9m.period1); rmse.temp9m.period1
rmse.temp9m.period1/sd(na.omit(mod2.match$obs.Temp9m[mod2.match$Year < 1975]))
NashSutcliffe.temp9m.period1 <- NSE(mod2.match$mod.Temp9m[mod2.match$Year < 1975], mod2.match$obs.Temp9m[mod2.match$Year < 1975]); NashSutcliffe.temp9m.period1

temp1m.regression.period2 <- lm(mod2.match$obs.Temp1m[mod2.match$Year > 1974 & mod2.match$Year < 1990] ~ mod2.match$mod.Temp1m[mod2.match$Year > 1974 & mod2.match$Year < 1990])
summary(temp1m.regression.period2)$adj.r.squared
mse.temp1m.period2 <- mean(residuals(temp1m.regression.period2)^2); rmse.temp1m.period2 <- sqrt(mse.temp1m.period2); rmse.temp1m.period2
rmse.temp1m.period2/sd(na.omit(mod2.match$obs.Temp1m[mod2.match$Year > 1974 & mod2.match$Year < 1990]))
NashSutcliffe.temp1m.period2 <- NSE(mod2.match$mod.Temp1m[mod2.match$Year > 1974 & mod2.match$Year < 1990], mod2.match$obs.Temp1m[mod2.match$Year > 1974 & mod2.match$Year < 1990]); NashSutcliffe.temp1m.period2

temp4m.regression.period2 <- lm(mod2.match$obs.Temp4m[mod2.match$Year > 1974 & mod2.match$Year < 1990] ~ mod2.match$mod.Temp4m[mod2.match$Year > 1974 & mod2.match$Year < 1990])
summary(temp4m.regression.period2)$adj.r.squared
mse.temp4m.period2 <- mean(residuals(temp4m.regression.period2)^2); rmse.temp4m.period2 <- sqrt(mse.temp4m.period2); rmse.temp4m.period2
rmse.temp4m.period2/sd(na.omit(mod2.match$obs.Temp4m[mod2.match$Year > 1974 & mod2.match$Year < 1990]))
NashSutcliffe.temp4m.period2 <- NSE(mod2.match$mod.Temp4m[mod2.match$Year > 1974 & mod2.match$Year < 1990], mod2.match$obs.Temp4m[mod2.match$Year > 1974 & mod2.match$Year < 1990]); NashSutcliffe.temp4m.period2

temp9m.regression.period2 <- lm(mod2.match$obs.Temp9m[mod2.match$Year > 1974 & mod2.match$Year < 1990] ~ mod2.match$mod.Temp9m[mod2.match$Year > 1974 & mod2.match$Year < 1990])
summary(temp9m.regression.period2)$adj.r.squared
mse.temp9m.period2 <- mean(residuals(temp9m.regression.period2)^2); rmse.temp9m.period2 <- sqrt(mse.temp9m.period2); rmse.temp9m.period2
rmse.temp9m.period2/sd(na.omit(mod2.match$obs.Temp9m[mod2.match$Year > 1974 & mod2.match$Year < 1990]))
NashSutcliffe.temp9m.period2 <- NSE(mod2.match$mod.Temp9m[mod2.match$Year > 1974 & mod2.match$Year < 1990], mod2.match$obs.Temp9m[mod2.match$Year > 1974 & mod2.match$Year < 1990]); NashSutcliffe.temp9m.period2

temp1m.regression.period3 <- lm(mod2.match$obs.Temp1m[mod2.match$Year > 1989] ~ mod2.match$mod.Temp1m[mod2.match$Year > 1989])
summary(temp1m.regression.period3)$adj.r.squared
mse.temp1m.period3 <- mean(residuals(temp1m.regression.period3)^2); rmse.temp1m.period3 <- sqrt(mse.temp1m.period3); rmse.temp1m.period3
rmse.temp1m.period3/sd(na.omit(mod2.match$obs.Temp1m[mod2.match$Year > 1989]))
NashSutcliffe.temp1m.period3 <- NSE(mod2.match$mod.Temp1m[mod2.match$Year > 1989], mod2.match$obs.Temp1m[mod2.match$Year > 1989]); NashSutcliffe.temp1m.period3

temp4m.regression.period3 <- lm(mod2.match$obs.Temp4m[mod2.match$Year > 1989] ~ mod2.match$mod.Temp4m[mod2.match$Year > 1989])
summary(temp4m.regression.period3)$adj.r.squared
mse.temp4m.period3 <- mean(residuals(temp4m.regression.period3)^2); rmse.temp4m.period3 <- sqrt(mse.temp4m.period3); rmse.temp4m.period3
rmse.temp4m.period3/sd(na.omit(mod2.match$obs.Temp4m[mod2.match$Year > 1989]))
NashSutcliffe.temp4m.period3 <- NSE(mod2.match$mod.Temp4m[mod2.match$Year > 1989], mod2.match$obs.Temp4m[mod2.match$Year > 1989]); NashSutcliffe.temp4m.period3

temp9m.regression.period3 <- lm(mod2.match$obs.Temp9m[mod2.match$Year > 1989] ~ mod2.match$mod.Temp9m[mod2.match$Year > 1989])
summary(temp9m.regression.period3)$adj.r.squared
mse.temp9m.period3 <- mean(residuals(temp9m.regression.period3)^2); rmse.temp9m.period3 <- sqrt(mse.temp9m.period3); rmse.temp9m.period3
rmse.temp9m.period3/sd(na.omit(mod2.match$obs.Temp9m[mod2.match$Year > 1989]))
NashSutcliffe.temp9m.period3 <- NSE(mod2.match$mod.Temp9m[mod2.match$Year > 1989], mod2.match$obs.Temp9m[mod2.match$Year > 1989]); NashSutcliffe.temp9m.period3

#### Temperature Plot ====
tempplot <- ggplot(mod2) +
  geom_rect(xmin = -Inf, xmax = as.numeric(as.Date("1975-01-01")), ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_rect(xmin = as.numeric(as.Date("1990-01-01")), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_line(data = mod2, aes(x = date, y = mod.Temp1m, col = "1 m"), size = 0.25) +
  geom_line(data = mod2, aes(x = date, y = mod.Temp4m, col = "4 m"), size = 0.25) +
  geom_line(data = mod2, aes(x = date, y = mod.Temp9m, col = "9 m"), size = 0.25) +  
  geom_point(data = mod2.match, aes(x = date, y = obs.Temp1m, col = "1 m"), pch = 19, size = 0.5) +
  geom_point(data = mod2.match, aes(x = date, y = obs.Temp4m, col = "4 m"), pch = 19, size = 0.5) +
  geom_point(data = mod2.match, aes(x = date, y = obs.Temp9m, col = "9 m"), pch = 19, size = 0.5) +
  ylab(expression("Temp " ( degree*C))) +
  xlab(" ") +
  scale_colour_manual("", breaks = c("1 m", "4 m", "9 m"), values = c("#f99d15ff", "#d14a42ff", "#240c4cff")) +
  theme(legend.position = "right")#, axis.text.x = element_blank())
print(tempplot)

#ggsave("Temp.jpg", tempplot, dpi = 300)
#### PP Fit Metrics ====
PP.regression.period1 <- lm(mod.match$obs.PP[mod.match$Year < 1975] ~ mod.match$mod.PP[mod.match$Year < 1975])
summary(PP.regression.period1)$adj.r.squared
mse.PP.period1 <- mean(residuals(PP.regression.period1)^2); rmse.PP.period1 <- sqrt(mse.PP.period1); rmse.PP.period1
rmse.PP.period1/sd(na.omit(mod.match$obs.PP[mod.match$Year < 1975]))
NashSutcliffe.PP.period1 <-NSE(mod.match$mod.PP[mod.match$Year < 1975], mod.match$obs.PP[mod.match$Year < 1975]); NashSutcliffe.PP.period1

PP.regression.period2 <- lm(mod.match$obs.PP[mod.match$Year > 1974 & mod.match$Year < 1990] ~ mod.match$mod.PP[mod.match$Year > 1974 & mod.match$Year < 1990])
summary(PP.regression.period2)$adj.r.squared
mse.PP.period2 <- mean(residuals(PP.regression.period2)^2); rmse.PP.period2 <- sqrt(mse.PP.period2); rmse.PP.period2
rmse.PP.period2/sd(na.omit(mod.match$obs.PP[mod.match$Year > 1974 & mod.match$Year < 1990]))
NashSutcliffe.PP.period2 <-NSE(mod.match$mod.PP[mod.match$Year > 1974 & mod.match$Year < 1990], mod.match$obs.PP[mod.match$Year > 1974 & mod.match$Year < 1990]); NashSutcliffe.PP.period2

PP.regression.period3 <- lm(mod.match$obs.PP[mod.match$Year > 1989] ~ mod.match$mod.PP[mod.match$Year > 1989])
summary(PP.regression.period3)$adj.r.squared
mse.PP.period3 <- mean(residuals(PP.regression.period3)^2); rmse.PP.period3 <- sqrt(mse.PP.period3); rmse.PP.period3
rmse.PP.period3/sd(na.omit(mod.match$obs.PP[mod.match$Year > 1989]))
NashSutcliffe.PP.period3 <-NSE(mod.match$mod.PP[mod.match$Year > 1989], mod.match$obs.PP[mod.match$Year > 1989]); NashSutcliffe.PP.period3

#### PP Plot ====
PPmodelplot <- ggplot(mod) +
  geom_rect(xmin = -Inf, xmax = as.numeric(as.Date("1975-01-01")), ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_rect(xmin = as.numeric(as.Date("1990-01-01")), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_line(data = mod, aes(x = date, y = mod.PP, col = "Modeled"), size = 0.5) +
  geom_point(data = mod.match, aes(x = date, y = obs.PP, col = "Observed"), pch = 19, size = 0.75) +
  ylab(expression(PP ~ (mu*g / L))) +
  xlab(" ") +
  ylim(c(0, 100)) +
  scale_colour_manual("", breaks = c("Observed", "Modeled"), values = c("#d14a42ff", "#240c4cff")) +
  theme(legend.position = "right")#, axis.text.x = element_blank()) 
print(PPmodelplot)

#ggsave("PP.jpg", PPmodelplot, dpi = 300)
#### Cumulative PP Fit Metrics ====
# Cumulative Sum of PP
mod.match.daily <- right_join(obs, mod, by = "date")
mod.match.daily$Month <- as.numeric(format(mod.match.daily$date, "%m"))
mod.match.daily$Year <- format(mod.match.daily$date, "%Y")
mod.match.daily <- mod.match.daily %>% 
  mutate(Day = as.numeric(format(mod.match.daily$date, "%d")))  %>% 
  filter(Year != 1969) %>% #remove spinup year
  filter(Year != 1996) %>% #remove pike year
  filter(Month > 4 & Month < 11)  %>% 
  filter(Month == 5 & Day > 15 | Month > 5 )
mod.match.daily$obs.PP[mod.match.daily$Year == 2016 & mod.match.daily$Month == 10 & mod.match.daily$Day == 31] <- 0
mod.match.daily$obs.PP[mod.match.daily$Month == 5 & mod.match.daily$Day == 16] <- 0
obs.PP.interpolated <- na.approx(mod.match.daily$obs.PP)
mod.match.daily$obs.PP <- obs.PP.interpolated
mod.match.cumulative <- ddply(mod.match.daily, .(Year), transform, Cum.Sum.mod.PP = cumsum(mod.PP), Cum.Sum.obs.PP = cumsum(obs.PP))

# calculate end-of season comparison between cumulative modeled vs. observed PP
end.of.season.cumulative <- as.data.frame(mod.match.cumulative) %>%
  filter(Day == 31 & Month == 10) %>%
  select(date, Month, Year, Day, Cum.Sum.mod.PP, Cum.Sum.obs.PP) %>%
  mutate(mod.minus.obs.PP = Cum.Sum.mod.PP -Cum.Sum.obs.PP) %>%
  mutate(mod.divided.obs.PP = Cum.Sum.mod.PP/Cum.Sum.obs.PP)

summary(end.of.season.cumulative$mod.divided.obs.PP[end.of.season.cumulative$Year < 1975])
sd(end.of.season.cumulative$mod.divided.obs.PP[end.of.season.cumulative$Year < 1975])

summary(end.of.season.cumulative$mod.divided.obs.PP[end.of.season.cumulative$Year > 1974 & end.of.season.cumulative$Year < 1990])
sd(end.of.season.cumulative$mod.divided.obs.PP[end.of.season.cumulative$Year > 1974 & end.of.season.cumulative$Year < 1990])

summary(end.of.season.cumulative$mod.divided.obs.PP[end.of.season.cumulative$Year > 1989])
sd(end.of.season.cumulative$mod.divided.obs.PP[end.of.season.cumulative$Year > 1989])

#### Cumulative PP Plot ====
PPcumulativeplot <- ggplot(data = mod.match.cumulative) +
  geom_rect(xmin = -Inf, xmax = as.numeric(as.Date("1975-01-01")), ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_rect(xmin = as.numeric(as.Date("1990-01-01")), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_point(aes(x = date, y = Cum.Sum.mod.PP, col = "Modeled"), size = 1) +
  geom_point(aes(x = date, y = Cum.Sum.obs.PP, col = "Observed"), size = 0.5) +
  ylab(expression(Cumulative ~ PP ~ (mu*g / L))) +
  xlab(" ") +
  scale_colour_manual("", breaks = c("Observed", "Modeled"), values = c("#d14a42ff", "#240c4cff")) +
  theme(legend.position = "right")#, axis.text.x = element_blank()) 
print(PPcumulativeplot)

#ggsave("PPcumulative.jpg", PPcumulativeplot, dpi = 300)

#### Model PP residuals Fit Metrics ====
PPresiduals <- mod.match %>%
  select(date, Year) %>%
  mutate(residuals.PP = mod.match$obs.PP - mod.match$mod.PP) %>%
  na.omit()

summary(PPresiduals$residuals.PP[PPresiduals$Year < 1975])
sd(PPresiduals$residuals.PP[PPresiduals$Year < 1975])

summary(PPresiduals$residuals.PP[PPresiduals$Year > 1974 & PPresiduals$Year < 1990])
sd(PPresiduals$residuals.PP[PPresiduals$Year > 1974 & PPresiduals$Year < 1990])

summary(PPresiduals$residuals.PP[PPresiduals$Year > 1989])
sd(PPresiduals$residuals.PP[PPresiduals$Year > 1989])

var.test(PPresiduals$residuals.PP[PPresiduals$Year < 1975], PPresiduals$residuals.PP[PPresiduals$Year > 1974 & PPresiduals$Year < 1990])
var.test(PPresiduals$residuals.PP[PPresiduals$Year < 1975], PPresiduals$residuals.PP[PPresiduals$Year > 1989])
var.test(PPresiduals$residuals.PP[PPresiduals$Year > 1974 & PPresiduals$Year < 1990], PPresiduals$residuals.PP[PPresiduals$Year > 1989])

# Variance of PP residuals is significantly lower in period 3 than in period 1 or 2.
# Variance of PP in periods 1 and 2 does not differ significantly.

t.test(PPresiduals$residuals.PP[PPresiduals$Year < 1975], PPresiduals$residuals.PP[PPresiduals$Year > 1974 & PPresiduals$Year < 1990])
t.test(PPresiduals$residuals.PP[PPresiduals$Year < 1975], PPresiduals$residuals.PP[PPresiduals$Year > 1989])
t.test(PPresiduals$residuals.PP[PPresiduals$Year > 1974 & PPresiduals$Year < 1990], PPresiduals$residuals.PP[PPresiduals$Year > 1989])

# PP residuals are significantly lower in period 3 than in period 1 or 2.
# PP residuals do not differ significantly in periods 1 and 2.

pettitt.test(PPresiduals$residuals.PP)
    # data:  PPresiduals$residuals.PP
    # U* = 30907, p-value = 4.319e-08
    # alternative hypothesis: two.sided
    # sample estimates:
    #   probable change point at time K 
    # 412 
mk.test(PPresiduals$residuals.PP[1:411])
    # data:  PPresiduals$residuals.PP[1:411]
    # z = 2.1524, n = 411, p-value = 0.03137
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S         varS          tau 
    # 5.990000e+03 7.742097e+06 7.109413e-02 

mk.test(PPresiduals$residuals.PP[412:687])
    # data:  PPresiduals$residuals.PP[412:687]
    # z = 0.90569, n = 276, p-value = 0.3651
    # alternative hypothesis: true S is not equal to 0
    # sample estimates:
    #   S         varS          tau 
    # 1.389000e+03 2.348676e+06 3.660417e-02 

PPresiduals <- mutate(PPresiduals, abs.residuals.PP = abs(residuals.PP))
tapply(PPresiduals$abs.residuals.PP, PPresiduals$Year, sum)
min(tapply(PPresiduals$abs.residuals.PP, PPresiduals$Year, sum))

#### Model PP residuals Plot ====
PPresidualsearly <- PPresiduals[1:411,]
PPresidualslate <- PPresiduals[412:687,]
PPresidualsplot <- ggplot(PPresiduals) +
  geom_rect(xmin = -Inf, xmax = as.numeric(as.Date("1975-01-01")), ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_rect(xmin = as.numeric(as.Date("1990-01-01")), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_point(data = PPresidualsearly, aes(x = date, y = residuals.PP), size = 0.5, color = "#f99d15ff") + 
  geom_point(data = PPresidualslate, aes(x = date, y = residuals.PP), size = 0.5, color = "gray40") + 
  ylab(expression(PP ~ residuals ~ (mu*g / L))) +
  xlab(" ") +
  scale_y_continuous(limits = c(-75, 75), breaks = c(-75, -50, -25, 0, 25, 50, 75)) +
  geom_vline(xintercept = as.numeric(PPresiduals$date[412]), lty = 5) +
  theme(legend.position = "none") 
print(PPresidualsplot)

#PPresidualsplot <- ggplot(PPresiduals) +
  # geom_rect(xmin = -Inf, xmax = as.numeric(as.Date("1975-01-01")), ymin = -Inf, ymax = Inf, fill = "gray90") +
  # geom_rect(xmin = as.numeric(as.Date("1990-01-01")), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  # geom_point(data = PPresiduals, aes(x = date, y = residuals.PP), size = 0.5, color = "#d14a42ff") + 
  # ylab(expression(PP ~ residuals ~ (mu*g / L))) +
  # xlab(" ") +
  # scale_y_continuous(limits = c(-75, 75), breaks = c(-75, -50, -25, 0, 25, 50, 75)) +
  # theme(legend.position = "none") 
#print(PPresidualsplot)

PPresiduals.bestyear <- filter(PPresiduals, Year == 1992)
PPresidualsplot.bestyear <- ggplot(PPresiduals.bestyear) +
    geom_point(aes(x = date, y = residuals.PP), size = 0.5, color = "#f99d15ff") + 
    geom_point(aes(x = date, y = residuals.PP), size = 0.5, color = "gray40") + 
    ylab(expression(PP ~ residuals ~ (mu*g / L))) +
    xlab(" ") +
    scale_y_continuous(limits = c(-75, 75), breaks = c(-75, -50, -25, 0, 25, 50, 75)) +
    theme(legend.position = "none") 
print(PPresidualsplot.bestyear)

mod.bestyear <- filter(mod, date > "1999-12-31" & date < "2002-01-01")
mod.match.bestyear <- filter(mod.match, Year == 2000 | Year == 2001)

  PPmodelplot.bestyear <- ggplot() +
    geom_line(data = mod.bestyear, aes(x = date, y = mod.PP, col = "PP modeled"), size = 0.5) +
    geom_point(data = mod.match.bestyear, aes(x = date, y = obs.PP, col = "PP observed"), pch = 19, size = 2) +
    #geom_line(data = mod.bestyear, aes(x = date, y = mod.TDP, col = "TDP modeled"), size = 0.5) +
    #geom_point(data = mod.match.bestyear, aes(x = date, y = obs.TDP, col = "TDP observed"), pch = 19, size = 2) +    
    ylab(expression(PP ~ (mu*g / L))) +
    xlab(" ") +
    scale_colour_manual("", breaks = c("PP modeled", "PP observed", "TDP modeled", "TDP observed"), values = c("#240c4cff", "#240c4cff", "#f27d16ff", "#f27d16ff")) +
    theme(legend.position = "none") 
print(PPmodelplot.bestyear)

#ggsave("PPbestyear.jpg", PPmodelplot.bestyear, dpi = 300)
#### TDP Fit Metrics ====
TDP.regression.period1 <- lm(mod.match$obs.TDP[mod.match$Year < 1975] ~ mod.match$mod.TDP[mod.match$Year < 1975])
summary(TDP.regression.period1)$adj.r.squared
mse.TDP.period1 <- mean(residuals(TDP.regression.period1)^2); rmse.TDP.period1 <- sqrt(mse.TDP.period1); rmse.TDP.period1
rmse.TDP.period1/sd(na.omit(mod.match$obs.TDP[mod.match$Year < 1975]))
NashSutcliffe.TDP.period1 <-NSE(mod.match$mod.TDP[mod.match$Year < 1975], mod.match$obs.TDP[mod.match$Year < 1975]); NashSutcliffe.TDP.period1

TDP.regression.period2 <- lm(mod.match$obs.TDP[mod.match$Year > 1974 & mod.match$Year < 1990] ~ mod.match$mod.TDP[mod.match$Year > 1974 & mod.match$Year < 1990])
summary(TDP.regression.period2)$adj.r.squared
mse.TDP.period2 <- mean(residuals(TDP.regression.period2)^2); rmse.TDP.period2 <- sqrt(mse.TDP.period2); rmse.TDP.period2
rmse.TDP.period2/sd(na.omit(mod.match$obs.TDP[mod.match$Year > 1974 & mod.match$Year < 1990]))
NashSutcliffe.TDP.period2 <-NSE(mod.match$mod.TDP[mod.match$Year > 1974 & mod.match$Year < 1990], mod.match$obs.TDP[mod.match$Year > 1974 & mod.match$Year < 1990]); NashSutcliffe.TDP.period2

TDP.regression.period3 <- lm(mod.match$obs.TDP[mod.match$Year > 1989] ~ mod.match$mod.TDP[mod.match$Year > 1989])
summary(TDP.regression.period3)$adj.r.squared
mse.TDP.period3 <- mean(residuals(TDP.regression.period3)^2); rmse.TDP.period3 <- sqrt(mse.TDP.period3); rmse.TDP.period3
rmse.TDP.period3/sd(na.omit(mod.match$obs.TDP[mod.match$Year > 1989]))
NashSutcliffe.TDP.period3 <-NSE(mod.match$mod.TDP[mod.match$Year > 1989], mod.match$obs.TDP[mod.match$Year > 1989]); NashSutcliffe.TDP.period3

#### TDP Plot ====
TDPplot <- ggplot(mod) +
  geom_rect(xmin = -Inf, xmax = as.numeric(as.Date("1975-01-01")), ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_rect(xmin = as.numeric(as.Date("1990-01-01")), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_line(data = mod, aes(x = date, y = mod.TDP, col = "Modeled"), size = 0.25) +
  geom_point(data = mod.match, aes(x = date, y = obs.TDP, col = "Observed"), pch = 19, size = 0.5) +
  ylab(expression(TDP ~ (mu*g / L))) +
  xlab(" ") +
  scale_colour_manual("", breaks = c("Observed", "Modeled"), values = c("#d14a42ff", "#240c4cff")) +
  theme(legend.position = "top") 
print(TDPplot)

#### O2 Fit Metrics ====
O2.regression.period1 <- lm(mod3.match$obs.O2.4m[mod3.match$Year < 1975] ~ mod3.match$mod.Oxy4m[mod3.match$Year < 1975])
summary(O2.regression.period1)$adj.r.squared
mse.O2.period1 <- mean(residuals(O2.regression.period1)^2); rmse.O2.period1 <- sqrt(mse.O2.period1); rmse.O2.period1
rmse.O2.period1/sd(na.omit(mod3.match$obs.O2.4m[mod3.match$Year < 1975]))
NashSutcliffe.O24m.period1 <- NSE(mod3.match$mod.Oxy4m[mod3.match$Year < 1975], mod3.match$obs.O2.4m[mod3.match$Year < 1975]); NashSutcliffe.O24m.period1 

O2.regression.period2 <- lm(mod3.match$obs.O2.4m[mod3.match$Year > 1974 & mod3.match$Year < 1990] ~ mod3.match$mod.Oxy4m[mod3.match$Year > 1974 & mod3.match$Year < 1990])
summary(O2.regression.period2)$adj.r.squared
mse.O2.period2 <- mean(residuals(O2.regression.period2)^2); rmse.O2.period2 <- sqrt(mse.O2.period2); rmse.O2.period2
rmse.O2.period2/sd(na.omit(mod3.match$obs.O2.4m[mod3.match$Year > 1974 & mod3.match$Year < 1990]))
NashSutcliffe.O24m.period2 <- NSE(mod3.match$mod.Oxy4m[mod3.match$Year > 1974 & mod3.match$Year < 1990], mod3.match$obs.O2.4m[mod3.match$Year > 1974 & mod3.match$Year < 1990]); NashSutcliffe.O24m.period2

O2.regression.period3 <- lm(mod3.match$obs.O2.4m[mod3.match$Year > 1989] ~ mod3.match$mod.Oxy4m[mod3.match$Year > 1989])
summary(O2.regression.period3)$adj.r.squared
mse.O2.period3 <- mean(residuals(O2.regression.period3)^2); rmse.O2.period3 <- sqrt(mse.O2.period3); rmse.O2.period3
rmse.O2.period3/sd(na.omit(mod3.match$obs.O2.4m[mod3.match$Year > 1989]))
NashSutcliffe.O24m.period3 <- NSE(mod3.match$mod.Oxy4m[mod3.match$Year > 1989], mod3.match$obs.O2.4m[mod3.match$Year > 1989]); NashSutcliffe.O24m.period3 

#### O2 Plot ====
O2plot <- ggplot(mod2) +
  geom_rect(xmin = -Inf, xmax = as.numeric(as.Date("1975-01-01")), ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_rect(xmin = as.numeric(as.Date("1990-01-01")), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_line(data = mod2, aes(x = date, y = mod.Oxy4m, col = "Modeled"), size = 0.25) +
  geom_point(data = mod3.match, aes(x = date, y = obs.O2.4m, col = "Observed"), pch = 19, size = 0.5) +
  ylab(expression(DO ~ (mg / L))) +
  xlab(" ") +
  scale_colour_manual("", breaks = c("Observed", "Modeled"), values = c("#d14a42ff", "#240c4cff")) +
  theme(legend.position = "top") 
print(O2plot)

#### DOC Fit Metrics ====
DOC.regression.period1 <- lm(mod.match$obs.DOC[mod.match$Year < 1975] ~ mod.match$mod.DOC[mod.match$Year < 1975])
summary(DOC.regression.period1)$adj.r.squared
mse.DOC.period1 <- mean(residuals(DOC.regression.period1)^2); rmse.DOC.period1 <- sqrt(mse.DOC.period1); rmse.DOC.period1
rmse.DOC.period1/sd(na.omit(mod.match$obs.DOC[mod.match$Year < 1975]))
NashSutcliffe.DOC.period1 <-NSE(mod.match$mod.DOC[mod.match$Year < 1975], mod.match$obs.DOC[mod.match$Year < 1975]); NashSutcliffe.DOC.period1

DOC.regression.period2 <- lm(mod.match$obs.DOC[mod.match$Year > 1974 & mod.match$Year < 1990] ~ mod.match$mod.DOC[mod.match$Year > 1974 & mod.match$Year < 1990])
summary(DOC.regression.period2)$adj.r.squared
mse.DOC.period2 <- mean(residuals(DOC.regression.period2)^2); rmse.DOC.period2 <- sqrt(mse.DOC.period2); rmse.DOC.period2
rmse.DOC.period2/sd(na.omit(mod.match$obs.DOC[mod.match$Year > 1974 & mod.match$Year < 1990]))
NashSutcliffe.DOC.period2 <-NSE(mod.match$mod.DOC[mod.match$Year > 1974 & mod.match$Year < 1990], mod.match$obs.DOC[mod.match$Year > 1974 & mod.match$Year < 1990]); NashSutcliffe.DOC.period2

DOC.regression.period3 <- lm(mod.match$obs.DOC[mod.match$Year > 1989] ~ mod.match$mod.DOC[mod.match$Year > 1989])
summary(DOC.regression.period3)$adj.r.squared
mse.DOC.period3 <- mean(residuals(DOC.regression.period3)^2); rmse.DOC.period3 <- sqrt(mse.DOC.period3); rmse.DOC.period3
rmse.DOC.period3/sd(na.omit(mod.match$obs.DOC[mod.match$Year > 1989]))
NashSutcliffe.DOC.period3 <-NSE(mod.match$mod.DOC[mod.match$Year > 1989], mod.match$obs.DOC[mod.match$Year > 1989]); NashSutcliffe.DOC.period3

#### DOC Plot ====
DOCplot <- ggplot(mod) +
  geom_rect(xmin = -Inf, xmax = as.numeric(as.Date("1975-01-01")), ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_rect(xmin = as.numeric(as.Date("1990-01-01")), xmax = Inf, ymin = -Inf, ymax = Inf, fill = "gray90") +
  geom_line(data = mod, aes(x = date, y = mod.DOC, col = "Modeled"), size = 0.25) +
  geom_point(data = mod.match, aes(x = date, y = obs.DOC, col = "Observed"), pch = 19, size = 0.5) +
  ylab(expression(DOC ~ (mg / L))) +
  xlab(" ") +
  scale_colour_manual("", breaks = c("Observed", "Modeled"), values = c("#d14a42ff", "#240c4cff")) +
  theme(legend.position = "top") 
print(DOCplot)

#### Combined Plot ====
grid.arrange(icedateplot, tempplot, PPmodelplot, PPcumulativeplot, PPresidualsplot,  ncol = 1)
grid.arrange(O2plot, DOCplot, TDPplot, ncol = 1)
#### Target Diagram ====
#### Temperature 1 m ####
model.mean.Temp1m.period1 <- mean(mod2.match$mod.Temp1m[mod2.match$Year < 1975])
obs.mean.Temp1m.period1 <- mean(mod2.match$obs.Temp1m[mod2.match$Year < 1975])
model.sd.Temp1m.period1 <- sd(mod2.match$mod.Temp1m[mod2.match$Year < 1975])
obs.sd.Temp1m.period1 <- sd(mod2.match$obs.Temp1m[mod2.match$Year < 1975])
model.residuals.Temp1m.period1 <- mod2.match$mod.Temp1m[mod2.match$Year < 1975] - model.mean.Temp1m.period1
obs.residuals.Temp1m.period1 <- mod2.match$obs.Temp1m[mod2.match$Year < 1975] - obs.mean.Temp1m.period1
unbiased.RMSD.Temp1m.period1 <- sqrt(mean((model.residuals.Temp1m.period1 - obs.residuals.Temp1m.period1)^2))
normalized.bias.Temp1m.period1 <- (model.mean.Temp1m.period1 - obs.mean.Temp1m.period1)/obs.sd.Temp1m.period1
normalized.unbiased.RMSD.Temp1m.period1 <- ((model.sd.Temp1m.period1 - obs.sd.Temp1m.period1)/obs.sd.Temp1m.period1) * unbiased.RMSD.Temp1m.period1

model.mean.Temp1m.period2 <- mean(mod2.match$mod.Temp1m[mod2.match$Year > 1974 & mod2.match$Year < 1990])
obs.mean.Temp1m.period2 <- mean(mod2.match$obs.Temp1m[mod2.match$Year > 1974 & mod2.match$Year < 1990])
model.sd.Temp1m.period2 <- sd(mod2.match$mod.Temp1m[mod2.match$Year > 1974 & mod2.match$Year < 1990])
obs.sd.Temp1m.period2 <- sd(mod2.match$obs.Temp1m[mod2.match$Year > 1974 & mod2.match$Year < 1990])
model.residuals.Temp1m.period2 <- mod2.match$mod.Temp1m[mod2.match$Year > 1974 & mod2.match$Year < 1990] - model.mean.Temp1m.period2
obs.residuals.Temp1m.period2 <- mod2.match$obs.Temp1m[mod2.match$Year > 1974 & mod2.match$Year < 1990] - obs.mean.Temp1m.period2
unbiased.RMSD.Temp1m.period2 <- sqrt(mean((model.residuals.Temp1m.period2 - obs.residuals.Temp1m.period2)^2))
normalized.bias.Temp1m.period2 <- (model.mean.Temp1m.period2 - obs.mean.Temp1m.period2)/obs.sd.Temp1m.period2
normalized.unbiased.RMSD.Temp1m.period2 <- ((model.sd.Temp1m.period2 - obs.sd.Temp1m.period2)/obs.sd.Temp1m.period2) * unbiased.RMSD.Temp1m.period2

model.mean.Temp1m.period3 <- mean(mod2.match$mod.Temp1m[mod2.match$Year > 1989])
obs.mean.Temp1m.period3 <- mean(mod2.match$obs.Temp1m[mod2.match$Year > 1989])
model.sd.Temp1m.period3 <- sd(mod2.match$mod.Temp1m[mod2.match$Year > 1989])
obs.sd.Temp1m.period3 <- sd(mod2.match$obs.Temp1m[mod2.match$Year > 1989])
model.residuals.Temp1m.period3 <- mod2.match$mod.Temp1m[mod2.match$Year > 1989] - model.mean.Temp1m.period3
obs.residuals.Temp1m.period3 <- mod2.match$obs.Temp1m[mod2.match$Year > 1989] - obs.mean.Temp1m.period3
unbiased.RMSD.Temp1m.period3 <- sqrt(mean((model.residuals.Temp1m.period3 - obs.residuals.Temp1m.period3)^2))
normalized.bias.Temp1m.period3 <- (model.mean.Temp1m.period3 - obs.mean.Temp1m.period3)/obs.sd.Temp1m.period3
normalized.unbiased.RMSD.Temp1m.period3 <- ((model.sd.Temp1m.period3 - obs.sd.Temp1m.period3)/obs.sd.Temp1m.period3) * unbiased.RMSD.Temp1m.period3

#### Temp 4 m ####
model.mean.Temp4m.period1 <- mean(mod2.match$mod.Temp4m[mod2.match$Year < 1975])
obs.mean.Temp4m.period1 <- mean(mod2.match$obs.Temp4m[mod2.match$Year < 1975])
model.sd.Temp4m.period1 <- sd(mod2.match$mod.Temp4m[mod2.match$Year < 1975])
obs.sd.Temp4m.period1 <- sd(mod2.match$obs.Temp4m[mod2.match$Year < 1975])
model.residuals.Temp4m.period1 <- mod2.match$mod.Temp4m[mod2.match$Year < 1975] - model.mean.Temp4m.period1
obs.residuals.Temp4m.period1 <- mod2.match$obs.Temp4m[mod2.match$Year < 1975] - obs.mean.Temp4m.period1
unbiased.RMSD.Temp4m.period1 <- sqrt(mean((model.residuals.Temp4m.period1 - obs.residuals.Temp4m.period1)^2))
normalized.bias.Temp4m.period1 <- (model.mean.Temp4m.period1 - obs.mean.Temp4m.period1)/obs.sd.Temp4m.period1
normalized.unbiased.RMSD.Temp4m.period1 <- ((model.sd.Temp4m.period1 - obs.sd.Temp4m.period1)/obs.sd.Temp4m.period1) * unbiased.RMSD.Temp4m.period1

model.mean.Temp4m.period2 <- mean(mod2.match$mod.Temp4m[mod2.match$Year > 1974 & mod2.match$Year < 1990])
obs.mean.Temp4m.period2 <- mean(mod2.match$obs.Temp4m[mod2.match$Year > 1974 & mod2.match$Year < 1990])
model.sd.Temp4m.period2 <- sd(mod2.match$mod.Temp4m[mod2.match$Year > 1974 & mod2.match$Year < 1990])
obs.sd.Temp4m.period2 <- sd(mod2.match$obs.Temp4m[mod2.match$Year > 1974 & mod2.match$Year < 1990])
model.residuals.Temp4m.period2 <- mod2.match$mod.Temp4m[mod2.match$Year > 1974 & mod2.match$Year < 1990] - model.mean.Temp4m.period2
obs.residuals.Temp4m.period2 <- mod2.match$obs.Temp4m[mod2.match$Year > 1974 & mod2.match$Year < 1990] - obs.mean.Temp4m.period2
unbiased.RMSD.Temp4m.period2 <- sqrt(mean((model.residuals.Temp4m.period2 - obs.residuals.Temp4m.period2)^2))
normalized.bias.Temp4m.period2 <- (model.mean.Temp4m.period2 - obs.mean.Temp4m.period2)/obs.sd.Temp4m.period2
normalized.unbiased.RMSD.Temp4m.period2 <- ((model.sd.Temp4m.period2 - obs.sd.Temp4m.period2)/obs.sd.Temp4m.period2) * unbiased.RMSD.Temp4m.period2

model.mean.Temp4m.period3 <- mean(mod2.match$mod.Temp4m[mod2.match$Year > 1989])
obs.mean.Temp4m.period3 <- mean(mod2.match$obs.Temp4m[mod2.match$Year > 1989])
model.sd.Temp4m.period3 <- sd(mod2.match$mod.Temp4m[mod2.match$Year > 1989])
obs.sd.Temp4m.period3 <- sd(mod2.match$obs.Temp4m[mod2.match$Year > 1989])
model.residuals.Temp4m.period3 <- mod2.match$mod.Temp4m[mod2.match$Year > 1989] - model.mean.Temp4m.period3
obs.residuals.Temp4m.period3 <- mod2.match$obs.Temp4m[mod2.match$Year > 1989] - obs.mean.Temp4m.period3
unbiased.RMSD.Temp4m.period3 <- sqrt(mean((model.residuals.Temp4m.period3 - obs.residuals.Temp4m.period3)^2))
normalized.bias.Temp4m.period3 <- (model.mean.Temp4m.period3 - obs.mean.Temp4m.period3)/obs.sd.Temp4m.period3
normalized.unbiased.RMSD.Temp4m.period3 <- ((model.sd.Temp4m.period3 - obs.sd.Temp4m.period3)/obs.sd.Temp4m.period3) * unbiased.RMSD.Temp4m.period3

#### Temp 9 m ####
Temp9m.match <- select(mod2.match, obs.Temp9m, mod.Temp9m, Year)
Temp9m.match <- na.omit(Temp9m.match)

model.mean.Temp9m.period1 <- mean(Temp9m.match$mod.Temp9m[Temp9m.match$Year < 1975])
obs.mean.Temp9m.period1 <- mean(Temp9m.match$obs.Temp9m[Temp9m.match$Year < 1975])
model.sd.Temp9m.period1 <- sd(Temp9m.match$mod.Temp9m[Temp9m.match$Year < 1975])
obs.sd.Temp9m.period1 <- sd(Temp9m.match$obs.Temp9m[Temp9m.match$Year < 1975])
model.residuals.Temp9m.period1 <- Temp9m.match$mod.Temp9m[Temp9m.match$Year < 1975] - model.mean.Temp9m.period1
obs.residuals.Temp9m.period1 <- Temp9m.match$obs.Temp9m[Temp9m.match$Year < 1975] - obs.mean.Temp9m.period1
unbiased.RMSD.Temp9m.period1 <- sqrt(mean((model.residuals.Temp9m.period1 - obs.residuals.Temp9m.period1)^2))
normalized.bias.Temp9m.period1 <- (model.mean.Temp9m.period1 - obs.mean.Temp9m.period1)/obs.sd.Temp9m.period1
normalized.unbiased.RMSD.Temp9m.period1 <- ((model.sd.Temp9m.period1 - obs.sd.Temp9m.period1)/obs.sd.Temp9m.period1) * unbiased.RMSD.Temp9m.period1

model.mean.Temp9m.period2 <- mean(Temp9m.match$mod.Temp9m[Temp9m.match$Year > 1974 & Temp9m.match$Year < 1990])
obs.mean.Temp9m.period2 <- mean(Temp9m.match$obs.Temp9m[Temp9m.match$Year > 1974 & Temp9m.match$Year < 1990])
model.sd.Temp9m.period2 <- sd(Temp9m.match$mod.Temp9m[Temp9m.match$Year > 1974 & Temp9m.match$Year < 1990])
obs.sd.Temp9m.period2 <- sd(Temp9m.match$obs.Temp9m[Temp9m.match$Year > 1974 & Temp9m.match$Year < 1990])
model.residuals.Temp9m.period2 <- Temp9m.match$mod.Temp9m[Temp9m.match$Year > 1974 & Temp9m.match$Year < 1990] - model.mean.Temp9m.period2
obs.residuals.Temp9m.period2 <- Temp9m.match$obs.Temp9m[Temp9m.match$Year > 1974 & Temp9m.match$Year < 1990] - obs.mean.Temp9m.period2
unbiased.RMSD.Temp9m.period2 <- sqrt(mean((model.residuals.Temp9m.period2 - obs.residuals.Temp9m.period2)^2))
normalized.bias.Temp9m.period2 <- (model.mean.Temp9m.period2 - obs.mean.Temp9m.period2)/obs.sd.Temp9m.period2
normalized.unbiased.RMSD.Temp9m.period2 <- ((model.sd.Temp9m.period2 - obs.sd.Temp9m.period2)/obs.sd.Temp9m.period2) * unbiased.RMSD.Temp9m.period2

model.mean.Temp9m.period3 <- mean(Temp9m.match$mod.Temp9m[Temp9m.match$Year > 1989])
obs.mean.Temp9m.period3 <- mean(Temp9m.match$obs.Temp9m[Temp9m.match$Year > 1989])
model.sd.Temp9m.period3 <- sd(Temp9m.match$mod.Temp9m[Temp9m.match$Year > 1989])
obs.sd.Temp9m.period3 <- sd(Temp9m.match$obs.Temp9m[Temp9m.match$Year > 1989])
model.residuals.Temp9m.period3 <- Temp9m.match$mod.Temp9m[Temp9m.match$Year > 1989] - model.mean.Temp9m.period3
obs.residuals.Temp9m.period3 <- Temp9m.match$obs.Temp9m[Temp9m.match$Year > 1989] - obs.mean.Temp9m.period3
unbiased.RMSD.Temp9m.period3 <- sqrt(mean((model.residuals.Temp9m.period3 - obs.residuals.Temp9m.period3)^2))
normalized.bias.Temp9m.period3 <- (model.mean.Temp9m.period3 - obs.mean.Temp9m.period3)/obs.sd.Temp9m.period3
normalized.unbiased.RMSD.Temp9m.period3 <- ((model.sd.Temp9m.period3 - obs.sd.Temp9m.period3)/obs.sd.Temp9m.period3) * unbiased.RMSD.Temp9m.period3

#### O2 4 m ####
mod3.match <- select(mod3.match, obs.O2.4m, mod.Oxy4m, Year)
mod3.match <- na.omit(mod3.match)

model.mean.O24m.period1 <- mean(mod3.match$mod.Oxy4m[mod3.match$Year < 1975])
obs.mean.O24m.period1 <- mean(mod3.match$obs.O2.4m[mod3.match$Year < 1975])
model.sd.O24m.period1 <- sd(mod3.match$mod.Oxy4m[mod3.match$Year < 1975])
obs.sd.O24m.period1 <- sd(mod3.match$obs.O2.4m[mod3.match$Year < 1975])
model.residuals.O24m.period1 <- mod3.match$mod.Oxy4m[mod3.match$Year < 1975] - model.mean.O24m.period1
obs.residuals.O24m.period1 <- mod3.match$obs.O2.4m[mod3.match$Year < 1975] - obs.mean.O24m.period1
unbiased.RMSD.O24m.period1 <- sqrt(mean((model.residuals.O24m.period1  - obs.residuals.O24m.period1)^2))
normalized.bias.O24m.period1 <- (model.mean.O24m.period1 - obs.mean.O24m.period1)/obs.sd.O24m.period1
normalized.unbiased.RMSD.O24m.period1 <- ((model.sd.O24m.period1 - obs.sd.O24m.period1)/obs.sd.O24m.period1) * unbiased.RMSD.O24m.period1

model.mean.O24m.period2 <- mean(mod3.match$mod.Oxy4m[mod3.match$Year > 1974 & mod3.match$Year < 1990])
obs.mean.O24m.period2 <- mean(mod3.match$obs.O2.4m[mod3.match$Year > 1974 & mod3.match$Year < 1990])
model.sd.O24m.period2 <- sd(mod3.match$mod.Oxy4m[mod3.match$Year > 1974 & mod3.match$Year < 1990])
obs.sd.O24m.period2 <- sd(mod3.match$obs.O2.4m[mod3.match$Year > 1974 & mod3.match$Year < 1990])
model.residuals.O24m.period2 <- mod3.match$mod.Oxy4m[mod3.match$Year > 1974 & mod3.match$Year < 1990] - model.mean.O24m.period2
obs.residuals.O24m.period2 <- mod3.match$obs.O2.4m[mod3.match$Year > 1974 & mod3.match$Year < 1990] - obs.mean.O24m.period2
unbiased.RMSD.O24m.period2 <- sqrt(mean((model.residuals.O24m.period2  - obs.residuals.O24m.period2)^2))
normalized.bias.O24m.period2 <- (model.mean.O24m.period2 - obs.mean.O24m.period2)/obs.sd.O24m.period2
normalized.unbiased.RMSD.O24m.period2 <- ((model.sd.O24m.period2 - obs.sd.O24m.period2)/obs.sd.O24m.period2) * unbiased.RMSD.O24m.period2

model.mean.O24m.period3 <- mean(mod3.match$mod.Oxy4m[mod3.match$Year > 1989])
obs.mean.O24m.period3 <- mean(mod3.match$obs.O2.4m[mod3.match$Year > 1989])
model.sd.O24m.period3 <- sd(mod3.match$mod.Oxy4m[mod3.match$Year > 1989])
obs.sd.O24m.period3 <- sd(mod3.match$obs.O2.4m[mod3.match$Year > 1989])
model.residuals.O24m.period3 <- mod3.match$mod.Oxy4m[mod3.match$Year > 1989] - model.mean.O24m.period3
obs.residuals.O24m.period3 <- mod3.match$obs.O2.4m[mod3.match$Year > 1989] - obs.mean.O24m.period3
unbiased.RMSD.O24m.period3 <- sqrt(mean((model.residuals.O24m.period3  - obs.residuals.O24m.period3)^2))
normalized.bias.O24m.period3 <- (model.mean.O24m.period3 - obs.mean.O24m.period3)/obs.sd.O24m.period3
normalized.unbiased.RMSD.O24m.period3 <- ((model.sd.O24m.period3 - obs.sd.O24m.period3)/obs.sd.O24m.period3) * unbiased.RMSD.O24m.period3

#### DOC integrated epi ####
DOC.match <- select(mod.match, obs.DOC, mod.DOC, Year)
DOC.match <- na.omit(DOC.match)

model.mean.DOC.period1 <- mean(DOC.match$mod.DOC[DOC.match$Year < 1975])
obs.mean.DOC.period1 <- mean(DOC.match$obs.DOC[DOC.match$Year < 1975])
model.sd.DOC.period1 <- sd(DOC.match$mod.DOC[DOC.match$Year < 1975])
obs.sd.DOC.period1 <- sd(DOC.match$obs.DOC[DOC.match$Year < 1975])
model.residuals.DOC.period1 <- DOC.match$mod.DOC[DOC.match$Year < 1975] - model.mean.DOC.period1
obs.residuals.DOC.period1 <- DOC.match$obs.DOC[DOC.match$Year < 1975] - obs.mean.DOC.period1
unbiased.RMSD.DOC.period1 <- sqrt(mean((model.residuals.DOC.period1 - obs.residuals.DOC.period1)^2))
normalized.bias.DOC.period1 <- (model.mean.DOC.period1 - obs.mean.DOC.period1)/obs.sd.DOC.period1
normalized.unbiased.RMSD.DOC.period1 <- ((model.sd.DOC.period1 - obs.sd.DOC.period1)/obs.sd.DOC.period1) * unbiased.RMSD.DOC.period1

model.mean.DOC.period2 <- mean(DOC.match$mod.DOC[DOC.match$Year > 1974 & DOC.match$Year < 1990])
obs.mean.DOC.period2 <- mean(DOC.match$obs.DOC[DOC.match$Year > 1974 & DOC.match$Year < 1990])
model.sd.DOC.period2 <- sd(DOC.match$mod.DOC[DOC.match$Year > 1974 & DOC.match$Year < 1990])
obs.sd.DOC.period2 <- sd(DOC.match$obs.DOC[DOC.match$Year > 1974 & DOC.match$Year < 1990])
model.residuals.DOC.period2 <- DOC.match$mod.DOC[DOC.match$Year > 1974 & DOC.match$Year < 1990] - model.mean.DOC.period2
obs.residuals.DOC.period2 <- DOC.match$obs.DOC[DOC.match$Year > 1974 & DOC.match$Year < 1990] - obs.mean.DOC.period2
unbiased.RMSD.DOC.period2 <- sqrt(mean((model.residuals.DOC.period2 - obs.residuals.DOC.period2)^2))
normalized.bias.DOC.period2 <- (model.mean.DOC.period2 - obs.mean.DOC.period2)/obs.sd.DOC.period2
normalized.unbiased.RMSD.DOC.period2 <- ((model.sd.DOC.period2 - obs.sd.DOC.period2)/obs.sd.DOC.period2) * unbiased.RMSD.DOC.period2

model.mean.DOC.period3 <- mean(DOC.match$mod.DOC[DOC.match$Year > 1989])
obs.mean.DOC.period3 <- mean(DOC.match$obs.DOC[DOC.match$Year > 1989])
model.sd.DOC.period3 <- sd(DOC.match$mod.DOC[DOC.match$Year > 1989])
obs.sd.DOC.period3 <- sd(DOC.match$obs.DOC[DOC.match$Year > 1989])
model.residuals.DOC.period3 <- DOC.match$mod.DOC[DOC.match$Year > 1989] - model.mean.DOC.period3
obs.residuals.DOC.period3 <- DOC.match$obs.DOC[DOC.match$Year > 1989] - obs.mean.DOC.period3
unbiased.RMSD.DOC.period3 <- sqrt(mean((model.residuals.DOC.period3 - obs.residuals.DOC.period3)^2))
normalized.bias.DOC.period3 <- (model.mean.DOC.period3 - obs.mean.DOC.period3)/obs.sd.DOC.period3
normalized.unbiased.RMSD.DOC.period3 <- ((model.sd.DOC.period3 - obs.sd.DOC.period3)/obs.sd.DOC.period3) * unbiased.RMSD.DOC.period3

#### TDP integrated epi ####
TDP.match <- select(mod.match, obs.TDP, mod.TDP, Year)
TDP.match <- na.omit(TDP.match)

model.mean.TDP.period1 <- mean(TDP.match$mod.TDP[TDP.match$Year < 1975])
obs.mean.TDP.period1 <- mean(TDP.match$obs.TDP[TDP.match$Year < 1975])
model.sd.TDP.period1 <- sd(TDP.match$mod.TDP[TDP.match$Year < 1975])
obs.sd.TDP.period1 <- sd(TDP.match$obs.TDP[TDP.match$Year < 1975])
model.residuals.TDP.period1 <- TDP.match$mod.TDP[TDP.match$Year < 1975] - model.mean.TDP.period1
obs.residuals.TDP.period1 <- TDP.match$obs.TDP[TDP.match$Year < 1975] - obs.mean.TDP.period1
unbiased.RMSD.TDP.period1 <- sqrt(mean((model.residuals.TDP.period1 - obs.residuals.TDP.period1)^2))
normalized.bias.TDP.period1 <- (model.mean.TDP.period1 - obs.mean.TDP.period1)/obs.sd.TDP.period1
normalized.unbiased.RMSD.TDP.period1 <- ((model.sd.TDP.period1 - obs.sd.TDP.period1)/obs.sd.TDP.period1) * unbiased.RMSD.TDP.period1

model.mean.TDP.period2 <- mean(TDP.match$mod.TDP[TDP.match$Year > 1974 & TDP.match$Year < 1990])
obs.mean.TDP.period2 <- mean(TDP.match$obs.TDP[TDP.match$Year > 1974 & TDP.match$Year < 1990])
model.sd.TDP.period2 <- sd(TDP.match$mod.TDP[TDP.match$Year > 1974 & TDP.match$Year < 1990])
obs.sd.TDP.period2 <- sd(TDP.match$obs.TDP[TDP.match$Year > 1974 & TDP.match$Year < 1990])
model.residuals.TDP.period2 <- TDP.match$mod.TDP[TDP.match$Year > 1974 & TDP.match$Year < 1990] - model.mean.TDP.period2
obs.residuals.TDP.period2 <- TDP.match$obs.TDP[TDP.match$Year > 1974 & TDP.match$Year < 1990] - obs.mean.TDP.period2
unbiased.RMSD.TDP.period2 <- sqrt(mean((model.residuals.TDP.period2 - obs.residuals.TDP.period2)^2))
normalized.bias.TDP.period2 <- (model.mean.TDP.period2 - obs.mean.TDP.period2)/obs.sd.TDP.period2
normalized.unbiased.RMSD.TDP.period2 <- ((model.sd.TDP.period2 - obs.sd.TDP.period2)/obs.sd.TDP.period2) * unbiased.RMSD.TDP.period2

model.mean.TDP.period3 <- mean(TDP.match$mod.TDP[TDP.match$Year > 1989])
obs.mean.TDP.period3 <- mean(TDP.match$obs.TDP[TDP.match$Year > 1989])
model.sd.TDP.period3 <- sd(TDP.match$mod.TDP[TDP.match$Year > 1989])
obs.sd.TDP.period3 <- sd(TDP.match$obs.TDP[TDP.match$Year > 1989])
model.residuals.TDP.period3 <- TDP.match$mod.TDP[TDP.match$Year > 1989] - model.mean.TDP.period3
obs.residuals.TDP.period3 <- TDP.match$obs.TDP[TDP.match$Year > 1989] - obs.mean.TDP.period3
unbiased.RMSD.TDP.period3 <- sqrt(mean((model.residuals.TDP.period3 - obs.residuals.TDP.period3)^2))
normalized.bias.TDP.period3 <- (model.mean.TDP.period3 - obs.mean.TDP.period3)/obs.sd.TDP.period3
normalized.unbiased.RMSD.TDP.period3 <- ((model.sd.TDP.period3 - obs.sd.TDP.period3)/obs.sd.TDP.period3) * unbiased.RMSD.TDP.period3

#### PP integrated epi ####
PP.match <- select(mod.match, obs.PP, mod.PP, Year)
PP.match <- na.omit(PP.match)

model.mean.PP.period1 <- mean(PP.match$mod.PP[PP.match$Year < 1975])
obs.mean.PP.period1 <- mean(PP.match$obs.PP[PP.match$Year < 1975])
model.sd.PP.period1 <- sd(PP.match$mod.PP[PP.match$Year < 1975])
obs.sd.PP.period1 <- sd(PP.match$obs.PP[PP.match$Year < 1975])
model.residuals.PP.period1 <- PP.match$mod.PP[PP.match$Year < 1975] - model.mean.PP.period1
obs.residuals.PP.period1 <- PP.match$obs.PP[PP.match$Year < 1975] - obs.mean.PP.period1
unbiased.RMSD.PP.period1 <- sqrt(mean((model.residuals.PP.period1 - obs.residuals.PP.period1)^2))
normalized.bias.PP.period1 <- (model.mean.PP.period1 - obs.mean.PP.period1)/obs.sd.PP.period1
normalized.unbiased.RMSD.PP.period1 <- ((model.sd.PP.period1 - obs.sd.PP.period1)/obs.sd.PP.period1) * unbiased.RMSD.PP.period1

model.mean.PP.period2 <- mean(PP.match$mod.PP[PP.match$Year > 1974 & PP.match$Year < 1990])
obs.mean.PP.period2 <- mean(PP.match$obs.PP[PP.match$Year > 1974 & PP.match$Year < 1990])
model.sd.PP.period2 <- sd(PP.match$mod.PP[PP.match$Year > 1974 & PP.match$Year < 1990])
obs.sd.PP.period2 <- sd(PP.match$obs.PP[PP.match$Year > 1974 & PP.match$Year < 1990])
model.residuals.PP.period2 <- PP.match$mod.PP[PP.match$Year > 1974 & PP.match$Year < 1990] - model.mean.PP.period2
obs.residuals.PP.period2 <- PP.match$obs.PP[PP.match$Year > 1974 & PP.match$Year < 1990] - obs.mean.PP.period2
unbiased.RMSD.PP.period2 <- sqrt(mean((model.residuals.PP.period2 - obs.residuals.PP.period2)^2))
normalized.bias.PP.period2 <- (model.mean.PP.period2 - obs.mean.PP.period2)/obs.sd.PP.period2
normalized.unbiased.RMSD.PP.period2 <- ((model.sd.PP.period2 - obs.sd.PP.period2)/obs.sd.PP.period2) * unbiased.RMSD.PP.period2

model.mean.PP.period3 <- mean(PP.match$mod.PP[PP.match$Year > 1989])
obs.mean.PP.period3 <- mean(PP.match$obs.PP[PP.match$Year > 1989])
model.sd.PP.period3 <- sd(PP.match$mod.PP[PP.match$Year > 1989])
obs.sd.PP.period3 <- sd(PP.match$obs.PP[PP.match$Year > 1989])
model.residuals.PP.period3 <- PP.match$mod.PP[PP.match$Year > 1989] - model.mean.PP.period3
obs.residuals.PP.period3 <- PP.match$obs.PP[PP.match$Year > 1989] - obs.mean.PP.period3
unbiased.RMSD.PP.period3 <- sqrt(mean((model.residuals.PP.period3 - obs.residuals.PP.period3)^2))
normalized.bias.PP.period3 <- (model.mean.PP.period3 - obs.mean.PP.period3)/obs.sd.PP.period3
normalized.unbiased.RMSD.PP.period3 <- ((model.sd.PP.period3 - obs.sd.PP.period3)/obs.sd.PP.period3) * unbiased.RMSD.PP.period3

#### Combined dataframe ####
TargetDiagramData <- data.frame(Variable = c("Temp 1 m ", "Temp 4 m", "Temp 9 m", "DO 4 m", "DOC", "TDP", "PP",
                                             "Temp 1 m ", "Temp 4 m", "Temp 9 m", "DO 4 m", "DOC", "TDP", "PP",
                                             "Temp 1 m ", "Temp 4 m", "Temp 9 m", "DO 4 m", "DOC", "TDP", "PP"), 
                                Period = c("High N:P","High N:P", "High N:P", "High N:P", "High N:P", "High N:P", "High N:P", 
                                           "Low N:P",  "Low N:P",  "Low N:P",  "Low N:P",  "Low N:P",  "Low N:P",  "Low N:P",
                                           "P only", "P only", "P only", "P only", "P only", "P only", "P only"), 
                                Normalized.Bias = c(normalized.bias.Temp1m.period1, normalized.bias.Temp4m.period1, normalized.bias.Temp9m.period1, 
                                                    normalized.bias.O24m.period1, normalized.bias.DOC.period1, normalized.bias.TDP.period1, normalized.bias.PP.period1,
                                                    normalized.bias.Temp1m.period2, normalized.bias.Temp4m.period2, normalized.bias.Temp9m.period2, 
                                                    normalized.bias.O24m.period2, normalized.bias.DOC.period2, normalized.bias.TDP.period2, normalized.bias.PP.period2,
                                                    normalized.bias.Temp1m.period3, normalized.bias.Temp4m.period3, normalized.bias.Temp9m.period3, 
                                                    normalized.bias.O24m.period3, normalized.bias.DOC.period3, normalized.bias.TDP.period3, normalized.bias.PP.period3), 
                                Normalized.Unbiased.RMSD = c(normalized.unbiased.RMSD.Temp1m.period1, normalized.unbiased.RMSD.Temp4m.period1, normalized.unbiased.RMSD.Temp9m.period1,
                                                             normalized.unbiased.RMSD.O24m.period1, normalized.unbiased.RMSD.DOC.period1, normalized.unbiased.RMSD.TDP.period1, normalized.unbiased.RMSD.PP.period1,
                                                             normalized.unbiased.RMSD.Temp1m.period2, normalized.unbiased.RMSD.Temp4m.period2, normalized.unbiased.RMSD.Temp9m.period2,
                                                             normalized.unbiased.RMSD.O24m.period2, normalized.unbiased.RMSD.DOC.period2, normalized.unbiased.RMSD.TDP.period2, normalized.unbiased.RMSD.PP.period2,
                                                             normalized.unbiased.RMSD.Temp1m.period3, normalized.unbiased.RMSD.Temp4m.period3, normalized.unbiased.RMSD.Temp9m.period3,
                                                             normalized.unbiased.RMSD.O24m.period3, normalized.unbiased.RMSD.DOC.period3, normalized.unbiased.RMSD.TDP.period3, normalized.unbiased.RMSD.PP.period3)) 

TargetDiagramData$Variable <- factor(TargetDiagramData$Variable, levels = c("Temp 1 m ", "Temp 4 m", "Temp 9 m", "DO 4 m", "DOC", "TDP", "PP"))

TargetDiagramData2 <- data.frame(Variable = c("Temp 1 m ", "Temp 4 m", "Temp 9 m", "PP",
                                             "Temp 1 m ", "Temp 4 m", "Temp 9 m", "PP",
                                             "Temp 1 m ", "Temp 4 m", "Temp 9 m", "PP"), 
                                Period = c("High N:P","High N:P", "High N:P", "High N:P", 
                                           "Low N:P",  "Low N:P",  "Low N:P",  "Low N:P", 
                                           "P only", "P only", "P only", "P only"), 
                                Normalized.Bias = c(normalized.bias.Temp1m.period1, normalized.bias.Temp4m.period1, normalized.bias.Temp9m.period1, normalized.bias.PP.period1,
                                                    normalized.bias.Temp1m.period2, normalized.bias.Temp4m.period2, normalized.bias.Temp9m.period2,  normalized.bias.PP.period2,
                                                    normalized.bias.Temp1m.period3, normalized.bias.Temp4m.period3, normalized.bias.Temp9m.period3, normalized.bias.PP.period3), 
                                Normalized.Unbiased.RMSD = c(normalized.unbiased.RMSD.Temp1m.period1, normalized.unbiased.RMSD.Temp4m.period1, normalized.unbiased.RMSD.Temp9m.period1,normalized.unbiased.RMSD.PP.period1,
                                                             normalized.unbiased.RMSD.Temp1m.period2, normalized.unbiased.RMSD.Temp4m.period2, normalized.unbiased.RMSD.Temp9m.period2,normalized.unbiased.RMSD.PP.period2,
                                                             normalized.unbiased.RMSD.Temp1m.period3, normalized.unbiased.RMSD.Temp4m.period3, normalized.unbiased.RMSD.Temp9m.period3,normalized.unbiased.RMSD.PP.period3)) 

TargetDiagramData$Variable <- factor(TargetDiagramData$Variable, levels = c("Temp 1 m ", "Temp 4 m", "Temp 9 m", "PP"))

#### Target Plot ####
TargetPlot <- 
ggplot(TargetDiagramData, aes(x = Normalized.Unbiased.RMSD, y = Normalized.Bias, shape = Variable, color = Period, fill = Period)) + 
  geom_point(size = 3) + 
  annotate("path", x=0+1*cos(seq(0,2*pi,length.out=100)), y=0+1*sin(seq(0,2*pi,length.out=100))) +
  annotate("path", x=0+0.75*cos(seq(0,2*pi,length.out=100)), y=0+0.75*sin(seq(0,2*pi,length.out=100))) +
  xlim(-15, 15) +
  ylim(-2, 2) +
  scale_color_manual(values = c("#f99d15ff", "#d14a42ff", "#240c4cff")) +
  scale_shape_manual(values = c(0, 1, 2, 18, 15, 19, 17)) + 
  ylab(expression(Normalized ~ Bias)) +
  xlab(expression(Normalized ~ Unbiased ~ RMSD)) +
  theme(legend.title = element_blank())
print(TargetPlot)

TargetPlot2 <- 
  ggplot(TargetDiagramData2, aes(x = Normalized.Unbiased.RMSD, y = Normalized.Bias, shape = Variable, color = Period, fill = Period)) + 
  geom_point(size = 3) + 
  annotate("path", x=0+1*cos(seq(0,2*pi,length.out=100)), y=0+1*sin(seq(0,2*pi,length.out=100))) +
  annotate("path", x=0+0.75*cos(seq(0,2*pi,length.out=100)), y=0+0.75*sin(seq(0,2*pi,length.out=100))) +
  xlim(-12, 1.5) +
  ylim(-1.5, 1.5) +
  scale_color_manual(values = c("#f99d15ff", "#d14a42ff", "#240c4cff")) +
  scale_shape_manual(values = c(19, 1, 0, 2)) + 
  ylab(expression(Normalized ~ Bias)) +
  xlab(expression(Normalized ~ Unbiased ~ RMSD)) +
  geom_vline(xintercept = 0, lty = 5) +
  geom_hline(yintercept = 0, lty = 5) +
  theme(legend.title = element_blank())
print(TargetPlot2)

#ggsave("Target.jpg", TargetPlot2, dpi = 300)
