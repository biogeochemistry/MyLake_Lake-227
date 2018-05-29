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
          plot.margin=unit(c(10,10,10,10),"pt"))}
theme_set(theme_std())

scales::show_col(colormap(colormap = colormaps$plasma, nshades=15))
scales::show_col(colormap(colormap = colormaps$inferno, nshades=15))
scales::show_col(colormap(colormap = colormaps$viridis, nshades=15))

#### Historical Data Analysis ----
# Figures 1 and S1 for manuscript
# Code taken from 
#### Stoichiometry ==== 
#### Cyano % of biomass ====
#### Phytoplankton ====
#### Long-term drivers ====

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


#### Ice ====
#### Fit metrics ####
icebreakregression.period1 <- lm (match.ice$obs.daybreak[match.ice$Year < 1975] ~ match.ice$out.daybreak[match.ice$Year < 1975])
summary(icebreakregression.period1)$adj.r.squared
msebreak.period1 <- mean(residuals(icebreakregression.period1)^2); rmsebreak.period1 <- sqrt(msebreak.period1); rmsebreak.period1

icefreezeregression.period1 <- lm (match.ice$obs.dayfreeze[match.ice$Year < 1975] ~ match.ice$out.dayfreeze[match.ice$Year < 1975])
summary(icefreezeregression.period1)$adj.r.squared
msefreeze.period1 <- mean(residuals(icefreezeregression.period1)^2); rmsefreeze.period1 <- sqrt(msefreeze.period1); rmsefreeze.period1

icebreakregression.period2 <- lm (match.ice$obs.daybreak[match.ice$Year > 1974 & match.ice$Year < 1990] ~ match.ice$out.daybreak[match.ice$Year > 1974 & match.ice$Year < 1990])
summary(icebreakregression.period2)$adj.r.squared
msebreak.period2 <- mean(residuals(icebreakregression.period2)^2); rmsebreak.period2 <- sqrt(msebreak.period2); rmsebreak.period2

icefreezeregression.period2 <- lm (match.ice$obs.dayfreeze[match.ice$Year > 1974 & match.ice$Year < 1990] ~ match.ice$out.dayfreeze[match.ice$Year > 1974 & match.ice$Year < 1990])
summary(icefreezeregression.period2)$adj.r.squared
msefreeze.period2 <- mean(residuals(icefreezeregression.period2)^2); rmsefreeze.period2 <- sqrt(msefreeze.period2); rmsefreeze.period2

icebreakregression.period3 <- lm (match.ice$obs.daybreak[match.ice$Year > 1989] ~ match.ice$out.daybreak[match.ice$Year > 1989])
summary(icebreakregression.period3)$adj.r.squared
msebreak.period3 <- mean(residuals(icebreakregression.period3)^2); rmsebreak.period3 <- sqrt(msebreak.period3); rmsebreak.period3

icefreezeregression.period3 <- lm (match.ice$obs.dayfreeze[match.ice$Year > 1989] ~ match.ice$out.dayfreeze[match.ice$Year > 1989])
summary(icefreezeregression.period3)$adj.r.squared
msefreeze.period3 <- mean(residuals(icefreezeregression.period3)^2); rmsefreeze.period3 <- sqrt(msefreeze.period3); rmsefreeze.period3

#### Plot ####
icedateplot <- 
  ggplot(data = match.ice, aes(x = Year, group = 1)) +
  geom_line(aes(y = out.daybreak, col = "Ice Break"), size = 0.5) +
  geom_point(aes(y = obs.daybreak, col = "Ice Break"), size = 0.5) +
  geom_line(aes(y = out.dayfreeze, col = "Ice Freeze"), size = 0.5) +
  geom_point(aes(y = obs.dayfreeze, col = "Ice Freeze"), size = 0.5) + 
  scale_y_reverse() + 
  ylab("Julian Day") +
  scale_colour_manual("", breaks = c("Ice Break", "Ice Freeze"), values = c("#f99d15ff", "#240c4cff")) +
  theme(legend.position = "top")
print(icedateplot)

#### Temperature ====
#### Fit metrics ####
temp1m.regression.period1 <- lm(mod2.match$obs.Temp1m[mod2.match$Year < 1975] ~ mod2.match$mod.Temp1m[mod2.match$Year < 1975])
summary(temp1m.regression.period1)$adj.r.squared
mse.temp1m.period1 <- mean(residuals(temp1m.regression.period1)^2); rmse.temp1m.period1 <- sqrt(mse.temp1m.period1); rmse.temp1m.period1

temp4m.regression.period1 <- lm(mod2.match$obs.Temp4m[mod2.match$Year < 1975] ~ mod2.match$mod.Temp4m[mod2.match$Year < 1975])
summary(temp4m.regression.period1)$adj.r.squared
mse.temp4m.period1 <- mean(residuals(temp4m.regression.period1)^2); rmse.temp4m.period1 <- sqrt(mse.temp4m.period1); rmse.temp4m.period1

temp9m.regression.period1 <- lm(mod2.match$obs.Temp9m[mod2.match$Year < 1975] ~ mod2.match$mod.Temp9m[mod2.match$Year < 1975])
summary(temp9m.regression.period1)$adj.r.squared
mse.temp9m.period1 <- mean(residuals(temp9m.regression.period1)^2); rmse.temp9m.period1 <- sqrt(mse.temp9m.period1); rmse.temp9m.period1

temp1m.regression.period2 <- lm(mod2.match$obs.Temp1m[mod2.match$Year > 1974 & mod2.match$Year < 1990] ~ mod2.match$mod.Temp1m[mod2.match$Year > 1974 & mod2.match$Year < 1990])
summary(temp1m.regression.period2)$adj.r.squared
mse.temp1m.period2 <- mean(residuals(temp1m.regression.period2)^2); rmse.temp1m.period2 <- sqrt(mse.temp1m.period2); rmse.temp1m.period2

temp4m.regression.period2 <- lm(mod2.match$obs.Temp4m[mod2.match$Year > 1974 & mod2.match$Year < 1990] ~ mod2.match$mod.Temp4m[mod2.match$Year > 1974 & mod2.match$Year < 1990])
summary(temp4m.regression.period2)$adj.r.squared
mse.temp4m.period2 <- mean(residuals(temp4m.regression.period2)^2); rmse.temp4m.period2 <- sqrt(mse.temp4m.period2); rmse.temp4m.period2

temp9m.regression.period2 <- lm(mod2.match$obs.Temp9m[mod2.match$Year > 1974 & mod2.match$Year < 1990] ~ mod2.match$mod.Temp9m[mod2.match$Year > 1974 & mod2.match$Year < 1990])
summary(temp9m.regression.period2)$adj.r.squared
mse.temp9m.period2 <- mean(residuals(temp9m.regression.period2)^2); rmse.temp9m.period2 <- sqrt(mse.temp9m.period2); rmse.temp9m.period2

temp1m.regression.period3 <- lm(mod2.match$obs.Temp1m[mod2.match$Year > 1989] ~ mod2.match$mod.Temp1m[mod2.match$Year > 1989])
summary(temp1m.regression.period3)$adj.r.squared
mse.temp1m.period3 <- mean(residuals(temp1m.regression.period3)^2); rmse.temp1m.period3 <- sqrt(mse.temp1m.period3); rmse.temp1m.period3

temp4m.regression.period3 <- lm(mod2.match$obs.Temp4m[mod2.match$Year > 1989] ~ mod2.match$mod.Temp4m[mod2.match$Year > 1989])
summary(temp4m.regression.period3)$adj.r.squared
mse.temp4m.period3 <- mean(residuals(temp4m.regression.period3)^2); rmse.temp4m.period3 <- sqrt(mse.temp4m.period3); rmse.temp4m.period3

temp9m.regression.period3 <- lm(mod2.match$obs.Temp9m[mod2.match$Year > 1989] ~ mod2.match$mod.Temp9m[mod2.match$Year > 1989])
summary(temp9m.regression.period3)$adj.r.squared
mse.temp9m.period3 <- mean(residuals(temp9m.regression.period3)^2); rmse.temp9m.period3 <- sqrt(mse.temp9m.period3); rmse.temp9m.period3

#### Plot ####
tempplot <- ggplot() +
  geom_line(data = mod2, aes(x = date, y = mod.Temp1m, col = "1 m"), size = 0.25) +
  geom_line(data = mod2, aes(x = date, y = mod.Temp4m, col = "4 m"), size = 0.25) +
  geom_line(data = mod2, aes(x = date, y = mod.Temp9m, col = "9 m"), size = 0.25) +  
  geom_point(data = mod2.match, aes(x = date, y = obs.Temp1m, col = "1 m"), pch = 19, size = 0.5) +
  geom_point(data = mod2.match, aes(x = date, y = obs.Temp4m, col = "4 m"), pch = 19, size = 0.5) +
  geom_point(data = mod2.match, aes(x = date, y = obs.Temp9m, col = "9 m"), pch = 19, size = 0.5) +
  ylab(expression("Temp " ( degree*C))) +
  xlab(" ") +
  scale_colour_manual("", breaks = c("1 m", "4 m", "9 m"), values = c("#f99d15ff", "#d14a42ff", "#240c4cff")) +
  #geom_vline(xintercept = 1975-01-01, lty = 5) +
  geom_vline(xintercept = 1990-01-01, lty = 5) +
  theme(legend.position = "top")
print(tempplot)

#### PP ====
# Fit Metrics

# Plot 

#### Cumulative PP ====
# Fit Metrics

# Plot 

#### Model PP residuals ====
# Fit Metrics

# Plot 

#### TDP ====
# Fit Metrics

# Plot 

#### O2 ====
# Fit Metrics

# Plot 

#### DOC ====
# Fit Metrics

# Plot 
DOCplot <- ggplot() +
  geom_line(data = mod, aes(x = date, y = mod.DOC, col = "Modeled"), size = 0.25) +
  geom_point(data = mod.match, aes(x = date, y = obs.DOC, col = "Observed"), pch = 19, size = 0.5) +
  ylab(expression(DOC ~ (mg / L))) +
  xlab(" ") +
  scale_colour_manual("", breaks = c("Observed", "Modeled"), values = c("#d14a42ff", "#240c4cff")) +
  theme(legend.position = "top") 
print(DOCplot)

#### Combined Plot ====
grid.arrange(icedateplot, tempplot, DOCplot, ncol = 1)

#### Target Diagram ====


