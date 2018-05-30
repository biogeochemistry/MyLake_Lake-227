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
  xlab(" ") +
  scale_colour_manual("", breaks = c("Ice Break", "Ice Freeze"), values = c("#f99d15ff", "#240c4cff")) +
  geom_vline(xintercept = 1975, lty = 5) +
  geom_vline(xintercept = 1990, lty = 5) +
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
  geom_vline(xintercept = as.numeric(as.Date("1975-01-01")), lty = 5) +
  geom_vline(xintercept = as.numeric(as.Date("1990-01-01")), lty = 5) +
  theme(legend.position = "top")
print(tempplot)

#### PP ====
#### Fit Metrics ####
PP.regression.period1 <- lm(mod.match$obs.PP[mod.match$Year < 1975] ~ mod.match$mod.PP[mod.match$Year < 1975])
summary(PP.regression.period1)$adj.r.squared
mse.PP.period1 <- mean(residuals(PP.regression.period1)^2); rmse.PP.period1 <- sqrt(mse.PP.period1); rmse.PP.period1

PP.regression.period2 <- lm(mod.match$obs.PP[mod.match$Year > 1974 & mod.match$Year < 1990] ~ mod.match$mod.PP[mod.match$Year > 1974 & mod.match$Year < 1990])
summary(PP.regression.period2)$adj.r.squared
mse.PP.period2 <- mean(residuals(PP.regression.period2)^2); rmse.PP.period2 <- sqrt(mse.PP.period2); rmse.PP.period2

PP.regression.period3 <- lm(mod.match$obs.PP[mod.match$Year > 1989] ~ mod.match$mod.PP[mod.match$Year > 1989])
summary(PP.regression.period3)$adj.r.squared
mse.PP.period3 <- mean(residuals(PP.regression.period3)^2); rmse.PP.period3 <- sqrt(mse.PP.period3); rmse.PP.period3

#### Plot ####
PPplot <- ggplot() +
  geom_line(data = mod, aes(x = date, y = mod.PP, col = "Modeled"), size = 0.25) +
  geom_point(data = mod.match, aes(x = date, y = obs.PP, col = "Observed"), pch = 19, size = 0.5) +
  ylab(expression(PP ~ (mu*g / L))) +
  xlab(" ") +
  scale_colour_manual("", breaks = c("Observed", "Modeled"), values = c("#d14a42ff", "#240c4cff")) +
  geom_vline(xintercept = as.numeric(as.Date("1975-01-01")), lty = 5) +
  geom_vline(xintercept = as.numeric(as.Date("1990-01-01")), lty = 5) +
  theme(legend.position = "top") 
print(PPplot)

#### Cumulative PP ====
#### Fit Metrics ####
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

#### Plot ####
PPcumulativeplot <- ggplot(data = mod.match.cumulative) +
  geom_point(aes(x = date, y = Cum.Sum.mod.PP, col = "Modeled"), size = 0.5) +
  geom_point(aes(x = date, y = Cum.Sum.obs.PP, col = "Observed"), size = 0.1, alpha = 0.5) +
  ylab(expression(Cumulative ~ PP ~ (mu*g / L))) +
  xlab(" ") +
  scale_colour_manual("", breaks = c("Observed", "Modeled"), values = c("#d14a42ff", "#240c4cff")) +
  geom_vline(xintercept = as.numeric(as.Date("1975-01-01")), lty = 5) +
  geom_vline(xintercept = as.numeric(as.Date("1990-01-01")), lty = 5) +
  theme(legend.position = "top") 
print(PPcumulativeplot)

#### Model PP residuals ====
#### Fit Metrics ####
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

#### Plot ####
PPresidualsplot <- ggplot() +
  geom_point(data = PPresiduals, aes(x = date, y = residuals.PP), size = 0.5, color = "#d14a42ff") + 
  ylab(expression(PP ~ residuals ~ (mu*g / L))) +
  xlab(" ") +
  geom_vline(xintercept = as.numeric(as.Date("1975-01-01")), lty = 5) +
  geom_vline(xintercept = as.numeric(as.Date("1990-01-01")), lty = 5) +
  theme(legend.position = "none") 
print(PPresidualsplot)

#### TDP ====
#### Fit Metrics ####
TDP.regression.period1 <- lm(mod.match$obs.TDP[mod.match$Year < 1975] ~ mod.match$mod.TDP[mod.match$Year < 1975])
summary(TDP.regression.period1)$adj.r.squared
mse.TDP.period1 <- mean(residuals(TDP.regression.period1)^2); rmse.TDP.period1 <- sqrt(mse.TDP.period1); rmse.TDP.period1

TDP.regression.period2 <- lm(mod.match$obs.TDP[mod.match$Year > 1974 & mod.match$Year < 1990] ~ mod.match$mod.TDP[mod.match$Year > 1974 & mod.match$Year < 1990])
summary(TDP.regression.period2)$adj.r.squared
mse.TDP.period2 <- mean(residuals(TDP.regression.period2)^2); rmse.TDP.period2 <- sqrt(mse.TDP.period2); rmse.TDP.period2

TDP.regression.period3 <- lm(mod.match$obs.TDP[mod.match$Year > 1989] ~ mod.match$mod.TDP[mod.match$Year > 1989])
summary(TDP.regression.period3)$adj.r.squared
mse.TDP.period3 <- mean(residuals(TDP.regression.period3)^2); rmse.TDP.period3 <- sqrt(mse.TDP.period3); rmse.TDP.period3

#### Plot ####
TDPplot <- ggplot() +
  geom_line(data = mod, aes(x = date, y = mod.TDP, col = "Modeled"), size = 0.25) +
  geom_point(data = mod.match, aes(x = date, y = obs.TDP, col = "Observed"), pch = 19, size = 0.5) +
  ylab(expression(TDP ~ (mu*g / L))) +
  xlab(" ") +
  scale_colour_manual("", breaks = c("Observed", "Modeled"), values = c("#d14a42ff", "#240c4cff")) +
  geom_vline(xintercept = as.numeric(as.Date("1975-01-01")), lty = 5) +
  geom_vline(xintercept = as.numeric(as.Date("1990-01-01")), lty = 5) +
  theme(legend.position = "top") 
print(TDPplot)

#### O2 ====
#### Fit Metrics ####
O2.regression.period1 <- lm(mod3.match$obs.O2.4m[mod3.match$Year < 1975] ~ mod3.match$mod.Oxy4m[mod3.match$Year < 1975])
summary(O2.regression.period1)$adj.r.squared
mse.O2.period1 <- mean(residuals(O2.regression.period1)^2); rmse.O2.period1 <- sqrt(mse.O2.period1); rmse.O2.period1

O2.regression.period2 <- lm(mod3.match$obs.O2.4m[mod3.match$Year > 1974 & mod3.match$Year < 1990] ~ mod3.match$mod.Oxy4m[mod3.match$Year > 1974 & mod3.match$Year < 1990])
summary(O2.regression.period2)$adj.r.squared
mse.O2.period2 <- mean(residuals(O2.regression.period2)^2); rmse.O2.period2 <- sqrt(mse.O2.period2); rmse.O2.period2

O2.regression.period3 <- lm(mod3.match$obs.O2.4m[mod3.match$Year > 1989] ~ mod3.match$mod.Oxy4m[mod3.match$Year > 1989])
summary(O2.regression.period3)$adj.r.squared
mse.O2.period3 <- mean(residuals(O2.regression.period3)^2); rmse.O2.period3 <- sqrt(mse.O2.period3); rmse.O2.period3

#### Plot ####
O2plot <- ggplot() +
  geom_line(data = mod2, aes(x = date, y = mod.Oxy4m, col = "Modeled"), size = 0.25) +
  geom_point(data = mod3.match, aes(x = date, y = obs.O2.4m, col = "Observed"), pch = 19, size = 0.5) +
  ylab(expression(DO ~ (mg / L))) +
  xlab(" ") +
  scale_colour_manual("", breaks = c("Observed", "Modeled"), values = c("#d14a42ff", "#240c4cff")) +
  geom_vline(xintercept = as.numeric(as.Date("1975-01-01")), lty = 5) +
  geom_vline(xintercept = as.numeric(as.Date("1990-01-01")), lty = 5) +
  theme(legend.position = "top") 
print(O2plot)

#### DOC ====
#### Fit Metrics ####
DOC.regression.period1 <- lm(mod.match$obs.DOC[mod.match$Year < 1975] ~ mod.match$mod.DOC[mod.match$Year < 1975])
summary(DOC.regression.period1)$adj.r.squared
mse.DOC.period1 <- mean(residuals(DOC.regression.period1)^2); rmse.DOC.period1 <- sqrt(mse.DOC.period1); rmse.DOC.period1

DOC.regression.period2 <- lm(mod.match$obs.DOC[mod.match$Year > 1974 & mod.match$Year < 1990] ~ mod.match$mod.DOC[mod.match$Year > 1974 & mod.match$Year < 1990])
summary(DOC.regression.period2)$adj.r.squared
mse.DOC.period2 <- mean(residuals(DOC.regression.period2)^2); rmse.DOC.period2 <- sqrt(mse.DOC.period2); rmse.DOC.period2

DOC.regression.period3 <- lm(mod.match$obs.DOC[mod.match$Year > 1989] ~ mod.match$mod.DOC[mod.match$Year > 1989])
summary(DOC.regression.period3)$adj.r.squared
mse.DOC.period3 <- mean(residuals(DOC.regression.period3)^2); rmse.DOC.period3 <- sqrt(mse.DOC.period3); rmse.DOC.period3

#### Plot ####
DOCplot <- ggplot() +
  geom_line(data = mod, aes(x = date, y = mod.DOC, col = "Modeled"), size = 0.25) +
  geom_point(data = mod.match, aes(x = date, y = obs.DOC, col = "Observed"), pch = 19, size = 0.5) +
  ylab(expression(DOC ~ (mg / L))) +
  xlab(" ") +
  scale_colour_manual("", breaks = c("Observed", "Modeled"), values = c("#d14a42ff", "#240c4cff")) +
  geom_vline(xintercept = as.numeric(as.Date("1975-01-01")), lty = 5) +
  geom_vline(xintercept = as.numeric(as.Date("1990-01-01")), lty = 5) +
  theme(legend.position = "top") 
print(DOCplot)

#### Combined Plot ====
grid.arrange(icedateplot, tempplot, TDPplot, PPplot, PPcumulativeplot, PPresidualsplot,  ncol = 2)

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

#### Target Plot ####
TargetPlot <- 
ggplot(TargetDiagramData, aes(x = Normalized.Unbiased.RMSD, y = Normalized.Bias, shape = Variable, color = Period)) + 
  geom_point() + 
  annotate("path", x=0+1*cos(seq(0,2*pi,length.out=100)), y=0+1*sin(seq(0,2*pi,length.out=100))) +
  annotate("path", x=0+0.75*cos(seq(0,2*pi,length.out=100)), y=0+0.75*sin(seq(0,2*pi,length.out=100))) +
  xlim(-15, 15) +
  ylim(-2, 2) +
  scale_color_manual(values = c("#f99d15ff", "#d14a42ff", "#240c4cff")) +
  scale_shape_manual(values = c(0, 1, 2, 5, 15, 16, 17))
print(TargetPlot)
