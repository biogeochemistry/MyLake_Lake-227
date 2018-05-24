#Testing Ice-on, Ice-off

setwd("/Users/krsalkgu/Documents/SourceTree/Lake227/Postproc_code/L227")

library(ggplot2)
library(knitr)
library(tidyverse)
knitr::opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE)

theme_std <- function (base_size = 16, base_family = "") {
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

#### Modify ice observation spreadsheet ----
obs.ice <- read.csv("Lake_239_ice-on_ice-off.csv", header = T)
attach(obs.ice)
obs.ice$Ice.Off.Date <- as.Date(obs.ice$Ice.Off.Date, format = "%m/%d/%y")
obs.ice$Ice.On.Date <- as.Date(obs.ice$Ice.On.Date, format = "%m/%d/%y")
obs.daybreak <- strftime(obs.ice$Ice.Off.Date, format = "%j")
obs.dayfreeze <- strftime(obs.ice$Ice.On.Date, format = "%j")
obs.year <- format(obs.ice$Ice.On.Date, "%Y")
obs.ice.edit <- cbind(obs.year, obs.ice[4], obs.ice[5], obs.daybreak, obs.dayfreeze)
colnames(obs.ice.edit) <- c("Year", "Ice.Off.Date", "Ice.On.Date", "obs.daybreak", "obs.dayfreeze")

#### Modify ice model output spreadsheet ----
out.ice <- read.csv("Output_Ice.csv", header = F)
colnames(out.ice) <- c("Year.Break", "Month.Break", "Day.Break", "Year.Freeze", "Month.Freeze", "Day.Freeze")
out.ice.edit <- out.ice %>% unite(Ice.Off.Date, Year.Break, Month.Break, Day.Break, sep = '-')
out.ice.edit <- out.ice.edit %>% unite(Ice.On.Date, Year.Freeze, Month.Freeze, Day.Freeze, sep = '-')
out.ice.edit$Ice.Off.Date <- as.Date(out.ice.edit$Ice.Off.Date, format = "%Y-%m-%d")
out.ice.edit$Ice.On.Date <- as.Date(out.ice.edit$Ice.On.Date, format = "%Y-%m-%d")
out.daybreak <- strftime(out.ice.edit$Ice.Off.Date, format = "%j")
out.dayfreeze <- strftime(out.ice.edit$Ice.On.Date, format = "%j")
out.year <- as.character(out.ice$Year.Freeze)
out.ice.edit2 <- cbind(out.year, out.ice.edit, out.daybreak, out.dayfreeze)
colnames(out.ice.edit2) <- c("Year", "Ice.Off.Date", "Ice.On.Date", "out.daybreak", "out.dayfreeze")

#### Combine observations and model output ----
match.ice <- inner_join(obs.ice.edit, out.ice.edit2, by = "Year") 
match.ice$obs.daybreak <- as.numeric(as.character(match.ice$obs.daybreak))
match.ice$out.daybreak <- as.numeric(as.character(match.ice$out.daybreak))
match.ice$obs.dayfreeze <- as.numeric(as.character(match.ice$obs.dayfreeze))
match.ice$out.dayfreeze <- as.numeric(as.character(match.ice$out.dayfreeze))
match.ice$Year <- as.numeric(match.ice$Year)

#### Calculate fit metrics
icebreakregression <- lm (match.ice$obs.daybreak ~ match.ice$out.daybreak)
summary(icebreakregression)
msebreak <- mean(residuals(icebreakregression)^2)
rmsebreak <- sqrt(msebreak)
rmsebreak
rmsebreak <- round(rmsebreak, digits = 2)

icefreezeregression <- lm (match.ice$obs.dayfreeze ~ match.ice$out.dayfreeze)
summary(icefreezeregression)
msefreeze <- mean(residuals(icefreezeregression)^2)
rmsefreeze <- sqrt(msefreeze)
rmsefreeze
rmsefreeze <- round(rmsefreeze, digits = 2)

#### Plots ----
icedateplot <- 
  ggplot(data = match.ice, aes(x = Year, group = 1)) +
    geom_point(aes(y = obs.daybreak), color = "red") +
    geom_line(aes(y = out.daybreak), color = "red") +
    geom_point(aes(y = obs.dayfreeze), color = "blue") +
    geom_line(aes(y = out.dayfreeze), color = "blue") + 
    scale_y_reverse() + 
    ylab("Julian Day") + 
    annotate("text", x = 1969, y = 280, label = "Ice Freeze", color = "blue", size = 5, hjust = 0) +
    annotate("text", x = 1969, y = 150, label = "Ice Break", color = "red", size = 5, hjust = 0) +
    annotate("text", x = 1974, y = 170, label = rmsebreak, color = "red", size = 5, hjust = 0) +
    annotate("text", x = 1969, y = 170, label = "RMSE =", color = "red", size = 5, hjust = 0) +
    annotate("text", x = 1974, y = 300, label = rmsefreeze, color = "blue", size = 5, hjust = 0) +
    annotate("text", x = 1969, y = 300, label = "RMSE =", color = "blue", size = 5, hjust = 0)
  
icebreakregressionplot <- 
  ggplot(data = match.ice) + 
    geom_point(aes(x = obs.daybreak, y = out.daybreak), color = "red") +
    geom_smooth(method = "lm", aes(x = obs.daybreak, y = out.daybreak), se = FALSE, color = "red")

icefreezeregressionplot <- 
  ggplot(data = match.ice) + 
    geom_point(aes(x = obs.dayfreeze, y = out.dayfreeze), color = "blue") +
    geom_smooth(method = "lm", aes(x = obs.dayfreeze, y = out.dayfreeze), se = FALSE, color = "blue") 
