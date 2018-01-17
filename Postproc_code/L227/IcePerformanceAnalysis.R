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


obs.ice <- read.csv("Lake_239_ice-on_ice-off.csv", header = T)
attach(obs.ice)
obs.ice$Ice.Off.Date <- as.Date(obs.ice$Ice.Off.Date, format = "%Y-%m-%d")
obs.ice$Ice.On.Date <- as.Date(obs.ice$Ice.On.Date, format = "%Y-%m-%d")
obs.daybreak <- strftime(obs.ice$Ice.Off.Date, format = "%j")
obs.dayfreeze <- strftime(obs.ice$Ice.On.Date, format = "%j")
obs.year <- format(obs.ice$Ice.On.Date, "%Y")
obs.ice.edit <- cbind(obs.year, obs.ice[4], obs.ice[5], obs.daybreak, obs.dayfreeze)

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


@iceplot <- 
  
  ggplot() +
  geom_point(data = obs.ice.edit, aes(x = obs.year, y = obs.daybreak)) +
  #geom_line(data = obs.ice.edit, aes(x = obs.year, y = obs.daybreak)) + 
  geom_point(data = out.ice.edit2, aes(x = out.year, y = out.daybreak), color = "red") + 
  geom_line(aes(x = out.year, y = out.daybreak))
  

