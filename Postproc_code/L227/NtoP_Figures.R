setwd("/Users/krsalkgu/Documents/SourceTree/L227/Postproc_code/L227")

NtoP <- read.csv("NP_Stoichiometry.csv")
attach(NtoP)
head(NtoP)

#Change Ratios from mass to molar
TNtoTP <- TN.TP.1*30.97/14.01
PNtoPP <- PN.PP.1*30.97/14.01
DINtoTDP <- DIN.TDP.1*30.97/14.01
TDNtoTDP <- TDN.TDP.1*30.97/14.01

NtoPmolar <- data.frame(Start.Date, TNtoTP, PNtoPP, DINtoTDP, TDNtoTDP)

library(ggplot2)
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
          plot.margin=unit(c(10,10,10,10),"pt"))}
theme_set(theme_std())

ggplot(NtoPmolar, aes(x = Start.Date)) +
  geom_point(aes(y = TNtoTP, col="Total"), size = 0.5) +
  geom_point(aes(y = TDNtoTDP, col="Dissolved"), size = 0.5) +
  ylim(0,200) +
  ylab(expression(TN:TP)) +
  xlab(" ") +
  theme(legend.position = c(0.9,0.9)) +
  scale_colour_manual("", breaks = c("Total", "Dissolved"), values = c("black", "red"))
