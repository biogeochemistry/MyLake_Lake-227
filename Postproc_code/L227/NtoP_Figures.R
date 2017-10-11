setwd("/Users/krsalkgu/Documents/SourceTree/L227/Postproc_code/L227")

NtoPLake <- read.csv("NP_Stoichiometry_L227.csv")
attach(NtoPLake)
head(NtoPLake)
NtoPInflow <- read.csv("NP_Stoichiometry_Inflow.csv")
attach(NtoPInflow)
head(NtoPInflow)

#Change concentrations from mass (ug/L) to molar (umol/L)
Fert_TP_molar <- Fert_TP/30.97
Inflow_TP_molar <- Inflow_TP/30.97
Fert_TN_molar <- Fert_TN/14.01
Inflow_TN_molar <- Inflow_TN/14.01
NO3_molar <- NO3/14.01
NH4_molar <- NH4/14.01
DIN_molar <- NO3_molar + NH4_molar
PN_molar <- PN/14.01
TDN_molar <- TDN/14.01
PP_molar <- PP/30.97
TDP_molar <- TDP/30.97
TN_molar <- TDN_molar + PN_molar
TP_molar <- TDP_molar + PP_molar

#Combine inflows + fertilization
Input_TN_molar <- Fert_TN_molar + Inflow_TN_molar
Input_TP_molar <- Fert_TP_molar + Inflow_TP_molar

#Create N:P stoichiometric ratios in lake and inputs
Fert_NtoP <- Fert_TN_molar/Fert_TP_molar
Inflow_NtoP <- Inflow_TN_molar/Inflow_TP_molar
Input_NtoP <- Input_TN_molar/Input_TP_molar
DINtoTDP <- DIN_molar/TDP_molar
TDNtoTDP <- TDN_molar/TDP_molar
PNtoPP <- PN_molar/PP_molar
TNtoTP <- TN_molar/TP_molar

#Convert dates from factor to date format
Datelake <- as.Date(NtoPLake$Date, "%m/%d/%y")
Dateinflow <- as.Date(NtoPInflow$Date, "%m/%d/%Y")

#Create data frames for N:P stoichiometry in lake and inflows
NtoPinsitu <- data.frame(Datelake, TNtoTP, PNtoPP, DINtoTDP, TDNtoTDP)
NtoPinput <- data.frame(Dateinflow, Fert_NtoP, Inflow_NtoP, Input_NtoP)

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

library(strucchange)
a = breakpoints (Input_NtoP ~ Dateinflow, data = NtoPinput)
breakdates(a, format.times=T)
TNtoTPbydate <- lm(TNtoTP ~ Datelake)



ggplot(NtoPinsitu, aes(x = Datelake)) +
  geom_point(aes(y = TNtoTP), size = 0.5) +
  ylim(0,200) +
  ylab(expression(TN:TP)) +
  xlab(" ") +
  theme(legend.position = c(0.9,0.9))

ggplot(NtoPinsitu, aes(x = Datelake)) +
  geom_point(aes(y = TDNtoTDP), size = 0.5) +
  ylim(0,500) +
  ylab(expression(TDN:TDP)) +
  xlab(" ") +
  theme(legend.position = c(0.9,0.9)) 

ggplot(NtoPinsitu, aes(x = Datelake)) +
  geom_point(aes(y = DINtoTDP), size = 0.5) +
  ylim(0,500) +
  ylab(expression(DIN:TDP)) +
  xlab(" ") +
  theme(legend.position = c(0.9,0.9)) 

ggplot(NtoPinsitu, aes(x = Datelake)) +
  geom_point(aes(y = PNtoPP), size = 0.5) +
  ylim(0,200) +
  ylab(expression(PN:PP)) +
  xlab(" ") +
  theme(legend.position = c(0.9,0.9)) 

ggplot(NtoPinput, aes(x = Dateinflow)) +
  geom_point(aes(y = Fert_NtoP), size = 0.5) +
  ylim(0,50) +
  ylab(expression(FertN:FertP)) +
  xlab(" ") +
  theme(legend.position = c(0.9,0.9)) 

ggplot(NtoPinput, aes(x = Dateinflow)) +
  geom_point(aes(y = Inflow_NtoP), size = 0.5) +
  ylab(expression(InflowN:InflowP)) +
  xlab(" ") +
  theme(legend.position = c(0.9,0.9)) 

ggplot(NtoPinput, aes(x = Dateinflow)) +
  geom_point(aes(y = Input_NtoP), size = 0.5) +
  ylim(0,50) +
  ylab(expression(InputN:InputP)) +
  xlab(" ") +
  theme(legend.position = c(0.9,0.9)) 
