
library(tidyverse)
InflowNP <- read.csv("InflowNP.csv")
TPbyyear <- aggregate(InflowNP$Inflow_TP_mgd, by = list(Year = InflowNP$Year), sum)
NO3byyear <- aggregate(InflowNP$Inflow_NO3_mgd, by = list(Year = InflowNP$Year), sum)
NPbyyear <- inner_join(TPbyyear, NO3byyear, by = "Year")
colnames(NPbyyear) <- c("Year", "TP", "NO3")

TPloading <- lm(NPbyyear$TP ~ NPbyyear$Year)
summary(TPloading)

NO3loading <- lm(NPbyyear$NO3 ~ NPbyyear$Year)
summary(NO3loading)

write.csv(NPbyyear, file = "InflowNPbyyear.csv", row.names = F)
