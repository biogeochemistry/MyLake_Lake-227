


integratedepi = read.csv("Observed_IntegratedEpi_icefree_lateMay.csv")
attach(integratedepi)
head(integratedepi)
integratedepi$date = as.Date(integratedepi$date, "%m/%d/%y")
O2 = read.csv("Observed_Oxygen.csv")
attach(O2)
head(O2)
O2$date = as.Date(O2$date, "%m/%d/%y")

phytopercent <- read.csv("../PhytoCommunity/Phytoplankton/L227_BiomassNonDiazo.csv")
phytopercent$start.date <- as.Date(phytopercent$start.date, "%Y-%m-%d")
phytopercent <- select(phytopercent, "start.date", "PropBiomass")
colnames(phytopercent)[1] <- "date"

library(tidyverse)
library(fuzzyjoin)

OutputForOptimization = merge(integratedepi, O2, by = "date", all=T)
OutputForOptimization_Simple = select(OutputForOptimization, date, Obs_TP:Obs_PP, obs.O2.4m:obs.O2.10m)
OutputForOptimization_Dates = mutate(OutputForOptimization_Simple, 
                                     Year = as.numeric(format(date, format="%Y")), 
                                     Month = as.numeric(format(date, format="%m")), 
                                     Day = as.numeric(format(date, format="%d")))

OutputForOptimization_Dates <- interval_left_join(OutputForOptimization_Dates, phytopercent, by = c("date", "date"), maxgap = 6)
colnames(OutputForOptimization_Dates)[1] <- "date"
colnames(OutputForOptimization_Dates)[13] <- "phyto.date"

OutputForOptimization_Dates$PropBiomass[OutputForOptimization_Dates$date < as.Date("1975-01-01")] <- 1

OutputForOptimization_Dates <- OutputForOptimization_Dates %>%
  mutate(diazoPP = Obs_PP * PropBiomass,
         nondiazoPP = Obs_PP * (1-PropBiomass))

write.csv(OutputForOptimization_Dates, file = "OutputForOptimization.csv", row.names = F, na = "")
