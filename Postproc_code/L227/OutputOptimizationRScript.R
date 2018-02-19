
setwd("/Users/krsalkgu/Documents/SourceTree/Lake227/Postproc_code/L227")

integratedepi = read.csv("Observed_IntegratedEpi_icefree.csv")
attach(integratedepi)
head(integratedepi)
integratedepi$date = as.Date(integratedepi$date, "%m/%d/%y")
O2 = read.csv("Observed_Oxygen.csv")
attach(O2)
head(O2)
O2$date = as.Date(O2$date, "%m/%d/%y")

library(tidyverse)

OutputForOptimization = merge(integratedepi, O2, by = "date", all=T)
OutputForOptimization_Simple = select(OutputForOptimization, date, Obs_TP:Obs_PP, obs.O2.4m:obs.O2.10m)
OutputForOptimization_Dates = mutate(OutputForOptimization_Simple, 
                                     Year = as.numeric(format(date, format="%Y")), 
                                     Month = as.numeric(format(date, format="%m")), 
                                     Day = as.numeric(format(date, format="%d")))

write.csv(OutputForOptimization_Dates, file = "OutputForOptimization.csv", row.names = F, na = "")
