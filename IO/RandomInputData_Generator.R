
InputData <- read.csv("Inflow_Interpolated.csv")
library(tidyverse)

InputData <- InputData %>% unite(date, Year, Month, Day, sep = '-')
InputData$date <- as.Date(InputData$date, format = "%Y-%m-%d")
InputData$date <- format(InputData$date, "%j")

#need to account for Julian day = 60 (February 29 or March 1)

MeanInputData <- aggregate(InputData, list(JulianDay = InputData$date), mean)
SDInputData <- aggregate(InputData, list(JulianDay = InputData$date), sd)

# RandomInputData <- matrix(NA, 366, 18)
# for(i in 1:366){
#   for(j in 1:18){
#     #mu <- MeanInputData[i, j]
#     #sigma <- SDInputData[i, j]
#     RandomInputData[i, j] <- rnorm((366*18, 5, 1), 366, 18)
#   }
# }
# 
# mu<-MeanInputData
# sigma<-SDInputData
# sample.size<-100
# norm.mat<-mapply(function(x,y){rnorm(x,y,n=sample.size)},x=mu,y=sigma)

MeanInputDataCloud <- as.numeric(MeanInputData$CloudI)
SDInputDataCloud <- as.numeric(SDInputData$CloudI)
RandoCloud <- as.data.frame(rnorm(3660, mean = MeanInputDataCloud, sd = SDInputDataCloud))

MeanInputData <- data.matrix(MeanInputData)
SDInputData <- data.matrix(SDInputData)
Rando <- as.data.frame(rnorm(MeanInputData, mean = MeanInputData, sd = SDInputData))
