
InputData <- read.csv("Inflow_Interpolated.csv")
library(tidyverse)

InputData <- InputData %>% unite(date, Year, Month, Day, sep = '-')
InputData$date <- as.Date(InputData$date, format = "%Y-%m-%d")
InputData$date <- format(InputData$date, "%j")

MeanInputData <- aggregate(InputData, list(JulianDay = InputData$date), mean)
SDInputData <- aggregate(InputData, list(JulianDay = InputData$date), sd)

JulianDay <- as.data.frame(c(1:366, 1:366, 1:366, 1:366, 1:366, 1:366, 1:366, 1:366, 1:366, 1:366))

MeanInputDataCloud <- as.numeric(MeanInputData$CloudI)
SDInputDataCloud <- as.numeric(SDInputData$CloudI)
RandoCloud <- as.data.frame(rnorm(3660, mean = MeanInputDataCloud, sd = SDInputDataCloud))

MeanInputDataAirTemp <- as.numeric(MeanInputData$AirTempI)
SDInputDataAirTemp <- as.numeric(SDInputData$AirTempI)
RandoAirTemp <- as.data.frame(rnorm(3660, mean = MeanInputDataAirTemp, sd = SDInputDataAirTemp))

MeanInputDataHumidity <- as.numeric(MeanInputData$HumidityI)
SDInputDataHumidity <- as.numeric(SDInputData$HumidityI)
RandoHumidity <- as.data.frame(rnorm(3660, mean = MeanInputDataHumidity, sd = SDInputDataHumidity))

MeanInputDataAirPressure <- as.numeric(MeanInputData$AirPressureI)
SDInputDataAirPressure <- as.numeric(SDInputData$AirPressureI)
RandoAirPressure <- as.data.frame(rnorm(3660, mean = MeanInputDataAirPressure, sd = SDInputDataAirPressure))

MeanInputDataWindSpeed<- as.numeric(MeanInputData$WindSpeedI)
SDInputDataWindSpeed <- as.numeric(SDInputData$WindSpeedI)
RandoWindSpeed <- as.data.frame(rnorm(3660, mean = MeanInputDataWindSpeed, sd = SDInputDataWindSpeed))

MeanInputDataPrecipitation <- as.numeric(MeanInputData$Precipitation)
SDInputDataPrecipitation <- as.numeric(SDInputData$Precipitation)
RandoPrecipitation <- as.data.frame(rnorm(3660, mean = MeanInputDataPrecipitation, sd = SDInputDataPrecipitation))

MeanInputDataInflowTemp<- as.numeric(MeanInputData$InflowTempI)
SDInputDataInflowTemp <- as.numeric(SDInputData$InflowTempI)
RandoInflowTemp <- as.data.frame(rnorm(3660, mean = MeanInputDataInflowTemp, sd = SDInputDataInflowTemp))

MeanInputDataTP<- as.numeric(MeanInputData$TPI)
SDInputDataTP <- as.numeric(SDInputData$TPI)
RandoTP<- as.data.frame(rnorm(3660, mean = MeanInputDataTP, sd = SDInputDataTP))

MeanInputDataDOC<- as.numeric(MeanInputData$DOCI)
SDInputDataDOC <- as.numeric(SDInputData$DOCI)
RandoDOC <- as.data.frame(rnorm(3660, mean = MeanInputDataDOC, sd = SDInputDataDOC))

MeanInputDataNO3<- as.numeric(MeanInputData$NO3I)
SDInputDataNO3 <- as.numeric(SDInputData$NO3I)
RandoNO3 <- as.data.frame(rnorm(3660, mean = MeanInputDataNO3, sd = SDInputDataNO3))

MeanInputDataNH4<- as.numeric(MeanInputData$NH4I)
SDInputDataNH4 <- as.numeric(SDInputData$NH4I)
RandoNH4 <- as.data.frame(rnorm(3660, mean = MeanInputDataNH4, sd = SDInputDataNH4))

MeanInputDataSO4<- as.numeric(MeanInputData$SO4I)
SDInputDataSO4 <- as.numeric(SDInputData$SO4I)
RandoSO4 <- as.data.frame(rnorm(3660, mean = MeanInputDataSO4, sd = SDInputDataSO4))

MeanInputDataFe2<- as.numeric(MeanInputData$Fe2I)
SDInputDataFe2 <- as.numeric(SDInputData$Fe2I)
RandoFe2 <- as.data.frame(rnorm(3660, mean = MeanInputDataFe2, sd = SDInputDataFe2))

MeanInputDataCa2<- as.numeric(MeanInputData$Ca2I)
SDInputDataCa2 <- as.numeric(SDInputData$Ca2I)
RandoCa2 <- as.data.frame(rnorm(3660, mean = MeanInputDataCa2, sd = SDInputDataCa2))

MeanInputDatapH<- as.numeric(MeanInputData$pHI)
SDInputDatapH <- as.numeric(SDInputData$pHI)
RandopH <- as.data.frame(rnorm(3660, mean = MeanInputDatapH, sd = SDInputDatapH))

MeanInputDataFe3 <- as.numeric(MeanInputData$Fe3I)
SDInputDataFe3 <- as.numeric(SDInputData$Fe3I)
RandoFe3 <- as.data.frame(rnorm(3660, mean = MeanInputDataFe3, sd = SDInputDataFe3))

MeanInputDataSiO2<- as.numeric(MeanInputData$SiO2I)
SDInputDataSiO2 <- as.numeric(SDInputData$SiO2I)
RandoSiO2 <- as.data.frame(rnorm(3660, mean = MeanInputDataSiO2, sd = SDInputDataSiO2))

Rando <- as.data.frame(cbind(JulianDay, RandoCloud, RandoAirTemp, RandoHumidity, 
                             RandoAirPressure, RandoWindSpeed, RandoPrecipitation, RandoInflowTemp, 
                             RandoTP, RandoDOC, RandoNO3, RandoNH4, RandoSO4, 
                             RandoFe2, RandoCa2, RandopH, RandoFe3, RandoSiO2))
colnames(Rando) <- c("JulianDay", "Cloud", "AirTemp", "Humidity", 
                     "AirPressure", "WindSpeed", "Precipitation", "InflowTemp", 
                     "TP", "DOC", "NO3", "NH4", "SO4", "Fe2", "Ca2", "pH", "Fe3", "SiO2")

# Change negative values to zero (exception: air temperature can be negative)
Rando$Cloud <- ifelse(Rando$Cloud < 0, 0, Rando$Cloud)
Rando$InflowTemp <- ifelse(Rando$InflowTemp < 0, 0, Rando$InflowTemp)
Rando$WindSpeed <- ifelse(Rando$WindSpeed < 0, 0, Rando$WindSpeed)
Rando$Precipitation <- ifelse(Rando$Precipitation < 0, 0, Rando$Precipitation)
Rando$TP <- ifelse(Rando$TP < 0, 0, Rando$TP)
Rando$DOC <- ifelse(Rando$DOC < 0, 0, Rando$DOC)
Rando$NO3 <- ifelse(Rando$NO3 < 0, 0, Rando$NO3)
Rando$NH4 <- ifelse(Rando$NH4 < 0, 0, Rando$NH4)
Rando$SO4 <- ifelse(Rando$SO4 < 0, 0, Rando$SO4)
Rando$Fe2 <- ifelse(Rando$Fe2 < 0, 0, Rando$Fe2)
Rando$Ca2 <- ifelse(Rando$Ca2 < 0, 0, Rando$Ca2)
Rando$Fe3 <- ifelse(Rando$Fe3 < 0, 0, Rando$Fe3)

# Check distributions
ggplot(Rando, aes(x = JulianDay, y = Precipitation)) + 
  geom_point()

# Generate correct number of samples for 2017-2026 (keep day 366 for 2020 and 2024 for leap day)
Rando <- Rando[-c(366, 732, 1098, 1830, 2196, 2562, 3294, 3660),]

# Write csv for random input data
write.csv(Rando, file = "RandomFutureInputs.csv", row.names = F)

# Add P fertilization
Pfertilization <- as.data.frame(c(1:366), c(1:366))

# Write csv for 